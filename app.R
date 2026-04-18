# =============================================================================
# SCF Fire Calibration Shiny App
# Calibrates ignition and spread parameters for LANDIS-II Social-Climate-Fire
#
# Spatial design:
#   calibration_vect  = large regional boundary (Cascades, Coast Range, etc.)
#                       Used to clip ERA data, FPA FOD, and GeoMAC perimeters.
#   template_r        = LANDIS grid raster (smaller park-specific extent)
#                       Used only for ignition surface generation and
#                       perimeter rasterization.
#
# All spatial objects are normalized to a single working CRS (taken from the
# template raster) before any processing. Leaflet display always reprojects
# to EPSG:4326 on the fly.
#
# Data requirements:
#   ERA-Interim FWI indices (FFMC, DMC, DC, FWI):
#     Vitolo et al. (2019) -- https://doi.org/10.5281/zenodo.1065400
#     Download: ffmc.nc, dmc.nc, dc.nc, fwi.nc
#
#   ERA5 wind components (u10, v10):
#     Copernicus CDS -- https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
#     Variables: 10m_u_component_of_wind, 10m_v_component_of_wind
#     Time: daily at 12:00 UTC, 1980-2018, NetCDF format
#
#   Fire occurrence (ignitions):
#     FPA FOD (Short 2022) -- https://www.fs.usda.gov/rds/archive/Catalog/RDS-2013-0009.6
#     File: FPA_FOD_20221014.gdb
#
#   Fire perimeters (spread):
#     GeoMAC Historic Perimeters 2000-2018
#     File: Historic_Geomac_Perimeters_All_Years_2000_2018.gdb
#
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  shiny, bslib, terra, dplyr, lubridate, readr, tidyr,
  ggplot2, glue, DT, tmap, purrr, stringr, broom, sf, scales,
  climateR, elevatr
)

source("modules/era_helpers.R")
source("modules/ignition_helpers.R")
source("modules/spread_helpers.R")
source("modules/species_helpers.R")
source("modules/ui_helpers.R")
options(jsonlite.named_vectors_as_objects = FALSE)

# =============================================================================
# UI
# =============================================================================

ui <- page_navbar(
  title = "SCF Fire Calibration",
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2c5f2e",
    secondary  = "#a0522d"
  ),
  id = "main_nav",

  # ---------------------------------------------------------------------------
  # Sidebar (shared across tabs)
  # ---------------------------------------------------------------------------
  sidebar = sidebar(
    width = 320,
    title = "Landscape Setup",

    # -- Calibration boundary (large region)
    h6("1. Calibration Boundary Shapefile"),
    p(class = "text-muted small",
      "Large regional boundary used to clip ERA data, FPA FOD, and GeoMAC ",
      "perimeters (e.g. WA Cascades + Coast Range). The .dbf, .shx, and .prj ",
      "files must be in the same directory as the .shp."),
    textInput("shp_path", label = NULL,
              value       = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/S_USA.ECOSYS_ECOMAPPROVINCES_2025/ECOSYS_CoastalRange.shp",
              placeholder = "/path/to/calibration_region.shp"),
    actionButton("load_shp", "Load Shapefile", class = "btn-sm btn-secondary w-100"),
    uiOutput("shp_status"),

    hr(),

    # -- Template raster (LANDIS grid)
    h6("2. LANDIS Template Raster"),
    p(class = "text-muted small",
      "Park-specific LANDIS grid raster (e.g. climate regions .tif). Used for ",
      "ignition surface generation and perimeter rasterization. Can be a ",
      "different CRS and extent than the calibration boundary — the app will ",
      "reproject everything to match this grid."),
    textInput("tif_path", label = NULL,
              value       = "/Users/jlamping/Desktop/LANDIS/NOCA/CoolDry/Rep1/input_maps/NOCA_ClimateRegions_09122025.tif",
              placeholder = "/path/to/template.tif"),
    actionButton("load_tif", "Load Raster", class = "btn-sm btn-secondary w-100"),
    uiOutput("tif_status"),

    hr(),

    # -- ERA-Interim NetCDF files
    h6("3. ERA-Interim FWI Files"),
    p(class = "text-muted small",
      tags$a("Download from Zenodo", href = "https://doi.org/10.5281/zenodo.1065400",
             target = "_blank"),
      " — ffmc.nc, dmc.nc, dc.nc, fwi.nc"),
    textInput("era_fwi_path",  "fwi.nc",
              value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/fwi.nc"),
    textInput("era_ffmc_path", "ffmc.nc",
              value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/ffmc.nc"),
    textInput("era_dmc_path",  "dmc.nc",
              value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/dmc.nc"),
    textInput("era_dc_path",   "dc.nc",
              value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/dc.nc"),
    actionButton("load_era_fwi", "Load FWI Files", class = "btn-sm btn-secondary w-100"),
    uiOutput("era_fwi_status"),

    hr(),

    # -- GRIDMET wind (remote fetch, no local files needed)
    h6("4. GRIDMET Wind (Remote)"),
    p(class = "text-muted small",
      "Daily 10m wind speed and direction fetched remotely from ",
      tags$a("GRIDMET", href = "https://www.climatologylab.org/gridmet.html",
             target = "_blank"),
      " via OPeNDAP. CONUS only, ~4 km, 1979–present. No download required — ",
      "load the calibration boundary and set the calibration period first."),
    actionButton("fetch_gridmet_wind", "Fetch GRIDMET Wind",
                 class = "btn-sm btn-secondary w-100",
                 icon  = icon("cloud-arrow-down")),
    uiOutput("era_wind_status"),

    hr(),

    # -- Calibration date range
    h6("5. Calibration Period"),
    sliderInput("cal_years", label = NULL,
                min = 1980, max = 2018,
                value = c(1992, 2018), step = 1, sep = ""),

    hr(),
    uiOutput("landscape_preview_ui")
  ),

  # ---------------------------------------------------------------------------
  # Tab 1: Ignition Calibration
  # ---------------------------------------------------------------------------
  nav_panel(
    title = "Ignition Calibration",
    icon  = icon("fire"),

    layout_columns(
      col_widths = c(3, 9),

      card(
        card_header("Ignition Inputs"),

        h6("FPA FOD Geodatabase"),
        p(class = "text-muted small",
          tags$a("Download FPA FOD",
                 href = "https://www.fs.usda.gov/rds/archive/Catalog/RDS-2013-0009.6",
                 target = "_blank")),
        textInput("fpa_gdb_path", label = NULL,
                  value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Code/calibration/scf_calibration/data_raw/Short2022/Data/FPA_FOD_20221014.gdb"),

        hr(),
        h6("Model Settings"),
        radioButtons("ign_distribution", "Distribution",
                     choices  = c("Zero-Inflated Poisson" = "ZIP",
                                  "Poisson"               = "Poisson"),
                     selected = "ZIP"),

        hr(),
        h6("Ignition Surface (Kernel)"),
        p(class = "text-muted small",
          "A Gaussian kernel is applied to historical ignition point locations ",
          "to produce a smooth probability surface on the LANDIS template grid. ",
          "Higher-weight cells are more likely to be selected as ignition origins ",
          "by SCF each time step."),
        sliderInput("kernel_bw", "Bandwidth (m)", min = 1000, max = 30000,
                    value = 7000, step = 500),
        p(class = "text-muted small",
          tags$b("Bandwidth:"), " Controls the spatial spread (standard deviation) of ",
          "the Gaussian kernel around each historical ignition point. Smaller values ",
          "concentrate ignition probability tightly around observed locations; larger ",
          "values spread probability more broadly across the landscape. Typical range: ",
          "5\u201315 km for regional calibrations."),
        sliderInput("kernel_maxdist", "Max Distance (m)", min = 5000, max = 75000,
                    value = 25000, step = 1000),
        p(class = "text-muted small",
          tags$b("Max Distance:"), " Hard cutoff radius beyond which the kernel weight ",
          "drops to zero, regardless of bandwidth. Prevents ignition probability from ",
          "bleeding into areas with no nearby historical fires. Set to roughly 3\u20134\u00d7 ",
          "the bandwidth to retain the full shape of the Gaussian tail."),

        hr(),
        actionButton("run_ignition", "Run Ignition Calibration",
                     class = "btn-success w-100", icon = icon("play")),
        br(), br(),
        uiOutput("ign_run_status")
      ),

      card(
        card_header("Ignition Outputs"),
        full_screen = TRUE,
        navset_tab(
          nav_panel("Overview Map",
          tmapOutput("ign_map", height = "500px")
                    
          ),
          nav_panel("Coefficients",
            br(),
            p(class = "text-muted",
              "Fitted ignition model coefficients (b0, b1, bz0, bz1) by ignition type."),
            DTOutput("ign_coef_table"),
            br(),
            downloadButton("dl_ign_coef", "Download Coefficients CSV",
                           class = "btn-sm btn-outline-primary")
          ),

          nav_panel("Scale Adjustment",
            br(),
            p(class = "text-muted small",
              "Ignition model coefficients are calibrated using all fires across the ",
              "large calibration boundary. Before using them in LANDIS, the Poisson ",
              "intercept (b0) must be scaled down to match the size of your LANDIS ",
              "landscape. The slope (b1) and zero-inflation terms (bz0, bz1) transfer ",
              "without modification."),
            fluidRow(
              column(5,
                h6("Scaling Method"),
                radioButtons("scale_method", NULL,
                  choices = c(
                    "A \u2014 Area Offset (recommended)" = "area_offset",
                    "B \u2014 Park-Specific Intercept from FPA FOD" = "park_intercept"
                  ),
                  selected = "area_offset"
                ),
                p(class = "text-muted small",
                  strong("Option A:"), " adjusts b0 by log(landscape area / calibration area). ",
                  "Assumes uniform fire density across the calibration region. ",
                  "Calculated automatically when calibration runs."),
                p(class = "text-muted small",
                  strong("Option B:"), " re-estimates b0 from FPA FOD fires within the ",
                  "LANDIS template extent, holding b1 fixed from the regional fit. ",
                  "Requires \u226515 years of fire history within the template. ",
                  "Falls back to Option A if fewer than 3 fires are found."),
                hr(),
                h6("Landscape Areas"),
                uiOutput("area_summary_ui"),
                br(),
                conditionalPanel(
                  "input.scale_method == 'park_intercept'",
                  actionButton("fit_park_b0", "Fit Park-Scale Intercept",
                               class = "btn-sm btn-outline-primary w-100"),
                  br(), br(),
                  uiOutput("park_fit_status_ui")
                )
              ),
              column(7,
                h6("Adjusted Coefficients"),
                p(class = "text-muted small",
                  "Coefficients to enter into the SCF parameter file. ",
                  "b0 has been adjusted; all other values are unchanged."),
                DTOutput("adjusted_coef_table"),
                br(),
                downloadButton("dl_adjusted_coef",
                               "Download Adjusted Coefficients CSV",
                               class = "btn-sm btn-outline-primary")
              )
            ),
            hr(),
            h6("Expected Annual Ignitions at LANDIS Landscape Scale"),
            p(class = "text-muted small",
              "Expected ignitions per year computed by applying model coefficients to ",
              "the ERA FWI daily climatology. Shows regional (unadjusted) and ",
              "landscape-adjusted rates side by side so you can verify the scaling ",
              "before running LANDIS."),
            fluidRow(
              column(5,
                DTOutput("expected_ign_table")
              ),
              column(7,
                plotOutput("expected_ign_plot", height = "480px")
              )
            )
          ),

          nav_panel("Diagnostics",
            br(),
            p(class = "text-muted small",
              "Four diagnostic panels for evaluating the fitted ignition model. ",
              "Run Ignition Calibration first. ZIP model selected in the inputs ",
              "adds the zero-inflation and count-component panels."),
            navset_tab(
              nav_panel("Observed vs Fitted",
                plotOutput("ign_diag_plot", height = "420px")
              ),
              nav_panel("Zero-Inflation (\u03c0)",
                plotOutput("ign_zero_inf_plot", height = "420px")
              ),
              nav_panel("Count Component (\u03bb)",
                plotOutput("ign_lambda_plot", height = "420px")
              ),
              nav_panel("Residuals",
                plotOutput("ign_resid_plot", height = "420px")
              )
            )
          ),
          nav_panel("ERA Startup Codes",
            br(),
            p(class = "text-muted",
              "Recommended FFMC, DMC, DC, FWI startup values derived from ERA ",
              "climatologies over the calibration boundary."),
            uiOutput("startup_date_ui"),
            actionButton("run_startup", "Compute Startup Values",
                         class = "btn-sm btn-outline-primary"),
            br(), br(),
            DTOutput("startup_table"),
            br(),
            downloadButton("dl_startup", "Download Startup Values CSV",
                           class = "btn-sm btn-outline-primary"),
            br(), br(),
            plotOutput("startup_clim_plot", height = "400px")
          ),
          nav_panel("Landscape Maps",
            br(),

            h5("Ignition Allocation Surfaces"),
            p(class = "text-muted small",
              "Spatial weights for fire ignition locations on the LANDIS template ",
              "grid. Exported as 32-bit signed integers (INT4S), scaled 0\u2013100. ",
              "SCF truncates raster values to integer at load time, so a 0\u20131 float ",
              "map collapses to weight 0 or 1 only \u2014 use this integer-scaled output."),
            fluidRow(
              column(4, downloadButton("dl_lightning_map",  "Lightning_Ignition_Map.tif",
                                       class = "btn-sm btn-outline-secondary w-100")),
              column(4, downloadButton("dl_accidental_map", "Accidental_Ignition_Map.tif",
                                       class = "btn-sm btn-outline-secondary w-100")),
              column(4, downloadButton("dl_rx_map",         "Rx_Ignition_Map.tif",
                                       class = "btn-sm btn-outline-secondary w-100"))
            ),

            hr(),
            h5("Terrain Maps"),
            p(class = "text-muted small",
              "Ground slope (degrees) and uphill slope azimuth (0\u2013359\u00b0 from N) ",
              "derived from a SRTM/3DEP DEM fetched via elevatr and resampled to the ",
              "LANDIS template grid. Required by SCF for spread direction physics."),
            actionButton("gen_terrain_maps", "Fetch DEM & Generate Terrain Maps",
                         class = "btn-sm btn-outline-primary", icon = icon("mountain")),
            br(),
            uiOutput("terrain_maps_status"),
            br(),
            fluidRow(
              column(6,
                h6("Ground Slope (degrees)"),
                tmapOutput("ground_slope_preview", height = "260px")
              ),
              column(6,
                h6("Uphill Slope Azimuth (degrees from N)"),
                tmapOutput("uphill_slope_preview", height = "260px")
              )
            ),
            br(),
            fluidRow(
              column(6, downloadButton("dl_ground_slope", "GroundSlope.tif  (INT2S)",
                                       class = "btn-sm btn-outline-secondary w-100")),
              column(6, downloadButton("dl_uphill_slope", "UphillSlope.tif  (INT2S)",
                                       class = "btn-sm btn-outline-secondary w-100"))
            ),

            hr(),
            h5("Suppression Maps"),
            p(class = "text-muted small",
              "Blank suppression maps (0.0 everywhere) for initial calibration. ",
              "Replace with zone-based maps once fire behavior is validated. ",
              "Exported as single-precision float (FLT4S) as required by SCF. ",
              "Requires the LANDIS template to be loaded first."),
            fluidRow(
              column(4, downloadButton("dl_supp_accidental",
                                       "Suppression_Accidental_3Zones.tif",
                                       class = "btn-sm btn-outline-secondary w-100")),
              column(4, downloadButton("dl_supp_lightning",
                                       "Suppression_Lightning_3Zones.tif",
                                       class = "btn-sm btn-outline-secondary w-100")),
              column(4, downloadButton("dl_supp_rx",
                                       "Suppression_Rx_3Zones.tif",
                                       class = "btn-sm btn-outline-secondary w-100"))
            )
          )
        )
      )
    )
  ),

  # ---------------------------------------------------------------------------
  # Tab 2: Species Table
  # ---------------------------------------------------------------------------
  nav_panel(
    title = "Species Table",
    icon  = icon("leaf"),

    layout_columns(
      col_widths = c(4, 8),

      # -- Left: inputs -------------------------------------------------------
      card(
        card_header("Species Table Inputs"),

        # ---- Step 1: Lookup CSV -------------------------------------------
        h6("Step 1 — LANDIS \u2192 FIA Species Lookup CSV"),
        p(class = "text-muted small",
          "A CSV with three columns: ",
          tags$code("SpeciesName"), " (LANDIS species code), ",
          tags$code("FIA_SPCD"), " (integer FIA species code), and ",
          tags$code("SCIENTIFIC_NAME"), ". All species in this file will be ",
          "included in the table. Load it once; it persists for the session."),
        textInput("lookup_csv_path", label = NULL,
                  placeholder = "/path/to/LANDIS_FIA_species_lookup.csv",
                  value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/FIA/LANDIS_FIA_species_lookup.csv"),
        actionButton("load_lookup_csv", "Load Lookup CSV",
                     class = "btn-sm btn-secondary w-100",
                     icon  = icon("table")),
        uiOutput("lookup_load_status"),

        hr(),

        # ---- Step 2: Optional NECN filter or manual add -------------------
        h6("Step 2 — Filter or Add Species (optional)"),
        radioButtons("spp_filter_mode", label = NULL,
          choices = c(
            "Use all species in lookup CSV"     = "all",
            "Filter to a NECN species CSV"      = "necn",
            "Enter species codes manually"      = "manual"
          ),
          selected = "all"
        ),
        conditionalPanel(
          "input.spp_filter_mode == 'necn'",
          p(class = "text-muted small",
            "Filters the lookup to only the species codes found in the ",
            tags$code("SpeciesCode"), " column of your NECN parameter CSV. ",
            "Useful when the lookup has more species than a specific project."),
          textInput("necn_csv_path", label = NULL,
                    placeholder = "/path/to/NECN_sp.csv",
                    value = "/Users/jlamping/Desktop/LANDIS/OLYM/WarmWet/Rep2/extensions/NECN_sp_02032026.csv")
        ),
        conditionalPanel(
          "input.spp_filter_mode == 'manual'",
          p(class = "text-muted small",
            "Enter one LANDIS species code per line. These are matched against ",
            "the lookup CSV; unmatched codes receive default values."),
          textAreaInput("manual_species_codes",
                        label   = NULL,
                        rows    = 6,
                        value   = "",
                        placeholder = "PseudotsugaMenziesii\nAbiesAmabilis\nPinusContorta\n...")
        ),
        actionButton("build_spp_table", "Build Species Table",
                     class = "btn-sm btn-secondary w-100",
                     icon  = icon("list-check")),
        uiOutput("spp_load_status"),

        hr(),

        # ---- Step 3: FIA bark parameters ----------------------------------
        h6("Step 3 — Fetch FIA Bark Parameters"),
        p(class = "text-muted small",
          "Reads FIA TREE table data (BARK_THICK, DIA, TOTAGE) for the selected ",
          "states. Estimates ", tags$b("MaximumBarkThickness"), " (asymptotic ",
          "bark thickness in cm) and ", tags$b("AgeDBH"), " (Michaelis-Menten ",
          "half-saturation age) from SCF Eq. 10. Reads local CSV files if ",
          "available and non-empty; downloads via rFIA for any missing states."),

        h6("FIA States"),
        selectizeInput("fia_states", label = NULL,
          choices  = c(
            "AK"="AK","AZ"="AZ","CA"="CA","CO"="CO","HI"="HI",
            "ID"="ID","MT"="MT","NM"="NM","NV"="NV","OR"="OR",
            "UT"="UT","WA"="WA","WY"="WY",
            "IA"="IA","IL"="IL","IN"="IN","KS"="KS","MI"="MI",
            "MN"="MN","MO"="MO","ND"="ND","NE"="NE","OH"="OH",
            "SD"="SD","WI"="WI",
            "AL"="AL","AR"="AR","FL"="FL","GA"="GA","KY"="KY",
            "LA"="LA","MS"="MS","NC"="NC","OK"="OK","SC"="SC",
            "TN"="TN","TX"="TX","VA"="VA","WV"="WV",
            "CT"="CT","DE"="DE","MA"="MA","MD"="MD","ME"="ME",
            "NH"="NH","NJ"="NJ","NY"="NY","PA"="PA","RI"="RI","VT"="VT"
          ),
          selected = c("OR", "WA", "ID", "MT"),
          multiple = TRUE,
          options  = list(plugins = list("remove_button"), maxItems = 51)
        ),

        h6("Local FIA Data Directory"),
        p(class = "text-muted small",
          "Directory containing pre-downloaded FIA state CSVs ",
          "(e.g. ", tags$code("OR_TREE.csv"), "). Files must be non-empty. ",
          "Empty or missing files for a state will be downloaded via rFIA ",
          "and saved here for future use. Leave blank to use a session ",
          "temp folder (data lost on restart)."),
        textInput("fia_local_dir", label = NULL,
                  placeholder = "/path/to/FIA_RAW/",
                  value = "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/FIA/FIA_RAW"),

        actionButton("fetch_fia_bark", "Fetch FIA Bark Parameters",
                     class = "btn-success w-100",
                     icon  = icon("database")),
        br(), br(),
        uiOutput("fia_fetch_status")
      ),

      # -- Right: outputs ------------------------------------------------------
      card(
        card_header("Fire_Spp_Table.csv"),
        full_screen = TRUE,

        navset_tab(

          nav_panel("Species Parameters",
            br(),
            p(class = "text-muted",
              "Double-click any cell in the ",
              tags$b("AgeDBH"), " or ", tags$b("MaximumBarkThickness"),
              " columns to edit values directly. Changes are saved automatically."),
            p(class = "text-muted small",
              tags$b("AgeDBH:"), " Age (years) at which bark thickness equals half of ",
              "MaximumBarkThickness. Lower values = faster bark development. ",
              "Typical range: 30\u2013150 years depending on species.",
              br(),
              tags$b("MaximumBarkThickness:"), " Asymptotic bark thickness (cm, one side) ",
              "approached at very large DBH. Ponderosa pine \u224840 mm, lodgepole pine ",
              "\u22488 mm, Douglas-fir \u224825 mm."),
            br(),
            DTOutput("spp_table_dt"),
            br(),
            fluidRow(
              column(6,
                downloadButton("dl_fire_spp_table", "Download Fire_Spp_Table.csv",
                               class = "btn-primary w-100")
              ),
              column(6,
                downloadButton("dl_spp_table_full", "Download Full Table (with notes)",
                               class = "btn-sm btn-outline-secondary w-100")
              )
            )
          ),

          nav_panel("Bark Thickness Curves",
            br(),
            p(class = "text-muted",
              "SCF Eq. 10: BarkThickness = (MaxBarkThickness \u00d7 Age) / (Age + AgeDBH). ",
              "Curves show predicted bark thickness vs stand age for each species in ",
              "the table. Thicker-barked species at a given age are more resistant to ",
              "fire-induced cambium kill."),
            plotOutput("bark_curve_plot", height = "500px")
          ),

          nav_panel("Species Lookup Table",
            br(),
            p(class = "text-muted",
              "Contents of the loaded LANDIS \u2192 FIA species lookup CSV. ",
              "Load a different CSV in Step 1 to update this table. ",
              "Add rows to the CSV file to support additional species."),
            DTOutput("fia_lookup_dt")
          )
        )
      )
    )
  ),

  # ---------------------------------------------------------------------------
  # Tab 3: Spread Calibration
  # ---------------------------------------------------------------------------
  nav_panel(
    title = "Spread Calibration",
    icon  = icon("wind"),

    layout_columns(
      col_widths = c(3, 9),

      card(
        card_header("Spread Inputs"),

        h6("GeoMAC Perimeter Geodatabase"),
        p(class = "text-muted small",
          "Historic_Geomac_Perimeters_All_Years_2000_2018.gdb"),
        textInput("geomac_gdb_path", label = NULL,
                  placeholder = "/path/to/Historic_Geomac.gdb"),

        hr(),
        h6("Cleaning Thresholds"),
        p(class = "text-muted small",
          "These filters control which perimeter pairs are used for each model."),
        numericInput("max_gap_spread_prob", "Max Gap Days — Spread Probability",
                     value = 3, min = 1, max = 14),
        p(class = "text-muted small",
          "Only pairs where consecutive perimeters are \u2264 N days apart feed the ",
          "cell-to-cell spread probability model. Shorter gaps = cleaner signal ",
          "but fewer pairs. Raise this if too few pairs are found."),
        numericInput("max_gap_daily_area", "Max Gap Days — Daily Spread Area",
                     value = 7, min = 1, max = 30),
        p(class = "text-muted small",
          "Maximum day gap allowed when estimating max daily spread area. ",
          "Growth is divided by gap length to get ha/day, so wider gaps are ",
          "tolerable here than for spread probability."),
        numericInput("neg_growth_tol", "Negative Growth Tolerance (ha)",
                     value = 1.0, min = 0, step = 0.5),
        p(class = "text-muted small",
          "Pairs where the later perimeter is more than N ha ",
          "smaller than the earlier are dropped. Small negative values (\u2264 1 ha) ",
          "are GeoMAC remapping noise; large negatives indicate a data problem."),

        hr(),
        h6("Spread Rasterization"),
        numericInput("spread_res_m", "Cell Size (m)",
                     value = 90, min = 30, max = 2000, step = 10),
        p(class = "text-muted small",
          "Resolution used to rasterize GeoMAC perimeters across the calibration ",
          "boundary. Defaults to the LANDIS template cell size. Finer resolution ",
          "captures more fire detail but increases processing time."),

        hr(),
        h6("Sampling Control"),
        numericInput("failure_sample_ratio", "Failure : Success Cell Ratio",
                     value = 3, min = 1, max = 10),
        p(class = "text-muted small",
          "For every cell that newly burned (success), sample N cells that were ",
          "adjacent but did not burn (failure). Higher values give a more balanced ",
          "logistic regression but increase memory use."),
        numericInput("max_samples_per_pair", "Max Cells Per Pair",
                     value = 25000, min = 1000, step = 1000),
        p(class = "text-muted small",
          "Cap on total raster cells per perimeter pair. Prevents very large fires ",
          "from dominating the fit."),

        hr(),
        h6("Fine Fuels"),
        radioButtons("fine_fuels_source", label = NULL,
          choices = c(
            "None — placeholder (1.0 everywhere)"  = "none",
            "Fetch LANDFIRE FBFM40 (remote)"       = "landfire",
            "Local file (0-1 scaled GeoTIFF)"      = "local"
          ),
          selected = "none"
        ),
        conditionalPanel(
          "input.fine_fuels_source == 'landfire'",
          p(class = "text-muted small",
            "Pulls LANDFIRE FBFM40 (Scott & Burgan 40 fuel models) for the ",
            "calibration boundary via ", tags$a("FedData", target = "_blank",
            href = "https://cran.r-project.org/package=FedData"),
            ". Reclassified to a 0-1 fine fuel load index using 1-hr + 10-hr ",
            "dead fuel loading values from Scott & Burgan (2005). Covers CONUS. ",
            "First run downloads and caches locally; subsequent runs are fast.")
        ),
        conditionalPanel(
          "input.fine_fuels_source == 'local'",
          textInput("fine_fuels_path", label = NULL,
                    placeholder = "/path/to/fine_fuels.tif"),
          p(class = "text-muted small",
            "Any GeoTIFF scaled 0-1 (e.g. pre-downloaded LANDFIRE 1-hr fuel, ",
            "FCCS surface fuel load, or custom map). Will be reprojected and ",
            "resampled to the spread grid automatically.")
        ),

        hr(),
        actionButton("run_spread_fit", "Run Spread Fitting",
                     class = "btn-success w-100", icon = icon("play")),
        br(), br(),
        uiOutput("spread_fit_status")
      ),

      card(
        card_header("Spread Outputs"),
        full_screen = TRUE,
        navset_tab(

          nav_panel("Fit Coefficients",
            br(),
            h5("Cell-to-Cell Spread Probability (Logistic)"),
            p(class = "text-muted small",
              "Maps to SpreadProbabilityB0-B3 in SCF parameter file."),
            DTOutput("spread_prob_coef_table"),
            br(),
            h5("Maximum Daily Spread Area (Linear)"),
            p(class = "text-muted small",
              "Maps to MaximumSpreadAreaB0-B2 in SCF parameter file."),
            DTOutput("max_area_coef_table"),
            br(),
            downloadButton("dl_spread_coef", "Download All Spread Coefficients CSV",
                           class = "btn-sm btn-outline-primary")
          ),

          nav_panel("Maps",
            br(),
            fluidRow(
              column(6,
                h5("Fine Fuels Index (FBFM40-derived)"),
                p(class = "text-muted small",
                  "1-hr + 10-hr dead fuel load normalized 0\u20131.",
                  "Used as B2 predictor in the spread probability equation."),
                uiOutput("fine_fuels_map_ui")
              ),
              column(6,
                h5("GeoMAC Calibration Fires"),
                p(class = "text-muted small",
                  "Historic fire perimeters loaded within the calibration",
                  "boundary and date range. Colored by fire year."),
                uiOutput("geomac_fires_map_ui")
              )
            )
          ),

          nav_panel("Diagnostics",
            br(),
            fluidRow(
              column(6,
                h6("Spread Probability vs FWI"),
                plotOutput("spread_prob_fwi_plot",  height = "280px")
              ),
              column(6,
                h6("Spread Probability vs Wind Speed"),
                plotOutput("spread_prob_wind_plot", height = "280px")
              )
            ),
            fluidRow(
              column(6,
                h6("ROC Curve — Spread Probability Model"),
                plotOutput("spread_roc_plot",        height = "280px")
              ),
              column(6,
                h6("Max Daily Spread Area: Observed vs Predicted"),
                plotOutput("max_area_scatter_plot",  height = "280px")
              )
            )
          ),

          nav_panel("Terrain",
            br(),
            p(class = "text-muted small",
              "SRTM/3DEP elevation fetched via elevatr at zoom level 8 (~600 m). ",
              "Slope and aspect drive effective wind speed at each fire centroid ",
              "via the upslope wind component plus a Rothermel slope equivalent."),
            fluidRow(
              column(4,
                h6("Effective vs Raw Wind Speed"),
                p(class = "text-muted small",
                  "Each point is one perimeter pair. Colour = slope at centroid."),
                plotOutput("wind_comparison_plot", height = "260px")
              ),
              column(4,
                h6("Slope Distribution at Fire Centroids"),
                plotOutput("slope_hist_plot",      height = "260px")
              ),
              column(4,
                h6("Wind-to-Upslope Alignment"),
                p(class = "text-muted small",
                  "0\u00b0 = wind blowing straight upslope; \u00b1180\u00b0 = downslope."),
                plotOutput("wind_align_plot",      height = "260px")
              )
            ),
            br(),
            fluidRow(
              column(6,
                h6("Slope (degrees)"),
                tmapOutput("slope_map", height = "360px")
              ),
              column(6,
                h6("Aspect (degrees from N)"),
                tmapOutput("aspect_map", height = "360px")
              )
            )
          ),

          nav_panel("Perimeter Pairs",
            br(),
            p(class = "text-muted", "Cleaned perimeter progression pairs."),
            DTOutput("pairs_table"),
            br(),
            downloadButton("dl_pairs", "Download Pairs CSV",
                           class = "btn-sm btn-outline-primary")
          ),

          nav_panel("Candidate Grid",
            br(),
            p(class = "text-muted",
              "Define a grid of candidate SCF parameter sets. Export snippets ",
              "to run LANDIS externally, then return here to score."),
            fluidRow(
              column(6,
                h6("Spread Prob B0 range"),
                numericInput("cand_b0_min",  "Min",  value = -5),
                numericInput("cand_b0_max",  "Max",  value = 0),
                numericInput("cand_b0_step", "Step", value = 1)
              ),
              column(6,
                h6("Spread Prob B1 (FWI) range"),
                numericInput("cand_b1_min",  "Min",  value = 0),
                numericInput("cand_b1_max",  "Max",  value = 0.5),
                numericInput("cand_b1_step", "Step", value = 0.1)
              )
            ),
            actionButton("build_candidates", "Build Candidate Grid",
                         class = "btn-sm btn-outline-primary"),
            br(), br(),
            DTOutput("candidate_table"),
            br(),
            downloadButton("dl_candidates", "Download Candidate Table CSV",
                           class = "btn-sm btn-outline-primary"),
            br(),
            downloadButton("dl_snippets_zip", "Download SCF Snippets (.zip)",
                           class = "btn-sm btn-outline-secondary")
          ),

          nav_panel("Score Candidates",
            br(),
            p(class = "text-muted",
              "After running LANDIS externally for each candidate, provide the ",
              "root directory containing run outputs to score them."),
            textInput("scf_runs_root", "SCF Runs Root Directory",
                      placeholder = "/path/to/scf_candidate_runs/PARK"),
            fluidRow(
              column(6, numericInput("score_w_size", "Weight: Event Size",
                                     value = 0.5, min = 0, max = 1, step = 0.1)),
              column(6, numericInput("score_w_aab",  "Weight: Annual Area",
                                     value = 0.5, min = 0, max = 1, step = 0.1))
            ),
            actionButton("run_scoring", "Score Candidates",
                         class = "btn-sm btn-success"),
            br(), br(),
            DTOutput("scores_table"),
            br(),
            plotOutput("cdf_plot",    height = "300px"),
            plotOutput("annual_plot", height = "300px"),
            br(),
            downloadButton("dl_scores", "Download Scores CSV",
                           class = "btn-sm btn-outline-primary")
          )
        )
      )
    )
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  tmap_mode("view")

  rv <- reactiveValues(
    # Spatial inputs
    cal_vect       = NULL,   # calibration boundary, native CRS
    cal_vect_wgs84 = NULL,   # calibration boundary, WGS84 (for leaflet)
    cal_vect_proj  = NULL,   # calibration boundary, working CRS
    template_r     = NULL,   # LANDIS template raster
    working_crs    = NULL,   # CRS string from template_r
    template_mask  = NULL,   # binary mask on template grid
    cell_area_ha   = NULL,

    # ERA time series (extracted over calibration boundary)
    era_fwi_daily  = NULL,
    era_ffmc_daily = NULL,
    era_dmc_daily  = NULL,
    era_dc_daily   = NULL,
    wind_daily     = NULL,

    # Ignition outputs
    ign_df               = NULL,
    ign_model_data       = NULL,
    ign_coef             = NULL,
    surf_lightning       = NULL,
    surf_accidental      = NULL,
    surf_rx              = NULL,
    # Scale adjustment
    landscape_areas      = NULL,   # list(cal_ha, tmpl_ha, ratio, log_ratio)
    ign_coef_adjusted    = NULL,   # scale-adjusted coefficients
    park_ign_df          = NULL,   # FPA FOD within template (Option B)
    park_model_data      = NULL,   # daily park counts (Option B)
    expected_ignitions   = NULL,   # tibble(year, ignition_type, expected_annual, scenario)
    startup_df      = NULL,
    ground_slope_r  = NULL,   # slope in degrees on template grid (INT2S export)
    uphill_slope_r  = NULL,   # upslope azimuth 0-359 deg on template grid (INT2S export)

    # Spread outputs
    pairs_clean         = NULL,
    spread_prob_coef    = NULL,
    spread_prob_model   = NULL,
    spread_prob_samples = NULL,
    max_area_coef       = NULL,
    spread_da_model     = NULL,
    spread_da_data      = NULL,
    fine_fuels_r        = NULL,
    geomac_result       = NULL,
    dem_r               = NULL,
    slope_r             = NULL,
    aspect_r            = NULL,
    candidate_grid      = NULL,
    scores_df           = NULL,

    # Species table outputs
    spp_lookup          = NULL,   # tibble from user lookup CSV: SpeciesName/FIA_SPCD/SCIENTIFIC_NAME
    spp_codes           = NULL,   # character vector of LANDIS species codes
    spp_table           = NULL,   # editable tibble: SpeciesCode/AgeDBH/MaxBT/Note
    fia_trees           = NULL    # raw FIA TREE records (cached to avoid re-download)
  )

  # ---------------------------------------------------------------------------
  # Internal helper: safely reproject a SpatVector to target CRS
  # ---------------------------------------------------------------------------

  safe_project_vect <- function(v, target_crs) {
    src_crs <- crs(v)
    if (is.na(src_crs) || nchar(src_crs) == 0) {
      warning("Input vector has no CRS — skipping reprojection.")
      return(v)
    }
    project(v, target_crs)
  }

  # ---------------------------------------------------------------------------
  # 1) Load calibration boundary
  # ---------------------------------------------------------------------------

  observeEvent(input$load_shp, {
    req(nchar(trimws(input$shp_path)) > 0)
    withProgress(message = "Loading calibration boundary...", {
      tryCatch({
        v <- vect(input$shp_path)
        rv$cal_vect       <- v
        rv$cal_vect_wgs84 <- project(v, "EPSG:4326")

        # If template already loaded, build projected copy immediately
        if (!is.null(rv$working_crs)) {
          rv$cal_vect_proj <- safe_project_vect(v, rv$working_crs)
        }

        output$shp_status <- renderUI(
          tags$span(class = "text-success small",
                    icon("check"),
                    sprintf(" Loaded: %d feature(s) | CRS: %s",
                            nrow(v), crs(v, describe = TRUE)$name))
        )
      }, error = function(e) {
        output$shp_status <- renderUI(
          tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
        )
      })
    })
  })

  # ---------------------------------------------------------------------------
  # 2) Load LANDIS template raster
  # ---------------------------------------------------------------------------

  observeEvent(input$load_tif, {
    req(nchar(trimws(input$tif_path)) > 0)
    withProgress(message = "Loading LANDIS template raster...", {
      tryCatch({
        r <- rast(input$tif_path)

        if (is.na(crs(r)) || nchar(crs(r)) == 0) {
          stop("Template raster has no CRS defined.")
        }

        rv$template_r   <- r
        rv$working_crs  <- crs(r)
        rv$cell_area_ha <- prod(res(r)) / 10000

        # Binary mask: 1 = active cell, NA = inactive
        msk <- r; msk[!is.na(msk)] <- 1
        rv$template_mask <- msk

        # Pre-fill spread rasterization cell size from template resolution
        updateNumericInput(session, "spread_res_m",
                           value = round(mean(res(r))))

        # Reproject calibration boundary to match if already loaded
        if (!is.null(rv$cal_vect)) {
          rv$cal_vect_proj <- safe_project_vect(rv$cal_vect, rv$working_crs)
        }

        output$tif_status <- renderUI(
          tags$span(class = "text-success small",
                    icon("check"),
                    sprintf(" Loaded: %d x %d | res %.0f m | %.2f ha/cell | CRS: %s",
                            nrow(r), ncol(r), mean(res(r)),
                            rv$cell_area_ha, crs(r, describe = TRUE)$name))
        )
      }, error = function(e) {
        output$tif_status <- renderUI(
          tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
        )
      })
    })
  })

  # ---------------------------------------------------------------------------
  # 3) Load ERA-Interim FWI — uses calibration boundary for spatial extraction
  # ---------------------------------------------------------------------------

  observeEvent(input$load_era_fwi, {
    req(rv$cal_vect)
    req(nchar(trimws(input$era_fwi_path)) > 0)

    paths <- list(FWI  = input$era_fwi_path,
                  FFMC = input$era_ffmc_path,
                  DMC  = input$era_dmc_path,
                  DC   = input$era_dc_path)

    withProgress(message = "Extracting ERA-Interim time series...", value = 0, {
      tryCatch({
        cal_start <- as.Date(sprintf("%d-01-01", input$cal_years[1]))
        cal_end   <- as.Date(sprintf("%d-12-31", input$cal_years[2]))

        for (nm in names(paths)) {
          p <- trimws(paths[[nm]])
          if (nchar(p) == 0) next
          setProgress(message = sprintf("Processing %s...", nm))
          ts <- extract_era_timeseries(p, rv$cal_vect,
                                        cal_start = cal_start,
                                        cal_end   = cal_end)
          if (nm == "FWI")  rv$era_fwi_daily  <- ts
          if (nm == "FFMC") rv$era_ffmc_daily <- ts
          if (nm == "DMC")  rv$era_dmc_daily  <- ts
          if (nm == "DC")   rv$era_dc_daily   <- ts
        }

        output$era_fwi_status <- renderUI(
          tags$span(class = "text-success small",
                    icon("check"),
                    sprintf(" %d daily FWI values | range %.1f - %.1f",
                            nrow(rv$era_fwi_daily),
                            min(rv$era_fwi_daily$value, na.rm = TRUE),
                            max(rv$era_fwi_daily$value, na.rm = TRUE)))
        )
      }, error = function(e) {
        output$era_fwi_status <- renderUI(
          tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
        )
      })
    })
  })

  # ---------------------------------------------------------------------------
  # 4) Fetch GRIDMET wind — remote OPeNDAP, no local files required
  # ---------------------------------------------------------------------------

  observeEvent(input$fetch_gridmet_wind, {
    req(rv$cal_vect)

    output$era_wind_status <- renderUI(
      tags$span(class = "text-warning small",
                icon("spinner", class = "fa-spin"),
                " Fetching GRIDMET wind data (this may take a minute)...")
    )

    withProgress(message = "Fetching GRIDMET wind...", value = 0, {
      tryCatch({
        cal_start <- as.Date(sprintf("%d-01-01", input$cal_years[1]))
        cal_end   <- as.Date(sprintf("%d-12-31", input$cal_years[2]))

        wind_ts <- extract_gridmet_wind(
          cal_vect    = rv$cal_vect,
          cal_start   = cal_start,
          cal_end     = cal_end,
          progress_fn = function(i, n, yr) {
            setProgress(i / n,
                        detail = sprintf("Year %d (%d of %d)...", yr, i, n))
          }
        )
        rv$wind_daily <- wind_ts

        output$era_wind_status <- renderUI(
          tags$span(class = "text-success small",
                    icon("check"),
                    sprintf(" %d daily GRIDMET wind records | mean %.1f km/h",
                            nrow(wind_ts),
                            mean(wind_ts$WindSpeed_kmh, na.rm = TRUE)))
        )
      }, error = function(e) {
        output$era_wind_status <- renderUI(
          tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
        )
      })
    })
  })

  # ---------------------------------------------------------------------------
  # Sidebar preview map
  # Green = calibration boundary | Orange = LANDIS template extent
  # ---------------------------------------------------------------------------
  output$landscape_preview_ui <- renderUI({
    req(rv$cal_vect_wgs84)
    tagList(
      h6("Landscape Preview"),
      p(class = "text-muted small",
        tags$span(style = "color:#2c5f2e; font-weight:bold;", "\u25a0"),
        " Calibration boundary  ",
        if (!is.null(rv$template_r))
          tagList(tags$span(style = "color:#a0522d; font-weight:bold;", "\u25a0"),
                  " LANDIS template extent")
      ),
      tmapOutput("landscape_map", height = "220px")
    )
  })
  
  
  output$landscape_map <- renderTmap({
    req(rv$cal_vect_wgs84)
    
    m <- tm_shape(sf::st_as_sf(rv$cal_vect_wgs84)) +
      tm_borders(col = "#2c5f2e", lwd = 2) +
      tm_fill(alpha = 0.1, col = "#2c5f2e")
    
    if (!is.null(rv$template_r)) {
      tmpl_bb <- sf::st_as_sfc(sf::st_bbox(
        sf::st_transform(
          sf::st_as_sf(as.polygons(ext(rv$template_r),
                                   crs = crs(rv$template_r))),
          4326
        )
      ))
      tmpl_sf <- sf::st_sf(geometry = tmpl_bb)
      m <- m +
        tm_shape(tmpl_sf) +
        tm_borders(col = "#a0522d", lwd = 2) +
        tm_fill(alpha = 0.08, col = "#a0522d")
    }
    
    m
  })


  # ---------------------------------------------------------------------------
  # TAB 1: Run Ignition Calibration
  # ---------------------------------------------------------------------------

  observeEvent(input$run_ignition, {
    req(rv$cal_vect_proj, rv$template_r, rv$era_fwi_daily)
    req(nchar(trimws(input$fpa_gdb_path)) > 0)

    output$ign_run_status <- renderUI(
      tags$span(class = "text-warning small",
                icon("spinner", class = "fa-spin"), " Running ignition calibration...")
    )

    withProgress(message = "Ignition calibration...", value = 0, {

      # Step 1: Load FPA FOD clipped to calibration boundary
      setProgress(0.1, detail = "Loading FPA FOD...")
      ign_df <- tryCatch(
        load_fpa_fod(
          gdb_path    = input$fpa_gdb_path,
          cal_vect    = rv$cal_vect_proj,
          working_crs = rv$working_crs,
          year_min    = input$cal_years[1],
          year_max    = input$cal_years[2]
        ),
        error = function(e) {
          output$ign_run_status <- renderUI(
            tags$span(class = "text-danger small", icon("xmark"),
                      " FPA FOD error: ", conditionMessage(e))
          )
          NULL
        }
      )
      req(ign_df)
      rv$ign_df <- ign_df

      # Step 2: Daily count x FWI
      setProgress(0.3, detail = "Joining FWI...")
      rv$ign_model_data <- build_ignition_model_data(ign_df, rv$era_fwi_daily)

      # Step 3: Fit models
      setProgress(0.5, detail = "Fitting ignition models...")
      rv$ign_coef <- fit_ignition_models(rv$ign_model_data,
                                          distribution = input$ign_distribution)

      # Step 4: Ignition surfaces on LANDIS template grid
      setProgress(0.7, detail = "Generating ignition surfaces on LANDIS grid...")
      surfaces <- make_ignition_surfaces(
        ign_df      = ign_df,
        template    = rv$template_r,
        working_crs = rv$working_crs,
        bw_m        = input$kernel_bw,
        maxdist_m   = input$kernel_maxdist
      )
      rv$surf_lightning  <- surfaces$Lightning
      rv$surf_accidental <- surfaces$Accidental
      rv$surf_rx         <- surfaces$Rx

      setProgress(1.0)
      output$ign_run_status <- renderUI(
        tags$span(class = "text-success small", icon("check"),
                  sprintf(" Complete. %d ignitions (%d Lightning, %d Accidental).",
                          nrow(ign_df),
                          sum(ign_df$ignition_type == "Lightning"),
                          sum(ign_df$ignition_type == "Accidental")))
      )
    })
  })

  # -- Ignition map: lightning surface + both boundaries
  output$ign_map <- renderTmap({
    req(rv$surf_lightning, rv$cal_vect_wgs84)
    
    # Aggregate if very large for display
    surf <- rv$surf_lightning
    if (ncell(surf) > 500000) {
      fact <- ceiling(sqrt(ncell(surf) / 500000))
      surf <- aggregate(surf, fact = fact, fun = "mean", na.rm = TRUE)
    }
    
    # Renormalize after any aggregation; replace NaN (from all-NA aggregate cells)
    v  <- values(surf, mat = FALSE)
    mx <- max(v, na.rm = TRUE)
    if (is.finite(mx) && mx > 0) v <- v / mx
    v[!is.finite(v)] <- NA
    values(surf) <- v

    m <- tm_shape(surf) +
      tm_raster(
        col        = names(surf),
        col.scale  = tm_scale_continuous(
          values   = "brewer.yl_or_rd",
          value.na = NA
        ),
        col.legend = tm_legend(title = "Ignition Weight"),
        col_alpha  = 0.75
      ) +
      tm_shape(sf::st_as_sf(rv$cal_vect_wgs84)) +
      tm_borders(col = "#2c5f2e", lwd = 2)
    
    if (!is.null(rv$template_r)) {
      tmpl_sf <- sf::st_as_sf(
        as.polygons(ext(rv$template_r), crs = crs(rv$template_r))
      ) %>% sf::st_transform(4326)
      
      m <- m +
        tm_shape(tmpl_sf) +
        tm_borders(col = "#a0522d", lwd = 2)
    }
    
    m
  })

  output$ign_coef_table <- renderDT({
    req(rv$ign_coef)
    rv$ign_coef %>%
      mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
      datatable(options = list(dom = "t", paging = FALSE),
                rownames = FALSE, class = "table-sm table-striped")
  })

  output$ign_diag_plot <- renderPlot({
    req(rv$ign_model_data, rv$ign_coef)
    plot_ignition_diagnostics(rv$ign_model_data, rv$ign_coef)
  })

  output$ign_zero_inf_plot <- renderPlot({
    req(rv$ign_model_data, rv$ign_coef)
    plot_zero_inflation(rv$ign_model_data, rv$ign_coef)
  })

  output$ign_lambda_plot <- renderPlot({
    req(rv$ign_model_data, rv$ign_coef)
    plot_count_component(rv$ign_model_data, rv$ign_coef)
  })

  output$ign_resid_plot <- renderPlot({
    req(rv$ign_model_data, rv$ign_coef)
    plot_ignition_residuals(rv$ign_model_data, rv$ign_coef)
  })

  # ---------------------------------------------------------------------------
  # Scale Adjustment — auto-compute Option A whenever calibration completes
  # ---------------------------------------------------------------------------

  observeEvent(rv$ign_coef, {
    req(rv$cal_vect_proj, rv$template_r, rv$ign_coef, rv$era_fwi_daily)

    # Landscape areas (used by both options)
    rv$landscape_areas <- compute_landscape_areas(rv$cal_vect_proj, rv$template_r)

    # Default adjusted coef = Option A (area offset)
    rv$ign_coef_adjusted <- apply_area_offset(
      rv$ign_coef,
      rv$landscape_areas$cal_area_ha,
      rv$landscape_areas$tmpl_area_ha
    )

    # Expected ignitions: regional (unadjusted) + landscape-adjusted
    regional_annual  <- predict_annual_ignitions(
      rv$ign_coef, rv$era_fwi_daily, scenario_label = "Regional (unadjusted)"
    )
    adjusted_annual  <- predict_annual_ignitions(
      rv$ign_coef_adjusted, rv$era_fwi_daily, scenario_label = "Landscape-adjusted"
    )
    rv$expected_ignitions <- bind_rows(regional_annual, adjusted_annual)
  })

  # Re-apply Option A when user switches back to it
  observeEvent(input$scale_method, {
    req(rv$ign_coef, rv$landscape_areas, rv$era_fwi_daily)
    if (input$scale_method == "area_offset") {
      rv$ign_coef_adjusted <- apply_area_offset(
        rv$ign_coef,
        rv$landscape_areas$cal_area_ha,
        rv$landscape_areas$tmpl_area_ha
      )
      adjusted_annual <- predict_annual_ignitions(
        rv$ign_coef_adjusted, rv$era_fwi_daily, scenario_label = "Landscape-adjusted"
      )
      regional_annual <- predict_annual_ignitions(
        rv$ign_coef, rv$era_fwi_daily, scenario_label = "Regional (unadjusted)"
      )
      rv$expected_ignitions <- bind_rows(regional_annual, adjusted_annual)
    }
  }, ignoreInit = TRUE)

  # Option B: fit park-specific intercept from FPA FOD within template
  observeEvent(input$fit_park_b0, {
    req(rv$ign_df, rv$template_r, rv$era_fwi_daily,
        rv$ign_coef, rv$landscape_areas)

    output$park_fit_status_ui <- renderUI(
      tags$span(class = "text-warning small",
                icon("spinner", class = "fa-spin"), " Filtering FPA FOD to template...")
    )

    tryCatch({
      # Filter FPA FOD to template extent
      rv$park_ign_df    <- filter_fpa_to_template(rv$ign_df, rv$template_r)
      n_park            <- nrow(rv$park_ign_df)

      output$park_fit_status_ui <- renderUI(
        tags$span(class = "text-info small",
                  sprintf(" %d FPA FOD fires within template extent. Fitting...", n_park))
      )

      # Build daily park counts
      rv$park_model_data <- build_ignition_model_data(rv$park_ign_df, rv$era_fwi_daily)

      # Fit park intercept (b1 fixed from regional)
      rv$ign_coef_adjusted <- fit_park_intercept(
        rv$park_model_data, rv$ign_coef,
        rv$landscape_areas$cal_area_ha,
        rv$landscape_areas$tmpl_area_ha
      )

      # Update expected ignitions with park intercept results
      park_annual    <- predict_annual_ignitions(
        rv$ign_coef_adjusted, rv$era_fwi_daily, scenario_label = "Park intercept"
      )
      regional_annual <- predict_annual_ignitions(
        rv$ign_coef, rv$era_fwi_daily, scenario_label = "Regional (unadjusted)"
      )
      rv$expected_ignitions <- bind_rows(regional_annual, park_annual)

      output$park_fit_status_ui <- renderUI(
        tags$span(class = "text-success small",
                  icon("check"), sprintf(" Done. %d fires in template used.", n_park))
      )
    }, error = function(e) {
      output$park_fit_status_ui <- renderUI(
        tags$span(class = "text-danger small",
                  icon("xmark"), " Error: ", conditionMessage(e))
      )
    })
  })

  # -- Area summary UI ---------------------------------------------------------
  output$area_summary_ui <- renderUI({
    req(rv$landscape_areas)
    a <- rv$landscape_areas
    tags$table(class = "table table-sm table-bordered small",
      tags$tbody(
        tags$tr(tags$th("Calibration boundary"),
                tags$td(sprintf("%.2f M ha", a$cal_area_ha / 1e6))),
        tags$tr(tags$th("LANDIS template"),
                tags$td(sprintf("%.0f ha  (%.2f%%)",
                                a$tmpl_area_ha,
                                a$ratio * 100))),
        tags$tr(tags$th("log(area ratio)"),
                tags$td(sprintf("%.4f  \u2192 b0 adjusted by this amount",
                                a$log_ratio)))
      )
    )
  })

  # -- Adjusted coefficients table --------------------------------------------
  output$adjusted_coef_table <- renderDT({
    req(rv$ign_coef_adjusted)
    rv$ign_coef_adjusted %>%
      select(ignition_type, distribution, b0, b1, bz0, bz1) %>%
      mutate(across(where(is.numeric), ~round(.x, 6))) %>%
      datatable(options = list(dom = "t", paging = FALSE),
                rownames = FALSE, class = "table-sm table-striped")
  })

  output$dl_adjusted_coef <- downloadHandler(
    filename = "ignition_coefficients_adjusted.csv",
    content  = function(f) {
      req(rv$ign_coef_adjusted, rv$landscape_areas)
      a <- rv$landscape_areas
      write.csv(
        rv$ign_coef_adjusted %>%
          mutate(cal_area_ha  = a$cal_area_ha,
                 tmpl_area_ha = a$tmpl_area_ha,
                 log_offset   = a$log_ratio),
        f, row.names = FALSE
      )
    }
  )

  # -- Expected ignitions table + plot ----------------------------------------
  output$expected_ign_table <- renderDT({
    req(rv$expected_ignitions)
    summarise_expected_ignitions(rv$expected_ignitions) %>%
      rename(Type = ignition_type, Scenario = scenario,
             Mean = mean_annual, Median = median_annual,
             SD = sd_annual, P10 = p10_annual, P90 = p90_annual) %>%
      datatable(
        options = list(dom = "t", paging = FALSE, ordering = FALSE),
        rownames = FALSE, class = "table-sm table-striped"
      ) %>%
      formatStyle("Type",
        target = "row",
        fontWeight = styleEqual("TOTAL", "bold")
      )
  })

  output$expected_ign_plot <- renderPlot({
    req(rv$expected_ignitions)
    plot_expected_ignitions(rv$expected_ignitions)
  })

  # ---------------------------------------------------------------------------
  # ERA Startup Codes
  # ---------------------------------------------------------------------------

  output$startup_date_ui <- renderUI({
    req(rv$era_fwi_daily)
    dateInput("startup_ref_date",
              "Climate data start date",
              value = min(rv$era_fwi_daily$date),
              min   = min(rv$era_fwi_daily$date),
              max   = max(rv$era_fwi_daily$date))
  })

  observeEvent(input$run_startup, {
    req(rv$era_ffmc_daily, rv$era_dmc_daily, rv$era_dc_daily, rv$era_fwi_daily)
    rv$startup_df <- compute_startup_values(
      ffmc_ts  = rv$era_ffmc_daily,
      dmc_ts   = rv$era_dmc_daily,
      dc_ts    = rv$era_dc_daily,
      fwi_ts   = rv$era_fwi_daily,
      ref_date = input$startup_ref_date
    )
    output$startup_table <- renderDT({
      rv$startup_df %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        datatable(options = list(dom = "t", paging = FALSE),
                  rownames = FALSE, class = "table-sm table-striped")
    })
    output$startup_clim_plot <- renderPlot({
      plot_era_climatologies(rv$era_ffmc_daily, rv$era_dmc_daily,
                              rv$era_dc_daily, rv$era_fwi_daily,
                              input$startup_ref_date)
    })
  })

  # ---------------------------------------------------------------------------
  # TAB 2: Run Spread Fitting
  # ---------------------------------------------------------------------------

  observeEvent(input$run_spread_fit, {
    req(rv$cal_vect_proj, rv$template_r,
        rv$era_fwi_daily, rv$wind_daily)
    req(nchar(trimws(input$geomac_gdb_path)) > 0)

    output$spread_fit_status <- renderUI(
      tags$span(class = "text-warning small",
                icon("spinner", class = "fa-spin"),
                " Running spread fitting (this may take several minutes)...")
    )

    withProgress(message = "Spread fitting...", value = 0, {

      setProgress(0.1, detail = "Loading GeoMAC perimeters...")
      geomac_result <- tryCatch(
        load_geomac_perimeters(
          gdb_path    = input$geomac_gdb_path,
          cal_vect    = rv$cal_vect_proj,
          template_r  = rv$template_r,
          working_crs = rv$working_crs,
          cal_start   = as.Date(sprintf("%d-01-01", input$cal_years[1])),
          cal_end     = as.Date(sprintf("%d-12-31", input$cal_years[2]))
        ),
        error = function(e) {
          output$spread_fit_status <- renderUI(
            tags$span(class = "text-danger small", icon("xmark"),
                      " GeoMAC error: ", conditionMessage(e))
          )
          NULL
        }
      )
      req(geomac_result)
      rv$geomac_result <- geomac_result

      # ---- Terrain (DEM -> slope + aspect) ------------------------------------
      setProgress(0.18, detail = "Fetching terrain (DEM)...")
      terrain_result <- tryCatch(
        fetch_terrain(
          cal_vect    = rv$cal_vect_proj,
          working_crs = rv$working_crs,
          z           = 8,
          progress_fn = function(msg) setProgress(0.18, detail = msg)
        ),
        error = function(e) {
          warning("Terrain fetch failed, using raw wind speed: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(terrain_result)) {
        rv$dem_r    <- terrain_result$dem
        rv$slope_r  <- terrain_result$slope
        rv$aspect_r <- terrain_result$aspect
        message("Terrain loaded successfully.")
      } else {
        rv$dem_r <- rv$slope_r <- rv$aspect_r <- NULL
      }

      setProgress(0.2, detail = "Building perimeter pairs with effective wind...")
      climate_joined <- bind_climate_to_pairs(
        geomac_v   = geomac_result,
        fwi_daily  = rv$era_fwi_daily,
        wind_daily = rv$wind_daily,
        slope_r    = rv$slope_r,
        aspect_r   = rv$aspect_r,
        max_gap_sp = input$max_gap_spread_prob,
        max_gap_da = input$max_gap_daily_area,
        neg_tol_ha = input$neg_growth_tol
      )
      rv$pairs_clean <- climate_joined$pairs

      setProgress(0.5, detail = "Fitting max daily area model...")
      da_fit             <- fit_max_daily_area(climate_joined$pairs)
      rv$max_area_coef   <- da_fit$coef
      rv$spread_da_model <- da_fit$model
      rv$spread_da_data  <- da_fit$data

      # Build a spread rasterization grid from the calibration boundary at the
      # user-specified cell size (defaults to LANDIS template resolution).
      setProgress(0.55, detail = "Building spread rasterization grid...")
      spread_res_m <- max(as.numeric(input$spread_res_m), 10)  # floor at 10m
      spread_template <- rast(
        ext(rv$cal_vect_proj),
        resolution = spread_res_m,
        crs        = rv$working_crs
      )
      spread_mask <- rasterize(rv$cal_vect_proj, spread_template,
                               field = 1, background = NA)
      spread_mask[!is.na(spread_mask)] <- 1
      message(sprintf(
        "Spread raster grid: %d x %d cells at %.0f m resolution over calibration boundary.",
        nrow(spread_template), ncol(spread_template), spread_res_m
      ))

      # ---- Fine fuels raster --------------------------------------------------
      setProgress(0.58, detail = "Preparing fine fuels...")
      fine_fuels_r <- tryCatch({
        if (input$fine_fuels_source == "landfire") {
          fetch_landfire_fine_fuels(
            cal_vect        = rv$cal_vect_proj,
            working_crs     = rv$working_crs,
            spread_template = spread_template,
            progress_fn     = function(msg) setProgress(0.58, detail = msg)
          )
        } else if (input$fine_fuels_source == "local" &&
                   nchar(trimws(input$fine_fuels_path)) > 0) {
          load_local_fine_fuels(
            path            = input$fine_fuels_path,
            working_crs     = rv$working_crs,
            spread_template = spread_template
          )
        } else {
          NULL   # placeholder — fit_spread_probability uses 1.0 everywhere
        }
      }, error = function(e) {
        warning("Fine fuels load failed, falling back to placeholder: ",
                conditionMessage(e))
        NULL
      })

      rv$fine_fuels_r <- fine_fuels_r

      if (!is.null(fine_fuels_r)) {
        message(sprintf("Fine fuels loaded: range [%.3f, %.3f]",
                        min(values(fine_fuels_r), na.rm = TRUE),
                        max(values(fine_fuels_r), na.rm = TRUE)))
      } else {
        message("Fine fuels: using placeholder (1.0 everywhere).")
      }

      setProgress(0.6, detail = "Fitting spread probability (rasterizing perimeters)...")
      sp_fit <- fit_spread_probability(
        pairs2        = climate_joined$spread_pairs,
        geomac_v      = geomac_result,
        template_r    = spread_template,
        park_mask     = spread_mask,
        fine_fuels_r  = fine_fuels_r,
        failure_ratio = input$failure_sample_ratio,
        max_samples   = input$max_samples_per_pair,
        progress_fn   = function(i, n) setProgress(
          0.6 + 0.35 * (i / n), detail = sprintf("Pair %d of %d...", i, n))
      )
      rv$spread_prob_coef    <- sp_fit$coef
      rv$spread_prob_model   <- sp_fit$model
      rv$spread_prob_samples <- sp_fit$samples

      setProgress(1.0)
      output$spread_fit_status <- renderUI(
        tags$span(class = "text-success small", icon("check"),
                  sprintf(" Complete. %d perimeter pairs processed.",
                          nrow(climate_joined$pairs)))
      )
    })
  })

  output$spread_prob_coef_table <- renderDT({
    req(rv$spread_prob_coef)
    rv$spread_prob_coef %>%
      mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
      datatable(options = list(dom = "t", paging = FALSE),
                rownames = FALSE, class = "table-sm table-striped")
  })

  output$max_area_coef_table <- renderDT({
    req(rv$max_area_coef)
    rv$max_area_coef %>%
      mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
      datatable(options = list(dom = "t", paging = FALSE),
                rownames = FALSE, class = "table-sm table-striped")
  })

  # ---------------------------------------------------------------------------
  # Maps tab outputs
  # ---------------------------------------------------------------------------

  output$fine_fuels_map_ui <- renderUI({
    if (!is.null(rv$fine_fuels_r)) {
      tmapOutput("fine_fuels_map", height = "420px")
    } else {
      div(class = "text-muted p-3 border rounded mt-2",
          icon("circle-info"), " ",
          "Fine fuels map not available — spread was fitted with a placeholder ",
          "(1.0 everywhere). Select LANDFIRE or local file and re-run spread fitting.")
    }
  })

  output$fine_fuels_map <- renderTmap({
    req(rv$fine_fuels_r, rv$cal_vect_proj)
    ff  <- rv$fine_fuels_r

    # leafem addStarsImage hard limit is 4 MB (4,194,304 bytes).
    # RGBA = 4 bytes/pixel -> max ~1,000,000 pixels safely.
    # Aggregate to stay well under that ceiling.
    max_cells <- 800000L
    n_cells   <- nrow(ff) * ncol(ff)
    if (n_cells > max_cells) {
      fact <- ceiling(sqrt(n_cells / max_cells))
      ff   <- aggregate(ff, fact = fact, fun = "mean", na.rm = TRUE)
    }

    bnd <- project(rv$cal_vect_proj, crs(ff))
    tm_shape(ff) +
      tm_raster(
        col        = names(ff)[1],
        col.scale  = tm_scale_intervals(
          values   = "brewer.yl_gn",
          breaks   = c(0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0),
          value.na = NA
        ),
        col.legend = tm_legend(title = "Fine Fuel Load\n(0\u20131)"),
        col_alpha  = 0.85
      ) +
      tm_shape(bnd) +
      tm_borders(col = "gray30", lwd = 1.5)
  })

  output$geomac_fires_map_ui <- renderUI({
    if (!is.null(rv$geomac_result)) {
      tmapOutput("geomac_fires_map", height = "420px")
    } else {
      div(class = "text-muted p-3 border rounded mt-2",
          icon("circle-info"), " ",
          "GeoMAC fires will appear here after spread fitting is run.")
    }
  })

  output$geomac_fires_map <- renderTmap({
    req(rv$geomac_result, rv$cal_vect_proj)
    fires <- project(rv$geomac_result, "EPSG:4326")
    bnd   <- project(rv$cal_vect_proj, "EPSG:4326")
    yr_col <- grep("FIRE_YEAR|fire_year|year",
                   names(fires), ignore.case = TRUE, value = TRUE)[1]
    fill_col <- if (!is.na(yr_col)) yr_col else names(fires)[1]
    tm_shape(bnd) +
      tm_borders(col = "gray30", lwd = 2) +
      tm_shape(fires) +
      tm_polygons(
        col        = fill_col,
        col.scale  = tm_scale_ordinal(values = "brewer.set1"),
        col_alpha  = 0.55,
        col.legend = tm_legend(title = "Fire Year")
      )
  })

  # ---------------------------------------------------------------------------
  # Terrain tab outputs
  # ---------------------------------------------------------------------------

  output$wind_comparison_plot <- renderPlot({
    req(rv$pairs_clean)
    dat <- rv$pairs_clean %>%
      filter(is.finite(WindSpeed_kmh), is.finite(EffectiveWind))
    if (nrow(dat) == 0) return(NULL)
    has_slope <- "Slope_deg" %in% names(dat) && any(is.finite(dat$Slope_deg))
    p <- ggplot(dat, aes(WindSpeed_kmh, EffectiveWind))
    if (has_slope) {
      p <- p + geom_point(aes(colour = Slope_deg), size = 2.5, alpha = 0.8) +
        scale_colour_gradientn(
          colours = c("#a8d5a2","#f5a623","#c0392b"),
          name    = "Slope\n(deg)")
    } else {
      p <- p + geom_point(colour = "#c0392b", size = 2.5, alpha = 0.8)
    }
    p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
      labs(title   = "Effective vs Raw Wind Speed",
           subtitle = if (has_slope) "Terrain-adjusted" else "No terrain (placeholder)",
           x = "Raw Wind Speed (km/h)",
           y = "Effective Wind Speed (km/h)") +
      theme_bw(base_size = 11)
  })

  output$slope_hist_plot <- renderPlot({
    req(rv$pairs_clean)
    dat <- rv$pairs_clean %>% filter(is.finite(Slope_deg))
    if (nrow(dat) == 0)
      return(ggplot() + labs(title = "No slope data (terrain not fetched)") + theme_bw())
    ggplot(dat, aes(Slope_deg)) +
      geom_histogram(fill = "#c0392b", colour = "white", bins = 20, alpha = 0.8) +
      labs(title   = "Slope at Fire Centroids",
           x = "Slope (degrees)", y = "Number of pairs") +
      theme_bw(base_size = 11)
  })

  output$wind_align_plot <- renderPlot({
    req(rv$pairs_clean)
    dat <- rv$pairs_clean %>%
      filter(is.finite(WindDir_deg), is.finite(Aspect_deg))
    if (nrow(dat) == 0)
      return(ggplot() + labs(title = "No aspect data (terrain not fetched)") + theme_bw())
    dat <- dat %>% mutate(
      upslope   = (Aspect_deg + 180) %% 360,
      align_deg = ((WindDir_deg - upslope + 180) %% 360) - 180
    )
    ggplot(dat, aes(align_deg)) +
      geom_histogram(fill = "#2c5f2e", colour = "white", bins = 24, alpha = 0.85) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
      scale_x_continuous(breaks = seq(-180, 180, 60),
                         limits = c(-180, 180)) +
      labs(title   = "Wind-to-Upslope Alignment",
           x = "Angle from upslope direction (degrees)",
           y = "Number of pairs") +
      theme_bw(base_size = 11)
  })

  output$slope_map <- renderTmap({
    req(rv$slope_r, rv$cal_vect_proj)
    sl <- rv$slope_r
    n_cells <- nrow(sl) * ncol(sl)
    if (n_cells > 800000L) {
      fact <- ceiling(sqrt(n_cells / 800000L))
      sl   <- aggregate(sl, fact = fact, fun = "mean", na.rm = TRUE)
    }
    bnd <- project(rv$cal_vect_proj, crs(sl))
    tm_shape(sl) +
      tm_raster(
        col        = names(sl)[1],
        col.scale  = tm_scale_intervals(
          values   = "brewer.yl_or_rd",
          n        = 7,
          style    = "pretty",
          value.na = NA
        ),
        col.legend = tm_legend(title = "Slope (deg)"),
        col_alpha  = 0.85
      ) +
      tm_shape(bnd) + tm_borders(col = "gray30", lwd = 1.2)
  })

  output$aspect_map <- renderTmap({
    req(rv$aspect_r, rv$cal_vect_proj)
    asp <- rv$aspect_r
    n_cells <- nrow(asp) * ncol(asp)
    if (n_cells > 800000L) {
      fact <- ceiling(sqrt(n_cells / 800000L))
      asp  <- aggregate(asp, fact = fact, fun = "mean", na.rm = TRUE)
    }
    bnd <- project(rv$cal_vect_proj, crs(asp))
    tm_shape(asp) +
      tm_raster(
        col        = names(asp)[1],
        col.scale  = tm_scale_intervals(
          values   = "brewer.rd_yl_bu",
          breaks   = c(0, 45, 90, 135, 180, 225, 270, 315, 360),
          labels   = c("N","NE","E","SE","S","SW","W","NW"),
          value.na = NA
        ),
        col.legend = tm_legend(title = "Aspect"),
        col_alpha  = 0.85
      ) +
      tm_shape(bnd) + tm_borders(col = "gray30", lwd = 1.2)
  })

  # ---------------------------------------------------------------------------
  # Diagnostics tab outputs
  # ---------------------------------------------------------------------------

  output$spread_prob_fwi_plot <- renderPlot({
    req(rv$spread_prob_model, rv$spread_prob_samples)
    samp    <- rv$spread_prob_samples
    m       <- rv$spread_prob_model
    fwi_seq <- seq(0, max(samp$FWI, na.rm = TRUE) * 1.05, length.out = 200)
    nd <- data.frame(
      FWI           = fwi_seq,
      WindSpeed_kmh = mean(samp$WindSpeed_kmh, na.rm = TRUE),
      FineFuels     = mean(samp$FineFuels,     na.rm = TRUE)
    )
    pr  <- predict(m, nd, type = "response", se.fit = TRUE)
    nd$fit   <- pr$fit
    nd$lower <- pmax(0, pr$fit - 1.96 * pr$se.fit)
    nd$upper <- pmin(1, pr$fit + 1.96 * pr$se.fit)
    ggplot(nd, aes(FWI, fit)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2, fill = "#c0392b") +
      geom_line(color = "#c0392b", linewidth = 1) +
      geom_rug(data = samp, aes(x = FWI, y = NULL),
               sides = "b", alpha = 0.04) +
      scale_y_continuous(limits = c(0, 1),
                         labels = scales::percent_format(accuracy = 1)) +
      labs(x = "FWI", y = "P(spread to neighbor cell)",
           subtitle = sprintf("At mean wind %.1f km/h, mean fine fuels %.2f",
                              mean(samp$WindSpeed_kmh, na.rm = TRUE),
                              mean(samp$FineFuels,     na.rm = TRUE))) +
      theme_bw(base_size = 11)
  })

  output$spread_prob_wind_plot <- renderPlot({
    req(rv$spread_prob_model, rv$spread_prob_samples)
    samp     <- rv$spread_prob_samples
    m        <- rv$spread_prob_model
    wind_seq <- seq(0, max(samp$WindSpeed_kmh, na.rm = TRUE) * 1.05, length.out = 200)
    nd <- data.frame(
      FWI           = mean(samp$FWI,       na.rm = TRUE),
      WindSpeed_kmh = wind_seq,
      FineFuels     = mean(samp$FineFuels, na.rm = TRUE)
    )
    pr  <- predict(m, nd, type = "response", se.fit = TRUE)
    nd$fit   <- pr$fit
    nd$lower <- pmax(0, pr$fit - 1.96 * pr$se.fit)
    nd$upper <- pmin(1, pr$fit + 1.96 * pr$se.fit)
    ggplot(nd, aes(WindSpeed_kmh, fit)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2, fill = "#2c5f2e") +
      geom_line(color = "#2c5f2e", linewidth = 1) +
      geom_rug(data = samp, aes(x = WindSpeed_kmh, y = NULL),
               sides = "b", alpha = 0.04) +
      scale_y_continuous(limits = c(0, 1),
                         labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Wind Speed (km/h)", y = "P(spread to neighbor cell)",
           subtitle = sprintf("At mean FWI %.1f, mean fine fuels %.2f",
                              mean(samp$FWI,       na.rm = TRUE),
                              mean(samp$FineFuels, na.rm = TRUE))) +
      theme_bw(base_size = 11)
  })

  output$spread_roc_plot <- renderPlot({
    req(rv$spread_prob_model, rv$spread_prob_samples)
    samp  <- rv$spread_prob_samples
    pred  <- predict(rv$spread_prob_model, samp, type = "response")
    ord   <- order(pred, decreasing = TRUE)
    y     <- samp$y[ord]
    n_pos <- sum(y == 1);  n_neg <- sum(y == 0)
    tpr   <- cumsum(y == 1) / n_pos
    fpr   <- cumsum(y == 0) / n_neg
    auc   <- sum(diff(c(0, fpr)) * (c(0, tpr[-length(tpr)]) + tpr) / 2)
    roc_df <- data.frame(fpr = c(0, fpr, 1), tpr = c(0, tpr, 1))
    ggplot(roc_df, aes(fpr, tpr)) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray60") +
      geom_line(color = "#c0392b", linewidth = 1) +
      annotate("text", x = 0.65, y = 0.12,
               label = sprintf("AUC = %.3f", auc), size = 4) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "False Positive Rate", y = "True Positive Rate") +
      theme_bw(base_size = 11)
  })

  output$max_area_scatter_plot <- renderPlot({
    req(rv$spread_da_model, rv$spread_da_data)
    dat <- rv$spread_da_data
    dat$pred_ha <- expm1(predict(rv$spread_da_model, dat))
    ggplot(dat, aes(daily_area_ha, pred_ha, color = FWI)) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray60") +
      geom_point(alpha = 0.7, size = 2.5) +
      scale_color_gradient(low = "#ffffb2", high = "#d7191c", name = "FWI") +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      labs(x = "Observed (ha/day)", y = "Predicted (ha/day)",
           subtitle = "Log scale \u2014 dashed = perfect fit") +
      theme_bw(base_size = 11)
  })

  output$pairs_table <- renderDT({
    req(rv$pairs_clean)
    cols <- intersect(
      c("FIRE_ID", "date", "gap_days", "daily_area_ha", "FWI",
        "WindSpeed_kmh", "WindDir_deg", "Slope_deg", "Aspect_deg",
        "EffectiveWind", "use_for_spread_prob", "use_for_daily_area"),
      names(rv$pairs_clean)
    )
    rv$pairs_clean %>%
      select(all_of(cols)) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      datatable(options = list(pageLength = 15, scrollX = TRUE),
                rownames = FALSE, class = "table-sm table-striped")
  })

  # ---------------------------------------------------------------------------
  # Candidate grid
  # ---------------------------------------------------------------------------

  observeEvent(input$build_candidates, {
    b0_seq <- seq(input$cand_b0_min, input$cand_b0_max, by = input$cand_b0_step)
    b1_seq <- seq(input$cand_b1_min, input$cand_b1_max, by = input$cand_b1_step)

    b3_default <- if (!is.null(rv$spread_prob_coef))
      rv$spread_prob_coef$estimate[grepl("B3", rv$spread_prob_coef$term)][1]
    else 0.02

    grid <- expand.grid(B0 = b0_seq, B1_FWI = b1_seq) %>%
      mutate(
        B2_FineFuels   = 0,
        B3_WindSpeed   = b3_default,
        MaxAreaB0      = if (!is.null(rv$max_area_coef)) rv$max_area_coef$estimate[1] else 0,
        MaxAreaB1_FWI  = if (!is.null(rv$max_area_coef)) rv$max_area_coef$estimate[2] else 10,
        MaxAreaB2_Wind = if (!is.null(rv$max_area_coef)) rv$max_area_coef$estimate[3] else 5,
        candidate_id   = row_number()
      )
    rv$candidate_grid <- grid

    output$candidate_table <- renderDT({
      grid %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        datatable(options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE, class = "table-sm table-striped")
    })
  })

  # ---------------------------------------------------------------------------
  # Candidate scoring
  # ---------------------------------------------------------------------------

  observeEvent(input$run_scoring, {
    req(rv$candidate_grid, rv$pairs_clean)
    req(nchar(trimws(input$scf_runs_root)) > 0)

    withProgress(message = "Scoring candidates...", {
      rv$scores_df <- score_candidates(
        candidate_grid = rv$candidate_grid,
        runs_root      = input$scf_runs_root,
        cell_area_ha   = rv$cell_area_ha,
        cal_start_year = input$cal_years[1],
        cal_end_year   = input$cal_years[2],
        w_size         = input$score_w_size,
        w_aab          = input$score_w_aab,
        fpa_sizes      = if (!is.null(rv$ign_df)) rv$ign_df$FIRE_SIZE_HA else NULL
      )
    })

    output$scores_table <- renderDT({
      rv$scores_df %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        datatable(options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE, class = "table-sm table-striped")
    })
    output$cdf_plot <- renderPlot({
      req(rv$ign_df)
      plot_candidate_cdf(rv$scores_df, rv$ign_df$FIRE_SIZE_HA,
                         input$scf_runs_root, rv$cell_area_ha,
                         input$cal_years[1], input$cal_years[2])
    })
    output$annual_plot <- renderPlot({
      plot_candidate_annual(rv$scores_df, input$scf_runs_root, rv$cell_area_ha,
                            input$cal_years[1], input$cal_years[2])
    })
  })

  # ---------------------------------------------------------------------------
  # Downloads
  # ---------------------------------------------------------------------------

  output$dl_ign_coef <- downloadHandler(
    filename = "ignition_coefficients.csv",
    content  = function(f) write_csv(rv$ign_coef, f))

  output$dl_startup <- downloadHandler(
    filename = "era_startup_values.csv",
    content  = function(f) write_csv(rv$startup_df, f))

  output$dl_lightning_map <- downloadHandler(
    filename = "Lightning_Ignition_Map.tif",
    content  = function(f) writeRaster(rv$surf_lightning,  f, datatype = "INT4S", overwrite = TRUE))

  output$dl_accidental_map <- downloadHandler(
    filename = "Accidental_Ignition_Map.tif",
    content  = function(f) writeRaster(rv$surf_accidental, f, datatype = "INT4S", overwrite = TRUE))

  output$dl_rx_map <- downloadHandler(
    filename = "Rx_Ignition_Map.tif",
    content  = function(f) writeRaster(rv$surf_rx,         f, datatype = "INT4S", overwrite = TRUE))

  # ---- Terrain maps -------------------------------------------------------
  observeEvent(input$gen_terrain_maps, {
    req(rv$template_r, rv$working_crs, rv$template_mask)

    output$terrain_maps_status <- renderUI(
      tags$span(class = "text-muted small", icon("spinner", class = "fa-spin"),
                " Fetching DEM and computing terrain maps...")
    )

    withProgress(message = "Generating terrain maps...", value = 0, {
      setProgress(0.1, detail = "Fetching DEM via elevatr...")

      template_ext_v <- vect(ext(rv$template_r), crs = rv$working_crs)

      terrain_result <- tryCatch(
        fetch_terrain(
          cal_vect    = template_ext_v,
          working_crs = rv$working_crs,
          z           = 10,
          progress_fn = function(msg) setProgress(0.3, detail = msg)
        ),
        error = function(e) {
          output$terrain_maps_status <- renderUI(
            tags$span(class = "text-danger small", icon("xmark"),
                      " Terrain fetch failed: ", conditionMessage(e))
          )
          NULL
        }
      )
      req(terrain_result)

      setProgress(0.75, detail = "Resampling to LANDIS template grid...")

      slope_res  <- resample(terrain_result$slope,  rv$template_r, method = "bilinear")
      aspect_res <- resample(terrain_result$aspect, rv$template_r, method = "bilinear")

      # Mask to active template cells (NA outside park boundary)
      slope_res  <- mask(slope_res,  rv$template_mask)
      aspect_res <- mask(aspect_res, rv$template_mask)

      # Ground slope: degrees, rounded to nearest integer
      gs <- round(slope_res)
      names(gs) <- "slope_deg"
      rv$ground_slope_r <- gs

      # Uphill azimuth: (aspect + 180) %% 360, clamped 0-359
      us <- clamp(round((aspect_res + 180) %% 360), lower = 0, upper = 359)
      names(us) <- "uphill_azimuth_deg"
      rv$uphill_slope_r <- us

      setProgress(1.0)
      output$terrain_maps_status <- renderUI(
        tags$span(class = "text-success small", icon("check"),
                  " Terrain maps ready for download.")
      )
    })
  })

  output$ground_slope_preview <- renderTmap({
    req(rv$ground_slope_r)
    sl <- rv$ground_slope_r
    if (ncell(sl) > 800000L)
      sl <- aggregate(sl, ceiling(sqrt(ncell(sl) / 800000L)), fun = "mean", na.rm = TRUE)
    tm_shape(sl) +
      tm_raster(
        col        = names(sl)[1],
        col.scale  = tm_scale_intervals(
          values   = "brewer.yl_or_rd",
          n        = 7,
          style    = "pretty",
          value.na = NA
        ),
        col.legend = tm_legend(title = "Slope (deg)")
      )
  })

  output$uphill_slope_preview <- renderTmap({
    req(rv$uphill_slope_r)
    az <- rv$uphill_slope_r
    if (ncell(az) > 800000L)
      az <- aggregate(az, ceiling(sqrt(ncell(az) / 800000L)), fun = "mean", na.rm = TRUE)
    tm_shape(az) +
      tm_raster(
        col        = names(az)[1],
        col.scale  = tm_scale_intervals(
          values   = "brewer.rd_yl_bu",
          breaks   = c(0, 45, 90, 135, 180, 225, 270, 315, 360),
          labels   = c("N","NE","E","SE","S","SW","W","NW"),
          value.na = NA
        ),
        col.legend = tm_legend(title = "Uphill Azimuth")
      )
  })

  output$dl_ground_slope <- downloadHandler(
    filename = "GroundSlope.tif",
    content  = function(f) {
      req(rv$ground_slope_r)
      writeRaster(rv$ground_slope_r, f, datatype = "INT2S", overwrite = TRUE)
    })

  output$dl_uphill_slope <- downloadHandler(
    filename = "UphillSlope.tif",
    content  = function(f) {
      req(rv$uphill_slope_r)
      writeRaster(rv$uphill_slope_r, f, datatype = "INT2S", overwrite = TRUE)
    })

  # ---- Blank suppression maps (0.0 everywhere on template grid) -----------
  make_blank_suppression <- function() {
    req(rv$template_mask)
    r <- rv$template_mask * 0.0
    names(r) <- "suppression"
    r
  }

  output$dl_supp_accidental <- downloadHandler(
    filename = "Suppression_Accidental_3Zones.tif",
    content  = function(f) writeRaster(make_blank_suppression(), f,
                                       datatype = "FLT4S", overwrite = TRUE))

  output$dl_supp_lightning <- downloadHandler(
    filename = "Suppression_Lightning_3Zones.tif",
    content  = function(f) writeRaster(make_blank_suppression(), f,
                                       datatype = "FLT4S", overwrite = TRUE))

  output$dl_supp_rx <- downloadHandler(
    filename = "Suppression_Rx_3Zones.tif",
    content  = function(f) writeRaster(make_blank_suppression(), f,
                                       datatype = "FLT4S", overwrite = TRUE))

  output$dl_spread_coef <- downloadHandler(
    filename = "spread_coefficients.csv",
    content  = function(f) {
      write_csv(bind_rows(
        rv$spread_prob_coef %>% mutate(model = "SpreadProbability"),
        rv$max_area_coef    %>% mutate(model = "MaxDailySpreadArea")
      ), f)
    })

  output$dl_pairs <- downloadHandler(
    filename = "geomac_clean_pairs.csv",
    content  = function(f) write_csv(rv$pairs_clean, f))

  output$dl_candidates <- downloadHandler(
    filename = "candidate_grid.csv",
    content  = function(f) write_csv(rv$candidate_grid, f))

  output$dl_snippets_zip <- downloadHandler(
    filename = "scf_param_snippets.zip",
    content  = function(f) {
      req(rv$candidate_grid)
      snip_dir <- file.path(tempdir(), "snippets")
      dir.create(snip_dir, showWarnings = FALSE)
      walk(seq_len(nrow(rv$candidate_grid)), function(i) {
        row <- rv$candidate_grid[i, ]
        writeLines(write_scf_snippet(row),
                   file.path(snip_dir, sprintf("candidate_%04d.txt", row$candidate_id)))
      })
      zip(f, files = list.files(snip_dir, full.names = TRUE), flags = "-j")
    })

  output$dl_scores <- downloadHandler(
    filename = "candidate_scores.csv",
    content  = function(f) write_csv(rv$scores_df, f))

  # ===========================================================================
  # Species Table tab
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Step 1: Load species lookup CSV
  # ---------------------------------------------------------------------------
  observeEvent(input$load_lookup_csv, {
    req(nchar(trimws(input$lookup_csv_path)) > 0)
    warn_msgs <- character(0)
    tryCatch({
      lkup <- withCallingHandlers(
        load_species_lookup(input$lookup_csv_path),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      rv$spp_lookup <- lkup

      # Build status UI: success banner + optional SPCD mismatch alert
      status_ui <- tagList(
        tags$span(class = "text-success small", icon("check"),
                  sprintf(" Loaded: %d species, %d with FIA SPCDs.",
                          nrow(lkup), sum(!is.na(lkup$FIA_SPCD))))
      )

      if (length(warn_msgs) > 0) {
        # Parse mismatch lines out of the warning message
        mismatch_lines <- warn_msgs[1] %>%
          strsplit("\n") %>% `[[`(1) %>%
          grep("•", ., value = TRUE, fixed = TRUE) %>%
          trimws()

        status_ui <- tagList(
          status_ui,
          br(),
          tags$div(
            class = "alert alert-warning p-2 mt-2",
            style = "font-size:0.82rem;",
            tags$b(icon("triangle-exclamation"),
                   " Possible SPCD mismatches — verify before fetching FIA data:"),
            tags$ul(style = "margin-bottom:0; padding-left:1.2em;",
                    lapply(mismatch_lines, tags$li)),
            tags$small(
              "Check against: ",
              tags$a("FIA REF_SPECIES.csv",
                     href = "https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES.csv",
                     target = "_blank")
            )
          )
        )
      }

      output$lookup_load_status <- renderUI(status_ui)

    }, error = function(e) {
      output$lookup_load_status <- renderUI(
        tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
      )
    })
  })

  # ---------------------------------------------------------------------------
  # Step 2: Build species table from lookup + optional filter
  # ---------------------------------------------------------------------------
  observeEvent(input$build_spp_table, {
    req(rv$spp_lookup)

    tryCatch({
      lkup <- rv$spp_lookup

      codes <- switch(input$spp_filter_mode,

        "all" = lkup$SpeciesName,

        "necn" = {
          req(nchar(trimws(input$necn_csv_path)) > 0)
          necn_codes <- parse_necn_species(input$necn_csv_path)
          # Keep only codes that appear in the lookup; warn about unmatched
          matched   <- intersect(necn_codes, lkup$SpeciesName)
          unmatched <- setdiff(necn_codes, lkup$SpeciesName)
          if (length(unmatched) > 0) {
            showNotification(
              paste0("NECN codes not found in lookup (will use defaults): ",
                     paste(unmatched, collapse = ", ")),
              type = "warning", duration = 8
            )
            # Still include unmatched so user can edit manually
            c(matched, unmatched)
          } else {
            matched
          }
        },

        "manual" = {
          raw <- input$manual_species_codes
          req(nchar(trimws(raw)) > 0)
          codes <- unique(trimws(strsplit(raw, "[\n,;]+")[[1]]))
          codes[nchar(codes) > 0]
        }
      )

      if (length(codes) == 0) {
        output$spp_load_status <- renderUI(
          tags$span(class = "text-warning small", "No species codes found.")
        )
        return()
      }

      rv$spp_codes <- codes
      rv$spp_table <- build_default_species_table(codes, lookup = lkup)

      output$spp_load_status <- renderUI(
        tags$span(class = "text-success small", icon("check"),
                  sprintf(" Table built: %d species (%d mapped to FIA SPCDs).",
                          length(codes),
                          sum(!is.na(rv$spp_table$FIA_SPCD))))
      )
    }, error = function(e) {
      output$spp_load_status <- renderUI(
        tags$span(class = "text-danger small", icon("xmark"), " ", conditionMessage(e))
      )
    })
  })

  # ---------------------------------------------------------------------------
  # Step 3: Fetch FIA bark parameters
  # ---------------------------------------------------------------------------
  observeEvent(input$fetch_fia_bark, {
    req(rv$spp_codes, rv$spp_lookup)
    req(length(input$fia_states) > 0)

    withProgress(message = "Loading FIA data...", value = 0, {
      tryCatch({

        # Load / download FIA TREE data
        # rv$fia_trees is invalidated when states change (see observer below)
        fia_dir_val <- trimws(input$fia_local_dir)
        dir_arg     <- if (nchar(fia_dir_val) > 0) fia_dir_val else NULL

        incProgress(0.1, detail = paste("States:",
                                        paste(input$fia_states, collapse = ", ")))

        # Re-use cached data if states haven't changed, otherwise reload
        if (is.null(rv$fia_trees) ||
            !identical(sort(attr(rv$fia_trees, "states")), sort(input$fia_states))) {
          rv$fia_trees <- load_fia_trees(states = input$fia_states,
                                          fia_dir = dir_arg)
          attr(rv$fia_trees, "states") <- input$fia_states
        }

        incProgress(0.7, detail = "Estimating bark parameters...")
        params <- estimate_bark_params_from_fia(
          species_codes = rv$spp_codes,
          fia_trees     = rv$fia_trees,
          lookup        = rv$spp_lookup
        )
        rv$spp_table <- params

        incProgress(0.2, detail = "Done.")
        n_fia <- sum(!is.na(params$FIA_SPCD) &
                       !grepl("defaults used", params$Note))
        output$fia_fetch_status <- renderUI(
          tags$span(class = "text-success small", icon("check"),
                    sprintf(" Done: FIA parameters estimated for %d of %d species.",
                            n_fia, nrow(params)))
        )
      }, error = function(e) {
        output$fia_fetch_status <- renderUI(
          tags$span(class = "text-danger small", icon("xmark"),
                    " FIA fetch failed: ", conditionMessage(e))
        )
      })
    })
  })

  # ---------------------------------------------------------------------------
  # Render editable DT
  # ---------------------------------------------------------------------------
  output$spp_table_dt <- renderDT({
    req(rv$spp_table)
    dt <- rv$spp_table %>%
      select(SpeciesCode, AgeDBH, MaximumBarkThickness, ScientificName, FIA_SPCD, Note) %>%
      mutate(AgeDBH = as.integer(AgeDBH),
             MaximumBarkThickness = round(MaximumBarkThickness, 2))

    datatable(
      dt,
      # Allow editing only AgeDBH (col 1) and MaximumBarkThickness (col 2)
      editable  = list(target = "cell",
                       disable = list(columns = c(0L, 3L, 4L, 5L))),
      rownames  = FALSE,
      selection = "none",
      options   = list(
        pageLength = 30,
        scrollX    = TRUE,
        columnDefs = list(
          list(width = "160px", targets = 0),  # SpeciesCode
          list(width =  "80px", targets = 1),  # AgeDBH
          list(width = "140px", targets = 2),  # MaxBT
          list(width = "180px", targets = 3),  # ScientificName
          list(width =  "70px", targets = 4),  # FIA_SPCD
          list(width = "360px", targets = 5)   # Note
        )
      ),
      colnames = c("SpeciesCode", "AgeDBH", "MaxBarkThickness (cm)",
                   "Scientific Name", "FIA SPCD", "Note")
    ) %>%
      formatStyle(c("AgeDBH", "MaximumBarkThickness"),
                  backgroundColor = "#fffde7", cursor = "pointer") %>%
      formatRound("MaximumBarkThickness", digits = 2)
  })

  # Write cell edits back to rv$spp_table
  observeEvent(input$spp_table_dt_cell_edit, {
    info     <- input$spp_table_dt_cell_edit
    # DT col index (0-based): 0=SpeciesCode,1=AgeDBH,2=MaxBT,3=Sci,4=SPCD,5=Note
    col_name <- c("SpeciesCode", "AgeDBH", "MaximumBarkThickness",
                  "ScientificName", "FIA_SPCD", "Note")[info$col + 1L]
    tbl      <- rv$spp_table
    if (col_name == "AgeDBH") {
      tbl[info$row, col_name] <- as.integer(round(as.numeric(info$value)))
    } else if (col_name == "MaximumBarkThickness") {
      tbl[info$row, col_name] <- round(as.numeric(info$value), 3)
    }
    rv$spp_table <- tbl
  })

  # ---------------------------------------------------------------------------
  # Bark thickness curve plot
  # ---------------------------------------------------------------------------
  output$bark_curve_plot <- renderPlot({
    req(rv$spp_table, nrow(rv$spp_table) > 0)
    curve_df <- bark_curve_data(rv$spp_table, max_age = 600)

    ggplot(curve_df, aes(x = Age, y = BarkThickness_cm,
                         colour = SpeciesCode, group = SpeciesCode)) +
      geom_line(linewidth = 0.8) +
      geom_hline(aes(yintercept = MaximumBarkThickness, colour = SpeciesCode),
                 linetype = "dashed", linewidth = 0.3, alpha = 0.5) +
      labs(
        x       = "Stand Age (years)",
        y       = "Bark Thickness (cm, one side)",
        colour  = "Species",
        title   = "SCF Eq. 10 — Bark Thickness vs Age",
        caption = "Dashed lines = MaximumBarkThickness asymptote. AgeDBH = age at half-maximum bark thickness."
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right",
            plot.caption = element_text(size = 8, colour = "grey50"))
  })

  # ---------------------------------------------------------------------------
  # Lookup table viewer
  # ---------------------------------------------------------------------------
  output$fia_lookup_dt <- renderDT({
    req(rv$spp_lookup)
    datatable(
      rv$spp_lookup,
      rownames = FALSE,
      options  = list(pageLength = 30, scrollX = TRUE),
      colnames = c("LANDIS Code", "FIA SPCD", "Scientific Name")
    )
  })

  # ---------------------------------------------------------------------------
  # Download handlers
  # ---------------------------------------------------------------------------
  output$dl_fire_spp_table <- downloadHandler(
    filename = "Fire_Spp_Table.csv",
    content  = function(f) {
      req(rv$spp_table)
      rv$spp_table %>%
        select(SpeciesCode, AgeDBH, MaximumBarkThickness) %>%
        write_csv(f)
    }
  )

  output$dl_spp_table_full <- downloadHandler(
    filename = "Fire_Spp_Table_annotated.csv",
    content  = function(f) {
      req(rv$spp_table)
      write_csv(rv$spp_table, f)
    }
  )
}

shinyApp(ui, server)
