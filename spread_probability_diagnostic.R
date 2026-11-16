# =============================================================================
# spread_probability_diagnostic.R
#
# Step-by-step diagnostic script for the SCF spread probability fitting pipeline.
# Run each section interactively in RStudio (Ctrl+Enter / by section) to inspect
# intermediate outputs and identify what is going wrong with the fit.
#
# Climate data processing exactly mirrors the app:
#   - ERA-Interim FWI:  extract_era_timeseries() on fwi.nc
#   - Wind:             extract_gridmet_wind()  (remote, no local files needed)
#                  OR   extract_era5_wind()     (local ERA5 u10/v10 NetCDF files)
#   - Terrain (EWS):    fetch_terrain()         (elevatr DEM -> slope/aspect)
#   - Fine fuels:       fetch_landfire_fine_fuels()  (LANDFIRE FBFM40, remote)
#                  OR   load_local_fine_fuels()      (user-provided 0-1 TIF)
#
# NOTE: spread rasterization uses a spread_template built from the calibration
# boundary extent (not the LANDIS template_r), exactly as the app does.
#
# All helper functions are loaded directly from the app's modules/ folder.
# =============================================================================


# =============================================================================
# SECTION 0 — Packages
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, ggplot2, tidyr, purrr, stringr,
               lubridate, glue, broom, scales, sf, future, furrr, progressr,
               quantreg)   # for quantile regression in MaxSpreadArea section

for (.pkg in c("climateR", "elevatr")) {
  if (!requireNamespace(.pkg, quietly = TRUE))
    message(sprintf("Optional: install.packages('%s') if using GRIDMET or terrain.", .pkg))
}
rm(.pkg)

# tigris: US Census TIGER/Line state boundaries (auto-downloaded when
# CAL_VECT_SHP = NULL and CAL_VECT_STATES is set).
if (!requireNamespace("tigris", quietly = TRUE)) {
  message("Installing tigris for automatic state boundary download...")
  install.packages("tigris")
}
options(tigris_use_cache = TRUE)   # cache downloaded shapefiles locally


# =============================================================================
# SECTION 1 — *** EDIT THESE PATHS AND SETTINGS ***
# =============================================================================

# ---- Spatial inputs ---------------------------------------------------------
# ---- Perimeter data source --------------------------------------------------
# Choose ONE:
#   "geomac" — Historic GeoMAC fire perimeters (2000–2018 GDB)
#   "fired"  — FIRED daily MODIS-based fire progression (recommended)
#
# FIRED download (one-time, ~500 MB):
#   https://scholar.colorado.edu/concern/datasets/r207tp32d
#   Balch et al. (2020) Remote Sensing 12(21): 3498
#   Download the CONUS daily GeoPackage (.gpkg) file and set FIRED_GPKG below.
PERIMETER_SOURCE <- "fired"    # "geomac" or "fired"

GEOMAC_GDB   <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Code/calibration/scf_calibration/data_raw/GEOMAC/Historic_Geomac_Perimeters_All_Years_2000_2018.gdb"

# TWO separate FIRED files exist — make sure you have the right one:
#
#   DAILY product  (needed for spread probability — one row per fire-day):
#     Download: https://scholar.colorado.edu/concern/datasets/765372341
#     Filename: fired_daily_conus_nov2001-jan2019.gpkg  (or similar)
#
#   EVENTS product (one row per fire, final perimeter — useful for MaxSpreadArea):
#     Already downloaded:
#     fired_events_conus_nov2001-jan2019.gpkg
#
# If PERIMETER_SOURCE = "fired" and FIRED_GPKG points to the events file,
# load_fired_perimeters() will detect this and stop with a helpful message.

# Path to the FIRED CONUS DAILY GeoPackage (needs to be downloaded — see above)
FIRED_GPKG        <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/FIRED_CONUS/daily/FIRED_CONUS_Daily_w_README/fired_conus_daily_nov2001-jan2019.gpkg"

# Path to the FIRED CONUS EVENTS GeoPackage (already downloaded — for MaxSpreadArea)
FIRED_EVENTS_GPKG <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/FIRED_CONUS/events/fired_events_conus_nov2001-jan2019.gpkg"

# FIRED layer name within the daily GeoPackage, or NULL to auto-detect
FIRED_LAYER       <- NULL   # auto-detect (looks for layer name containing "daily")

# Minimum final fire size (ha) for spread probability pairs.
# 100 ha removes MODIS artefacts while retaining calibration-relevant fires.
FIRED_MIN_FIRE_HA <- 100

# Number of CPU cores for parallel steps (FIRED cumulative union + Section 12 pair loop).
# NULL = auto (all physical cores minus 1).
# Set to 1 to run sequentially (e.g. for debugging).
#parallel::detectCores()
N_CORES <- 20
TEMPLATE_TIF <- "/Users/jlamping/Desktop/LANDIS/MORA/CoolDry/Rep1/input_maps/MORA_ClimateRegions_09122025.tif"
# Calibration boundary — TWO options (set exactly one):
#
# OPTION A: Load from a local shapefile / GeoPackage / GeoJSON.
#   Set CAL_VECT_SHP to the full file path and CAL_VECT_STATES to NULL.
#
# OPTION B: Automatically download US state boundaries via the tigris package
#   (Census TIGER/Line cartographic files, cached locally after first download).
#   Set CAL_VECT_SHP to NULL and list state names in CAL_VECT_STATES.
#   Multiple states are dissolved into a single polygon.
#   Example values: "Oregon", "Washington", "Idaho", "Montana", "California"
#
# Active option B — OR + WA:
CAL_VECT_SHP    <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/S_USA.ECOSYS_ECOMAPPROVINCES_2025/ECOSYS_CoastalRange.shp" # set to a file path to use Option A
CAL_VECT_STATES <- NULL #c("Oregon", "Washington")     # ignored when CAL_VECT_SHP is not NULL

# ---- ERA-Interim FWI NetCDF -------------------------------------------------
# Vitolo et al. 2019: https://doi.org/10.5281/zenodo.1065400
# Download: fwi.nc  (daily FWI, 1980-2018, global, ERA-Interim grid)
ERA_FWI_NC  <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/fwi.nc"

# ---- Wind source — choose ONE option ----------------------------------------

# OPTION A (recommended): GRIDMET remote fetch — no local files needed.
# Works for CONUS only. Fetches from THREDDS server via the climateR package.
WIND_SOURCE <- "gridmet"   # set to "gridmet" or "era5"

# OPTION B: Local ERA5 u10/v10 NetCDF files from Copernicus CDS.
# Leave as NULL if using GRIDMET.
ERA5_U10_NC <- NULL  # e.g. "/path/to/adaptor.mars.internal-u10.nc"
ERA5_V10_NC <- NULL  # e.g. "/path/to/adaptor.mars.internal-v10.nc"
# If both u10 and v10 are in one combined file, set ERA5_U10_NC to that file
# and leave ERA5_V10_NC as NULL — extract_era5_wind() will auto-detect both vars.

# ---- Calibration date range -------------------------------------------------
CAL_START <- as.Date("2000-01-01")
CAL_END   <- as.Date("2018-12-31")

# ---- Terrain (for Effective Wind Speed) -------------------------------------
# Set USE_TERRAIN = TRUE to fetch a DEM via elevatr and compute slope/aspect.
# Set to FALSE to use raw wind speed as EWS (same as app when terrain not loaded).
USE_TERRAIN         <- TRUE
TERRAIN_ZOOM        <- 9L    # elevatr zoom: 8 ≈ 600 m resolution, 9 ≈ 300 m

# Combustion buoyancy: 10 = low (SCF default), 25 = moderate, 50 = high.
# Only matters when USE_TERRAIN = TRUE.
COMBUSTION_BUOYANCY <- 10.0

# ---- Pair-building settings (match app defaults) ----------------------------
MAX_GAP_SP  <- 5L    # max days between perimeters for spread-prob pairs
MAX_GAP_DA  <- 5L    # max days between perimeters for daily-area pairs
NEG_TOL_HA  <- 10.0   # tolerate up to this many ha of shrinkage (digitizing noise)

# ---- Sample settings --------------------------------------------------------
FAILURE_RATIO <- 3L       # no-spread cells sampled per spread-success cell
MAX_SAMPLES   <- 25000L   # total cell cap across the GLM training set

# ---- Dilation radius (frontier ring width) ----------------------------------
# ROOK_W (the default 1-cell dilation) only captures cells directly adjacent to
# the t-perimeter.  When fires grow rapidly between perimeter observations — as
# is common in the PNW — the t+1 perimeter can expand far beyond this thin ring
# and no frontier cells register as "burned", giving n_succ=0.
#
# DILATION_CELLS controls how many cells wide the frontier ring is:
#   1  = original behaviour (1-cell rook neighbourhood, ~1 cell width)
#   3  = moderate expansion (~3-cell ring; good starting point for rapid-growth fires)
#   5  = wide ring (useful for coarse grids or fires with large daily increments)
#
# The dilation kernel is built as a disk (cells within DILATION_CELLS of center
# are included) so the ring is roughly circular rather than square.
DILATION_CELLS <- 5L

# ---- Spread rasterization grid ----------------------------------------------
# The app builds a separate spread_template from the calibration boundary at a
# user-specified resolution (default = LANDIS template resolution).  Perimeters
# are rasterized and fine fuels are resampled onto this grid — NOT template_r.
# Set to NULL to auto-detect from template_r resolution.
SPREAD_RES_M <- NULL   # e.g. 240 for 240 m cells; NULL = use template_r res

# ---- Fine fuels source ------------------------------------------------------
# Mirrors the radio button in the app's Spread tab:
#   "landfire" — fetch LANDFIRE FBFM40 remotely via the USGS ImageServer REST API
#                (reclassified using Scott & Burgan 2005 fuel load table)
#   "local"    — load a user-supplied raster
#   "none"     — placeholder FF=1.0 everywhere (B2_FineFuels will be unidentifiable)
FINE_FUELS_SOURCE <- "landfire"   # "landfire", "local", or "none"

# Path to local fine fuels raster — only used when FINE_FUELS_SOURCE = "local"
FINE_FUELS_TIF <- NULL  # e.g. "/path/to/fbfm40.tif"

# ---- Fine fuels normalization ------------------------------------------------
# TRUE  (app default): values divided by the maximum possible load so the raster
#        is scaled 0-1.  B2 coefficient is dimensionless ("per unit FF index").
#        For LANDFIRE FBFM40 the divisor is 9.0 t/ac (GR9 fuel model).
#        For "local" the divisor is the raster's own maximum value.
#
# FALSE: keep raw fuel load values in their native units.
#        For LANDFIRE FBFM40 this is tons/acre (1-hr + 10-hr dead fuel,
#        Scott & Burgan 2005).  B2 will be in units of "per ton/acre" — a
#        smaller number than the normalized coefficient, but directly comparable
#        to other fuel load data if SCF is changed to accept physical units.
#
# Toggle to FALSE to explore whether un-normalized FF changes the model fit
# or coefficient sign, without re-downloading the raster.
NORMALIZE_FINE_FUELS <- FALSE

# ---- Spatial coherence filter ------------------------------------------------
# pct_proximal is the fraction of newly-burned cells (y=1) that were adjacent to
# the frontier ring of the previous perimeter.  100% = contiguous day-to-day
# growth; 0% = the two perimeters are spatially disconnected (likely a GeoMAC
# mapping artefact: different survey methods on consecutive days).
#
# MIN_PCT_PROXIMAL controls which pairs enter the GLM:
#   0   = use ALL pairs regardless of spatial coherence (default, matches Robbins)
#   5   = drop pairs where <5% of new area overlaps the frontier ring
#   25  = moderate filter; keeps only roughly contiguous growth pairs
#   50  = strict; only clearly contiguous pairs
#   100 = only pairs where all new area was adjacent to the old frontier
#
# Raising this threshold REDUCES sample size but IMPROVES training data quality.
# The trade-off: with GeoMAC data, strict filtering often leaves too few pairs
# for a stable GLM.  Inspect the pair_diag table (Section 12) before tightening.
MIN_PCT_PROXIMAL <- 0L   # set >0 to filter spatially incoherent pairs

# ---- Seasonal and weather filters -------------------------------------------
# CAL_MONTHS: restrict pairs to specific calendar months.
#   NULL = use all months (default)
#   6:9  = June–September only (peak fire season; removes fall/winter wind fires)
#   This is the most important filter when calibrating in maritime regions
#   (e.g. Coast Range) where fall east-wind fires spread rapidly at low FWI and
#   would otherwise produce wrong-sign coefficients.
CAL_MONTHS <- 6:9    # NULL = all months; 6:9 = June–Sep only

# FWI_MIN: exclude pairs with daily FWI below this threshold.
#   NULL or 0 = no minimum (default)
#   5         = exclude near-zero FWI days (fall/winter smoldering)
#   10        = stricter; fire-season only in most western regions
#   Works alongside CAL_MONTHS — either or both can be set.
FWI_MIN <- 5         # NULL or 0 = no minimum


# =============================================================================
# SECTION 2 — Load helper functions from app modules
# =============================================================================

# Detect app root automatically when running in RStudio; otherwise set manually.
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  APP_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  APP_DIR <- "/Users/jlamping/Documents/Anthropic_Projects/scf_calibration_app"
}

source(file.path(APP_DIR, "modules/era_helpers.R"))
source(file.path(APP_DIR, "modules/ignition_helpers.R"))   # pick_first_col(), parse_any_date()
source(file.path(APP_DIR, "modules/spread_helpers.R"))

message("App modules loaded from: ", APP_DIR)


# =============================================================================
# SECTION 3 — Load spatial inputs
# =============================================================================

message("\n===== SECTION 3: Load spatial inputs =====")

template_r  <- rast(TEMPLATE_TIF)
working_crs <- crs(template_r)

message(sprintf("Template raster : %d rows x %d cols | res %.1f m",
                nrow(template_r), ncol(template_r), res(template_r)[1]))
message(sprintf("  CRS           : %s", crs(template_r, describe = TRUE)$name))
message(sprintf("  Extent        : %s", as.character(ext(template_r))))

# --- Load calibration boundary (shapefile OR tigris state download) -----------
if (!is.null(CAL_VECT_SHP)) {

  # Option A: user-supplied file
  cal_vect <- vect(CAL_VECT_SHP)
  message(sprintf("Cal boundary (file): %s", basename(CAL_VECT_SHP)))

} else if (!is.null(CAL_VECT_STATES) && length(CAL_VECT_STATES) > 0) {

  # Option B: tigris Census state boundaries
  message(sprintf("Fetching state boundaries via tigris: %s ...",
                  paste(CAL_VECT_STATES, collapse = ", ")))
  states_sf <- tigris::states(cb = TRUE, year = 2020, progress_bar = FALSE)
  states_sf <- states_sf[states_sf$NAME %in% CAL_VECT_STATES, ]

  missing_st <- setdiff(CAL_VECT_STATES, states_sf$NAME)
  if (length(missing_st) > 0)
    warning("State(s) not found in tigris: ", paste(missing_st, collapse = ", "),
            "\nCheck spelling — tigris uses full names, e.g. 'Oregon' not 'OR'.")

  # Dissolve to a single polygon so cal_vect is one feature
  cal_vect <- terra::vect(sf::st_union(states_sf))
  message(sprintf("  -> %d state(s) dissolved to 1 polygon.", nrow(states_sf)))

} else {
  stop("Set either CAL_VECT_SHP (file path) or CAL_VECT_STATES (character vector) in Section 1.")
}

cal_vect <- project(cal_vect, working_crs)
message(sprintf("Cal boundary    : %d feature(s) | extent %s",
                nrow(cal_vect), as.character(ext(cal_vect))))

# Sanity: template should overlap calibration boundary
if (!relate(ext(template_r), ext(cal_vect), "intersects")[1, 1]) {
  warning("Template extent does NOT intersect calibration boundary — check CRS/paths!")
} else {
  message("  -> Template and cal boundary extents overlap. OK.")
}


# =============================================================================
# SECTION 4 — Extract ERA-Interim FWI (same as app's "Load FWI Files" button)
# =============================================================================

message("\n===== SECTION 4: Extract ERA-Interim FWI from NetCDF =====")
message("Calling extract_era_timeseries() on: ", ERA_FWI_NC)

# The app uses cal_vect (the large regional boundary) to spatially average FWI,
# NOT the smaller template_r extent.  This matches observeEvent(input$load_era_fwi).
fwi_daily <- extract_era_timeseries(
  nc_path   = ERA_FWI_NC,
  park_vect = cal_vect,
  cal_start = CAL_START,
  cal_end   = CAL_END
)

message(sprintf("FWI extracted   : %d days | %s to %s",
                nrow(fwi_daily), min(fwi_daily$date), max(fwi_daily$date)))
message(sprintf("  min=%.1f  mean=%.1f  max=%.1f  NAs=%d",
                min(fwi_daily$value,  na.rm = TRUE),
                mean(fwi_daily$value, na.rm = TRUE),
                max(fwi_daily$value,  na.rm = TRUE),
                sum(is.na(fwi_daily$value))))

# Check: FWI should be zero in winter and peak in summer
fwi_daily %>%
  mutate(month = month(date, label = TRUE)) %>%
  group_by(month) %>%
  summarise(mean_fwi = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(month, mean_fwi, group = 1)) +
  geom_col(fill = "steelblue") +
  labs(title = "ERA-Interim FWI — monthly mean over calibration boundary",
       x = NULL, y = "Mean FWI") +
  theme_bw()


# =============================================================================
# SECTION 5 — Extract wind (GRIDMET or ERA5, same as app)
# =============================================================================

message("\n===== SECTION 5: Extract wind =====")

if (WIND_SOURCE == "gridmet") {
  # ---- OPTION A: GRIDMET remote fetch ----------------------------------------
  # Mirrors observeEvent(input$fetch_gridmet_wind) in app.R.
  # Requires internet access and the climateR package.
  if (!requireNamespace("climateR", quietly = TRUE))
    stop("Install climateR: install.packages('climateR')")

  message("Fetching GRIDMET wind (", year(CAL_START), "-", year(CAL_END), ") ...")
  message("  This can take several minutes for a long date range.")

  wind_daily <- extract_gridmet_wind(
    cal_vect  = cal_vect,
    cal_start = CAL_START,
    cal_end   = CAL_END,
    progress_fn = function(i, n, yr)
      message(sprintf("  GRIDMET: year %d of %d (%d)", i, n, yr))
  )

} else {
  # ---- OPTION B: Local ERA5 u10/v10 NetCDF files -----------------------------
  # Mirrors extract_era5_wind() in era_helpers.R.
  if (is.null(ERA5_U10_NC))
    stop("WIND_SOURCE='era5' but ERA5_U10_NC is NULL. Set the path or switch to 'gridmet'.")

  message("Extracting ERA5 wind from local NetCDF: ", ERA5_U10_NC)

  wind_daily <- extract_era5_wind(
    u10_path  = ERA5_U10_NC,
    v10_path  = ERA5_V10_NC,
    park_vect = cal_vect,
    cal_start = CAL_START,
    cal_end   = CAL_END
  )
}

message(sprintf("Wind extracted  : %d days | %s to %s",
                nrow(wind_daily), min(wind_daily$date), max(wind_daily$date)))
message(sprintf("  WindSpeed_kmh : min=%.1f  mean=%.1f  max=%.1f  NAs=%d",
                min(wind_daily$WindSpeed_kmh,  na.rm = TRUE),
                mean(wind_daily$WindSpeed_kmh, na.rm = TRUE),
                max(wind_daily$WindSpeed_kmh,  na.rm = TRUE),
                sum(is.na(wind_daily$WindSpeed_kmh))))

# SANITY: wind should be in km/h — typical range 0–80 km/h over western US.
# If mean < 5, data may still be in m/s (multiply by 3.6).
# If mean > 150, something is wrong.
if (mean(wind_daily$WindSpeed_kmh, na.rm = TRUE) < 5)
  warning("Mean WindSpeed_kmh < 5 — may still be in m/s. Check extract function.")
if (mean(wind_daily$WindSpeed_kmh, na.rm = TRUE) > 150)
  warning("Mean WindSpeed_kmh > 150 — unexpectedly high. Check wind units.")

# Plot FWI and wind together
bind_rows(
  fwi_daily  %>% transmute(date, variable = "FWI (ERA-Interim)",     value),
  wind_daily %>% transmute(date, variable = "Wind speed (km/h)", value = WindSpeed_kmh)
) %>%
  ggplot(aes(date, value)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  labs(title = "Climate time series (calibration boundary spatial mean)",
       x = NULL, y = NULL) +
  theme_bw()


# =============================================================================
# SECTION 6 — Fetch terrain for Effective Wind Speed (optional)
# Same as app's "Fetch DEM & Generate Terrain Maps" button in the Spread tab.
# =============================================================================

message("\n===== SECTION 6: Terrain (slope / aspect for EWS) =====")

if (USE_TERRAIN) {
  if (!requireNamespace("elevatr", quietly = TRUE))
    stop("Install elevatr: install.packages('elevatr')")

  message(sprintf("Fetching DEM via elevatr (zoom=%d) ...", TERRAIN_ZOOM))
  terrain_data <- fetch_terrain(
    cal_vect    = cal_vect,
    working_crs = working_crs,
    z           = TERRAIN_ZOOM
  )

  slope_r  <- terrain_data$slope
  aspect_r <- terrain_data$aspect

  message(sprintf("Slope  range: %.1f – %.1f degrees",
                  minmax(slope_r)[1], minmax(slope_r)[2]))
  message(sprintf("Aspect range: %.1f – %.1f degrees",
                  minmax(aspect_r)[1], minmax(aspect_r)[2]))

  # Warn about Cb sensitivity on steep terrain
  mean_slope <- mean(as.numeric(values(slope_r, mat = FALSE)), na.rm = TRUE)
  message(sprintf("Mean slope across cal boundary: %.1f degrees", mean_slope))
  if (mean_slope > 20 && COMBUSTION_BUOYANCY >= 25) {
    message(sprintf(paste0(
      "⚠  HIGH SLOPE (%.1f deg) + HIGH Cb (%.0f): slope contribution to EWS ",
      "is QUADRATIC in Cb. On steep terrain EWS will dominate, potentially ",
      "making FWI irrelevant and producing a NEGATIVE B1 in the GLM. ",
      "Consider Cb=10 (app default)."
    ), mean_slope, COMBUSTION_BUOYANCY))
  }

} else {
  slope_r  <- NULL
  aspect_r <- NULL
  message("Terrain skipped (USE_TERRAIN = FALSE) — raw wind speed used as EWS.")
}


# =============================================================================
# SECTION 7 — Build spread rasterization grid + load fine fuels
#
# The app builds a spread_template from the calibration boundary extent at a
# user-specified resolution, then loads fine fuels (LANDFIRE FBFM40 or local
# TIF) resampled onto that same grid.  All perimeter rasterization and fine
# fuel extraction in the GLM loop uses spread_template, NOT template_r.
#
# Mirrors app.R lines ~3900-3950 inside observeEvent(input$run_spread_fit).
# =============================================================================

message("\n===== SECTION 7: Build spread grid + load fine fuels =====")

# ---- Spread rasterization grid ----------------------------------------------
spread_res_m <- if (!is.null(SPREAD_RES_M) && is.finite(SPREAD_RES_M)) {
  max(as.numeric(SPREAD_RES_M), 10)
} else {
  res(template_r)[1]   # default: match LANDIS template resolution
}

spread_template <- rast(
  ext(cal_vect),
  resolution = spread_res_m,
  crs        = working_crs
)

spread_mask <- rasterize(cal_vect, spread_template, field = 1, background = NA)
spread_mask[!is.na(spread_mask)] <- 1L

message(sprintf("Spread grid     : %d x %d cells at %.0f m | extent %s",
                nrow(spread_template), ncol(spread_template),
                spread_res_m, as.character(ext(spread_template))))
message(sprintf("  Active cells  : %d of %d (%.1f%% of bounding box)",
                sum(values(spread_mask) == 1L, na.rm = TRUE),
                ncell(spread_template),
                100 * mean(values(spread_mask) == 1L, na.rm = TRUE)))

# ---- Fine fuels -------------------------------------------------------------
# Matches app.R tryCatch block: fetch LANDFIRE FBFM40 remote, load local, or NULL.
message(sprintf("\nFine fuels source: %s", FINE_FUELS_SOURCE))

# ---- Units label and scale factor (must be defined BEFORE the tryCatch) -----
# ff_units referenced inside the fetch message — defining here avoids a
# forward-reference error that silently causes the tryCatch to return NULL.
ff_units <- if (!NORMALIZE_FINE_FUELS && FINE_FUELS_SOURCE == "landfire") {
  "tons/acre (1-hr + 10-hr dead, Scott & Burgan 2005)"
} else if (!NORMALIZE_FINE_FUELS && FINE_FUELS_SOURCE == "local") {
  "original units (local raster, not normalized)"
} else {
  "0-1 index"
}
ff_scale <- if (!NORMALIZE_FINE_FUELS && FINE_FUELS_SOURCE == "landfire") {
  9.0   # GR9 fuel model maximum (Scott & Burgan 2005)
} else {
  1.0
}
message(sprintf("Fine fuels units: %s  |  ff_scale=%.4f  (NORMALIZE=%s)",
                ff_units, ff_scale, NORMALIZE_FINE_FUELS))

fine_fuels_r <- tryCatch({

  if (FINE_FUELS_SOURCE == "landfire") {
    # fetch_landfire_fine_fuels() in spread_helpers.R:
    #   - Requests LANDFIRE LF2024 FBFM40 for the cal_vect bounding box via
    #     the USGS ImageServer REST API (no local download needed)
    #   - Reclassifies integer Scott & Burgan fuel model codes to a 0-1 fine
    #     fuel load index using the FBFM40_FINE_FUELS lookup table (Scott &
    #     Burgan 2005: 1-hr + 10-hr dead fuel load, normalized to 0-1)
    #   - Reprojects and resamples to spread_template
    message("  Fetching LANDFIRE FBFM40 via USGS ImageServer REST API ...")
    message(sprintf("  normalize=%s  -> output units: %s", NORMALIZE_FINE_FUELS, ff_units))
    message("  (requires internet; may take 30-120 s depending on bbox size)")
    fetch_landfire_fine_fuels(
      cal_vect        = cal_vect,
      working_crs     = working_crs,
      spread_template = spread_template,
      normalize       = NORMALIZE_FINE_FUELS,
      progress_fn     = function(msg) message("  ", msg)
    )

  } else if (FINE_FUELS_SOURCE == "local") {
    if (is.null(FINE_FUELS_TIF) || !file.exists(FINE_FUELS_TIF))
      stop("FINE_FUELS_SOURCE='local' but FINE_FUELS_TIF path is NULL or not found.")
    message("  Loading local fine fuels raster: ", FINE_FUELS_TIF)
    message(sprintf("  normalize=%s  -> output units: %s", NORMALIZE_FINE_FUELS, ff_units))
    load_local_fine_fuels(
      path            = FINE_FUELS_TIF,
      working_crs     = working_crs,
      spread_template = spread_template,
      normalize       = NORMALIZE_FINE_FUELS
    )

  } else {
    message("  Using placeholder FF=1.0 everywhere (B2_FineFuels will be unidentifiable).")
    NULL
  }

}, error = function(e) {
  warning("Fine fuels load failed — falling back to placeholder FF=1.0: ",
          conditionMessage(e))
  NULL
})

# ff_units and ff_scale are defined earlier, before the tryCatch block.

# ---- Fine fuels diagnostics -------------------------------------------------
if (!is.null(fine_fuels_r)) {
  ff_vals_all <- values(fine_fuels_r, mat = FALSE)
  ff_vals_all <- ff_vals_all[is.finite(ff_vals_all)]
  message(sprintf("Fine fuels range: [%.3f, %.3f]  mean=%.3f  NAs=%d",
                  min(ff_vals_all), max(ff_vals_all), mean(ff_vals_all),
                  sum(is.na(values(fine_fuels_r, mat = FALSE)))))
  message(sprintf("Fine fuels variance: %.6f  (need > 0.001 for B2 to be identifiable)",
                  var(ff_vals_all)))

  # Map — reprojects to WGS84 for display
  ff_disp <- project(fine_fuels_r, "EPSG:4326")
  plot(ff_disp,
       main = sprintf("Fine fuels (%s)", ff_units),
       col  = hcl.colors(50, "YlOrRd"))

  # Distribution
  hist(ff_vals_all, breaks = 40, col = "steelblue", border = "white",
       main = sprintf("Fine fuels distribution (%s)", ff_units),
       xlab = sprintf("Fine fuel load (%s)", ff_units))

  # ---- NA breakdown -----------------------------------------------------------
  # NAs come from two sources:
  #   (a) cells outside the calibration boundary polygon (ocean, outside mask)
  #   (b) non-burnable FBFM40 codes (91=urban, 92=snow/ice, 93=agriculture,
  #       98=open water, 99=bare ground) → reclassified to NA during fetch
  # Only (b) is meaningful for the model: these cells are true non-burnable and
  # should not be candidate spread cells.  Check how many fall inside the mask.
  n_total_cells <- ncell(fine_fuels_r)
  n_inside_mask <- sum(values(spread_mask) == 1L, na.rm = TRUE)
  n_na_total    <- sum(is.na(values(fine_fuels_r, mat = FALSE)))
  ff_inside     <- values(fine_fuels_r * spread_mask, mat = FALSE)
  n_na_inside   <- sum(is.na(ff_inside))
  n_nonburn_inside <- sum(ff_inside == 0, na.rm = TRUE)

  message(sprintf("\nNA breakdown:"))
  message(sprintf("  Total raster cells    : %d", n_total_cells))
  message(sprintf("  Inside cal boundary   : %d", n_inside_mask))
  message(sprintf("  NA total (all cells)  : %d  (%.1f%%)", n_na_total,
                  100 * n_na_total / n_total_cells))
  message(sprintf("  NA inside boundary    : %d  (%.1f%% of masked cells)",
                  n_na_inside, 100 * n_na_inside / max(1L, n_inside_mask)))
  message(sprintf("  Zero (non-burnable)   : %d  (%.1f%% of masked cells)",
                  n_nonburn_inside, 100 * n_nonburn_inside / max(1L, n_inside_mask)))

  if (n_na_inside > 0.3 * n_inside_mask)
    message(paste0(
      "⚠  >30% of cells inside the calibration boundary have NA fine fuels. ",
      "These will be imputed with the mean FF value during the GLM sample loop, ",
      "which reduces spatial signal. Consider checking the LANDFIRE fetch coverage."
    ))

  # ---- FF quantile breakdown --------------------------------------------------
  ff_q <- quantile(ff_vals_all[ff_vals_all > 0], c(0.10, 0.25, 0.50, 0.75, 0.90),
                   na.rm = TRUE)
  message(sprintf("\nFF quantiles (burnable cells only, units: %s):", ff_units))
  message(sprintf("  Q10=%.3f  Q25=%.3f  Median=%.3f  Q75=%.3f  Q90=%.3f",
                  ff_q[1], ff_q[2], ff_q[3], ff_q[4], ff_q[5]))

  # ---- Broad fuel category breakdown ------------------------------------------
  # Thresholds scale with ff_scale so they mean the same fuel types regardless
  # of whether values are normalized (0-1) or in tons/acre.
  #   non-burnable  0 t/ac  (0.00 normalized)  urban, water, bare ground
  #   grass/low     0-0.9   (0-0.10 norm)       GR1-GR3, sparse TU/TL
  #   shrub/mid     0.9-2.7 (0.10-0.30 norm)    SH, GS, most TL/TU
  #   timber/high   > 2.7   (> 0.30 norm)       GR6-GR9, dense SH, slash SB
  lo_grass  <- 1e-9  * ff_scale
  hi_grass  <- 0.10  * ff_scale
  hi_shrub  <- 0.30  * ff_scale

  ff_in <- ff_inside[!is.na(ff_inside)]
  pct <- function(lo, hi) 100 * mean(ff_in >= lo & ff_in < hi)
  message(sprintf("\nFF category breakdown (masked cells; units: %s):", ff_units))
  message(sprintf("  non-burnable (= 0)              : %5.1f%%", pct(0, 1e-9 * ff_scale)))
  message(sprintf("  grass/low   (0 – %.2f)          : %5.1f%%", hi_grass,  pct(lo_grass, hi_grass)))
  message(sprintf("  shrub/mid   (%.2f – %.2f)       : %5.1f%%", hi_grass, hi_shrub, pct(hi_grass, hi_shrub)))
  message(sprintf("  timber/high (> %.2f)            : %5.1f%%", hi_shrub,  pct(hi_shrub, ff_scale * 1.1)))

} else {
  message("Fine fuels: NULL — placeholder 1.0 will be used in the GLM loop.")
  message("  B2_FineFuels will be set to 0 (unidentifiable — no spatial variation).")
}

# ---- Build dilation kernel from DILATION_CELLS ------------------------------
# A disk kernel so the frontier ring is roughly circular (not square-cornered).
# For DILATION_CELLS=1 this reproduces the classic rook neighbourhood exactly.
dil_kernel <- local({
  r <- as.integer(DILATION_CELLS)
  sz <- 2L * r + 1L
  center <- r + 1L
  m <- matrix(0L, sz, sz)
  for (ri in seq_len(sz))
    for (ci in seq_len(sz))
      if (sqrt((ri - center)^2 + (ci - center)^2) <= r) m[ri, ci] <- 1L
  m
})

message(sprintf(
  "Dilation kernel : %dx%d disk (DILATION_CELLS=%d) | frontier ring ~%d cells wide",
  nrow(dil_kernel), ncol(dil_kernel), DILATION_CELLS, DILATION_CELLS
))


# =============================================================================
# SECTION 8 — Load fire perimeters (GeoMAC or FIRED)
# Controlled by PERIMETER_SOURCE in Section 1.
# =============================================================================

message(sprintf("\n===== SECTION 8: Load perimeters (source = '%s') =====",
                PERIMETER_SOURCE))

if (PERIMETER_SOURCE == "geomac") {

  geomac_v <- load_geomac_perimeters(
    gdb_path    = GEOMAC_GDB,
    cal_vect    = cal_vect,
    template_r  = template_r,
    working_crs = working_crs,
    cal_start   = CAL_START,
    cal_end     = CAL_END
  )

} else if (PERIMETER_SOURCE == "fired") {

  # load_fired_perimeters() returns a SpatVector with the same schema
  # (FIRE_ID, perimeter_date, AREA_HA, row_id) as load_geomac_perimeters().
  # Everything downstream (bind_climate_to_pairs, rasterize_perim, etc.) is
  # source-agnostic — no other sections need to change.
  geomac_v <- load_fired_perimeters(
    fired_gpkg  = FIRED_GPKG,
    cal_vect    = cal_vect,
    working_crs = working_crs,
    cal_start   = CAL_START,
    cal_end     = CAL_END,
    min_fire_ha = FIRED_MIN_FIRE_HA,
    layer       = FIRED_LAYER,
    n_cores     = N_CORES
  )

} else {
  stop("PERIMETER_SOURCE must be 'geomac' or 'fired'. Got: ", PERIMETER_SOURCE)
}

# ---- Common summary (same for both sources) ---------------------------------
gm_df <- as.data.frame(geomac_v)

message(sprintf("Perimeters loaded: %d for %d unique fires",
                nrow(gm_df), length(unique(gm_df$FIRE_ID))))
message(sprintf("Date range       : %s  to  %s",
                min(as.Date(gm_df$perimeter_date), na.rm = TRUE),
                max(as.Date(gm_df$perimeter_date), na.rm = TRUE)))
message(sprintf("Area (ha)        : min=%.1f  median=%.1f  max=%.1f",
                min(gm_df$AREA_HA, na.rm = TRUE),
                median(gm_df$AREA_HA, na.rm = TRUE),
                max(gm_df$AREA_HA, na.rm = TRUE)))

perim_per_fire <- gm_df %>%
  group_by(FIRE_ID) %>%
  summarise(n_perims = n(), .groups = "drop")

message(sprintf("Perimeters/fire  : min=%d  median=%.0f  max=%d",
                min(perim_per_fire$n_perims),
                median(perim_per_fire$n_perims),
                max(perim_per_fire$n_perims)))
message(sprintf("  Single-perimeter fires (no pairs): %d of %d",
                sum(perim_per_fire$n_perims == 1),
                nrow(perim_per_fire)))

# Plots
ggplot(gm_df, aes(AREA_HA)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::comma) +
  labs(title = sprintf("%s perimeter area distribution", toupper(PERIMETER_SOURCE)),
       x = "Area (ha, log scale)", y = "Count") +
  theme_bw()

gm_df %>%
  mutate(year = year(as.Date(perimeter_date))) %>%
  count(year) %>%
  ggplot(aes(year, n)) +
  geom_col(fill = "steelblue") +
  labs(title = "GeoMAC perimeters per year in calibration boundary",
       x = NULL, y = "Count") +
  theme_bw()

# ---- Fine fuels at fire locations --------------------------------------------
# Only runs if fine fuels were loaded.  Extracts FF values at each GeoMAC fire
# centroid and compares them to the background FF distribution.
# Key question: do fires occur in higher-FF areas than the landscape average?
# If fires cluster in low-FF cells (e.g. grass/shrub), B2 will be negative or
# near-zero because most spread-success cells are low-FF.
if (!is.null(fine_fuels_r)) {
  message("\n--- Fine fuels at fire locations ---")

  fire_centroids_v <- centroids(geomac_v)
  ff_at_fires <- terra::extract(fine_fuels_r, fire_centroids_v)[, 2]

  fire_ff_df <- tibble(
    FIRE_ID       = as.character(geomac_v$FIRE_ID),
    AREA_HA       = as.numeric(geomac_v$AREA_HA),
    perimeter_date = as.Date(geomac_v$perimeter_date),
    FF_centroid   = ff_at_fires
  ) %>%
    filter(is.finite(FF_centroid))

  message(sprintf("FF at fire centroids: n=%d  mean=%.3f  median=%.3f",
                  nrow(fire_ff_df),
                  mean(fire_ff_df$FF_centroid, na.rm = TRUE),
                  median(fire_ff_df$FF_centroid, na.rm = TRUE)))
  message(sprintf("Background FF mean  : %.3f  (all masked cells)",
                  mean(values(fine_fuels_r * spread_mask, mat = FALSE), na.rm = TRUE)))

  # Are fires in higher-FF areas than background?
  bg_ff  <- values(fine_fuels_r * spread_mask, mat = FALSE)
  bg_ff  <- bg_ff[is.finite(bg_ff)]
  t_test <- tryCatch(
    t.test(fire_ff_df$FF_centroid, mu = mean(bg_ff)),
    error = function(e) NULL
  )
  if (!is.null(t_test)) {
    message(sprintf("  t-test (fire FF vs background mean): t=%.2f  p=%.4f  [%s]",
                    t_test$statistic, t_test$p.value,
                    if (t_test$p.value < 0.05) {
                      if (mean(fire_ff_df$FF_centroid, na.rm=TRUE) > mean(bg_ff))
                        "fires in HIGHER FF than background — B2 could be positive"
                      else
                        "fires in LOWER FF than background — B2 may be near zero/negative"
                    } else "no significant difference from background"))
  }

  # Density: fire centroid FF vs background FF
  ggplot() +
    geom_density(data = tibble(FF = bg_ff[sample(length(bg_ff),
                               min(length(bg_ff), 100000))]),
                 aes(FF, fill = "Background landscape"), alpha = 0.4) +
    geom_density(data = fire_ff_df, aes(FF_centroid, fill = "Fire centroids"),
                 alpha = 0.4) +
    scale_fill_manual(values = c("Background landscape" = "steelblue",
                                 "Fire centroids"       = "tomato")) +
    labs(title = "Fine fuels: fire centroid locations vs background landscape",
         subtitle = "Rightward shift = fires burn in higher-FF areas (supports positive B2)",
         x = sprintf("Fine fuel load (%s)", ff_units), y = "Density", fill = NULL) +
    theme_bw()

  # FF vs fire size — do large fires occur in high-FF areas?
  ggplot(fire_ff_df %>% filter(AREA_HA > 0), aes(FF_centroid, AREA_HA)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Fine fuels vs fire size",
         subtitle = sprintf("r = %.3f",
                            cor(fire_ff_df$FF_centroid, log(fire_ff_df$AREA_HA + 1),
                                use = "complete.obs")),
         x = sprintf("FF at fire centroid (%s)", ff_units),
         y = "Fire area (ha, log scale)") +
    theme_bw()
}


# =============================================================================
# SECTION 9 — Build perimeter pairs + join climate
# Mirrors bind_climate_to_pairs() called in app server around line 3881.
# =============================================================================

message("\n===== SECTION 9: Build perimeter pairs + join climate =====")

pairs_result <- bind_climate_to_pairs(
  geomac_v            = geomac_v,
  fwi_daily           = fwi_daily,
  wind_daily          = wind_daily,
  slope_r             = slope_r,
  aspect_r            = aspect_r,
  max_gap_sp          = MAX_GAP_SP,
  max_gap_da          = MAX_GAP_DA,
  neg_tol_ha          = NEG_TOL_HA,
  combustion_buoyancy = COMBUSTION_BUOYANCY
)

all_pairs    <- pairs_result$pairs
spread_pairs <- pairs_result$spread_pairs   # pre-filtered: use_for_spread_prob == TRUE

message(sprintf("All pairs       : %d", nrow(all_pairs)))
message(sprintf("Spread pairs    : %d  (gap <= %d day, delta_ha > 0)",
                nrow(spread_pairs), MAX_GAP_SP))
message(sprintf("Daily-area pairs: %d  (gap <= %d days, delta_ha > 0)",
                sum(all_pairs$use_for_daily_area, na.rm = TRUE), MAX_GAP_DA))
message(sprintf("Dropped (shrink): %d  (delta_ha < -%.1f ha)",
                sum(all_pairs$drop_neg, na.rm = TRUE), NEG_TOL_HA))

# ---- Climate join quality ---------------------------------------------------
fwi_join_rate  <- mean(is.finite(all_pairs$FWI))
wind_join_rate <- mean(is.finite(all_pairs$WindSpeed_kmh))

message(sprintf("\nFWI join rate   : %.0f%% of pairs have FWI  (%d of %d)",
                100 * fwi_join_rate, sum(is.finite(all_pairs$FWI)), nrow(all_pairs)))
message(sprintf("Wind join rate  : %.0f%% of pairs have wind (%d of %d)",
                100 * wind_join_rate, sum(is.finite(all_pairs$WindSpeed_kmh)), nrow(all_pairs)))

if (fwi_join_rate < 0.5) {
  message(paste0(
    "⚠  Less than 50% of pairs matched to FWI. Likely causes:\n",
    "   - Perimeter dates fall outside FWI NetCDF date range\n",
    "   - Date class mismatch between GeoMAC and ERA NetCDF\n",
    "   - ERA FWI file does not cover the calibration date range\n",
    "   Sample pair dates : ", paste(head(all_pairs$date, 5), collapse = ", "), "\n",
    "   Sample FWI dates  : ", paste(head(fwi_daily$date, 5), collapse = ", ")
  ))
}

# Plot: gap distribution
ggplot(all_pairs, aes(gap_days)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "gray") +
  labs(title = "Days between consecutive perimeters (all pairs)",
       x = "Gap (days)", y = "Count") +
  theme_bw()

# Plot: delta area (non-shrinking pairs)
ggplot(all_pairs %>% filter(!drop_neg, delta_ha > 0), aes(delta_ha)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "grey") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Area increment (delta_ha) — non-shrinking pairs",
       x = "delta_ha (log scale, ha)", y = "Count") +
  theme_bw()


# =============================================================================
# SECTION 10 — Effective Wind Speed diagnostics
# =============================================================================

message("\n===== SECTION 10: Effective wind speed diagnostics =====")

wind_col <- if ("EffectiveWind" %in% names(all_pairs) &&
                any(is.finite(all_pairs$EffectiveWind))) "EffectiveWind" else "WindSpeed_kmh"

message(sprintf("Wind predictor  : %s", wind_col))

if (wind_col == "EffectiveWind") {

  ews_df <- all_pairs %>%
    filter(is.finite(EffectiveWind), is.finite(WindSpeed_kmh))

  message(sprintf("Raw wind (km/h) : mean=%.1f  sd=%.1f  range [%.1f, %.1f]",
                  mean(ews_df$WindSpeed_kmh), sd(ews_df$WindSpeed_kmh),
                  min(ews_df$WindSpeed_kmh),  max(ews_df$WindSpeed_kmh)))
  message(sprintf("EWS  (km/h,Cb=%.0f): mean=%.1f  sd=%.1f  range [%.1f, %.1f]",
                  COMBUSTION_BUOYANCY,
                  mean(ews_df$EffectiveWind), sd(ews_df$EffectiveWind),
                  min(ews_df$EffectiveWind),  max(ews_df$EffectiveWind)))
  message(sprintf("Slope add (km/h): mean=%.1f  sd=%.1f",
                  mean(ews_df$EffectiveWind - ews_df$WindSpeed_kmh),
                  sd(ews_df$EffectiveWind   - ews_df$WindSpeed_kmh)))

  # Cb sensitivity — show what EWS would look like at Cb 10/25/50
  message("\nCb sensitivity check for the pair with largest slope contribution:")
  ex <- ews_df %>%
    filter(!is.na(Slope_deg), !is.na(Aspect_deg)) %>%
    slice_max(abs(EffectiveWind - WindSpeed_kmh), n = 1)
  if (nrow(ex) > 0) {
    for (cb_try in c(10, 25, 50)) {
      ews_try <- compute_effective_wind(
        wind_speed = ex$WindSpeed_kmh, wind_dir = ex$WindDir_deg,
        slope_deg  = ex$Slope_deg,    aspect_deg = ex$Aspect_deg,
        combustion_buoyancy = cb_try
      )
      message(sprintf("  Cb=%2d: wind=%.1f km/h  slope=%.1f°  -> EWS=%.1f km/h",
                      cb_try, ex$WindSpeed_kmh, ex$Slope_deg, ews_try))
    }
  }

  # FWI / EWS correlation
  r_fwi_ews <- cor(ews_df$FWI, ews_df$EffectiveWind, use = "complete.obs")
  message(sprintf("\nCor(FWI, EWS)      : %.3f", r_fwi_ews))
  message(sprintf("Cor(FWI, RawWind)  : %.3f",
                  cor(ews_df$FWI, ews_df$WindSpeed_kmh, use = "complete.obs")))
  if (abs(r_fwi_ews) > 0.8)
    message("⚠  High FWI/EWS correlation — multicollinearity will inflate SEs.")

  ggplot(ews_df, aes(WindSpeed_kmh, EffectiveWind)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = sprintf("Raw wind vs Effective Wind Speed (Cb=%.0f)", COMBUSTION_BUOYANCY),
         x = "Raw wind (km/h)", y = "EWS (km/h)",
         caption = "Points on the diagonal = no slope contribution") +
    theme_bw()

  ggplot(ews_df, aes(FWI, EffectiveWind)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(title = "FWI vs EWS — high correlation causes multicollinearity",
         subtitle = sprintf("Pearson r = %.3f", r_fwi_ews),
         x = "FWI", y = "EWS (km/h)") +
    theme_bw()

} else {
  message("No EWS (terrain not loaded). Using raw wind speed (km/h).")
  message(sprintf("Wind (km/h): mean=%.1f  sd=%.1f  range [%.1f, %.1f]",
                  mean(all_pairs$WindSpeed_kmh, na.rm = TRUE),
                  sd(all_pairs$WindSpeed_kmh,   na.rm = TRUE),
                  min(all_pairs$WindSpeed_kmh,   na.rm = TRUE),
                  max(all_pairs$WindSpeed_kmh,   na.rm = TRUE)))
}


# =============================================================================
# SECTION 11 — Single-pair rasterization example
# Visualises the binary spread/no-spread classification for one pair.
# =============================================================================

message("\n===== SECTION 11: Single-pair rasterization example =====")
message(sprintf("  Using dilation kernel: %dx%d disk (DILATION_CELLS=%d)",
                nrow(dil_kernel), ncol(dil_kernel), DILATION_CELLS))
message("  Perimeters are rasterized WITHOUT the cal_vect mask so fires that")
message("  extend beyond the boundary are NOT clipped (avoids false no-successes).")

example_pairs <- spread_pairs %>%
  filter(is.finite(FWI)) %>%
  arrange(desc(delta_ha))

if (nrow(example_pairs) == 0) {
  message("⚠  No eligible spread pairs — cannot show single-pair example.")
} else {
  ex_row <- example_pairs[1, ]
  message(sprintf("Example pair: fire=%-20s  date=%s  delta_ha=%.1f ha  FWI=%.1f",
                  ex_row$FIRE_ID, ex_row$date, ex_row$delta_ha, ex_row$FWI))

  # NULL mask: rasterize the full fire polygon without clipping to cal_vect.
  # Fires were selected because they intersect cal_vect, so they are spatially
  # relevant — but clipping their geometry would zero out cells on the frontier
  # ring that actually burned, producing spurious n_succ=0 results.
  r_t  <- rasterize_perim(geomac_v, ex_row$row_id,      spread_template, NULL)
  r_t1 <- rasterize_perim(geomac_v, ex_row$row_id_next, spread_template, NULL)

  if (is.null(r_t) || is.null(r_t1)) {
    message("⚠  rasterize_perim() returned NULL.")
    message(sprintf("  row_id=%s  row_id_next=%s", ex_row$row_id, ex_row$row_id_next))
    message(sprintf("  Available row_ids (first 20): %s",
                    paste(head(sort(as.integer(geomac_v$row_id)), 20), collapse = " ")))
    message("  Likely cause: CRS mismatch between geomac_v and spread_template.")
    message(sprintf("  geomac_v CRS     : %s", crs(geomac_v, describe = TRUE)$name))
    message(sprintf("  spread_template CRS: %s", crs(spread_template, describe = TRUE)$name))
  } else {
    cell_ha <- prod(res(spread_template)) / 10000
    n_t  <- sum(values(r_t)  == 1, na.rm = TRUE)
    n_t1 <- sum(values(r_t1) == 1, na.rm = TRUE)
    message(sprintf("r_t  cells=%d (%.1f ha)  |  r_t1 cells=%d (%.1f ha)",
                    n_t, n_t * cell_ha, n_t1, n_t1 * cell_ha))

    # SUCCESS: cells newly burned between t and t+1 (Robbins 2021 approach)
    new_burned <- (r_t1 == 1) & (r_t == 0)
    n_succ     <- sum(values(new_burned) == 1, na.rm = TRUE)

    # FAILURE: frontier ring of r_t that did NOT burn
    dil    <- focal(r_t, w = dil_kernel, fun = "max", na.policy = "omit", fillvalue = 0)
    cand   <- (dil == 1) & (r_t == 0)
    fail   <- cand & (r_t1 == 0)
    n_fail <- sum(values(fail) == 1, na.rm = TRUE)
    n_cand <- sum(values(cand) == 1, na.rm = TRUE)

    message(sprintf("New burned cells (success): %d  |  Frontier unburned (failure): %d  |  Frontier total: %d",
                    n_succ, n_fail, n_cand))

    if (n_succ == 0 && n_fail == 0) {
      message(paste0(
        "⚠  No new burned cells and empty frontier. Most likely causes:\n",
        "   1. Perimeter projects outside the spread_template extent\n",
        "      Spread grid extent : ", as.character(ext(spread_template)), "\n",
        "      GeoMAC extent      : ", as.character(ext(geomac_v)), "\n",
        "   2. Fire did not grow (delta_ha ≤ 0 or within rounding at this resolution)\n",
        "      Spread cell        : ", round(cell_ha, 2), " ha  (",
        round(spread_res_m), " m resolution)\n",
        "      Fire delta_ha      : ", round(ex_row$delta_ha, 1), " ha\n",
        "   Tip: reduce SPREAD_RES_M to better resolve small perimeters."
      ))
    } else {
      # Build classification raster:
      #   1 = already burned at t (grey)
      #   2 = frontier ring of r_t, didn't burn (failure, tomato)
      #   3 = newly burned at t+1, outside frontier (distant spread / spotting, gold)
      #   4 = newly burned at t+1, inside frontier (proximal spread, royalblue)
      cls_r <- r_t * 0L
      cls_r[r_t         == 1] <- 1L
      cls_r[cand        == 1] <- 2L   # frontier unburned (failure)
      cls_r[(new_burned == 1) & (cand == 0)] <- 3L   # new area outside frontier
      cls_r[(new_burned == 1) & (cand == 1)] <- 4L   # new area inside frontier
      cls_r[cls_r == 0]   <- NA

      # Zoom to fire bounding box with margin
      bbox_t  <- ext(geomac_v[geomac_v$row_id == ex_row$row_id, ])
      bbox_t1 <- ext(geomac_v[geomac_v$row_id == ex_row$row_id_next, ])
      # Union both bboxes so both t and t+1 perimeters are visible
      union_bbox <- ext(
        min(xmin(bbox_t), xmin(bbox_t1)), max(xmax(bbox_t), xmax(bbox_t1)),
        min(ymin(bbox_t), ymin(bbox_t1)), max(ymax(bbox_t), ymax(bbox_t1))
      )
      margin   <- max(diff(c(xmin(union_bbox), xmax(union_bbox))),
                      diff(c(ymin(union_bbox), ymax(union_bbox)))) * 0.15
      disp_ext <- ext(xmin(union_bbox) - margin, xmax(union_bbox) + margin,
                      ymin(union_bbox) - margin, ymax(union_bbox) + margin)
      cls_crop <- crop(cls_r, disp_ext)

      plot(cls_crop,
           col    = c("grey80", "tomato", "gold2", "royalblue"),
           main   = sprintf("Pair: %s  %s\nn_new=%d  n_frontier_fail=%d",
                            ex_row$FIRE_ID, ex_row$date, n_succ, n_fail),
           legend = FALSE)
      legend("topright",
             legend = c("Already burned (t)", "Frontier unburned (failure)",
                        "New area — distant spread", "New area — proximal spread"),
             fill   = c("grey80", "tomato", "gold2", "royalblue"),
             cex    = 0.7, bg = "white")

      # Spatial offset diagnostic: if most new cells are "distant" (gold, not blue),
      # it confirms the perimeters are spatially inconsistent.
      n_proximal <- sum(values(cls_crop) == 4L, na.rm = TRUE)
      n_distant  <- sum(values(cls_crop) == 3L, na.rm = TRUE)
      if (n_proximal + n_distant > 0) {
        pct_proximal <- 100 * n_proximal / (n_proximal + n_distant)
        message(sprintf(
          "New area breakdown: %.0f%% proximal (inside frontier, blue) | %.0f%% distant (outside frontier, gold)",
          pct_proximal, 100 - pct_proximal
        ))
        if (pct_proximal < 20)
          message(paste0(
            "⚠  >80% of new burned area is OUTSIDE the frontier ring.\n",
            "   This confirms spatial inconsistency between the t and t+1 GeoMAC perimeters.\n",
            "   The fire grew, but the t+1 polygon is offset or disconnected from t.\n",
            "   These samples are included in the GLM (success cells are real),\n",
            "   but FWI/EWS conditions are assigned from the fire centroid, not cell-level."
          ))
      }
    }
  }
}


# =============================================================================
# SECTION 12 — Build full sample dataset (all eligible pairs)
# Mirrors the loop inside fit_spread_probability() in spread_helpers.R.
# =============================================================================

message("\n===== SECTION 12: Build full sample dataset =====")

eligible <- spread_pairs %>%
  filter(use_for_spread_prob, is.finite(.data[[wind_col]]), is.finite(FWI))

# ---- Seasonal filter (CAL_MONTHS) -------------------------------------------
if (!is.null(CAL_MONTHS) && length(CAL_MONTHS) > 0) {
  n_before <- nrow(eligible)
  eligible <- eligible %>% filter(lubridate::month(date) %in% CAL_MONTHS)
  message(sprintf("CAL_MONTHS filter [%s]: kept %d of %d pairs.",
                  paste(CAL_MONTHS, collapse=","), nrow(eligible), n_before))
}

# ---- FWI minimum filter (FWI_MIN) -------------------------------------------
if (!is.null(FWI_MIN) && FWI_MIN > 0) {
  n_before <- nrow(eligible)
  eligible <- eligible %>% filter(FWI >= FWI_MIN)
  message(sprintf("FWI_MIN filter (>= %.1f): kept %d of %d pairs.",
                  FWI_MIN, nrow(eligible), n_before))
}

message(sprintf("Eligible pairs  : %d", nrow(eligible)))

if (nrow(eligible) == 0)
  stop("No eligible pairs. Check FWI join (Section 9) and wind data (Section 5).")

# fine_fuels_r was already loaded and inspected in Section 7 — reuse it here.
# (If it is NULL, the loop uses FF=1.0 placeholder and B2 will be set to 0.)
message(sprintf("Fine fuels      : %s",
                if (is.null(fine_fuels_r)) "NULL (placeholder FF=1.0)"
                else sprintf("loaded, range [%.3f, %.3f]",
                             min(values(fine_fuels_r), na.rm = TRUE),
                             max(values(fine_fuels_r), na.rm = TRUE))))

# ---- Sample loop (identical to fit_spread_probability) ----------------------
# NOTE: rasterize_perim() is called with mask_r = NULL (no cal_vect clipping).
# Reason: fires were selected because they intersect cal_vect, so they are
# spatially relevant.  Clipping r_t and r_t1 to the cal_vect polygon would zero
# out frontier cells that extend past the administrative boundary, producing
# n_succ=0 for fires that genuinely spread outside that edge — a false signal.
# Fine fuels are extracted from fine_fuels_r which is already bounded by the
# spread_template extent; cells outside cal_vect will have NA FF and are
# imputed with the pair's mean FF value below.
message(sprintf(
  "\nSection 12 loop: mask=NULL (no boundary clip)  |  dilation kernel=%dx%d (DILATION_CELLS=%d)",
  nrow(dil_kernel), ncol(dil_kernel), DILATION_CELLS
))

# ---- Resolve worker count (reuse N_CORES from Section 1) ---------------------
.n_workers <- {
  avail <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  n     <- if (is.null(N_CORES)) avail else max(1L, min(as.integer(N_CORES), avail))
  min(n, nrow(eligible))   # never more workers than pairs
}
message(sprintf("Section 12 parallel: %d pairs across %d workers", nrow(eligible), .n_workers))

# ---- Prepare terra objects for worker access ---------------------------------
# multisession workers are fresh R sessions — terra SpatRaster/SpatVector objects
# cannot be serialized across process boundaries (XPtr environments cause ls()
# errors inside future). Strategy:
#   - perimeters (SpatVector): write to temp GPKG, workers read with vect()
#   - spread_template: workers only need its geometry for rasterize(), not values.
#                      Pass a plain-list of metadata; reconstruct inside worker.
#                      (Writing a lazy-backed raster to disk can silently fail.)
#   - fine_fuels_r: actual cell values are needed, so write to temp TIF.
eligible$.orig_idx <- seq_len(nrow(eligible))   # preserve pair index after chunking

.tmp_geomac <- tempfile(fileext = ".gpkg")
.tmp_ff     <- if (!is.null(fine_fuels_r)) tempfile(fileext = ".tif") else NULL

message("  Writing perimeter vector and fine fuels raster to temp files...")
terra::writeVector(geomac_v, .tmp_geomac, filetype = "GPKG", overwrite = TRUE)
if (!is.null(.tmp_ff)) terra::writeRaster(fine_fuels_r, .tmp_ff, overwrite = TRUE)

# Template metadata — plain R list, no terra objects, no file needed
.template_meta <- list(
  nrows = terra::nrow(spread_template),
  ncols = terra::ncol(spread_template),
  xmin  = terra::xmin(spread_template),
  xmax  = terra::xmax(spread_template),
  ymin  = terra::ymin(spread_template),
  ymax  = terra::ymax(spread_template),
  crs   = terra::crs(spread_template)
)

# ---- Split eligible rows into one chunk per worker ---------------------------
# include.lowest=TRUE ensures row 1 is not dropped (cut() default is open-left).
.chunk_ids <- cut(seq_len(nrow(eligible)), breaks = .n_workers,
                  labels = FALSE, include.lowest = TRUE)
.chunks    <- split(eligible, .chunk_ids)

# ---- Worker function ---------------------------------------------------------
# Receives file paths as plain strings and loads terra objects fresh in each
# worker session. p is a progressr::progressor — captured as a closure variable
# inside the lambda in the dispatch block (NOT passed as a named arg, which
# would cause furrr to call ls() on its C-level environment).
.process_pairs_chunk <- function(chunk_rows,
                                  geomac_path, template_meta, ff_path,
                                  dil_kernel, wind_col,
                                  MAX_SAMPLES, FAILURE_RATIO,
                                  p = NULL) {
  # Load / reconstruct terra objects fresh in this worker session
  gv   <- terra::vect(geomac_path)
  tmpl <- terra::rast(
    nrows = template_meta$nrows, ncols = template_meta$ncols,
    xmin  = template_meta$xmin,  xmax  = template_meta$xmax,
    ymin  = template_meta$ymin,  ymax  = template_meta$ymax,
    crs   = template_meta$crs
  )
  ff   <- if (!is.null(ff_path)) terra::rast(ff_path) else NULL

  # rasterize_perim is small enough to define inline (avoids closure detection)
  .rasterize_perim <- function(v, row_id_val, template, mask_r = NULL) {
    vv <- v[v$row_id == row_id_val, ]
    if (nrow(vv) == 0) return(NULL)
    r <- terra::rasterize(vv, template, field = 1, touches = TRUE)
    r[!is.na(r)] <- 1
    r[is.na(r)]  <- 0
    if (!is.null(mask_r)) { r <- r * mask_r; r[is.na(r)] <- 0 }
    r
  }

  n            <- nrow(chunk_rows)
  pair_stats_l <- vector("list", n)
  samples_l    <- vector("list", n)

  for (j in seq_len(n)) {
    row <- chunk_rows[j, ]
    i   <- row$.orig_idx
    fwi <- row$FWI
    ews <- row[[wind_col]]

    r_t  <- .rasterize_perim(gv, row$row_id,      tmpl, NULL)
    r_t1 <- .rasterize_perim(gv, row$row_id_next, tmpl, NULL)

    if (is.null(r_t) || is.null(r_t1)) {
      pair_stats_l[[j]] <- dplyr::tibble(
        pair_idx = i, status = "null_raster",
        n_new_burned = 0L, n_fail = 0L, pct_proximal = NA_real_,
        FWI = fwi, EWS = ews, delta_ha = row$delta_ha)
      if (!is.null(p)) p()
      next
    }

    new_burned <- (r_t1 == 1) & (r_t == 0)
    succ_cells <- which(terra::values(new_burned) == 1)

    dil        <- terra::focal(r_t, w = dil_kernel, fun = "max",
                               na.policy = "omit", fillvalue = 0)
    cand       <- (dil == 1) & (r_t == 0)
    fail_cells <- which(terra::values(cand & (r_t1 == 0)) == 1)

    n_succ <- length(succ_cells)
    n_fail <- length(fail_cells)
    n_new  <- n_succ

    n_proximal <- sum(terra::values((new_burned == 1) & (cand == 1)) == 1, na.rm = TRUE)
    pct_prox   <- if (n_new > 0) round(100 * n_proximal / n_new, 1) else NA_real_

    pair_stats_l[[j]] <- dplyr::tibble(
      pair_idx     = i, FIRE_ID = row$FIRE_ID, date = row$date,
      FWI = fwi, EWS = ews, delta_ha = row$delta_ha,
      n_new_burned = n_new, n_fail = n_fail, pct_proximal = pct_prox,
      status       = if (n_succ == 0) "no_success" else "ok"
    )

    if (n_succ == 0) { if (!is.null(p)) p(); next }

    succ_cap   <- as.integer(MAX_SAMPLES / (1L + FAILURE_RATIO))
    succ_cells <- if (succ_cap < n_succ) sample(succ_cells, succ_cap) else succ_cells
    n_succ     <- length(succ_cells)

    n_fail_keep <- min(n_fail, n_succ * FAILURE_RATIO)
    fail_keep   <- if (n_fail_keep > 0) sample(fail_cells, n_fail_keep) else integer(0)
    n_f         <- length(fail_keep)
    if (n_succ + n_f == 0) { if (!is.null(p)) p(); next }

    ff_vals <- if (!is.null(ff)) {
      v <- as.numeric(terra::values(ff)[c(succ_cells, fail_keep)])
      v[!is.finite(v)] <- mean(v, na.rm = TRUE)
      v
    } else {
      rep(1.0, n_succ + n_f)
    }

    samples_l[[j]] <- dplyr::tibble(
      pair_idx      = rep(i,           n_succ + n_f),
      FIRE_ID       = rep(row$FIRE_ID, n_succ + n_f),
      date          = rep(row$date,    n_succ + n_f),
      y             = c(rep(1L, n_succ), rep(0L, n_f)),
      FWI           = rep(fwi, n_succ + n_f),
      EffectiveWind = rep(ews, n_succ + n_f),
      FineFuels     = ff_vals
    )

    if (!is.null(p)) p()
  }

  list(
    pair_stats = dplyr::bind_rows(pair_stats_l),
    samples    = dplyr::bind_rows(samples_l)
  )
}

# ---- Dispatch with progress bar ----------------------------------------------
if (.n_workers > 1L) {
  future::plan(future::multisession, workers = .n_workers)
} else {
  future::plan(future::sequential)
}

t0 <- proc.time()[["elapsed"]]
progressr::handlers("cli")
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(eligible))
  .chunk_results <- furrr::future_map(
    .chunks,
    function(chunk) .process_pairs_chunk(
      chunk_rows    = chunk,
      geomac_path   = .tmp_geomac,
      template_meta = .template_meta,
      ff_path       = .tmp_ff,
      dil_kernel    = dil_kernel,
      wind_col      = wind_col,
      MAX_SAMPLES   = MAX_SAMPLES,
      FAILURE_RATIO = FAILURE_RATIO,
      p             = p
    ),
    .options = furrr::furrr_options(seed = TRUE)
  )
})

future::plan(future::sequential)
unlink(c(.tmp_geomac, .tmp_ff))   # clean up temp files
message(sprintf("  Section 12 complete in %.0f s (~%.1f min).",
                proc.time()[["elapsed"]] - t0,
                (proc.time()[["elapsed"]] - t0) / 60))

# ---- Reassemble results in original pair order -------------------------------
pair_diag   <- dplyr::bind_rows(lapply(.chunk_results, `[[`, "pair_stats")) %>%
  dplyr::arrange(pair_idx)
all_samples <- dplyr::bind_rows(lapply(.chunk_results, `[[`, "samples")) %>%
  dplyr::arrange(pair_idx)

# ---- Pair-level diagnostics -------------------------------------------------
# Success cell definition (Robbins 2021):
#   n_new_burned = total cells newly burned between t and t+1 (used as y=1 in GLM)
#   n_fail       = frontier ring cells of r_t that did NOT burn (y=0 in GLM)
#   pct_proximal = % of new burned cells that were adjacent to r_t (100% = fire
#                  grew contiguously; low% = jump spread or GeoMAC spatial offset)
message("\nPair-level summary:")
message(sprintf("  null rasterization : %d", sum(pair_diag$status == "null_raster", na.rm = TRUE)))
message(sprintf("  no new burned cells: %d", sum(pair_diag$status == "no_success",  na.rm = TRUE)))
message(sprintf("  usable pairs       : %d", sum(pair_diag$status == "ok",          na.rm = TRUE)))

pd_ok <- pair_diag %>% filter(status == "ok")
if (nrow(pd_ok) > 0) {
  message(sprintf("  New burned cells (median/max): %d / %d",
                  as.integer(median(pd_ok$n_new_burned)), max(pd_ok$n_new_burned)))
  message(sprintf("  Frontier unburned (median/max): %d / %d",
                  as.integer(median(pd_ok$n_fail)), max(pd_ok$n_fail)))
  message(sprintf("  Pct proximal spread (median): %.1f%%  (100%% = contiguous growth)",
                  median(pd_ok$pct_proximal, na.rm = TRUE)))
  # Spatial offset warning
  low_prox <- pd_ok %>% filter(!is.na(pct_proximal) & pct_proximal < 20)
  if (nrow(low_prox) > 0)
    message(sprintf(
      "  ⚠  %d pair(s) have <20%% proximal spread — new burned area mostly outside frontier ring.",
      nrow(low_prox)
    ))
}

# ---- MIN_PCT_PROXIMAL filter --------------------------------------------------
# Drop pairs whose new burned area is spatially disconnected from the frontier ring.
# Controlled by MIN_PCT_PROXIMAL in Section 1 (0 = no filter, default).
if (MIN_PCT_PROXIMAL > 0 && nrow(all_samples) > 0) {
  keep_idx <- pair_diag %>%
    filter(status == "ok") %>%
    filter(!is.na(pct_proximal) & pct_proximal >= MIN_PCT_PROXIMAL) %>%
    pull(pair_idx)
  n_before <- nrow(all_samples)
  p_before <- sum(pair_diag$status == "ok")
  all_samples <- all_samples %>% filter(pair_idx %in% keep_idx)
  pd_ok       <- pair_diag %>% filter(status == "ok", pair_idx %in% keep_idx)
  message(sprintf(
    "\nMIN_PCT_PROXIMAL=%d%%: kept %d / %d ok pairs (%d / %d samples)",
    MIN_PCT_PROXIMAL, length(keep_idx), p_before, nrow(all_samples), n_before
  ))
  if (length(keep_idx) < 5)
    message("  ⚠  Very few pairs after proximal filter — GLM will be unstable.",
            "\n     Lower MIN_PCT_PROXIMAL or expand cal_vect to include more fires.")
} else {
  message(sprintf("\nMIN_PCT_PROXIMAL=0: using all %d ok pairs (no proximal filter)", nrow(pd_ok)))
}

message(sprintf("\nTotal GLM sample: %d rows  (%d successes, %d failures)",
                nrow(all_samples), sum(all_samples$y == 1), sum(all_samples$y == 0)))

# Pair-level signal check: log(n_new_burned) vs FWI and EWS.
# With the Robbins approach, n_new_burned is a proxy for spread magnitude per pair.
# If these correlations are near zero BEFORE fitting, the GLM has no signal.
if (nrow(pd_ok) > 2) {
  pd_ok_sig <- pd_ok %>% filter(n_new_burned > 0)
  log_new   <- log(pd_ok_sig$n_new_burned)
  cor_fwi   <- cor(pd_ok_sig$FWI, log_new, use = "complete.obs")
  cor_ews   <- cor(pd_ok_sig$EWS, log_new, use = "complete.obs")
  message(sprintf("\nPair-level Cor(FWI, log(n_new_burned)) = %.3f  <- should be > 0", cor_fwi))
  message(sprintf("Pair-level Cor(EWS, log(n_new_burned)) = %.3f  <- should be > 0", cor_ews))
  if (cor_fwi < 0.05)
    message(paste0(
      "⚠  Near-zero FWI-spread correlation at the pair level. Possible causes:\n",
      "   - FWI dates don't match perimeter dates (check join)\n",
      "   - Regional FWI is a poor proxy for fire-specific FWI in this area\n",
      "   - Too few pairs (noisy estimate)\n",
      "   - Spread is wind-dominated rather than FWI-dominated in this landscape"
    ))

  p1 <- ggplot(pd_ok_sig, aes(FWI, n_new_burned)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    scale_y_log10(labels = scales::comma) +
    labs(title = sprintf("New burned cells vs FWI  (r = %.2f, log scale)", cor_fwi),
         x = "FWI", y = "n new burned cells (log)") + theme_bw()

  p2 <- ggplot(pd_ok_sig, aes(EWS, n_new_burned)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    scale_y_log10(labels = scales::comma) +
    labs(title = sprintf("New burned cells vs EWS  (r = %.2f, log scale)", cor_ews),
         x = paste0(wind_col, " (km/h)"), y = "n new burned cells (log)") + theme_bw()

  print(p1); print(p2)
}

# ---- Fine fuels signal diagnostics (per-pair) --------------------------------
# This block answers the core question: does FF actually differ between spread
# success and failure cells WITHIN each perimeter pair?
#
# Because FWI and EWS are pair-level constants (every cell in a pair shares the
# same value), FF is the ONLY predictor that varies at the cell level within a
# pair.  B2 will only be identified if success cells consistently have higher FF
# than failure cells across multiple pairs.
if (!is.null(fine_fuels_r) && nrow(all_samples) > 0 &&
    "FineFuels" %in% names(all_samples) && var(all_samples$FineFuels) > 1e-9) {

  message("\n--- Fine fuels cell-level signal ---")

  # Per-pair: mean FF for success vs failure cells
  ff_pair <- all_samples %>%
    group_by(FIRE_ID, date, FWI) %>%
    summarise(
      FF_success = mean(FineFuels[y == 1], na.rm = TRUE),
      FF_failure = mean(FineFuels[y == 0], na.rm = TRUE),
      n_succ     = sum(y == 1),
      n_fail     = sum(y == 0),
      .groups    = "drop"
    ) %>%
    filter(is.finite(FF_success), is.finite(FF_failure))

  if (nrow(ff_pair) > 0) {
    ff_diff    <- ff_pair$FF_success - ff_pair$FF_failure
    # Pooled sd from the full sample
    ff_s1 <- all_samples$FineFuels[all_samples$y == 1]
    ff_s0 <- all_samples$FineFuels[all_samples$y == 0]
    ff_pooled_sd <- sqrt((var(ff_s1, na.rm=TRUE) + var(ff_s0, na.rm=TRUE)) / 2)
    ff_cohen_d   <- (mean(ff_s1, na.rm=TRUE) - mean(ff_s0, na.rm=TRUE)) /
                     max(ff_pooled_sd, 1e-9)

    message(sprintf("  mean FF in success cells (y=1): %.4f", mean(ff_s1, na.rm=TRUE)))
    message(sprintf("  mean FF in failure cells (y=0): %.4f", mean(ff_s0, na.rm=TRUE)))
    message(sprintf("  difference (success - failure): %+.4f", mean(ff_s1, na.rm=TRUE) - mean(ff_s0, na.rm=TRUE)))
    message(sprintf("  Cohen's d (cell level)        : %+.3f%s", ff_cohen_d,
                    if (abs(ff_cohen_d) < 0.1) "  ⚠ VERY WEAK — B2 may be near zero"
                    else if (ff_cohen_d > 0) "  (success > failure → expect B2 > 0)"
                    else "  (success < failure → B2 may be negative or near zero)"))

    # Per-pair comparison
    n_pairs_pos <- sum(ff_pair$FF_success > ff_pair$FF_failure, na.rm = TRUE)
    n_pairs_neg <- sum(ff_pair$FF_success < ff_pair$FF_failure, na.rm = TRUE)
    message(sprintf("  Pairs where FF_success > FF_failure: %d of %d (%.0f%%)",
                    n_pairs_pos, nrow(ff_pair),
                    100 * n_pairs_pos / nrow(ff_pair)))
    message(sprintf("  Pairs where FF_success < FF_failure: %d of %d (%.0f%%)",
                    n_pairs_neg, nrow(ff_pair),
                    100 * n_pairs_neg / nrow(ff_pair)))
    if (n_pairs_pos < 0.6 * nrow(ff_pair))
      message(paste0(
        "  ⚠ FF_success > FF_failure in fewer than 60% of pairs.\n",
        "     B2_FineFuels will likely be near zero or negative in the GLM.\n",
        "     This landscape may spread into low-FF areas (open areas, grass)\n",
        "     just as readily as high-FF areas (dense timber)."
      ))

    # Scatter: per-pair mean FF success vs failure
    ggplot(ff_pair, aes(FF_failure, FF_success)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = "Per-pair mean FF: spread-success cells vs spread-failure cells",
           subtitle = sprintf("Points above the diagonal = B2 > 0 signal  |  %d of %d pairs above",
                              n_pairs_pos, nrow(ff_pair)),
           x = "Mean FF in failure cells (y=0)",
           y = "Mean FF in success cells (y=1)") +
      theme_bw()
  }

  # Density overlay: FF for spread success vs failure (cell level)
  ggplot(all_samples, aes(FineFuels,
                          fill = factor(y, labels = c("Failure (y=0)", "Success (y=1)")))) +
    geom_density(alpha = 0.45) +
    scale_fill_manual(values = c("Failure (y=0)" = "steelblue",
                                 "Success (y=1)" = "tomato")) +
    labs(title = "Fine fuels distribution: spread success vs failure (all sample cells)",
         subtitle = "Rightward shift for y=1 = fires spread into higher-FF areas (B2 > 0)",
         x = sprintf("Fine fuel load (%s)", ff_units), y = "Density", fill = NULL) +
    theme_bw()

  # FF vs FWI at the pair level — is there a spurious correlation?
  # If high-FWI days happen to be in higher-FF fire areas, the GLM will struggle
  # to separate FWI and FF effects.
  if (nrow(ff_pair) > 5) {
    cor_fwi_ff <- cor(ff_pair$FWI,
                      (ff_pair$FF_success * ff_pair$n_succ +
                       ff_pair$FF_failure * ff_pair$n_fail) /
                      (ff_pair$n_succ + ff_pair$n_fail),
                      use = "complete.obs")
    message(sprintf("\n  Cor(FWI, mean pair FF): %.3f", cor_fwi_ff))
    if (abs(cor_fwi_ff) > 0.4)
      message(paste0(
        "  ⚠ FWI and mean FF are correlated at the pair level (|r|>0.4).\n",
        "     This suggests fires at different FWI levels also differ in fuel type,\n",
        "     making it hard for the GLM to cleanly separate B1 (FWI) and B2 (FF).\n",
        "     Check: do high-FWI fires occur in different vegetation than low-FWI fires?"
      ))

    ggplot(ff_pair,
           aes(FWI,
               (FF_success * n_succ + FF_failure * n_fail) / (n_succ + n_fail))) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(title = "Pair-level mean FF vs FWI",
           subtitle = sprintf("r = %.3f — strong correlation complicates joint B1/B2 estimation",
                              cor_fwi_ff),
           x = "FWI (day of pair)", y = "Mean FF in candidate cells") +
      theme_bw()
  }
}


# =============================================================================
# SECTION 13 — Predictor distributions and collinearity
# =============================================================================

message("\n===== SECTION 13: Predictor distributions and collinearity =====")

if (nrow(all_samples) == 0) stop("No samples — cannot proceed to GLM.")

# Variance check for fine fuels — determines whether B2 is identifiable
ff_var  <- var(all_samples$FineFuels, na.rm = TRUE)
drop_b2 <- is.null(fine_fuels_r) || ff_var < 0.001
message(sprintf("FineFuels variance : %.6f  -> B2 %s",
                ff_var, if (drop_b2) "DROPPED (unidentifiable)" else "identifiable"))

# Predictor summary
message("\nPredictor summary across all GLM samples:")
for (col in c("FWI", "EffectiveWind", "FineFuels")) {
  if (!col %in% names(all_samples)) next
  v <- all_samples[[col]]
  message(sprintf("  %-15s  mean=%7.2f  sd=%6.2f  range [%.2f, %.2f]  NAs=%d",
                  col, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE),
                  min(v, na.rm=TRUE), max(v, na.rm=TRUE), sum(is.na(v))))
}

# Cohen's d — do y=1 and y=0 samples differ on each predictor?
message("\nMean difference between spread (y=1) and no-spread (y=0) cells:")
for (col in c("FWI", "EffectiveWind")) {
  if (!col %in% names(all_samples)) next
  s1 <- all_samples[[col]][all_samples$y == 1]
  s0 <- all_samples[[col]][all_samples$y == 0]
  d  <- mean(s1, na.rm=TRUE) - mean(s0, na.rm=TRUE)
  pd_sd <- sqrt((var(s1, na.rm=TRUE) + var(s0, na.rm=TRUE)) / 2)
  cd <- d / max(pd_sd, 1e-9)
  message(sprintf("  %-15s  mean(y=1)=%6.2f  mean(y=0)=%6.2f  diff=%+.2f  Cohen's d=%+.2f%s",
                  col, mean(s1, na.rm=TRUE), mean(s0, na.rm=TRUE), d, cd,
                  if (abs(cd) < 0.1) "  ⚠ VERY WEAK SIGNAL" else ""))
}

# Collinearity
X_pred <- as.matrix(all_samples[, c("FWI", "EffectiveWind")])
X_pred <- X_pred[complete.cases(X_pred), ]
vif_vals <- tryCatch(diag(solve(cor(X_pred))), error = function(e) c(NA_real_, NA_real_))
sv       <- svd(cbind(1, X_pred))$d
kappa    <- max(sv) / min(sv[sv > .Machine$double.eps * max(sv)])
r_fwi_ew <- cor(X_pred[,"FWI"], X_pred[,"EffectiveWind"])
message(sprintf("\nVIF: FWI=%.2f  EWS=%.2f  (>10 = problem)",   vif_vals[1], vif_vals[2]))
message(sprintf("Condition number: %.1f  (>100 = problem)",   kappa))
message(sprintf("Cor(FWI, EWS)  : %.3f  (>0.8 = collinearity)", r_fwi_ew))

# Density overlay by outcome
all_samples %>%
  pivot_longer(c(FWI, EffectiveWind, FineFuels), names_to = "predictor", values_to = "value") %>%
  mutate(outcome = factor(y, levels = c(0L, 1L),
                          labels = c("No spread (y=0)", "Spread (y=1)"))) %>%
  ggplot(aes(value, fill = outcome)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~predictor, scales = "free") +
  labs(title = "Predictor distributions by spread outcome",
       subtitle = "Overlapping distributions = weak predictor signal",
       fill = NULL) +
  theme_bw()


# =============================================================================
# SECTION 14 — Fit and compare candidate GLMs
# =============================================================================
#
# Four models are always fit and compared side-by-side:
#
#   M1: y ~ FWI                          (FWI only — baseline signal check)
#   M2: y ~ FWI + FineFuels              (FWI + FF, no wind — Robbins "Par_logit")
#   M3: y ~ FWI + EffectiveWind          (FWI + wind, no FF)
#   M4: y ~ FWI + FineFuels + EWS        (full model — SCF default)
#
# If drop_b2=TRUE (FF has near-zero variance), M2 and M4 use FWI only / FWI+EWS.
#
# The primary output for SCF is m_spread (M4, or M2 if EWS is counterintuitive).
# Use the sign column to diagnose which predictors have the correct direction.
# =============================================================================

message("\n===== SECTION 14: Fit and compare spread probability GLMs =====")

.auc <- function(fit) {
  pr <- predict(fit, type = "response")
  y  <- fit$y
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  # Wilcoxon / Mann-Whitney AUC
  mean(sample(pr[y == 1], min(n1, 5000), replace = n1 < 5000) >
       sample(pr[y == 0], min(n0, 5000), replace = n0 < 5000))
}

.r2 <- function(fit) 1 - fit$deviance / fit$null.deviance

.fit_model <- function(formula, label) {
  message(sprintf("  Fitting %-50s ...", deparse(formula)))
  tryCatch(
    glm(formula, data = all_samples, family = binomial(link = "logit")),
    error   = function(e) { message("    ERROR: ", e$message); NULL },
    warning = function(w) {
      fit <- withCallingHandlers(
        glm(formula, data = all_samples, family = binomial(link = "logit")),
        warning = function(w) invokeRestart("muffleWarning"))
      message("    WARNING: ", w$message)
      fit
    })
}

# Build model list
model_specs <- if (drop_b2) {
  list(
    M1_FWI      = y ~ FWI,
    M3_FWI_EWS  = y ~ FWI + EffectiveWind
  )
} else {
  list(
    M1_FWI          = y ~ FWI,
    M2_FWI_FF       = y ~ FWI + FineFuels,
    M3_FWI_EWS      = y ~ FWI + EffectiveWind,
    M4_FWI_FF_EWS   = y ~ FWI + FineFuels + EffectiveWind
  )
}

model_fits <- lapply(model_specs, .fit_model)
model_fits <- Filter(Negate(is.null), model_fits)

# Primary model used for SCF output (prefer full if available)
m_spread <- if (!is.null(model_fits[["M4_FWI_FF_EWS"]])) {
  model_fits[["M4_FWI_FF_EWS"]]
} else if (!is.null(model_fits[["M3_FWI_EWS"]])) {
  model_fits[["M3_FWI_EWS"]]
} else model_fits[[1]]

# ---- Comparison table -------------------------------------------------------
message("\n--- Model comparison ---\n")
comp_rows <- lapply(names(model_fits), function(nm) {
  fit <- model_fits[[nm]]
  cf  <- coef(fit)
  tibble(
    Model     = nm,
    Formula   = deparse(formula(fit)),
    n         = nobs(fit),
    McF_R2    = round(.r2(fit), 4),
    AUC       = round(.auc(fit), 4),
    B0_sign   = if ("(Intercept)"   %in% names(cf)) sign(cf["(Intercept)"])   else NA_integer_,
    B1_FWI    = if ("FWI"           %in% names(cf)) round(cf["FWI"],         5) else NA_real_,
    B2_FF     = if ("FineFuels"     %in% names(cf)) round(cf["FineFuels"],   4) else NA_real_,
    B3_EWS    = if ("EffectiveWind" %in% names(cf)) round(cf["EffectiveWind"], 5) else NA_real_,
    FWI_sign  = if ("FWI" %in% names(cf)) ifelse(cf["FWI"] > 0, "✓ POS", "✗ NEG") else "—",
    FF_sign   = if ("FineFuels" %in% names(cf)) ifelse(cf["FineFuels"] > 0, "✓ POS", "✗ NEG") else "—",
    EWS_sign  = if ("EffectiveWind" %in% names(cf)) ifelse(cf["EffectiveWind"] > 0, "✓ POS", "✗ NEG") else "—"
  )
})
comp_tbl <- bind_rows(comp_rows)
print(as.data.frame(comp_tbl[, c("Model","McF_R2","AUC","B1_FWI","B2_FF","B3_EWS",
                                  "FWI_sign","FF_sign","EWS_sign")]),
      digits = 5, row.names = FALSE)

message(sprintf(
  "\nPrimary model (m_spread): %s  —  McFadden R2=%.4f  AUC=%.4f",
  deparse(formula(m_spread)), .r2(m_spread), .auc(m_spread)
))

coef_sum <- broom::tidy(m_spread) %>%
  rename(std_error = std.error, p_value = p.value) %>%
  mutate(
    expected_sign = case_when(
      term == "FWI"           ~ "+",
      term == "EffectiveWind" ~ "+",
      term == "FineFuels"     ~ "+",
      TRUE ~ NA_character_
    ),
    sign_ok = case_when(
      term == "(Intercept)"   ~ TRUE,
      TRUE ~ (estimate > 0)
    ),
    flag = dplyr::case_when(
      !sign_ok & !is.na(expected_sign) ~ "⚠ WRONG SIGN",
      p_value > 0.05                   ~ "⚠ NOT SIGNIFICANT",
      TRUE                              ~ "OK"
    )
  )

message("\nGLM coefficient table:")
for (i in seq_len(nrow(coef_sum))) {
  message(sprintf("  %-35s  est=%+8.4f  se=%.4f  p=%.4f  [%s]",
                  coef_sum$term[i],
                  coef_sum$estimate[i],
                  coef_sum$std_error[i],
                  coef_sum$p_value[i],
                  coef_sum$flag[i]))
}

mcfadden_r2 <- 1 - m_spread$deviance / m_spread$null.deviance
message(sprintf("\nNull deviance    : %.1f (df %d)", m_spread$null.deviance, m_spread$df.null))
message(sprintf("Residual deviance: %.1f (df %d)", m_spread$deviance, m_spread$df.residual))
message(sprintf("McFadden R²      : %.3f  (>0.2 reasonable, <0.05 poor)", mcfadden_r2))

# AUC — sampling-based (wilcox.test fails on large samples)
probs_fit <- predict(m_spread, type = "response")
y_true    <- all_samples$y
auc_est   <- tryCatch({
  n1 <- sum(y_true == 1); n0 <- sum(y_true == 0)
  if (n1 == 0 || n0 == 0) NA_real_ else
    mean(sample(probs_fit[y_true == 1], min(n1, 10000), replace = n1 < 10000) >
         sample(probs_fit[y_true == 0], min(n0, 10000), replace = n0 < 10000))
}, error = function(e) NA_real_)
message(sprintf("AUC              : %.3f  (0.5=chance, >0.7=good)", auc_est))

# Predicted probability curves
ews_q   <- quantile(all_samples$EffectiveWind, c(0.25, 0.5, 0.75), na.rm = TRUE)
fwi_seq <- seq(0, max(all_samples$FWI, na.rm = TRUE), length.out = 200)

pred_df <- map_dfr(
  list(list(ews = ews_q[1], lbl = sprintf("EWS Q25 = %.1f km/h", ews_q[1])),
       list(ews = ews_q[2], lbl = sprintf("EWS Q50 = %.1f km/h", ews_q[2])),
       list(ews = ews_q[3], lbl = sprintf("EWS Q75 = %.1f km/h", ews_q[3]))),
  function(s) tibble(FWI = fwi_seq, EffectiveWind = s$ews,
                     FineFuels = 1.0, EWS_label = s$lbl)
) %>%
  mutate(pred_prob = predict(m_spread, newdata = ., type = "response"))

ggplot(pred_df, aes(FWI, pred_prob, colour = EWS_label)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Fitted spread probability vs FWI",
       subtitle = "Curves should increase with FWI — if flat or decreasing, see diagnostics",
       x = "FWI", y = "P(spread)", colour = NULL) +
  theme_bw()

# Calibration plot
tibble(y = y_true, pred = probs_fit) %>%
  mutate(decile = ntile(pred, 10)) %>%
  group_by(decile) %>%
  summarise(mean_pred = mean(pred), obs_rate = mean(y), n = n(), .groups = "drop") %>%
  ggplot(aes(mean_pred, obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), color = "steelblue") +
  geom_line(color = "steelblue") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Calibration plot (by decile of predicted probability)",
       subtitle = "Points on the diagonal = perfect calibration",
       x = "Mean predicted probability", y = "Observed spread rate", size = "n cells") +
  theme_bw()


# =============================================================================
# SECTION 15 — Residual diagnostics
# =============================================================================

message("\n===== SECTION 15: Residual diagnostics =====")

resid_df <- tibble(
  FWI           = all_samples$FWI,
  EffectiveWind = all_samples$EffectiveWind,
  y             = y_true,
  fitted        = probs_fit,
  pearson_resid = residuals(m_spread, type = "pearson")
)

# ggplot(resid_df, aes(fitted, pearson_resid)) +
#   geom_point(alpha = 0.08, size = 0.4) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_smooth(method = "loess", color = "blue", se = FALSE) +
#   labs(title = "Pearson residuals vs fitted probability",
#        x = "Fitted P(spread)", y = "Pearson residual") +
#   theme_bw()
# 
# ggplot(resid_df, aes(FWI, pearson_resid)) +
#   geom_point(alpha = 0.08, size = 0.4) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_smooth(method = "loess", color = "blue", se = FALSE) +
#   labs(title = "Pearson residuals vs FWI",
#        x = "FWI", y = "Pearson residual") +
#   theme_bw()

# Overdispersion (binary logistic: expected phi ≈ 1; large phi = clustering)
phi <- sum(residuals(m_spread, type = "pearson")^2) / m_spread$df.residual
message(sprintf("Pearson dispersion phi = %.2f", phi))
if (phi > 1.5)
  message(paste0(
    "⚠  Overdispersion (phi=", round(phi, 2), "). All cells within a pair share the ",
    "same FWI and EWS values, so they are not independent observations. ",
    "This inflates effective sample size and deflates standard errors, making ",
    "coefficients appear more significant than they are. ",
    "Consider: fit on pair-level observed spread rates (binomial with n=n_cand) ",
    "instead of cell-level binary outcomes."
  ))


# =============================================================================
# SECTION 16 — Maximum Daily Spread Area calibration (FIRED Events)
#
# Uses the FIRED *events* product (fired_events_conus_*.gpkg — one row per
# fire event) to calibrate the SCF MaxSpreadArea model:
#
#   MaxSpreadArea = B0 + B1 * FWI + B2 * EWS   [units: hectares]
#
# max_growth_ha   = single largest daily burned area (ha) for each fire event
# max_growth_date = the calendar date on which maximum growth occurred
#
# Calibration target: 75th percentile of observed max daily areas.
# MaxSpreadArea is a CEILING in SCF — the fire advances to the next day when
# daily_burned > MaxSpreadArea.  Calibrating to the mean (OLS) systematically
# underpredicts spread.
#
# Two methods are compared:
#   QR75       : quantile regression at tau=0.75 (primary; requires quantreg)
#   OLS + Q75  : OLS fit with 75th-percentile residual intercept shift
#
# Runs independently of PERIMETER_SOURCE — FIRED events data is always used
# here regardless of whether spread probability used GeoMAC or FIRED daily.
# =============================================================================

message("\n===== SECTION 16: Maximum Daily Spread Area calibration (FIRED Events) =====")

# Initialise outputs so Section 17/18 don't fail if this section is skipped
fired_max      <- NULL
maxarea_coef   <- NULL
maxarea_models <- NULL

if (!file.exists(FIRED_EVENTS_GPKG)) {

  message("⚠  FIRED_EVENTS_GPKG not found:")
  message("     ", FIRED_EVENTS_GPKG)
  message("   Skipping MaxSpreadArea calibration. Set FIRED_EVENTS_GPKG in Section 1.")

} else {

  # ---- 16.1 Load FIRED events -------------------------------------------------
  message(sprintf("Loading FIRED events from: %s", basename(FIRED_EVENTS_GPKG)))

  fired_max <- tryCatch(
    load_fired_events_for_max_area(
      fired_events_gpkg = FIRED_EVENTS_GPKG,
      cal_vect          = cal_vect,
      working_crs       = working_crs,
      fwi_daily         = fwi_daily,
      wind_daily        = wind_daily,
      cal_start         = CAL_START,
      cal_end           = CAL_END,
      min_fire_ha       = FIRED_MIN_FIRE_HA
    ),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fired_max) || nrow(fired_max) == 0) {

    message("⚠  No FIRED events loaded. Check FIRED_EVENTS_GPKG path and date range.")

  } else {

    message(sprintf("\n%d events loaded (min_fire_ha=%.0f):", nrow(fired_max), FIRED_MIN_FIRE_HA))
    message(sprintf("  Max daily area range : %.0f – %.0f ha",
                    min(fired_max$daily_area_ha, na.rm=TRUE),
                    max(fired_max$daily_area_ha, na.rm=TRUE)))
    message(sprintf("  Date range           : %s to %s",
                    min(fired_max$date, na.rm=TRUE),
                    max(fired_max$date, na.rm=TRUE)))

    # ---- 16.2 Compute Effective Wind Speed ------------------------------------
    # load_fired_events_for_max_area() returns WindSpeed_kmh (raw).
    # Compute EWS using domain-mean slope/aspect (centroid geometry not available
    # in the events product — this is a reasonable calibration-domain approximation).
    ews_src_max <- "raw wind speed (WindSpeed_kmh)"   # default label

    if (USE_TERRAIN && !is.null(slope_r) && !is.null(aspect_r) &&
        "WindDir_deg" %in% names(wind_daily)) {

      # Join wind direction for each max-growth date
      fired_max <- fired_max %>%
        left_join(
          wind_daily %>%
            dplyr::select(date, WindDir_deg) %>%
            dplyr::filter(!is.na(WindDir_deg)),
          by = "date"
        )

      # Domain-mean slope and aspect across the calibration boundary
      maxarea_slope_mean  <- mean(values(slope_r,  mat = FALSE), na.rm = TRUE)
      maxarea_aspect_mean <- mean(values(aspect_r, mat = FALSE), na.rm = TRUE)

      fired_max <- fired_max %>%
        mutate(
          EWS = compute_effective_wind(
            wind_speed          = WindSpeed_kmh,
            wind_dir            = if ("WindDir_deg" %in% names(fired_max)) WindDir_deg else 0,
            slope_deg           = maxarea_slope_mean,
            aspect_deg          = maxarea_aspect_mean,
            combustion_buoyancy = COMBUSTION_BUOYANCY
          )
        )

      ews_src_max <- sprintf(
        "Nelson (2002), Cb=%.0f, domain-mean slope=%.1f deg",
        COMBUSTION_BUOYANCY, maxarea_slope_mean
      )
      message(sprintf("EWS computed (%s):", ews_src_max))
      message(sprintf("  EWS mean=%.1f km/h  raw wind mean=%.1f km/h",
                      mean(fired_max$EWS,          na.rm = TRUE),
                      mean(fired_max$WindSpeed_kmh, na.rm = TRUE)))

    } else {
      # Terrain not loaded or wind direction not available
      fired_max <- fired_max %>% mutate(EWS = WindSpeed_kmh)
      if (!USE_TERRAIN) {
        ews_src_max <- "raw wind speed (USE_TERRAIN=FALSE)"
      } else if (!("WindDir_deg" %in% names(wind_daily))) {
        ews_src_max <- "raw wind speed (WindDir_deg not in wind_daily)"
      } else {
        ews_src_max <- "raw wind speed (terrain rasters not available)"
      }
      message(sprintf("EWS = raw wind speed — %s.", ews_src_max))
    }

    fired_max$EWS <- as.numeric(fired_max$EWS)

    # ---- 16.3 Seasonal and FWI filters (mirror Section 12) -------------------
    if (!is.null(CAL_MONTHS) && length(CAL_MONTHS) > 0) {
      n_before  <- nrow(fired_max)
      fired_max <- fired_max %>% filter(lubridate::month(date) %in% CAL_MONTHS)
      message(sprintf("CAL_MONTHS filter [%s]: kept %d of %d events.",
                      paste(CAL_MONTHS, collapse = ","), nrow(fired_max), n_before))
    }
    if (!is.null(FWI_MIN) && FWI_MIN > 0) {
      n_before  <- nrow(fired_max)
      fired_max <- fired_max %>% filter(FWI >= FWI_MIN)
      message(sprintf("FWI_MIN filter (>= %.1f): kept %d of %d events.",
                      FWI_MIN, nrow(fired_max), n_before))
    }
    fired_max <- fired_max %>% filter(is.finite(FWI), is.finite(EWS))

    message(sprintf("Events for MaxSpreadArea fit: %d", nrow(fired_max)))

    if (nrow(fired_max) < 10) {
      message("⚠  Fewer than 10 events — model will be very unstable.")
      message("   Consider: lower FIRED_MIN_FIRE_HA, relax CAL_MONTHS or FWI_MIN.")
    }

    # ---- 16.4 Predictor summary and distribution plots -----------------------
    q_area <- quantile(fired_max$daily_area_ha, c(0.25, 0.50, 0.75, 0.90, 0.95),
                       na.rm = TRUE)
    message(sprintf("\nMax daily area distribution (ha):"))
    message(sprintf("  Q25=%.0f  Q50=%.0f  Q75=%.0f  Q90=%.0f  Q95=%.0f",
                    q_area[1], q_area[2], q_area[3], q_area[4], q_area[5]))
    message(sprintf("  Mean=%.0f  SD=%.0f",
                    mean(fired_max$daily_area_ha), sd(fired_max$daily_area_ha)))
    message(sprintf("FWI range  (max-growth days): %.1f to %.1f  mean=%.1f",
                    min(fired_max$FWI), max(fired_max$FWI), mean(fired_max$FWI)))
    message(sprintf("EWS range  (km/h)           : %.1f to %.1f  mean=%.1f",
                    min(fired_max$EWS), max(fired_max$EWS), mean(fired_max$EWS)))

    # Correlation diagnostics
    cor_fwi_area  <- cor(fired_max$FWI, fired_max$daily_area_ha, use = "complete.obs")
    cor_ews_area  <- cor(fired_max$EWS, fired_max$daily_area_ha, use = "complete.obs")
    cor_fwi_ews_m <- cor(fired_max$FWI, fired_max$EWS,           use = "complete.obs")
    message(sprintf("\nCor(FWI, max_area)  = %.3f  (should be > 0)", cor_fwi_area))
    message(sprintf("Cor(EWS, max_area)  = %.3f  (should be > 0)", cor_ews_area))
    message(sprintf("Cor(FWI, EWS)       = %.3f  (if > 0.8, multicollinearity risk)", cor_fwi_ews_m))

    if (cor_fwi_area < 0)
      message(paste0(
        "⚠  Negative FWI-area correlation. Likely causes:\n",
        "   - Fall wind-driven fires (high area, low FWI) dominate the sample.\n",
        "     Try CAL_MONTHS = 6:9 in Section 1.\n",
        "   - FWI range too narrow to show a monotone signal."
      ))

    # Distribution plot
    ggplot(fired_max, aes(daily_area_ha)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "white") +
      geom_vline(xintercept = q_area["75%"], color = "red",    linetype = "dashed",  linewidth = 1) +
      geom_vline(xintercept = q_area["50%"], color = "grey40", linetype = "dotted",  linewidth = 0.8) +
      scale_x_continuous(labels = scales::comma) +
      annotate("text", x = q_area["75%"], y = Inf,
               label = "Q75 (calibration target)", hjust = -0.1, vjust = 1.5,
               color = "red", size = 3) +
      labs(title = "FIRED events: max daily burned area distribution",
           subtitle = sprintf("n=%d  |  Q75=%.0f ha  |  Mean=%.0f ha",
                              nrow(fired_max), q_area["75%"], mean(fired_max$daily_area_ha)),
           x = "Max daily area (ha)", y = "Count") +
      theme_bw()

    # Scatter: FWI vs max daily area
    ggplot(fired_max, aes(FWI, daily_area_ha)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 1) +
      scale_y_continuous(labels = scales::comma) +
      labs(title = sprintf("Max daily area vs FWI  (Pearson r = %.3f)", cor_fwi_area),
           x = "FWI (max-growth day)", y = "Max daily area (ha)") +
      theme_bw()

    # Scatter: EWS vs max daily area
    ggplot(fired_max, aes(EWS, daily_area_ha)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", color = "blue", se = TRUE, linewidth = 1) +
      scale_y_continuous(labels = scales::comma) +
      labs(title = sprintf("Max daily area vs EWS  (Pearson r = %.3f)", cor_ews_area),
           x = "EWS (km/h)", y = "Max daily area (ha)") +
      theme_bw()

    # ---- 16.5 Fit models (OLS + Quantile regression at tau=0.75) -------------
    message("\n--- Fitting MaxSpreadArea models ---")

    # OLS: mean regression — reference / comparison only
    fit_ols_max <- tryCatch(
      lm(daily_area_ha ~ FWI + EWS, data = fired_max),
      error = function(e) { message("  OLS ERROR: ", e$message); NULL }
    )

    # Quantile regression at tau=0.75 — full model (FWI + EWS)
    fit_q75_max <- tryCatch(
      quantreg::rq(daily_area_ha ~ FWI + EWS, tau = 0.75, data = fired_max),
      error = function(e) { message("  Quantile regression (tau=0.75) ERROR: ", e$message); NULL }
    )

    # Quantile regression at tau=0.75 — FWI-only model
    # Useful when B0 in the full model is deeply negative (EWS/FWI anti-correlation
    # forces the intercept below zero to fit the conditional quantile in the data
    # range). FWI-only often yields a positive, interpretable intercept.
    fit_q75_fwi_max <- tryCatch(
      quantreg::rq(daily_area_ha ~ FWI, tau = 0.75, data = fired_max),
      error = function(e) { message("  QR (FWI-only) ERROR: ", e$message); NULL }
    )

    # OLS + Q75 residual shift: simple alternative if quantreg is unavailable
    # Intercept is shifted upward by the 75th percentile of the OLS residuals so
    # that ~75% of observations fall below the predicted line.
    fit_ols_shifted <- if (!is.null(fit_ols_max)) {
      q75_resid <- quantile(residuals(fit_ols_max), 0.75, na.rm = TRUE)
      cf        <- coef(fit_ols_max)
      list(
        B0     = cf[["(Intercept)"]] + q75_resid,
        B1     = cf[["FWI"]],
        B2     = cf[["EWS"]],
        method = sprintf("OLS + Q75 intercept shift (+%.1f ha)", q75_resid)
      )
    } else NULL

    # Print comparison
    if (!is.null(fit_ols_max)) {
      cf <- coef(fit_ols_max)
      message(sprintf("  OLS (mean):        B0=%+7.1f  B1=%+.4f  B2=%+.4f  R²=%.3f",
                      cf[["(Intercept)"]], cf[["FWI"]], cf[["EWS"]],
                      summary(fit_ols_max)$r.squared))
    }
    if (!is.null(fit_ols_shifted)) {
      message(sprintf("  OLS + Q75 shift:   B0=%+7.1f  B1=%+.4f  B2=%+.4f",
                      fit_ols_shifted$B0, fit_ols_shifted$B1, fit_ols_shifted$B2))
    }
    if (!is.null(fit_q75_fwi_max)) {
      cf <- coef(fit_q75_fwi_max)
      # Compute coverage for FWI-only model
      pred_fwi <- cf[["(Intercept)"]] + cf[["FWI"]] * fired_max$FWI
      pct_fwi  <- 100 * mean(fired_max$daily_area_ha < pred_fwi, na.rm = TRUE)
      message(sprintf("  QR tau=0.75 (FWI): B0=%+7.1f  B1=%+.4f  B2=— (FWI-only)  cov=%.1f%%",
                      cf[["(Intercept)"]], cf[["FWI"]], pct_fwi))
    }
    if (!is.null(fit_q75_max)) {
      cf <- coef(fit_q75_max)
      pred_full <- cf[["(Intercept)"]] + cf[["FWI"]] * fired_max$FWI + cf[["EWS"]] * fired_max$EWS
      pct_full  <- 100 * mean(fired_max$daily_area_ha < pred_full, na.rm = TRUE)
      message(sprintf("  QR tau=0.75 (full):B0=%+7.1f  B1=%+.4f  B2=%+.4f  cov=%.1f%%",
                      cf[["(Intercept)"]], cf[["FWI"]], cf[["EWS"]], pct_full))
    }

    # ---- 16.6 Select primary model and validate ------------------------------
    # Preference order:
    #   1. QR(FWI+EWS) if B0 > 0
    #   2. QR(FWI-only) — preferred when full model has negative B0, since removing
    #      the anti-correlated EWS predictor lets the intercept absorb baseline area
    #   3. OLS+Q75 shift
    #   4. OLS (mean) — not recommended for a ceiling parameter

    # Determine if full QR model has a valid (>0) intercept
    full_q75_b0_ok <- !is.null(fit_q75_max) && coef(fit_q75_max)[["(Intercept)"]] > 0

    maxarea_coef <- if (full_q75_b0_ok) {
      cf <- coef(fit_q75_max)
      list(
        B0     = cf[["(Intercept)"]],
        B1     = cf[["FWI"]],
        B2     = cf[["EWS"]],
        method = "quantile regression (tau=0.75), FWI + EWS"
      )
    } else if (!is.null(fit_q75_fwi_max) && coef(fit_q75_fwi_max)[["(Intercept)"]] > 0) {
      cf <- coef(fit_q75_fwi_max)
      if (!full_q75_b0_ok && !is.null(fit_q75_max))
        message(paste0(
          "  → Full QR model has B0 <= 0 — falling back to FWI-only QR.\n",
          "    This commonly occurs when FWI and EWS are anti-correlated (Cor < 0)\n",
          "    and high EWS/low FWI days drive the intercept below zero.\n",
          "    The FWI-only model trades EWS precision for a valid intercept."
        ))
      list(
        B0     = cf[["(Intercept)"]],
        B1     = cf[["FWI"]],
        B2     = 0.0,
        method = "quantile regression (tau=0.75), FWI-only (full model B0 < 0)"
      )
    } else if (!is.null(fit_ols_shifted)) {
      message("  → All QR models have B0 <= 0 — falling back to OLS + Q75 intercept shift.")
      fit_ols_shifted
    } else if (!is.null(fit_ols_max)) {
      message("  → WARNING: Using OLS mean regression as last resort — not recommended.")
      cf <- coef(fit_ols_max)
      list(
        B0     = cf[["(Intercept)"]],
        B1     = cf[["FWI"]],
        B2     = cf[["EWS"]],
        method = "OLS — mean regression (WARNING: not recommended for SCF ceiling parameter)"
      )
    } else NULL

    if (!is.null(maxarea_coef)) {

      message(sprintf("\n--- MaxSpreadArea PRIMARY model (%s) ---", maxarea_coef$method))
      message(sprintf("  B0 = %+.2f   (intercept; must be > 0 to avoid SCF warnings)", maxarea_coef$B0))
      message(sprintf("  B1 = %+.4f   (FWI; should be > 0)", maxarea_coef$B1))
      message(sprintf("  B2 = %+.4f   (EWS; should be > 0)", maxarea_coef$B2))

      # Warning flags
      if (maxarea_coef$B0 <= 0)
        message(paste0(
          "  ⚠  B0 <= 0. At FWI=0 and EWS=0, MaxSpreadArea will be <= 0.\n",
          "     SCF will generate 'MaxSpreadArea < 0' warnings.\n",
          "     Fix: clamp with  maxarea_coef$B0 <- max(maxarea_coef$B0, 0.1)"
        ))
      if (maxarea_coef$B1 < 0)
        message(paste0(
          "  ⚠  B1 < 0 (negative FWI coefficient — MaxSpreadArea decreases with FWI).\n",
          "     Likely cause: fall wind-driven fires dominate the sample at low FWI.\n",
          "     Try CAL_MONTHS = 6:9 in Section 1."
        ))
      if (maxarea_coef$B2 < 0)
        message("  ⚠  B2 < 0 (negative EWS coefficient) — check wind data quality.")

      # Coverage check
      pred_vals <- maxarea_coef$B0 + maxarea_coef$B1 * fired_max$FWI +
                   maxarea_coef$B2 * fired_max$EWS
      pct_below <- 100 * mean(fired_max$daily_area_ha < pred_vals, na.rm = TRUE)
      message(sprintf("\n  %% events below predicted MaxSpreadArea: %.1f%%  (target: 75%%)", pct_below))
      if (abs(pct_below - 75) > 10)
        message(sprintf(paste0(
          "  ⚠  Coverage (%.1f%%) deviates >10pp from target 75%%.\n",
          "     With quantile regression this may indicate extreme leverage from very\n",
          "     large fires.  Consider raising FIRED_MIN_FIRE_HA or trimming outliers."
        ), pct_below))

      # Spot-check predicted values at key FWI values
      ews_mean_max <- mean(fired_max$EWS, na.rm = TRUE)
      message(sprintf("\nSpot-check MaxSpreadArea at mean EWS = %.1f km/h:", ews_mean_max))
      message(sprintf("  %-8s  %s", "FWI", "MaxSpreadArea (ha/day)"))
      for (fwi_v in c(5, 10, 20, 30, 50, 70, 90)) {
        pred_ha <- maxarea_coef$B0 + maxarea_coef$B1 * fwi_v + maxarea_coef$B2 * ews_mean_max
        message(sprintf("  FWI = %3d   %6.0f ha%s",
                        fwi_v, pred_ha,
                        if (pred_ha <= 0) "  ⚠ NEGATIVE" else ""))
      }

      # ---- 16.7 Plot: observed vs Q75 model line across FWI ----------------
      fwi_seq_max <- seq(min(fired_max$FWI), max(fired_max$FWI), length.out = 100)
      pred_line_df <- data.frame(
        FWI           = fwi_seq_max,
        MaxSpreadArea = maxarea_coef$B0 + maxarea_coef$B1 * fwi_seq_max +
                        maxarea_coef$B2 * ews_mean_max
      )

      ggplot(fired_max, aes(FWI, daily_area_ha)) +
        geom_point(alpha = 0.5, color = "steelblue") +
        geom_line(data = pred_line_df, aes(y = MaxSpreadArea),
                  color = "red", linewidth = 1.2) +
        scale_y_continuous(labels = scales::comma) +
        labs(
          title    = "MaxSpreadArea model: observed vs Q75 fit",
          subtitle = sprintf(
            "%s\nB0=%+.1f  B1=%+.4f  B2=%+.4f  |  %.1f%% of events below line",
            maxarea_coef$method,
            maxarea_coef$B0, maxarea_coef$B1, maxarea_coef$B2, pct_below),
          x       = "FWI (max-growth day)",
          y       = "Max daily area (ha)",
          caption = sprintf("Red line = MaxSpreadArea at mean EWS = %.1f km/h", ews_mean_max)
        ) +
        theme_bw()

      # Also compare OLS mean line vs Q75 line on the same plot
      if (!is.null(fit_ols_max)) {
        ols_cf  <- coef(fit_ols_max)
        compare_df <- rbind(
          data.frame(
            FWI   = fwi_seq_max,
            Value = maxarea_coef$B0 + maxarea_coef$B1 * fwi_seq_max +
                    maxarea_coef$B2 * ews_mean_max,
            Model = sprintf("%s (PRIMARY)", maxarea_coef$method)
          ),
          data.frame(
            FWI   = fwi_seq_max,
            Value = ols_cf[["(Intercept)"]] + ols_cf[["FWI"]] * fwi_seq_max +
                    ols_cf[["EWS"]] * ews_mean_max,
            Model = "OLS mean (reference)"
          )
        )
        ggplot(fired_max, aes(FWI, daily_area_ha)) +
          geom_point(alpha = 0.4, color = "steelblue") +
          geom_line(data = compare_df, aes(y = Value, color = Model, linetype = Model),
                    linewidth = 1.0) +
          scale_color_manual(values = c("red", "grey40")) +
          scale_linetype_manual(values = c("solid", "dashed")) +
          scale_y_continuous(labels = scales::comma) +
          labs(title    = "MaxSpreadArea: OLS mean vs Q75 (primary) model",
               subtitle = sprintf("At mean EWS = %.1f km/h", ews_mean_max),
               x = "FWI (max-growth day)", y = "Max daily area (ha)",
               color = NULL, linetype = NULL) +
          theme_bw() +
          theme(legend.position = "bottom")
      }
    }

    maxarea_models <- list(
      ols         = fit_ols_max,
      q75         = fit_q75_max,
      q75_fwi     = fit_q75_fwi_max,
      ols_shifted = fit_ols_shifted,
      coef        = maxarea_coef,
      data        = fired_max,
      ews_src     = ews_src_max
    )
  }
}


# =============================================================================
# SECTION 17 — SCF parameter summary
# =============================================================================

message("\n===== SECTION 17: SCF parameter summary =====")

b_coef <- coef(m_spread)
b0 <- b_coef["(Intercept)"]
b1 <- b_coef["FWI"]
b2 <- if (!drop_b2 && "FineFuels" %in% names(b_coef)) b_coef["FineFuels"] else 0.0
b3 <- b_coef["EffectiveWind"]

message("\n--- SCF Spread Probability parameters (enter in SCF parameter file) ---")
message(sprintf("Fine fuels units used in calibration: %s  (NORMALIZE_FINE_FUELS=%s)",
                ff_units, NORMALIZE_FINE_FUELS))
if (!NORMALIZE_FINE_FUELS && !drop_b2)
  message(sprintf(paste0(
    "NOTE: B2 is in logit units per %s.\n",
    "  SCF must receive FF values in the same units at runtime, OR\n",
    "  convert B2 to normalized units by dividing by ff_scale (%.4f):\n",
    "  B2_normalized = B2_raw / ff_scale = %+.6f / %.4f = %+.6f"
  ), ff_units, ff_scale, b2, ff_scale, b2 / ff_scale))

message(sprintf("SpreadProbabilityB0 = %+.4f   (intercept)", b0))
message(sprintf("SpreadProbabilityB1 = %+.4f   (FWI; should be > 0)", b1))
message(sprintf("SpreadProbabilityB2 = %+.4f   (FineFuels [%s]; %s)",
                b2, ff_units,
                if (drop_b2) "set to 0 — FF unidentifiable"
                else if (NORMALIZE_FINE_FUELS) "fitted, 0-1 normalized"
                else sprintf("fitted, raw units — divide by %.4f for normalized equiv.", ff_scale)))
message(sprintf("SpreadProbabilityB3 = %+.4f   (EWS; should be > 0)", b3))
message("\nFormula: P(spread) = 1 / (1 + exp(-(B0 + B1*FWI + B2*FF + B3*EWS)))")

# Spot-check P(spread) at a few FWI values using mean FF in its actual units
ff_mean_val <- mean(all_samples$FineFuels, na.rm = TRUE)
ews_mean    <- mean(all_samples$EffectiveWind, na.rm = TRUE)
message(sprintf("\nSpot-check at mean EWS=%.1f km/h, mean FF=%.3f %s:",
                ews_mean, ff_mean_val, ff_units))
for (fwi_v in c(10, 30, 60, 90)) {
  p_check <- plogis(b0 + b1 * fwi_v + b2 * ff_mean_val + b3 * ews_mean)
  message(sprintf("  FWI=%3d -> P(spread) = %.1f%%", fwi_v, 100 * p_check))
}

# Also show the equivalent normalized B2 for easy comparison across runs
if (!drop_b2 && !NORMALIZE_FINE_FUELS && ff_scale != 1.0) {
  message(sprintf("\nFor comparison — equivalent normalized B2 (if SCF uses 0-1 FF):"))
  message(sprintf("  B2_normalized = %+.4f / %.4f = %+.6f", b2, ff_scale, b2 / ff_scale))
  message(sprintf("  (Run again with NORMALIZE_FINE_FUELS=TRUE to confirm this value)"))
}

# ---- MaxSpreadArea parameters (from Section 16) ------------------------------
message("\n--- SCF Maximum Spread Area parameters (enter in SCF parameter file) ---")
if (exists("maxarea_coef") && !is.null(maxarea_coef)) {
  message(sprintf("MaximumSpreadAreaB0 = %+.2f   (intercept; must be > 0)", maxarea_coef$B0))
  message(sprintf("MaximumSpreadAreaB1 = %+.4f   (FWI; should be > 0)", maxarea_coef$B1))
  message(sprintf("MaximumSpreadAreaB2 = %+.4f   (EWS; should be > 0)", maxarea_coef$B2))
  message(sprintf("Calibration method  : %s", maxarea_coef$method))
  message("\nFormula: MaxSpreadArea = (int)(B0 + B1*FWI + B2*EWS)  [units: hectares]")
  message("NOTE: MaxSpreadArea is a CEILING — SCF advances fire to next day when")
  message("      daily_burned > MaxSpreadArea. B0 must be > 0 to prevent SCF warnings.")
  # Spot-check at same FWI values as spread probability
  if (exists("fired_max") && !is.null(fired_max)) {
    ews_show_max <- mean(fired_max$EWS, na.rm = TRUE)
    message(sprintf("\nSpot-check MaxSpreadArea at mean EWS=%.1f km/h:", ews_show_max))
    for (fwi_v in c(10, 30, 60, 90)) {
      ha <- maxarea_coef$B0 + maxarea_coef$B1 * fwi_v + maxarea_coef$B2 * ews_show_max
      message(sprintf("  FWI=%3d -> MaxSpreadArea = %.0f ha/day", fwi_v, ha))
    }
  }
} else {
  message("MaxSpreadArea not fitted — run Section 16 with a valid FIRED_EVENTS_GPKG.")
  message("  Set FIRED_EVENTS_GPKG in Section 1 (fired_events_conus_*.gpkg).")
}

# ---- Issue summary ----------------------------------------------------------
message("\n--- Issue summary ---")
n_issues <- 0L
if (any(!coef_sum$sign_ok & !is.na(coef_sum$expected_sign))) {
  n_issues <- n_issues + 1L
  bad <- coef_sum %>% filter(!sign_ok, !is.na(expected_sign)) %>% pull(term)
  message(sprintf("  [%d] Wrong-sign coefficient(s): %s", n_issues, paste(bad, collapse = ", ")))
  message("       -> See Sections 10, 13: FWI/EWS collinearity or near-zero signal")
}
if (!is.na(auc_est) && auc_est < 0.60) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] Low AUC (%.3f) — weak discrimination", n_issues, auc_est))
  message("       -> See Section 12: pair-level spread rate vs FWI/EWS plots")
}
if (phi > 2.0) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] Overdispersion phi=%.1f — cell samples not independent", n_issues, phi))
  message("       -> Consider fitting on pair-level rates instead of cell-level binary")
}
if (!is.na(r_fwi_ew) && abs(r_fwi_ew) > 0.8) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] FWI/EWS correlation=%.3f — multicollinearity", n_issues, r_fwi_ew))
  message("       -> Try dropping EffectiveWind or using raw wind instead of EWS")
}
if (mcfadden_r2 < 0.05) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] McFadden R²=%.3f — very poor fit", n_issues, mcfadden_r2))
  message("       -> Check FWI date join (Section 9) and single-pair example (Section 11)")
}
if (fwi_join_rate < 0.5) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] FWI join rate=%.0f%% — most pairs missing climate data",
                  n_issues, 100 * fwi_join_rate))
  message("       -> Check ERA NetCDF date range vs GeoMAC perimeter dates (Section 9)")
}
ff_var_check <- if (!is.null(fine_fuels_r)) {
  var(values(fine_fuels_r, mat = FALSE), na.rm = TRUE)
} else 0.0
if (ff_var_check < 0.001) {
  n_issues <- n_issues + 1L
  message(sprintf("  [%d] Fine fuels variance=%.6f — B2_FineFuels unidentifiable",
                  n_issues, ff_var_check))
  if (is.null(fine_fuels_r))
    message("       -> Set FINE_FUELS_SOURCE='landfire' in Section 1 and rerun Section 7")
  else
    message("       -> Fine fuels loaded but spatially constant — check LANDFIRE fetch")
}

if (n_issues == 0L) {
  message("  No major issues detected. Review Section 14 plots for overall fit quality.")
} else {
  message(sprintf("\n%d issue(s) flagged above — work through sections interactively to investigate.", n_issues))
}


# =============================================================================
# SECTION 18 — Diagnostic report (paste back to Claude for interpretation)
#
# Run this section after completing Sections 8–17 (or as far as you got).
# It collects every key statistic into a single structured text report and
# writes it to a .txt file in your working directory.  Copy the full contents
# of that file and paste it into Claude for targeted advice.
#
# Safe to run even if some sections failed — all variables are checked with
# exists() / tryCatch() before use.
# =============================================================================

message("\n===== SECTION 18: Building diagnostic report =====")

# ---- Helper: safe scalar fetch -----------------------------------------------
.s <- function(expr, missing = "NA") {
  tryCatch(
    { v <- eval(expr); if (length(v) == 0 || all(is.na(v))) missing else as.character(v[1]) },
    error = function(e) missing
  )
}

.fmt <- function(x, digits = 4) {
  tryCatch(sprintf(paste0("%.", digits, "f"), as.numeric(x)), error = function(e) "NA")
}

# ---- Build report lines -------------------------------------------------------
rpt <- character(0)
rpt_add <- function(...) rpt <<- c(rpt, paste0(...))
hr  <- function(ch = "=", n = 70) rpt_add(paste(rep(ch, n), collapse = ""))
blk <- function(label, val)       rpt_add(sprintf("  %-38s %s", paste0(label, ":"), val))

hr()
rpt_add("SCF SPREAD PROBABILITY DIAGNOSTIC REPORT")
rpt_add(sprintf("Generated : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
hr()

# ---- [1] Settings ------------------------------------------------------------
hr("-")
rpt_add("[1] SETTINGS")
hr("-")
blk("CAL_START",         .s(CAL_START))
blk("CAL_END",           .s(CAL_END))
blk("SPREAD_RES_M",      .s(if (exists("spread_res_m")) spread_res_m else SPREAD_RES_M))
blk("DILATION_CELLS",    .s(DILATION_CELLS))
blk("MAX_GAP_SP",        .s(MAX_GAP_SP))
blk("MAX_GAP_DA",        .s(MAX_GAP_DA))
blk("NEG_TOL_HA",        .s(NEG_TOL_HA))
blk("FAILURE_RATIO",     .s(FAILURE_RATIO))
blk("MAX_SAMPLES",       .s(MAX_SAMPLES))
blk("MIN_PCT_PROXIMAL",  .s(MIN_PCT_PROXIMAL))
blk("CAL_MONTHS",       if (is.null(CAL_MONTHS)) "all" else paste(CAL_MONTHS, collapse=","))
blk("FWI_MIN",          if (is.null(FWI_MIN) || FWI_MIN <= 0) "none" else .s(FWI_MIN))
blk("FINE_FUELS_SOURCE", .s(FINE_FUELS_SOURCE))
blk("NORMALIZE_FF",      .s(NORMALIZE_FINE_FUELS))
blk("USE_TERRAIN",       .s(USE_TERRAIN))
blk("COMBUSTION_BUOYANCY", .s(COMBUSTION_BUOYANCY))
blk("WIND_SOURCE",       .s(WIND_SOURCE))
blk("CAL_VECT",          tryCatch(
  if (!is.null(CAL_VECT_SHP)) CAL_VECT_SHP
  else paste("tigris states:", paste(CAL_VECT_STATES, collapse = " + ")),
  error = function(e) "NA"))
blk("TEMPLATE_TIF",      tryCatch(basename(TEMPLATE_TIF), error = function(e) "NA"))
blk("PERIMETER_SOURCE",  tryCatch(.s(PERIMETER_SOURCE), error = function(e) "NA"))
if (exists("PERIMETER_SOURCE") && PERIMETER_SOURCE == "fired") {
  blk("FIRED_GPKG",       tryCatch(basename(FIRED_GPKG),       error = function(e) "NA"))
  blk("FIRED_MIN_FIRE_HA",tryCatch(.s(FIRED_MIN_FIRE_HA),      error = function(e) "NA"))
} else {
  blk("GEOMAC_GDB",       tryCatch(basename(GEOMAC_GDB),       error = function(e) "NA"))
}

# ---- [2] Spatial grid --------------------------------------------------------
hr("-")
rpt_add("[2] SPATIAL GRID")
hr("-")
if (exists("spread_template")) {
  blk("spread grid (rows x cols)", sprintf("%d x %d", nrow(spread_template), ncol(spread_template)))
  blk("spread grid resolution (m)",  .fmt(res(spread_template)[1], 1))
  blk("spread grid cell area (ha)",  .fmt(prod(res(spread_template)) / 10000, 2))
  blk("spread grid extent",          as.character(ext(spread_template)))
  blk("n active cells (in mask)",
      if (exists("spread_mask")) as.character(sum(values(spread_mask) == 1L, na.rm = TRUE))
      else "mask not built")
} else {
  rpt_add("  spread_template not found — run Section 7")
}
if (exists("template_r"))
  blk("template_r resolution (m)", .fmt(res(template_r)[1], 1))
if (exists("working_crs"))
  blk("working CRS", tryCatch(crs(template_r, describe = TRUE)$name, error = function(e) "unknown"))
if (exists("dil_kernel"))
  blk("dilation kernel",
      sprintf("%dx%d disk (DILATION_CELLS=%d)", nrow(dil_kernel), ncol(dil_kernel), DILATION_CELLS))

# ---- [3] FWI ----------------------------------------------------------------
hr("-")
rpt_add("[3] FWI (ERA-Interim)")
hr("-")
if (exists("fwi_daily") && nrow(fwi_daily) > 0) {
  blk("date range",   sprintf("%s  to  %s", min(fwi_daily$date), max(fwi_daily$date)))
  blk("n days",       nrow(fwi_daily))
  blk("min / mean / max", sprintf("%.1f / %.1f / %.1f",
                                   min(fwi_daily$value, na.rm=TRUE),
                                   mean(fwi_daily$value, na.rm=TRUE),
                                   max(fwi_daily$value, na.rm=TRUE)))
  blk("n NA",         sum(is.na(fwi_daily$value)))
} else {
  rpt_add("  fwi_daily not found — run Section 4")
}

# ---- [4] Wind ---------------------------------------------------------------
hr("-")
rpt_add("[4] WIND")
hr("-")
if (exists("wind_daily") && nrow(wind_daily) > 0) {
  blk("source",       .s(WIND_SOURCE))
  blk("date range",   sprintf("%s  to  %s", min(wind_daily$date), max(wind_daily$date)))
  blk("n days",       nrow(wind_daily))
  blk("WindSpeed_kmh min/mean/max",
      sprintf("%.1f / %.1f / %.1f",
              min(wind_daily$WindSpeed_kmh, na.rm=TRUE),
              mean(wind_daily$WindSpeed_kmh, na.rm=TRUE),
              max(wind_daily$WindSpeed_kmh, na.rm=TRUE)))
} else {
  rpt_add("  wind_daily not found — run Section 5")
}

# ---- [5] EWS / terrain -------------------------------------------------------
hr("-")
rpt_add("[5] EFFECTIVE WIND SPEED / TERRAIN")
hr("-")
if (exists("all_pairs") && "EffectiveWind" %in% names(all_pairs) &&
    any(is.finite(all_pairs$EffectiveWind))) {
  blk("terrain loaded",   "YES")
  blk("EWS mean (km/h)",  .fmt(mean(all_pairs$EffectiveWind, na.rm=TRUE), 2))
  blk("EWS sd   (km/h)",  .fmt(sd(all_pairs$EffectiveWind, na.rm=TRUE), 2))
  blk("EWS range (km/h)", sprintf("%.1f  to  %.1f",
                                    min(all_pairs$EffectiveWind, na.rm=TRUE),
                                    max(all_pairs$EffectiveWind, na.rm=TRUE)))
  blk("Cor(FWI, EWS)",
      .fmt(cor(all_pairs$FWI, all_pairs$EffectiveWind, use="complete.obs"), 3))
  blk("Cor(FWI, raw wind)",
      .fmt(cor(all_pairs$FWI, all_pairs$WindSpeed_kmh, use="complete.obs"), 3))
  blk("Cb (combustion buoyancy)", .s(COMBUSTION_BUOYANCY))
  if (exists("slope_r"))
    blk("mean slope (deg)",
        .fmt(mean(as.numeric(values(slope_r, mat=FALSE)), na.rm=TRUE), 1))
} else if (exists("all_pairs") && "WindSpeed_kmh" %in% names(all_pairs)) {
  blk("terrain loaded",   "NO (raw wind used as EWS)")
  blk("WindSpeed_kmh mean", .fmt(mean(all_pairs$WindSpeed_kmh, na.rm=TRUE), 2))
  blk("Cor(FWI, wind)",
      .fmt(cor(all_pairs$FWI, all_pairs$WindSpeed_kmh, use="complete.obs"), 3))
} else {
  rpt_add("  Pairs not built — run Sections 5–9")
}

# ---- [6] Fine fuels ----------------------------------------------------------
hr("-")
rpt_add("[6] FINE FUELS")
hr("-")
blk("source",       .s(FINE_FUELS_SOURCE))
blk("normalize",    .s(NORMALIZE_FINE_FUELS))
blk("units",        tryCatch(ff_units, error = function(e) "unknown"))
if (exists("fine_fuels_r") && !is.null(fine_fuels_r)) {
  ff_v <- values(fine_fuels_r, mat=FALSE)
  ff_v_fin <- ff_v[is.finite(ff_v)]
  blk("range",      sprintf("%.4f  to  %.4f", min(ff_v_fin), max(ff_v_fin)))
  blk("mean",       .fmt(mean(ff_v_fin), 4))
  blk("variance",   .fmt(var(ff_v_fin), 6))
  blk("n NA",       sum(is.na(ff_v)))
  blk("NA fraction (all cells)", sprintf("%.1f%%", 100 * mean(is.na(ff_v))))
  if (exists("spread_mask")) {
    ff_in <- values(fine_fuels_r * spread_mask, mat=FALSE)
    blk("NA fraction (inside mask)",
        sprintf("%.1f%%", 100 * mean(is.na(ff_in))))
  }
} else {
  blk("fine_fuels_r", "NULL (placeholder 1.0 used)")
}
blk("drop_b2 (FF unidentifiable)", tryCatch(as.character(drop_b2), error = function(e) "unknown"))

# ---- [7] GeoMAC --------------------------------------------------------------
hr("-")
rpt_add("[7] GEOMAC PERIMETERS")
hr("-")
if (exists("geomac_v")) {
  blk("n perimeters loaded",  nrow(geomac_v))
  blk("n unique fires",       length(unique(geomac_v$FIRE_ID)))
  blk("date range",           sprintf("%s  to  %s",
                                       min(as.Date(geomac_v$perimeter_date), na.rm=TRUE),
                                       max(as.Date(geomac_v$perimeter_date), na.rm=TRUE)))
  blk("area range (ha)",      sprintf("%.0f  to  %.0f",
                                       min(geomac_v$AREA_HA, na.rm=TRUE),
                                       max(geomac_v$AREA_HA, na.rm=TRUE)))
} else {
  rpt_add("  geomac_v not found — run Section 8")
}

# ---- [8] Pairs ---------------------------------------------------------------
hr("-")
rpt_add("[8] PERIMETER PAIRS")
hr("-")
if (exists("all_pairs")) {
  blk("total pairs",         nrow(all_pairs))
  blk("spread-prob pairs",   sum(all_pairs$use_for_spread_prob, na.rm=TRUE))
  blk("daily-area pairs",    sum(all_pairs$use_for_daily_area,  na.rm=TRUE))
  blk("dropped (shrinkage)", sum(all_pairs$drop_neg, na.rm=TRUE))
  blk("FWI join rate",       sprintf("%.0f%%  (%d of %d)",
                                      100 * mean(is.finite(all_pairs$FWI)),
                                      sum(is.finite(all_pairs$FWI)), nrow(all_pairs)))
  blk("wind join rate",      sprintf("%.0f%%  (%d of %d)",
                                      100 * mean(is.finite(all_pairs$WindSpeed_kmh)),
                                      sum(is.finite(all_pairs$WindSpeed_kmh)), nrow(all_pairs)))
  if (exists("spread_pairs"))
    blk("eligible (finite FWI+EWS)", sum(is.finite(spread_pairs$FWI) &
                                          is.finite(spread_pairs[[if ("EffectiveWind" %in% names(spread_pairs)) "EffectiveWind" else "WindSpeed_kmh"]]),
                                         na.rm=TRUE))
  # Gap distribution
  gap_tbl <- table(all_pairs$gap_days[!is.na(all_pairs$gap_days)])
  blk("gap_days distribution",
      paste(names(gap_tbl), gap_tbl, sep="=", collapse="  "))
  # delta_ha distribution (non-shrinking)
  da <- all_pairs$delta_ha[!all_pairs$drop_neg & all_pairs$delta_ha > 0]
  if (length(da) > 0)
    blk("delta_ha q10/med/q90",
        sprintf("%.0f / %.0f / %.0f",
                quantile(da, .10), quantile(da, .50), quantile(da, .90)))
} else {
  rpt_add("  all_pairs not found — run Section 9")
}

# ---- [9] pair_diag table -----------------------------------------------------
hr("-")
rpt_add("[9] PAIR_DIAG TABLE  (from Section 12 loop — Robbins 2021 cell definition)")
rpt_add("    pair_idx     = row in eligible spread_pairs")
rpt_add("    n_new_burned = cells newly burned between t and t+1  (success before subsampling)")
rpt_add("    n_fail       = frontier ring cells that did NOT burn  (y=0 in GLM)")
rpt_add("    pct_prox     = % of new cells adjacent to old frontier (100=contiguous, low=spatial offset)")
rpt_add("    status       = ok | no_success | null_raster")
hr("-")
if (exists("pair_diag") && nrow(pair_diag) > 0) {
  blk("null rasterization",  sum(pair_diag$status == "null_raster", na.rm=TRUE))
  blk("no new burned cells", sum(pair_diag$status == "no_success",  na.rm=TRUE))
  blk("ok (usable)",         sum(pair_diag$status == "ok",          na.rm=TRUE))
  rpt_add("")
  # Fixed-width table header
  rpt_add(sprintf("  %-4s  %-24s  %-10s  %8s  %8s  %7s  %7s  %5s  %-10s",
                  "idx", "FIRE_ID", "date", "delta_ha", "n_new", "n_fail",
                  "prox%", "FWI", "status"))
  rpt_add(sprintf("  %s  %s  %s  %s  %s  %s  %s  %s  %s",
                  strrep("-",4), strrep("-",24), strrep("-",10),
                  strrep("-",8), strrep("-",8), strrep("-",7),
                  strrep("-",7), strrep("-",5), strrep("-",10)))
  for (k in seq_len(nrow(pair_diag))) {
    pd <- pair_diag[k, ]
    rpt_add(sprintf("  %-4d  %-24s  %-10s  %8.0f  %8d  %7d  %7s  %5.1f  %-10s",
                    pd$pair_idx,
                    substr(as.character(pd$FIRE_ID), 1, 24),
                    as.character(pd$date),
                    pd$delta_ha,
                    if ("n_new_burned" %in% names(pd)) pd$n_new_burned else NA_integer_,
                    if ("n_fail"       %in% names(pd)) pd$n_fail       else NA_integer_,
                    if ("pct_proximal" %in% names(pd) && is.finite(pd$pct_proximal))
                      sprintf("%.0f%%", pd$pct_proximal) else "NA",
                    pd$FWI,
                    pd$status))
  }
  # CSV block for easy parsing
  rpt_add("")
  rpt_add("  [CSV]")
  rpt_add(paste(c("pair_idx","FIRE_ID","date","delta_ha","n_new_burned","n_fail",
                  "pct_proximal","FWI","EWS","status"), collapse=","))
  for (k in seq_len(nrow(pair_diag))) {
    pd <- pair_diag[k, ]
    rpt_add(paste(c(pd$pair_idx, pd$FIRE_ID, as.character(pd$date),
                    round(pd$delta_ha,1),
                    if ("n_new_burned" %in% names(pd)) pd$n_new_burned else NA,
                    if ("n_fail"       %in% names(pd)) pd$n_fail       else NA,
                    if ("pct_proximal" %in% names(pd)) round(pd$pct_proximal,1) else NA,
                    round(pd$FWI,2),
                    if ("EWS" %in% names(pd)) round(pd$EWS,2) else NA,
                    pd$status), collapse=","))
  }
} else {
  rpt_add("  pair_diag not found — run Section 12")
}

# ---- [10] GLM sample summary -------------------------------------------------
hr("-")
rpt_add("[10] GLM SAMPLE SUMMARY")
hr("-")
if (exists("all_samples") && nrow(all_samples) > 0) {
  blk("n rows total",     nrow(all_samples))
  blk("n successes (y=1)", sum(all_samples$y == 1))
  blk("n failures  (y=0)", sum(all_samples$y == 0))
  blk("n unique fires",   length(unique(all_samples$FIRE_ID)))
  blk("n unique dates",   length(unique(all_samples$date)))
  blk("FWI range",        sprintf("%.1f  to  %.1f",
                                   min(all_samples$FWI, na.rm=TRUE),
                                   max(all_samples$FWI, na.rm=TRUE)))
  blk("EWS range",        sprintf("%.1f  to  %.1f",
                                   min(all_samples$EffectiveWind, na.rm=TRUE),
                                   max(all_samples$EffectiveWind, na.rm=TRUE)))
  if ("FineFuels" %in% names(all_samples)) {
    blk("FF range",       sprintf("%.4f  to  %.4f",
                                   min(all_samples$FineFuels, na.rm=TRUE),
                                   max(all_samples$FineFuels, na.rm=TRUE)))
    blk("FF variance",    .fmt(var(all_samples$FineFuels, na.rm=TRUE), 6))
  }
  # Pair-level signal (Robbins approach: use log(n_new_burned) as spread magnitude)
  if (exists("pair_diag") && "n_new_burned" %in% names(pair_diag)) {
    pd_ok_r <- pair_diag[pair_diag$status == "ok" & pair_diag$n_new_burned > 0, ]
    if (nrow(pd_ok_r) > 1) {
      ew_col <- if ("EWS" %in% names(pd_ok_r)) pd_ok_r$EWS else pd_ok_r$WindSpeed_kmh
      blk("pair-level Cor(FWI, log(n_new))",
          .fmt(cor(pd_ok_r$FWI, log(pd_ok_r$n_new_burned), use="complete.obs"), 3))
      blk("pair-level Cor(EWS, log(n_new))",
          .fmt(cor(ew_col, log(pd_ok_r$n_new_burned), use="complete.obs"), 3))
      blk("median pct_proximal (contiguity)",
          sprintf("%.1f%%", median(pd_ok_r$pct_proximal, na.rm=TRUE)))
    }
  }
} else {
  rpt_add("  all_samples not found — run Section 12")
}

# ---- [11] Collinearity -------------------------------------------------------
hr("-")
rpt_add("[11] COLLINEARITY")
hr("-")
if (exists("all_samples") && nrow(all_samples) > 0) {
  X_rpt <- as.matrix(all_samples[, c("FWI","EffectiveWind")])
  X_rpt <- X_rpt[complete.cases(X_rpt),]
  if (nrow(X_rpt) > 2) {
    vif_rpt <- tryCatch(diag(solve(cor(X_rpt))), error = function(e) c(NA,NA))
    sv_rpt  <- svd(cbind(1, X_rpt))$d
    kappa_rpt <- max(sv_rpt) / min(sv_rpt[sv_rpt > .Machine$double.eps * max(sv_rpt)])
    r_rpt <- cor(X_rpt[,"FWI"], X_rpt[,"EffectiveWind"])
    blk("Cor(FWI, EWS) — cell level", .fmt(r_rpt, 3))
    blk("VIF FWI",      .fmt(vif_rpt[1], 2))
    blk("VIF EWS",      .fmt(vif_rpt[2], 2))
    blk("condition number", .fmt(kappa_rpt, 1))
    blk(">10 VIF or >100 condition",
        if (any(vif_rpt > 10, na.rm=TRUE) || kappa_rpt > 100) "YES — multicollinearity" else "no")
  }
}

# ---- [12] GLM COEFFICIENTS ---------------------------------------------------
hr("-")
rpt_add("[12] GLM COEFFICIENTS")
hr("-")
if (exists("m_spread") && !is.null(m_spread)) {
  blk("primary model",    deparse(formula(m_spread)))
  blk("n observations",   nobs(m_spread))
  blk("McFadden R2",      .fmt(1 - m_spread$deviance / m_spread$null.deviance, 4))
  blk("AUC",              tryCatch(.fmt(auc_est, 4), error = function(e) "not computed"))
  blk("Pearson phi",      tryCatch(.fmt(phi, 3), error = function(e) "not computed"))

  # Model comparison table (all candidate models from Section 14)
  if (exists("comp_tbl") && nrow(comp_tbl) > 0) {
    rpt_add("")
    rpt_add("  Model comparison:")
    rpt_add(sprintf("  %-20s  %8s  %6s  %10s  %10s  %10s  %8s  %6s  %6s",
                    "Model", "McF_R2", "AUC", "B1_FWI", "B2_FF", "B3_EWS",
                    "FWI", "FF", "EWS"))
    rpt_add(sprintf("  %s  %s  %s  %s  %s  %s  %s  %s  %s",
                    strrep("-",20), strrep("-",8), strrep("-",6),
                    strrep("-",10), strrep("-",10), strrep("-",10),
                    strrep("-",8), strrep("-",6), strrep("-",6)))
    for (k in seq_len(nrow(comp_tbl))) {
      ct <- comp_tbl[k,]
      rpt_add(sprintf("  %-20s  %8.4f  %6.4f  %10s  %10s  %10s  %8s  %6s  %6s",
                      ct$Model, ct$McF_R2, ct$AUC,
                      ifelse(is.na(ct$B1_FWI), "—", sprintf("%+.5f", ct$B1_FWI)),
                      ifelse(is.na(ct$B2_FF),  "—", sprintf("%+.4f", ct$B2_FF)),
                      ifelse(is.na(ct$B3_EWS), "—", sprintf("%+.5f", ct$B3_EWS)),
                      ct$FWI_sign, ct$FF_sign, ct$EWS_sign))
    }
  }

  # Primary model full coefficient table
  rpt_add("")
  rpt_add("  Primary model full coefficients:")
  rpt_add(sprintf("  %-20s  %10s  %10s  %8s  %8s  %-18s",
                  "term", "estimate", "std_error", "z_stat", "p_value", "flag"))
  rpt_add(sprintf("  %s  %s  %s  %s  %s  %s",
                  strrep("-",20), strrep("-",10), strrep("-",10),
                  strrep("-",8), strrep("-",8), strrep("-",18)))
  if (exists("coef_sum")) {
    for (k in seq_len(nrow(coef_sum))) {
      cs <- coef_sum[k,]
      rpt_add(sprintf("  %-20s  %+10.5f  %10.5f  %8.3f  %8.5f  %-18s",
                      cs$term, cs$estimate, cs$std_error,
                      cs$estimate / max(cs$std_error, 1e-12),
                      cs$p_value, cs$flag))
    }
  }

  # Spot-check table
  rpt_add("")
  rpt_add("  P(spread) spot-check at mean EWS and mean FF:")
  if (exists("b0") && exists("b1") && exists("b2") && exists("b3") &&
      exists("ews_mean") && exists("ff_mean_val")) {
    rpt_add(sprintf("  %-10s  %s", "FWI", "P(spread)"))
    for (fwi_v in c(10, 20, 30, 40, 60, 80, 90)) {
      pv <- plogis(b0 + b1*fwi_v + b2*ff_mean_val + b3*ews_mean)
      rpt_add(sprintf("  FWI = %3d   %.1f%%", fwi_v, 100*pv))
    }
  }
} else {
  rpt_add("  GLM not fitted — run Section 14")
}

# ---- [13] Issues flagged in Section 17 --------------------------------------
hr("-")
rpt_add("[13] ISSUES FLAGGED (from Section 17)")
hr("-")
if (exists("n_issues")) {
  blk("n issues", n_issues)
  if (exists("coef_sum")) {
    bad_sign <- coef_sum$term[!coef_sum$sign_ok & !is.na(coef_sum$expected_sign)]
    if (length(bad_sign) > 0)
      rpt_add(sprintf("  Wrong-sign coef(s): %s", paste(bad_sign, collapse=", ")))
  }
  if (exists("auc_est") && !is.na(auc_est) && auc_est < 0.60)
    rpt_add(sprintf("  Low AUC: %.3f", auc_est))
  if (exists("phi") && phi > 2.0)
    rpt_add(sprintf("  Overdispersion phi=%.2f", phi))
  if (exists("r_fwi_ew") && !is.na(r_fwi_ew) && abs(r_fwi_ew) > 0.8)
    rpt_add(sprintf("  FWI/EWS collinearity r=%.3f", r_fwi_ew))
  if (exists("mcfadden_r2") && mcfadden_r2 < 0.05)
    rpt_add(sprintf("  Poor fit McFadden R2=%.4f", mcfadden_r2))
  if (exists("fwi_join_rate") && fwi_join_rate < 0.5)
    rpt_add(sprintf("  Low FWI join rate: %.0f%%", 100*fwi_join_rate))
  if (n_issues == 0)
    rpt_add("  None")
} else {
  rpt_add("  Section 16 not run")
}

# ---- [14] MaxSpreadArea -------------------------------------------------------
hr("-")
rpt_add("[14] MAXSPREADAREA (FIRED Events — Section 16)")
hr("-")
if (exists("fired_max") && !is.null(fired_max) && nrow(fired_max) > 0) {
  blk("n events (after filters)", nrow(fired_max))
  blk("daily_area_ha range (ha)", sprintf("%.0f  to  %.0f",
      min(fired_max$daily_area_ha, na.rm=TRUE),
      max(fired_max$daily_area_ha, na.rm=TRUE)))
  q75_area_rpt <- quantile(fired_max$daily_area_ha, 0.75, na.rm=TRUE)
  blk("Q50 / Q75 max daily area (ha)", sprintf("%.0f / %.0f",
      median(fired_max$daily_area_ha, na.rm=TRUE), q75_area_rpt))
  blk("FWI range on max-growth days", sprintf("%.1f  to  %.1f",
      min(fired_max$FWI, na.rm=TRUE), max(fired_max$FWI, na.rm=TRUE)))
  blk("EWS range (km/h)", sprintf("%.1f  to  %.1f",
      min(fired_max$EWS, na.rm=TRUE), max(fired_max$EWS, na.rm=TRUE)))
  blk("Cor(FWI, max_area)",
      .fmt(cor(fired_max$FWI, fired_max$daily_area_ha, use="complete.obs"), 3))
  blk("Cor(EWS, max_area)",
      .fmt(cor(fired_max$EWS, fired_max$daily_area_ha, use="complete.obs"), 3))
  blk("Cor(FWI, EWS)",
      .fmt(cor(fired_max$FWI, fired_max$EWS, use="complete.obs"), 3))
  blk("EWS source", tryCatch(ews_src_max, error = function(e) "unknown"))

  if (exists("maxarea_coef") && !is.null(maxarea_coef)) {
    rpt_add("")
    rpt_add(sprintf("  Model: %s", maxarea_coef$method))
    rpt_add(sprintf("  MaximumSpreadAreaB0 = %+.2f   (intercept; must be > 0)", maxarea_coef$B0))
    rpt_add(sprintf("  MaximumSpreadAreaB1 = %+.4f   (FWI; should be > 0)", maxarea_coef$B1))
    rpt_add(sprintf("  MaximumSpreadAreaB2 = %+.4f   (EWS; should be > 0)", maxarea_coef$B2))

    # Sign / validity flags
    if (maxarea_coef$B0 <= 0)
      rpt_add("  ⚠ B0 <= 0 — clamp to > 0 before using in SCF (set to 0.1 minimum)")
    if (maxarea_coef$B1 < 0)
      rpt_add("  ⚠ B1 < 0 — negative FWI coefficient (wrong sign; try CAL_MONTHS=6:9)")
    if (maxarea_coef$B2 < 0)
      rpt_add("  ⚠ B2 < 0 — negative EWS coefficient (check wind data)")

    # Coverage
    pred_rpt <- maxarea_coef$B0 + maxarea_coef$B1 * fired_max$FWI +
                maxarea_coef$B2 * fired_max$EWS
    pct_rpt  <- 100 * mean(fired_max$daily_area_ha < pred_rpt, na.rm=TRUE)
    blk("% events below predicted line", sprintf("%.1f%%  (target: 75%%)", pct_rpt))

    # OLS R² if available
    if (exists("fit_ols_max") && !is.null(fit_ols_max))
      blk("OLS R² (reference)", .fmt(summary(fit_ols_max)$r.squared, 4))

    # Spot-check
    rpt_add("")
    ews_rpt_max <- mean(fired_max$EWS, na.rm=TRUE)
    rpt_add(sprintf("  Spot-check at mean EWS = %.1f km/h:", ews_rpt_max))
    rpt_add(sprintf("  %-8s  %s", "FWI", "MaxSpreadArea (ha/day)"))
    for (fwi_v in c(5, 10, 20, 30, 50, 70, 90)) {
      ha_rpt <- maxarea_coef$B0 + maxarea_coef$B1 * fwi_v + maxarea_coef$B2 * ews_rpt_max
      rpt_add(sprintf("  FWI = %3d   %6.0f ha%s",
                      fwi_v, ha_rpt,
                      if (ha_rpt <= 0) "  ⚠ NEGATIVE" else ""))
    }
  } else {
    rpt_add("  MaxSpreadArea model not fitted (see Section 16 messages).")
  }
} else {
  rpt_add("  fired_max not built — run Section 16 with FIRED_EVENTS_GPKG set.")
}

# ---- Footer ------------------------------------------------------------------
hr()
rpt_add("END OF REPORT")
hr()

# ---- Print to console --------------------------------------------------------
cat(paste(rpt, collapse = "\n"), "\n")

# ---- Write to file -----------------------------------------------------------
report_path <- file.path(
  tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
           error = function(e) getwd()),
  sprintf("spread_diag_report_%s.txt",
          format(Sys.time(), "%Y%m%d_%H%M%S"))
)
writeLines(rpt, report_path)
message(sprintf("\nReport written to: %s", report_path))
message("Paste the full contents of that file into Claude for interpretation.")
