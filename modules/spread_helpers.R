# =============================================================================
# spread_helpers.R
# Functions for SCF spread parameter calibration.
#
# Key change from v1:
#   load_geomac_perimeters() now accepts cal_vect (already in working_crs)
#   and working_crs separately. GeoMAC perimeters are projected to working_crs
#   before the spatial intersect filter, consistent with ignition_helpers.R.
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, ggplot2, tidyr, purrr, stringr,
               lubridate, glue, broom, scales, httr, elevatr)

ROOK_W <- matrix(c(0,1,0, 1,1,1, 0,1,0), nrow = 3, byrow = TRUE)

# Reuse utilities from ignition_helpers (sourced before this file in app.R)
# pick_first_col(), parse_any_date() are defined there.


# =============================================================================
# Fine fuels helpers
# =============================================================================

# -----------------------------------------------------------------------------
# FBFM40 -> fine fuel load lookup table
#
# Fine fuel load = 1-hr + 10-hr dead fuel (tons/acre) for each Scott & Burgan
# 40 fuel model, sourced from Scott & Burgan (2005).  LANDFIRE stores FBFM40
# as integer codes; this table maps those codes to physical fuel loads which
# are then normalized 0-1 for use in the SCF spread probability equation.
# -----------------------------------------------------------------------------

FBFM40_FINE_FUELS <- c(
  # Non-burnable
  "91" = 0.00, "92" = 0.00, "93" = 0.05,
  "98" = 0.00, "99" = 0.00,
  # Grass (GR)
  "101" = 0.10, "102" = 0.10, "103" = 0.50,
  "104" = 2.25, "105" = 3.45, "106" = 3.60,
  "107" = 5.50, "108" = 4.50, "109" = 9.00,
  # Grass-Shrub (GS)
  "121" = 0.35, "122" = 1.55, "123" = 1.80, "124" = 6.40,
  # Shrub (SH)
  "141" = 0.45, "142" = 1.35, "143" = 0.50, "144" = 1.80,
  "145" = 3.60, "146" = 2.90, "147" = 6.20, "149" = 8.35,
  # Timber-Understory (TU)
  "161" = 0.80, "162" = 1.60, "163" = 1.50,
  "164" = 0.70, "165" = 3.00,
  # Timber Litter (TL)
  "181" = 1.00, "182" = 1.40, "183" = 0.80, "184" = 0.80,
  "185" = 1.40, "186" = 2.40, "187" = 1.40,
  "188" = 3.50, "189" = 4.50,
  # Slash-Blowdown (SB)
  "201" = 1.55, "202" = 4.50, "203" = 5.50, "204" = 3.50
)


# -----------------------------------------------------------------------------
# fetch_landfire_fine_fuels
#
# Pulls LANDFIRE FBFM40 for the calibration boundary via FedData, reclassifies
# integer fuel model codes to a 0-1 fine fuel load index using FBFM40_FINE_FUELS,
# reprojects and resamples to the spread rasterization grid.
#
# Arguments:
#   cal_vect       SpatVector of calibration boundary (any CRS)
#   working_crs    CRS string for output
#   spread_template SpatRaster defining the output grid (extent + resolution)
#   progress_fn    optional function(msg) for status messages
#
# Returns a SpatRaster on spread_template grid, values 0-1 (NA = non-burnable)
# -----------------------------------------------------------------------------

fetch_landfire_fine_fuels <- function(cal_vect, working_crs, spread_template,
                                       progress_fn = NULL) {

  # ---- Bounding box in WGS84 -----------------------------------------------
  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  bb        <- as.vector(ext(cal_wgs84))   # xmin, xmax, ymin, ymax

  # ---- Request pixel dimensions (cap at 4096; LANDFIRE native = 30 m) -------
  deg_per_30m <- 30 / 111320          # ~0.000270 degrees per 30 m
  nx <- min(as.integer(ceiling((bb["xmax"] - bb["xmin"]) / deg_per_30m)), 4096L)
  ny <- min(as.integer(ceiling((bb["ymax"] - bb["ymin"]) / deg_per_30m)), 4096L)
  message(sprintf("Requesting LANDFIRE FBFM40: %d x %d pixels over bbox [%.4f %.4f %.4f %.4f]",
                  nx, ny, bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]))

  # ---- LANDFIRE LF2024 FBFM40 CONUS ImageServer ----------------------------
  # Hosted at lfps.usgs.gov; native CRS EPSG:5070, pixel type S16, range 91-204.
  # rasterFunction "None" returns raw integer fuel model codes.
  base_url <- paste0(
    "https://lfps.usgs.gov/arcgis/rest/services/",
    "Landfire_LF2024/LF2024_FBFM40_CONUS/ImageServer/exportImage"
  )

  if (!is.null(progress_fn)) progress_fn("Downloading LANDFIRE FBFM40 from USGS server...")

  tmp_tif <- tempfile(fileext = ".tif")

  resp <- tryCatch(
    httr::GET(
      base_url,
      query = list(
        bbox          = sprintf("%.6f,%.6f,%.6f,%.6f",
                                bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]),
        bboxSR        = "4326",
        imageSR       = "4326",
        size          = sprintf("%d,%d", nx, ny),
        format        = "tiff",
        pixelType     = "S16",
        renderingRule = '{"rasterFunction":"None"}',
        f             = "image"
      ),
      httr::write_disk(tmp_tif, overwrite = TRUE),
      httr::timeout(300)
    ),
    error = function(e)
      stop("LANDFIRE HTTP request failed: ", conditionMessage(e))
  )

  if (httr::http_error(resp))
    stop("LANDFIRE REST API returned HTTP ", httr::status_code(resp),
         ". Check internet connection or try again later.")

  # ---- Read and validate ----------------------------------------------------
  fbfm_r <- tryCatch(
    rast(tmp_tif),
    error = function(e)
      stop("Could not read LANDFIRE response as raster: ", conditionMessage(e))
  )

  raw_vals  <- values(fbfm_r, mat = FALSE)
  val_range <- range(raw_vals, na.rm = TRUE)
  message(sprintf("LANDFIRE raw pixel range: [%.0f, %.0f]  (expected 91-204)",
                  val_range[1], val_range[2]))

  # ---- Reclassify integer codes -> fine fuel load (tons/acre), normalise ----
  if (!is.null(progress_fn)) progress_fn("Reclassifying FBFM40 to fine fuel index...")

  code_strs <- as.character(as.integer(round(raw_vals)))
  loads     <- as.numeric(FBFM40_FINE_FUELS[code_strs])   # NA for unknown codes

  max_load <- max(FBFM40_FINE_FUELS)                       # GR9 = 9.0 t/ac
  loads    <- loads / max_load
  loads[loads < 0 | !is.finite(loads)] <- NA

  values(fbfm_r) <- loads

  matched_frac <- mean(!is.na(loads))
  message(sprintf("FBFM40 reclassification: %.1f%% of pixels matched a known fuel code",
                  matched_frac * 100))
  if (matched_frac < 0.1)
    warning("Fewer than 10% of pixels matched known FBFM40 codes. ",
            "The downloaded raster may not contain raw fuel model integers. ",
            "Fine fuels will fall back to placeholder (1.0) for most cells.")

  # ---- Reproject + resample to spread grid ---------------------------------
  if (!is.null(progress_fn)) progress_fn("Reprojecting and resampling fine fuels...")
  fbfm_proj <- project(fbfm_r, working_crs, method = "bilinear")
  resample(fbfm_proj, spread_template, method = "bilinear")
}


# -----------------------------------------------------------------------------
# load_local_fine_fuels
#
# Loads a user-provided raster, validates it is 0-1 scaled, reprojects and
# resamples to the spread grid.
# -----------------------------------------------------------------------------

load_local_fine_fuels <- function(path, working_crs, spread_template) {

  if (!file.exists(path))
    stop("Fine fuels raster not found: ", path)

  r <- rast(path)

  rng <- range(values(r, mat = FALSE), na.rm = TRUE)
  if (rng[2] > 10)
    warning(sprintf(
      "Fine fuels raster max value is %.1f — expected 0-1 scale. ",
      rng[2]),
      "Values will be normalized to 0-1.")

  # Normalize if not already 0-1
  v <- values(r, mat = FALSE)
  mx <- max(v, na.rm = TRUE)
  if (is.finite(mx) && mx > 1) v <- v / mx
  v[!is.finite(v)] <- NA
  values(r) <- v

  r_proj <- project(r, working_crs, method = "bilinear")
  resample(r_proj, spread_template, method = "bilinear")
}


# -----------------------------------------------------------------------------
# load_geomac_perimeters
#
# Loads GeoMAC historic perimeter geodatabase, projects to working_crs,
# filters to perimeters intersecting the calibration boundary, then pulls
# all perimeters for those fire IDs within the calibration date window.
#
# Arguments:
#   gdb_path    path to Historic_Geomac_*.gdb
#   cal_vect    calibration boundary SpatVector, already in working_crs
#   template_r  LANDIS template SpatRaster (defines rasterization grid)
#   working_crs CRS string from template_r
#   cal_start/end Date range
#
# Returns SpatVector with standardized fields:
#   FIRE_ID, perimeter_date, AREA_HA, row_id
# -----------------------------------------------------------------------------

load_geomac_perimeters <- function(gdb_path, cal_vect, template_r,
                                    working_crs, cal_start, cal_end) {

  if (!file.exists(gdb_path)) stop("GeoMAC GDB not found: ", gdb_path)

  # ---- Detect layer name --------------------------------------------------
  lyrs <- tryCatch(terra::vector_layers(gdb_path)$layer, error = function(e) NULL)
  lyr  <- if (!is.null(lyrs)) {
    hit <- lyrs[grepl("perim|fire", tolower(lyrs))]
    if (length(hit) > 0) hit[1] else lyrs[1]
  } else "US_HIST_FIRE_PERIMTRS_2000_2018_DD83"

  message("Loading GeoMAC layer: ", lyr)
  gm_v <- vect(gdb_path, layer = lyr)

  # ---- Reproject to working CRS BEFORE spatial operations ----------------
  if (is.na(crs(gm_v)) || nchar(crs(gm_v)) == 0)
    stop("GeoMAC layer has no CRS defined.")

  message("Reprojecting GeoMAC to working CRS...")
  gm_v <- project(gm_v, working_crs)

  # ---- Column detection ---------------------------------------------------
  nm       <- names(gm_v)
  id_col   <- pick_first_col(nm, c("uniquefireidentifier", "unique.*id", "irwinid", "fire.*id"))
  year_col <- pick_first_col(nm, c("^fireyear$", "fireyear", "fire_year"))
  name_col <- pick_first_col(nm, c("incidentname", "incident.*name", "fire_name"))
  dt_col   <- pick_first_col(nm, c("perimeterdatetime", "datecurrent",
                                    "perimeterdate", "datetime", "date"))
  area_col <- pick_first_col(nm, c("gisacres", "acres", "area_ha", "area"))

  if (is.na(dt_col)) stop("Could not find perimeter date field in GeoMAC layer.")

  # ---- Spatial filter: perimeters intersecting calibration boundary --------
  message("Filtering GeoMAC to calibration boundary...")
  intersects_cal <- relate(gm_v, cal_vect, "intersects")
  gm_cal <- gm_v[intersects_cal, ]

  message(sprintf("  -> %d perimeters intersect calibration boundary.", nrow(gm_cal)))

  if (nrow(gm_cal) == 0) {
    stop(paste0(
      "No GeoMAC perimeters intersect the calibration boundary.\n",
      "  GeoMAC extent: ", as.character(ext(gm_v)), "\n",
      "  Calibration boundary extent: ", as.character(ext(cal_vect))
    ))
  }

  # ---- Pull fire IDs from intersecting perimeters -------------------------
  # Goal: retrieve ALL perimeters for each fire that touches the boundary
  # (not just the boundary-intersecting subset) so pair-building gets the
  # full progression. Fall back to gm_cal if the ID column is absent or
  # all-NA for the intersecting set.
  all_fires <- gm_cal   # safe default

  if (!is.na(id_col)) {
    fire_ids <- unique(gm_cal[[id_col]])
    fire_ids <- fire_ids[!is.na(fire_ids) & nchar(as.character(fire_ids)) > 0]

    if (length(fire_ids) > 0) {
      candidate <- gm_v[as.character(gm_v[[id_col]]) %in% as.character(fire_ids), ]
      if (nrow(candidate) > 0) all_fires <- candidate
    }
  }

  message(sprintf("  -> %d total perimeters for %d unique fires passed to date filter.",
                  nrow(all_fires),
                  if (!is.na(id_col) && nrow(all_fires) > 0)
                    length(unique(all_fires[[id_col]]))
                  else nrow(all_fires)))

  # ---- Parse attributes + date filter -------------------------------------
  # Show a sample of raw date values to aid diagnosis if parsing fails
  raw_dates <- as.data.frame(all_fires)[[dt_col]]
  message(sprintf("GeoMAC date column '%s' — sample raw values: %s",
                  dt_col,
                  paste(head(unique(as.character(raw_dates)), 5), collapse = " | ")))

  parse_geomac_date <- function(x) {
    # Already a Date or POSIXct
    if (inherits(x, "Date"))                  return(as.Date(x))
    if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
    # Character / factor: try common formats in order
    x <- as.character(x)
    x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
    d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
    if (!all(is.na(d))) return(d)
    d <- suppressWarnings(as.Date(as.POSIXct(x, tz = "UTC")))
    if (!all(is.na(d))) return(d)
    d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))
    if (!all(is.na(d))) return(d)
    d <- suppressWarnings(as.Date(x, "%Y/%m/%d"))
    if (!all(is.na(d))) return(d)
    d <- suppressWarnings(as.Date(x, "%d-%b-%Y"))
    if (!all(is.na(d))) return(d)
    d  # return last attempt (all NA) so error fires below
  }

  full_df <- as.data.frame(all_fires) %>%
    mutate(
      perimeter_date = parse_geomac_date(.data[[dt_col]]),
      FIRE_YEAR      = if (!is.na(year_col)) as.integer(.data[[year_col]])
                       else year(perimeter_date),
      INCIDENT       = if (!is.na(name_col)) as.character(.data[[name_col]])
                       else NA_character_,
      FIRE_ID        = if (!is.na(id_col)) as.character(.data[[id_col]])
                       else paste0(INCIDENT, "_", FIRE_YEAR),
      GIS_ACRES      = if (!is.na(area_col))
                         suppressWarnings(as.numeric(.data[[area_col]])) else NA_real_
    ) %>%
    filter(!is.na(perimeter_date),
           perimeter_date >= cal_start,
           perimeter_date <= cal_end)

  if (nrow(full_df) == 0) {
    n_parsed <- sum(!is.na(parse_geomac_date(raw_dates)))
    stop(sprintf(
      "No GeoMAC perimeters fall within the calibration date range (%s to %s).\n",
      cal_start, cal_end),
      sprintf("  %d of %d perimeters had parseable dates.\n", n_parsed, length(raw_dates)),
      sprintf("  Sample raw values: %s",
              paste(head(unique(as.character(raw_dates)), 5), collapse = " | "))
    )
  }

  # ---- Area ---------------------------------------------------------------
  if (all(is.na(full_df$GIS_ACRES))) {
    area_ha <- expanse(all_fires, unit = "ha")
  } else {
    area_ha <- full_df$GIS_ACRES * 0.40468564224
  }

  # ---- Attach standardized fields to vector -------------------------------
  all_fires$FIRE_ID        <- full_df$FIRE_ID
  all_fires$perimeter_date <- full_df$perimeter_date
  all_fires$AREA_HA        <- area_ha
  all_fires$row_id         <- seq_len(nrow(all_fires))

  # Filter to date window (geometry subset must match full_df rows)
  date_keep <- !is.na(as.Date(all_fires$perimeter_date)) &
               as.Date(all_fires$perimeter_date) >= cal_start &
               as.Date(all_fires$perimeter_date) <= cal_end
  all_fires <- all_fires[date_keep, ]
  all_fires$row_id <- seq_len(nrow(all_fires))  # re-index after filter

  message(sprintf("GeoMAC: %d perimeters for %d fires touching calibration boundary.",
                  nrow(all_fires), length(unique(all_fires$FIRE_ID))))

  all_fires
}


# -----------------------------------------------------------------------------
# fetch_terrain
#
# Downloads a DEM for the calibration boundary via elevatr (SRTM/3DEP tiles),
# then computes slope (degrees) and aspect (degrees clockwise from N) with
# terra::terrain().
#
# Arguments:
#   cal_vect    SpatVector in working CRS defining the calibration area
#   working_crs CRS string (taken from template_r)
#   z           elevatr zoom level (default 8 ≈ 600 m; 9 ≈ 300 m)
#   progress_fn optional function(msg) for status messages
#
# Returns list(dem, slope, aspect) — all SpatRasters in working_crs
# -----------------------------------------------------------------------------

fetch_terrain <- function(cal_vect, working_crs, z = 8, progress_fn = NULL) {

  if (!requireNamespace("elevatr", quietly = TRUE))
    stop("elevatr package required. Install with: install.packages('elevatr')")

  if (!is.null(progress_fn)) progress_fn("Fetching DEM via elevatr (SRTM/3DEP)...")

  cal_sf <- sf::st_as_sf(project(cal_vect, "EPSG:4326"))

  dem_raster <- tryCatch(
    elevatr::get_elev_raster(cal_sf, z = z, clip = "bbox", neg_to_na = FALSE),
    error = function(e) stop("elevatr DEM fetch failed: ", conditionMessage(e))
  )

  dem_r    <- rast(dem_raster)
  dem_proj <- project(dem_r, working_crs, method = "bilinear")

  if (!is.null(progress_fn)) progress_fn("Computing slope and aspect...")

  slope_r  <- terrain(dem_proj, v = "slope",  unit = "degrees")
  aspect_r <- terrain(dem_proj, v = "aspect", unit = "degrees")

  message(sprintf(
    "Terrain: DEM range [%.0f, %.0f] m  |  slope [%.1f, %.1f] deg  |  %d x %d cells",
    minmax(dem_proj)[1], minmax(dem_proj)[2],
    minmax(slope_r)[1],  minmax(slope_r)[2],
    nrow(dem_proj), ncol(dem_proj)
  ))

  list(dem = dem_proj, slope = slope_r, aspect = aspect_r)
}


# -----------------------------------------------------------------------------
# compute_effective_wind
#
# Converts raw wind speed + direction + terrain (slope/aspect) into an
# "effective wind speed" using the Nelson (2002) vector-addition formula —
# EXACTLY the formula used by SCF's CalculateEffectiveWindSpeed().
#
# Reference: Nelson, R.M. (2002). An evaluation of the fire behavior model
#   BURNUP. International Journal of Wildland Fire, 11(1), 17-28. (Eq. 5)
#
# Formula:
#   UaUb  = wind_speed / combustion_buoyancy
#   EWS   = Cb × √(UaUb² + 2·UaUb·sin(θ)·cos(φ) + sin²(θ))
#
#   Expanding: EWS = √(wind_speed² + 2·wind_speed·Cb·sin(θ)·cos(φ) + Cb²·sin²(θ))
#
#   The Cb² factor means slope contribution scales QUADRATICALLY with Cb.
#   On a 20° slope at 5 m/s wind:
#     Cb=10: EWS ≈  7 m/s   (low intensity)
#     Cb=25: EWS ≈ 14 m/s   (moderate intensity)
#     Cb=50: EWS ≈ 22 m/s   (high intensity — 3× the Cb=10 value!)
#
#   where:
#     Cb  = combustion_buoyancy:
#             10  = low intensity  (SCF default when SiteVars.Intensity <= 3)
#             25  = moderate       (SCF uses this when 3 < Intensity <= 6)
#             50  = high intensity (SCF uses this when Intensity > 6)
#     θ   = slope angle (radians)
#     φ   = angle between wind direction and uphill direction (radians)
#
# *** CALIBRATION NOTE: Use Cb=10 as the starting point (SCF default for     ***
# *** low-intensity fires). Only raise Cb if your landscape reliably burns at ***
# *** intensity index >3. Cb=25/50 amplifies slope contribution QUADRATICALLY***
# *** — on steep terrain this can dominate EWS and create counter-intuitive  ***
# *** negative correlations with FWI, causing backwards GLM signs.           ***
#
# UNITS: wind_speed must be in km/h to match the LANDIS Climate Library.
#   SCF User Guide (p.5, §2.43, §2.47, §2.51):
#   "The climate library converts all wind speed units into kilometers / hour.
#    Be sure to convert your wind data into the correct units when inputting
#    into the climate library."
#   Pass ERA5/GRIDMET wind speed in km/h (multiply m/s by 3.6 before calling).
#   EWS output will also be in km/h, matching what SCF receives at runtime.
#
# terra::terrain(aspect) returns the DOWNHILL direction (clockwise from N);
# uphill direction = (aspect + 180) %% 360, which matches SCF's
# UphillSlopeAzimuth site variable.
#
# Special cases:
#   slope = 0  →  EWS = |wind_speed|  (flat terrain, wind only)
#   wind  = 0  →  EWS = Cb × sin(θ)  (slope contribution only)
#
# All inputs are vectors of equal length. Returns vector of same length (km/h).
# -----------------------------------------------------------------------------

compute_effective_wind <- function(wind_speed, wind_dir, slope_deg, aspect_deg,
                                    combustion_buoyancy = 10.0) {

  # Replace NAs with neutral values so vector stays full-length
  slope_deg  <- ifelse(is.finite(slope_deg),  slope_deg,  0)
  aspect_deg <- ifelse(is.finite(aspect_deg), aspect_deg, wind_dir)
  wind_speed <- ifelse(is.finite(wind_speed), wind_speed, 0)
  wind_dir   <- ifelse(is.finite(wind_dir),   wind_dir,   0)

  # SCF uses UphillSlopeAzimuth; terra aspect is downhill → add 180
  upslope_dir <- (aspect_deg + 180) %% 360

  # Relative wind direction (angle between wind and uphill slope)
  # SCF: relativeWindDirection = (windDirection - slopeAngle) / 180 * PI
  rwd_rad <- ((wind_dir - upslope_dir) / 180.0) * pi

  # Slope in radians
  slope_rad <- (slope_deg / 180.0) * pi

  # Nelson (2002) Eq. 5 — matches SCF CalculateEffectiveWindSpeed exactly
  Ua_Ub <- wind_speed / combustion_buoyancy
  ews   <- combustion_buoyancy * sqrt(
    Ua_Ub^2 +
    2.0 * Ua_Ub * sin(slope_rad) * cos(rwd_rad) +
    sin(slope_rad)^2
  )

  pmax(0.0, ews)
}


# -----------------------------------------------------------------------------
# bind_climate_to_pairs
# -----------------------------------------------------------------------------

bind_climate_to_pairs <- function(geomac_v, fwi_daily, wind_daily,
                                   slope_r             = NULL,
                                   aspect_r            = NULL,
                                   max_gap_sp          = 1,
                                   max_gap_da          = 3,
                                   neg_tol_ha          = 1.0,
                                   combustion_buoyancy = 10.0) {

  gm_df <- as.data.frame(geomac_v) %>%
    transmute(
      row_id  = as.integer(row_id),
      FIRE_ID = as.character(FIRE_ID),
      date    = as.Date(perimeter_date),
      area_ha = suppressWarnings(as.numeric(AREA_HA))
    ) %>%
    filter(!is.na(FIRE_ID), !is.na(date), is.finite(area_ha), area_ha > 0)

  # Keep largest perimeter per fire x date
  gm_clean <- gm_df %>%
    group_by(FIRE_ID, date) %>%
    slice_max(order_by = area_ha, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Build sequential pairs
  pairs <- gm_clean %>%
    arrange(FIRE_ID, date) %>%
    group_by(FIRE_ID) %>%
    mutate(
      row_id_next = lead(row_id),
      date_next   = lead(date),
      area_next   = lead(area_ha),
      gap_days    = as.integer(date_next - date),
      delta_ha    = area_next - area_ha
    ) %>%
    ungroup() %>%
    filter(!is.na(row_id_next), !is.na(date_next), gap_days >= 1) %>%
    mutate(
      daily_area_ha       = delta_ha / gap_days,
      drop_neg            = delta_ha < (-neg_tol_ha),
      use_for_spread_prob = !drop_neg & gap_days <= max_gap_sp & delta_ha > 0,
      use_for_daily_area  = !drop_neg & gap_days <= max_gap_da & delta_ha > 0
    )

  # Extract centroid coordinates for each pair's t perimeter (for terrain lookup)
  gm_row_ids  <- as.integer(geomac_v$row_id)
  pair_idx    <- match(pairs$row_id, gm_row_ids)
  valid_idx   <- !is.na(pair_idx)

  ctr_v       <- centroids(geomac_v[pair_idx[valid_idx], ])
  ctr_xy      <- geom(ctr_v)[, c("x", "y"), drop = FALSE]

  pairs$centroid_x <- NA_real_
  pairs$centroid_y <- NA_real_
  pairs$centroid_x[valid_idx] <- ctr_xy[, "x"]
  pairs$centroid_y[valid_idx] <- ctr_xy[, "y"]

  # Join climate
  climate <- fwi_daily %>%
    rename(FWI = value) %>%
    left_join(wind_daily, by = "date")

  pairs <- pairs %>%
    left_join(climate, by = "date") %>%
    mutate(
      FWI           = suppressWarnings(as.numeric(FWI)),
      WindSpeed_kmh = suppressWarnings(as.numeric(WindSpeed_kmh)),
      WindDir_deg   = suppressWarnings(as.numeric(WindDir_deg))
    )

  # Compute effective wind: topographic formula if terrain available, else raw speed
  if (!is.null(slope_r) && !is.null(aspect_r)) {
    message(sprintf(
      "Computing effective wind using slope/aspect at fire centroids (Cb=%.0f)...",
      combustion_buoyancy
    ))
    valid_pts <- !is.na(pairs$centroid_x) & !is.na(pairs$centroid_y) &
                 is.finite(pairs$WindSpeed_kmh) & is.finite(pairs$WindDir_deg)

    slope_vals  <- rep(NA_real_, nrow(pairs))
    aspect_vals <- rep(NA_real_, nrow(pairs))

    if (any(valid_pts)) {
      pts_v <- vect(
        data.frame(x = pairs$centroid_x[valid_pts], y = pairs$centroid_y[valid_pts]),
        geom = c("x", "y"), crs = crs(slope_r)
      )
      slope_vals[valid_pts]  <- as.numeric(terra::extract(slope_r,  pts_v)[, 2])
      aspect_vals[valid_pts] <- as.numeric(terra::extract(aspect_r, pts_v)[, 2])
    }

    pairs$Slope_deg   <- slope_vals
    pairs$Aspect_deg  <- aspect_vals
    pairs$EffectiveWind <- compute_effective_wind(
      wind_speed          = pairs$WindSpeed_kmh,
      wind_dir            = pairs$WindDir_deg,
      slope_deg           = slope_vals,
      aspect_deg          = aspect_vals,
      combustion_buoyancy = combustion_buoyancy
    )
    message(sprintf(
      "Effective wind: mean=%.2f km/h  raw wind mean=%.2f km/h  slope addition mean=%.2f km/h  (Cb=%.0f)",
      mean(pairs$EffectiveWind,  na.rm = TRUE),
      mean(pairs$WindSpeed_kmh,  na.rm = TRUE),
      mean(pairs$EffectiveWind - pmax(0, pairs$WindSpeed_kmh), na.rm = TRUE),
      combustion_buoyancy
    ))
  } else {
    message(sprintf(
      "No terrain rasters provided — using raw wind speed (km/h) as EffectiveWind (Cb=%.0f not applied).",
      combustion_buoyancy
    ))
    pairs$Slope_deg     <- NA_real_
    pairs$Aspect_deg    <- NA_real_
    pairs$EffectiveWind <- pairs$WindSpeed_kmh
  }

  spread_pairs <- pairs %>% filter(use_for_spread_prob)

  message(sprintf(
    "Pairs: %d total | %d spread-prob (gap<=%d) | %d daily-area (gap<=%d)",
    nrow(pairs), sum(pairs$use_for_spread_prob, na.rm = TRUE), max_gap_sp,
    sum(pairs$use_for_daily_area, na.rm = TRUE), max_gap_da
  ))

  list(pairs = pairs, spread_pairs = spread_pairs)
}


# -----------------------------------------------------------------------------
# fit_max_daily_area
#
# Calibrates MaximumSpreadAreaB0/B1/B2 for SCF's maximum daily spread area
# equation, which in FireEvent.cs is implemented as a strictly LINEAR model:
#
#   MaxSpreadArea = (int)(B0 + B1 * FWI + B2 * effectiveWindSpeed)
#
# The result is in HECTARES and is used directly as a ceiling on how much
# area a fire can burn in a single day before the simulation advances to the
# next day.  Coefficients must therefore be fit on the RAW hectare scale
# (no log transform) so they can be entered directly into the SCF parameter
# file without any back-transformation.
#
# Calibration target:
#   MaxSpreadArea is an UPPER LIMIT (ceiling), not an average.  Calibrating
#   to the mean daily increment means ~50 % of observed fire-days would be
#   capped, slowing simulated fires unrealistically.  We therefore provide
#   fits at the mean AND at the 75th and 90th percentile of observed daily
#   area increments.  The 75th-percentile fit ("high" target) is recommended
#   as the starting point for SCF — it allows most fire-days to burn near
#   the observed historical maximum while still providing a constraint.
#
# Returns: list(coef = tibble of coefficients at mean/q75/q90,
#               model_mean, model_q75, model_q90, data)
# -----------------------------------------------------------------------------

fit_max_daily_area <- function(pairs) {

  dat <- pairs %>%
    filter(use_for_daily_area,
           is.finite(daily_area_ha), daily_area_ha > 0,
           is.finite(FWI), is.finite(EffectiveWind))

  if (nrow(dat) < 5)
    stop("Fewer than 5 usable pairs for max daily area fit.")

  # ---- Helper: build tidy coef table from an lm object ----------------------
  tidy_max_area <- function(m, label) {
    broom::tidy(m) %>%
      rename(std_error = std.error, p_value = p.value) %>%
      mutate(
        term = case_when(
          term == "(Intercept)"   ~ "B0_Intercept",
          term == "FWI"           ~ "B1_FWI",
          term == "EffectiveWind" ~ "B2_EffectiveWind",
          TRUE ~ term
        ),
        fit_target = label,
        scf_parameter = case_when(
          grepl("B0", term) ~ "MaximumSpreadAreaB0",
          grepl("B1", term) ~ "MaximumSpreadAreaB1",
          grepl("B2", term) ~ "MaximumSpreadAreaB2",
          TRUE ~ NA_character_
        ),
        note = "Raw linear fit (ha/day) — enter directly in SCF parameter file"
      )
  }

  # ---- Weighted quantile helper (no extra package needed) -------------------
  wt_quantile <- function(y, p) {
    # Converts the target quantile to a weight adjustment so a standard
    # WLS with asymmetric L1 loss approximates the p-th quantile regression.
    # For large samples this is nearly identical to rq(); for small samples
    # it's a reasonable approximation.
    wts <- ifelse(y >= quantile(y, p, na.rm = TRUE), p, 1 - p) + 1e-6
    lm(daily_area_ha ~ FWI + EffectiveWind, data = dat, weights = wts)
  }

  # ---- Three fits -----------------------------------------------------------
  m_mean <- lm(daily_area_ha ~ FWI + EffectiveWind, data = dat)
  m_q75  <- wt_quantile(dat$daily_area_ha, 0.75)
  m_q90  <- wt_quantile(dat$daily_area_ha, 0.90)

  coef_df <- bind_rows(
    tidy_max_area(m_mean, "Mean (OLS)"),
    tidy_max_area(m_q75,  "75th percentile (recommended)"),
    tidy_max_area(m_q90,  "90th percentile")
  )

  # ---- Sanity check: warn if B0 < 0 -----------------------------------------
  b0_vals <- coef_df %>% filter(term == "B0_Intercept") %>% pull(estimate)
  if (any(b0_vals < 0)) {
    warning(
      "One or more B0_Intercept estimates are negative (",
      paste(round(b0_vals[b0_vals < 0], 2), collapse = ", "), " ha). ",
      "SCF will warn 'MaxSpreadArea < 0' when FWI and wind are both low. ",
      "Consider using the 75th or 90th percentile fit, or constraining B0 >= 1."
    )
  }

  message(sprintf(
    "Max daily area fit: n=%d pairs | daily area range [%.1f, %.1f] ha | mean=%.0f ha",
    nrow(dat), min(dat$daily_area_ha), max(dat$daily_area_ha),
    mean(dat$daily_area_ha)
  ))
  message("  Recommended target: 75th percentile fit (MaxSpreadArea is a ceiling, not a mean)")

  list(coef   = coef_df,
       model_mean = m_mean,
       model_q75  = m_q75,
       model_q90  = m_q90,
       data   = dat)
}


# -----------------------------------------------------------------------------
# fit_spread_probability
# Logistic model from rasterized perimeter pairs.
# -----------------------------------------------------------------------------

fit_spread_probability <- function(pairs2, geomac_v, template_r, park_mask,
                                    fine_fuels_r  = NULL,
                                    failure_ratio = 3, max_samples = 25000,
                                    progress_fn = NULL) {

  # Use EffectiveWind when available; fall back to WindSpeed_kmh if no terrain
  has_eff_wind <- "EffectiveWind" %in% names(pairs2) &&
                  any(is.finite(pairs2$EffectiveWind))

  wind_col <- if (has_eff_wind) "EffectiveWind" else "WindSpeed_kmh"
  message(sprintf(
    "Spread probability: using '%s' as wind predictor (matches SCF's effectiveWindSpeed).",
    wind_col
  ))
  if (!has_eff_wind)
    message("  NOTE: terrain rasters not loaded — raw wind speed used as placeholder.")

  eligible <- pairs2 %>%
    filter(use_for_spread_prob, is.finite(.data[[wind_col]]), is.finite(FWI))

  if (nrow(eligible) == 0)
    stop("No eligible perimeter pairs for spread probability fit.")

  samples_list <- vector("list", nrow(eligible))

  for (i in seq_len(nrow(eligible))) {
    if (!is.null(progress_fn)) progress_fn(i, nrow(eligible))

    row  <- eligible[i, ]
    fwi  <- row$FWI
    ews  <- row[[wind_col]]     # effective wind speed (or raw if no terrain)
    if (!is.finite(ews)) next

    r_t  <- rasterize_perim(geomac_v, row$row_id,      template_r, park_mask)
    r_t1 <- rasterize_perim(geomac_v, row$row_id_next, template_r, park_mask)

    if (is.null(r_t) || is.null(r_t1)) {
      message(sprintf("  Pair %d: skipped — rasterize_perim returned NULL (row_id %s / %s)",
                      i, row$row_id, row$row_id_next))
      next
    }

    n_burned_t  <- sum(values(r_t)  == 1, na.rm = TRUE)
    n_burned_t1 <- sum(values(r_t1) == 1, na.rm = TRUE)
    message(sprintf("  Pair %d (fire %s, %s): r_t cells=%d  r_t1 cells=%d",
                    i, row$FIRE_ID, row$date, n_burned_t, n_burned_t1))

    dil        <- focal(r_t, w = ROOK_W, fun = "max", na.policy = "omit", fillvalue = 0)
    cand       <- (dil == 1) & (r_t == 0)
    succ_cells <- which(values(cand & (r_t1 == 1)) == 1)
    fail_cells <- which(values(cand & (r_t1 == 0)) == 1)

    message(sprintf("    cand=%d  succ=%d  fail=%d",
                    sum(values(cand) == 1, na.rm = TRUE),
                    length(succ_cells), length(fail_cells)))

    n_succ <- length(succ_cells)
    if (n_succ == 0) next

    n_fail_keep <- min(length(fail_cells), n_succ * failure_ratio)
    fail_keep   <- if (n_fail_keep > 0) sample(fail_cells, n_fail_keep) else integer(0)

    total <- n_succ + length(fail_keep)
    if (total > max_samples) {
      cap       <- max(0L, max_samples - n_succ)
      fail_keep <- if (cap < length(fail_keep)) sample(fail_keep, cap) else fail_keep
    }

    n_f <- length(fail_keep)
    if (n_succ + n_f == 0) next

    # Extract fine fuel values at sample cell locations (or use 1.0 placeholder)
    if (!is.null(fine_fuels_r)) {
      all_cells    <- c(succ_cells, fail_keep)
      ff_vals      <- as.numeric(values(fine_fuels_r)[all_cells])
      ff_mean      <- mean(ff_vals, na.rm = TRUE)
      if (!is.finite(ff_mean)) ff_mean <- 1.0
      ff_vals[!is.finite(ff_vals)] <- ff_mean
    } else {
      ff_vals <- rep(1.0, n_succ + n_f)
    }

    samples_list[[i]] <- tibble(
      FIRE_ID       = rep(row$FIRE_ID, n_succ + n_f),
      date          = rep(row$date,    n_succ + n_f),
      y             = c(rep(1L, n_succ), rep(0L, n_f)),
      FWI           = rep(fwi, n_succ + n_f),
      EffectiveWind = rep(ews, n_succ + n_f),   # always named EWS in samples
      FineFuels     = ff_vals
    )
  }

  all_samples <- bind_rows(samples_list)
  if (nrow(all_samples) == 0)
    stop("No raster samples generated. Check that perimeters fall within template extent.")

  # GLM: logit link — coefficients go DIRECTLY into SCF parameter file.
  # SCF CanSpread():
  #   Pspread = 1 / (1 + exp(-(B0 + B1*FWI + B2*FineFuels + B3*EWS)))
  # which is exactly what binomial(logit) produces.
  # All coefficients should be positive (higher FWI/fuels/wind → more spread).
  m <- glm(y ~ FWI + FineFuels + EffectiveWind,
           data   = all_samples,
           family = binomial(link = "logit"))

  fine_fuels_label <- if (is.null(fine_fuels_r)) "B2_FineFuels (placeholder=1.0)" else "B2_FineFuels"
  wind_label       <- if (has_eff_wind) "B3_EffectiveWind" else "B3_EffectiveWind (raw wind, no terrain)"

  coef_df <- broom::tidy(m) %>%
    rename(std_error = std.error, p_value = p.value) %>%
    mutate(
      term = case_when(
        term == "(Intercept)"   ~ "B0_Intercept",
        term == "FWI"           ~ "B1_FWI",
        term == "FineFuels"     ~ fine_fuels_label,
        term == "EffectiveWind" ~ wind_label,
        TRUE ~ term
      ),
      scf_parameter = case_when(
        grepl("B0", term) ~ "SpreadProbabilityB0",
        grepl("B1", term) ~ "SpreadProbabilityB1",
        grepl("B2", term) ~ "SpreadProbabilityB2",
        grepl("B3", term) ~ "SpreadProbabilityB3",
        TRUE ~ NA_character_
      ),
      note = "logit-link coefficients — enter directly in SCF parameter file (no back-transform needed)"
    )

  message(sprintf("Spread prob fit: %d samples (%d successes, %d failures)",
                  nrow(all_samples),
                  sum(all_samples$y == 1),
                  sum(all_samples$y == 0)))

  list(coef = coef_df, model = m, samples = all_samples,
       wind_predictor = wind_col)
}


# -----------------------------------------------------------------------------
# rasterize_perim  (internal)
# -----------------------------------------------------------------------------

rasterize_perim <- function(v, row_id_val, template, mask_r = NULL) {
  vv <- v[v$row_id == row_id_val, ]
  if (nrow(vv) == 0) return(NULL)
  r <- rasterize(vv, template, field = 1, touches = TRUE)
  r[!is.na(r)] <- 1
  r[is.na(r)]  <- 0
  if (!is.null(mask_r)) { r <- r * mask_r; r[is.na(r)] <- 0 }
  r
}


# -----------------------------------------------------------------------------
# Diagnostic plots
# -----------------------------------------------------------------------------

plot_spread_prob_diagnostics <- function(pairs, spread_prob_fit) {

  if (is.null(spread_prob_fit))
    return(ggplot() + labs(title = "Spread probability not yet fitted.") + theme_bw())

  coef_df <- spread_prob_fit$coef

  # Match coefficients by name prefix (term may have suffix like "(placeholder)")
  get_b <- function(prefix) {
    idx <- grep(prefix, coef_df$term, fixed = FALSE)[1]
    if (is.na(idx)) NA_real_ else coef_df$estimate[idx]
  }
  b0 <- get_b("B0_Intercept")
  b1 <- get_b("B1_FWI")
  b2 <- get_b("B2_Fine")
  b3 <- get_b("B3_Effective")

  if (any(is.na(c(b0, b1))))
    return(ggplot() + labs(title = "Spread probability not yet fitted.") + theme_bw())

  # Use effective wind when available in pairs
  ews_col  <- if ("EffectiveWind" %in% names(pairs)) "EffectiveWind" else "WindSpeed_kmh"
  ews_lo   <- quantile(pairs[[ews_col]], 0.25, na.rm = TRUE)
  ews_mean <- mean(pairs[[ews_col]], na.rm = TRUE)
  ews_hi   <- quantile(pairs[[ews_col]], 0.75, na.rm = TRUE)

  fwi_seq <- seq(0, max(pairs$FWI, na.rm = TRUE), length.out = 100)

  make_line <- function(ews_val, label) {
    tibble(
      FWI      = fwi_seq,
      prob     = plogis(coalesce(b0, 0) + coalesce(b1, 0) * fwi_seq +
                        coalesce(b2, 0) * 1.0 +       # fine fuels = 1 (fully loaded)
                        coalesce(b3, 0) * ews_val),
      scenario = label
    )
  }

  pred_df <- bind_rows(
    make_line(ews_lo,   sprintf("EWS Q25 = %.1f km/h", ews_lo)),
    make_line(ews_mean, sprintf("EWS mean = %.1f km/h", ews_mean)),
    make_line(ews_hi,   sprintf("EWS Q75 = %.1f km/h", ews_hi))
  )

  ggplot(pred_df, aes(FWI, prob, colour = scenario)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_colour_manual(values = c("grey55", "#2c3e50", "#c0392b")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title    = "Fitted Spread Probability vs FWI",
      subtitle = paste0(
        "P(spread) = logistic(B0 + B1\u00d7FWI + B2\u00d7FineFuels + B3\u00d7EWS)\n",
        "Shown at fine fuels = 1.0 (fully loaded); three effective wind scenarios (km/h)."
      ),
      x = "FWI", y = "P(spread to neighbor cell)",
      colour = "Wind scenario"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}


plot_max_area_diagnostics <- function(pairs, max_area_fit) {

  dat <- pairs %>%
    filter(use_for_daily_area,
           is.finite(daily_area_ha), daily_area_ha > 0,
           is.finite(FWI), is.finite(EffectiveWind))

  if (nrow(dat) == 0 || is.null(max_area_fit))
    return(ggplot() + labs(title = "No max daily area data.") + theme_bw())

  coef_df <- max_area_fit$coef

  # Extract coefficients for each fit target
  make_pred <- function(target_label) {
    b <- coef_df %>%
      filter(fit_target == target_label) %>%
      { setNames(.$estimate, .$term) }
    b0 <- b["B0_Intercept"];   b1 <- b["B1_FWI"];  b2 <- b["B2_EffectiveWind"]
    if (any(is.na(c(b0, b1, b2)))) return(NULL)
    dat %>%
      mutate(pred_ha = b0 + b1 * FWI + b2 * EffectiveWind,
             fit_target = target_label)
  }

  pred_df <- bind_rows(
    make_pred("Mean (OLS)"),
    make_pred("75th percentile (recommended)"),
    make_pred("90th percentile")
  )

  if (nrow(pred_df) == 0)
    return(ggplot() + labs(title = "No predicted values.") + theme_bw())

  q75 <- quantile(dat$daily_area_ha, 0.75, na.rm = TRUE)
  q90 <- quantile(dat$daily_area_ha, 0.90, na.rm = TRUE)

  ggplot() +
    # Observed data as points
    geom_point(data = dat, aes(FWI, daily_area_ha),
               alpha = 0.35, color = "grey40", size = 1.2) +
    # Fitted lines for each target
    geom_line(data = pred_df,
              aes(FWI, pred_ha, color = fit_target, linetype = fit_target),
              linewidth = 1.1) +
    geom_hline(yintercept = q75, linetype = "dotted", colour = "#2c5f2e", linewidth = 0.6) +
    geom_hline(yintercept = q90, linetype = "dotted", colour = "#a0522d", linewidth = 0.6) +
    annotate("text", x = max(dat$FWI, na.rm = TRUE), y = q75,
             label = sprintf("obs Q75 = %.0f ha", q75),
             hjust = 1.05, vjust = -0.4, size = 3, colour = "#2c5f2e") +
    annotate("text", x = max(dat$FWI, na.rm = TRUE), y = q90,
             label = sprintf("obs Q90 = %.0f ha", q90),
             hjust = 1.05, vjust = -0.4, size = 3, colour = "#a0522d") +
    scale_color_manual(values = c(
      "Mean (OLS)"                      = "grey30",
      "75th percentile (recommended)"   = "#2c5f2e",
      "90th percentile"                 = "#a0522d"
    )) +
    scale_linetype_manual(values = c(
      "Mean (OLS)"                      = "dashed",
      "75th percentile (recommended)"   = "solid",
      "90th percentile"                 = "longdash"
    )) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
    labs(
      title    = "Max Daily Spread Area: Observed vs Fitted",
      subtitle = paste0(
        "Points = observed fire-day increments (ha).  Lines = SCF MaxSpreadArea = B0 + B1\u00d7FWI + B2\u00d7EWS.\n",
        "Recommended: 75th-percentile fit (MaxSpreadArea is a ceiling, not an average)."
      ),
      x = "FWI at ignition", y = "Daily area increment (ha)",
      color = "Fit target", linetype = "Fit target"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# score_candidates
# -----------------------------------------------------------------------------

score_candidates <- function(candidate_grid, runs_root, cell_area_ha,
                              cal_start_year, cal_end_year,
                              w_size = 0.5, w_aab = 0.5,
                              fpa_sizes = NULL) {

  map_dfr(seq_len(nrow(candidate_grid)), function(i) {
    cand    <- candidate_grid[i, ]
    run_dir <- file.path(runs_root,
                         sprintf("candidate_%04d", cand$candidate_id),
                         "output")

    if (!dir.exists(run_dir)) {
      return(cand %>% mutate(run_dir = run_dir,
                             score_size = NA_real_,
                             score_aab  = NA_real_,
                             score_total = NA_real_))
    }

    ev <- read_scf_events(run_dir)
    sim_sizes <- ev %>%
      mutate(event_area_ha = sites_burned * cell_area_ha) %>%
      filter(FIRE_YEAR >= cal_start_year, FIRE_YEAR <= cal_end_year,
             is.finite(event_area_ha), event_area_ha > 0)

    sim_annual <- sim_sizes %>%
      group_by(FIRE_YEAR) %>%
      summarise(area_burned_ha = sum(event_area_ha), .groups = "drop") %>%
      complete(FIRE_YEAR = seq(cal_start_year, cal_end_year),
               fill = list(area_burned_ha = 0))

    s_size <- if (!is.null(fpa_sizes) && nrow(sim_sizes) > 0)
      score_size_quantiles(fpa_sizes, sim_sizes$event_area_ha) else NA_real_

    s_aab <- score_annual_area_rmse(sim_annual)

    cand %>% mutate(run_dir = run_dir,
                    score_size  = s_size,
                    score_aab   = s_aab,
                    score_total = w_size * replace_na(s_size, 1) +
                                  w_aab  * replace_na(s_aab,  1))
  }) %>% arrange(score_total)
}


score_size_quantiles <- function(obs, sim, probs = seq(0.1, 0.9, 0.1)) {
  if (length(obs) == 0 || length(sim) == 0) return(NA_real_)
  mean(abs(log1p(quantile(obs, probs, na.rm = TRUE)) -
           log1p(quantile(sim, probs, na.rm = TRUE))))
}

score_annual_area_rmse <- function(sim_annual) {
  if (nrow(sim_annual) == 0) return(NA_real_)
  sqrt(mean(sim_annual$area_burned_ha^2, na.rm = TRUE))
}

read_scf_events <- function(run_dir) {
  candidates <- c(
    file.path(run_dir, "FireEventLog.csv"),
    file.path(run_dir, "fire-event-log.csv"),
    file.path(run_dir, "SCF-event-log.csv")
  )
  f <- candidates[file.exists(candidates)]
  if (length(f) == 0)
    return(tibble(FIRE_YEAR = integer(), sites_burned = integer()))

  ev  <- read_csv(f[1], show_col_types = FALSE)
  nm  <- tolower(names(ev))
  yr  <- names(ev)[which(nm %in% c("year","fire_year","fireyear"))[1]]
  st  <- names(ev)[which(nm %in% c("sites_burned","sitesburned","total_sites"))[1]]

  if (is.na(yr) || is.na(st)) {
    warning("Could not parse FireEventLog columns in: ", f[1])
    return(tibble(FIRE_YEAR = integer(), sites_burned = integer()))
  }

  ev %>%
    transmute(FIRE_YEAR    = as.integer(.data[[yr]]),
              sites_burned = as.numeric(.data[[st]])) %>%
    filter(!is.na(FIRE_YEAR), is.finite(sites_burned))
}


# -----------------------------------------------------------------------------
# write_scf_snippet
# -----------------------------------------------------------------------------

write_scf_snippet <- function(spread_coef_df, max_area_fit,
                               max_area_target = "75th percentile (recommended)") {
  # ---- Spread probability coefficients (logit-link, enter directly) ----------
  get_sp <- function(prefix) {
    idx <- grep(prefix, spread_coef_df$term)[1]
    if (is.na(idx)) NA_real_ else spread_coef_df$estimate[idx]
  }
  sp_b0 <- get_sp("B0_Intercept")
  sp_b1 <- get_sp("B1_FWI")
  sp_b2 <- get_sp("B2_Fine")
  sp_b3 <- get_sp("B3_Effective")

  # ---- Max daily area coefficients (linear, ha — enter directly) -------------
  ma_coef <- max_area_fit$coef %>%
    filter(fit_target == max_area_target)
  get_ma <- function(prefix) {
    idx <- grep(prefix, ma_coef$term)[1]
    if (is.na(idx)) NA_real_ else ma_coef$estimate[idx]
  }
  ma_b0 <- get_ma("B0_Intercept")
  ma_b1 <- get_ma("B1_FWI")
  ma_b2 <- get_ma("B2_Effective")

  paste(c(
    ">> =====================================================================",
    ">> SCF Social Climate Fire — Calibrated Spread Parameters",
    ">> Generated by scf_calibration_app",
    ">> =====================================================================",
    "",
    ">> IMPORTANT UNIT NOTES:",
    ">>   EWS units:        km/h  (SCF User Guide §2.43/2.47/2.51: climate library",
    ">>                     converts all wind speed units to km/h; ERA5 m/s * 3.6)",
    ">>   Fine fuels scale: 0-1  (SiteVars.FineFuels / MaxFineFuels in SCF)",
    ">>   MaxFineFuels:     set in your SCF parameter file; calibration normalised",
    ">>                     to max FBFM40 load = 9.0 t/acre (GR9 fuel model).",
    ">>                     If your SCF MaxFineFuels differs, B2 needs rescaling.",
    ">>   CombustionBuoyancy used in calibration (Cb) is shown in the EWS note.",
    ">>   SCF uses Cb=10 (intensity<=3), Cb=25 (4-6), Cb=50 (>6); calibrate at",
    ">>   the Cb matching your landscape's typical fire intensity.",
    "",
    ">> Cell-to-cell spread probability (logistic / binomial GLM)",
    ">>   P(spread) = 1 / (1 + exp(-(B0 + B1*FWI + B2*FineFuels + B3*EWS)))",
    ">>   FineFuels = SiteVars.FineFuels / MaxFineFuels  (0-1 scale)",
    ">>   EWS       = effective wind speed (Nelson 2002 Eq.5, km/h)",
    ">>   Coefficients are logit-scale — enter directly, no back-transform needed.",
    sprintf("SpreadProbabilityB0    %9.6f", coalesce(sp_b0, NA_real_)),
    sprintf("SpreadProbabilityB1    %9.6f", coalesce(sp_b1, NA_real_)),
    sprintf("SpreadProbabilityB2    %9.6f", coalesce(sp_b2, NA_real_)),
    sprintf("SpreadProbabilityB3    %9.6f", coalesce(sp_b3, NA_real_)),
    "",
    sprintf(">> Maximum daily spread area (%s fit)", max_area_target),
    ">>   MaxSpreadArea = (int)(B0 + B1*FWI + B2*EWS)   [units: HECTARES]",
    ">>   EWS in km/h (same as above). B0 must be > 0 to avoid MaxSpreadArea<=0.",
    ">>   When fire burns > MaxSpreadArea ha in a day, simulation advances day+1.",
    ">>   Coefficients are on the RAW ha scale — enter directly.",
    sprintf("MaximumSpreadAreaB0    %9.4f", coalesce(ma_b0, NA_real_)),
    sprintf("MaximumSpreadAreaB1    %9.4f", coalesce(ma_b1, NA_real_)),
    sprintf("MaximumSpreadAreaB2    %9.4f", coalesce(ma_b2, NA_real_)),
    ""
  ), collapse = "\n")
}


# -----------------------------------------------------------------------------
# Scoring diagnostic plots
# -----------------------------------------------------------------------------

plot_candidate_cdf <- function(scores_df, obs_sizes, runs_root, cell_area_ha,
                                cal_start_year, cal_end_year) {

  top <- scores_df %>% filter(is.finite(score_total)) %>% slice_head(n = 3)
  if (nrow(top) == 0 || is.null(obs_sizes))
    return(ggplot() + labs(title = "No scored candidates yet.") + theme_bw())

  obs_cdf <- tibble(size_ha = sort(obs_sizes),
                    cdf     = seq_along(obs_sizes) / length(obs_sizes),
                    source  = "Observed (FPA FOD)")

  sim_list <- map(seq_len(nrow(top)), function(i) {
    cand    <- top[i, ]
    run_dir <- file.path(runs_root,
                         sprintf("candidate_%04d", cand$candidate_id), "output")
    ev  <- read_scf_events(run_dir)
    sz  <- ev %>%
      mutate(size_ha = sites_burned * cell_area_ha) %>%
      filter(FIRE_YEAR >= cal_start_year, FIRE_YEAR <= cal_end_year,
             is.finite(size_ha), size_ha > 0) %>%
      pull(size_ha)
    if (length(sz) == 0) return(NULL)
    tibble(size_ha = sort(sz),
           cdf     = seq_along(sz) / length(sz),
           source  = sprintf("Candidate %d (score=%.3f)",
                             cand$candidate_id, cand$score_total))
  }) %>% compact()

  bind_rows(obs_cdf, sim_list) %>%
    ggplot(aes(size_ha, cdf, color = source)) +
    geom_line(linewidth = 1) +
    scale_x_log10(labels = comma) +
    labs(title = "Event Size CDF: Observed vs Top Candidates",
         x = "Fire Event Size (ha, log scale)",
         y = "Cumulative probability", color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}


plot_candidate_annual <- function(scores_df, runs_root, cell_area_ha,
                                   cal_start_year, cal_end_year) {

  top <- scores_df %>% filter(is.finite(score_total)) %>% slice_head(n = 1)
  if (nrow(top) == 0)
    return(ggplot() + labs(title = "No scored candidates yet.") + theme_bw())

  run_dir <- file.path(runs_root,
                       sprintf("candidate_%04d", top$candidate_id[1]), "output")
  ev <- read_scf_events(run_dir)

  ev %>%
    mutate(area_ha = sites_burned * cell_area_ha) %>%
    filter(FIRE_YEAR >= cal_start_year, FIRE_YEAR <= cal_end_year) %>%
    group_by(FIRE_YEAR) %>%
    summarise(area_burned_ha = sum(area_ha), .groups = "drop") %>%
    complete(FIRE_YEAR = cal_start_year:cal_end_year,
             fill = list(area_burned_ha = 0)) %>%
    ggplot(aes(FIRE_YEAR, area_burned_ha)) +
    geom_bar(stat = "identity", fill = "#a0522d", alpha = 0.7) +
    labs(title = sprintf("Annual Area Burned — Candidate %d", top$candidate_id[1]),
         x = "Year", y = "Area Burned (ha)") +
    theme_bw(base_size = 12)
}


# =============================================================================
# Threshold Grid Search Functions
# =============================================================================

# -----------------------------------------------------------------------------
# compute_auc
#
# Simple AUC via Wilcoxon rank-sum statistic. No external packages required.
# Returns NA if only one class is present.
# -----------------------------------------------------------------------------

compute_auc <- function(y_true, y_prob) {
  n1 <- sum(y_true == 1, na.rm = TRUE)
  n0 <- sum(y_true == 0, na.rm = TRUE)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(y_prob, ties.method = "average")
  (sum(r[y_true == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


# -----------------------------------------------------------------------------
# build_spread_samples
#
# Rasterizes ALL perimeter pairs and returns a list of per-pair sample tibbles.
# This is the one-time expensive step that builds the cache for threshold search.
#
# Arguments:
#   pairs_all            Full pairs tibble (already filtered to widest thresholds,
#                        with is.finite(EffectiveWind) & is.finite(FWI))
#   geomac_v             SpatVector of GeoMAC perimeters (from load_geomac_perimeters)
#   template_r           SpatRaster defining the rasterization grid
#   park_mask            Binary SpatRaster mask (1 = inside boundary, NA = outside)
#   fine_fuels_r         Optional SpatRaster of fine fuels (0-1). NULL -> use 1.0.
#   max_cache_fail_ratio Maximum failure:success ratio to store per pair (pre-subsampling)
#   progress_fn          Optional function(i, n_total) for progress reporting
#
# Returns a list of non-NULL tibbles (one per eligible pair), each with columns:
#   FIRE_ID, date, gap_days, delta_ha, y, FWI, EffectiveWind, FineFuels
# -----------------------------------------------------------------------------

build_spread_samples <- function(pairs_all, geomac_v, template_r, park_mask,
                                  fine_fuels_r = NULL,
                                  max_cache_fail_ratio = 10,
                                  progress_fn = NULL) {

  n_total <- nrow(pairs_all)
  message(sprintf("build_spread_samples: rasterizing %d pairs for threshold search cache",
                  n_total))

  result_list <- vector("list", n_total)

  for (i in seq_len(n_total)) {
    if (!is.null(progress_fn)) progress_fn(i, n_total)

    row <- pairs_all[i, ]

    r_t  <- rasterize_perim(geomac_v, row$row_id,      template_r, park_mask)
    r_t1 <- rasterize_perim(geomac_v, row$row_id_next, template_r, park_mask)

    if (is.null(r_t) || is.null(r_t1)) next

    dil        <- focal(r_t, w = ROOK_W, fun = "max", na.policy = "omit", fillvalue = 0)
    cand       <- (dil == 1) & (r_t == 0)
    succ_cells <- which(values(cand & (r_t1 == 1)) == 1)
    fail_cells <- which(values(cand & (r_t1 == 0)) == 1)

    n_succ <- length(succ_cells)
    if (n_succ == 0) next

    n_fail_keep <- min(length(fail_cells), n_succ * max_cache_fail_ratio)
    fail_keep   <- if (n_fail_keep > 0) sample(fail_cells, n_fail_keep) else integer(0)

    all_cells <- c(succ_cells, fail_keep)
    n_f       <- length(fail_keep)

    if (!is.null(fine_fuels_r)) {
      ff_vals  <- as.numeric(values(fine_fuels_r)[all_cells])
      ff_mean  <- mean(ff_vals, na.rm = TRUE)
      if (!is.finite(ff_mean)) ff_mean <- 1.0
      ff_vals[!is.finite(ff_vals)] <- ff_mean
    } else {
      ff_vals <- rep(1.0, n_succ + n_f)
    }

    result_list[[i]] <- tibble(
      FIRE_ID       = rep(as.character(row$FIRE_ID), n_succ + n_f),
      date          = rep(as.Date(row$date),         n_succ + n_f),
      gap_days      = rep(as.integer(row$gap_days),  n_succ + n_f),
      delta_ha      = rep(as.numeric(row$delta_ha),  n_succ + n_f),
      y             = c(rep(1L, n_succ), rep(0L, n_f)),
      FWI           = rep(as.numeric(row$FWI),           n_succ + n_f),
      EffectiveWind = rep(as.numeric(row$EffectiveWind),  n_succ + n_f),
      FineFuels     = ff_vals
    )
  }

  Filter(Negate(is.null), result_list)
}


# -----------------------------------------------------------------------------
# fit_spread_prob_from_cache
#
# Applies threshold filters and subsampling to cached per-pair samples,
# then fits the binomial GLM.
#
# Arguments:
#   per_pair_samples  List of tibbles from build_spread_samples()
#   max_gap_days      Maximum gap_days to include
#   neg_tol_ha        Pairs with delta_ha < -neg_tol_ha are excluded
#   failure_ratio     Max failures kept per success per pair
#   max_samples_pair  Cap on total samples per pair
#
# Returns list(auc, n_pairs, n_samples, model) or NULL on failure.
# -----------------------------------------------------------------------------

fit_spread_prob_from_cache <- function(per_pair_samples,
                                        max_gap_days    = 3,
                                        neg_tol_ha      = 1.0,
                                        failure_ratio   = 3,
                                        max_samples_pair = 25000) {

  filtered <- lapply(per_pair_samples, function(tbl) {
    tbl <- tbl[tbl$gap_days <= max_gap_days & tbl$delta_ha >= -neg_tol_ha, ]
    if (nrow(tbl) == 0) return(NULL)

    succ_rows <- tbl[tbl$y == 1, ]
    fail_rows <- tbl[tbl$y == 0, ]
    n_succ    <- nrow(succ_rows)
    if (n_succ == 0) return(NULL)

    n_fail_keep <- min(nrow(fail_rows), n_succ * failure_ratio)
    if (n_fail_keep > 0) {
      fail_rows <- fail_rows[sample(nrow(fail_rows), n_fail_keep), ]
    } else {
      fail_rows <- fail_rows[integer(0), ]
    }

    combined <- rbind(succ_rows, fail_rows)
    if (nrow(combined) > max_samples_pair) {
      combined <- combined[sample(nrow(combined), max_samples_pair), ]
    }
    combined
  })

  filtered <- Filter(Negate(is.null), filtered)
  if (length(filtered) == 0) return(NULL)

  all_samples <- bind_rows(filtered)
  if (nrow(all_samples) == 0) return(NULL)

  m <- tryCatch(
    glm(y ~ FWI + FineFuels + EffectiveWind,
        data   = all_samples,
        family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  auc <- compute_auc(all_samples$y, predict(m, type = "response"))

  list(
    auc       = auc,
    n_pairs   = length(filtered),
    n_samples = nrow(all_samples),
    model     = m
  )
}


# -----------------------------------------------------------------------------
# score_max_daily_area_variant
#
# Scores the max daily area model for a given threshold combination.
# Fast — no rasterization needed.
#
# Arguments:
#   pairs_all   Full pairs tibble (all gaps/tolerances present)
#   max_gap_da  Maximum gap_days to include for daily area model
#   neg_tol_ha  Exclude pairs with delta_ha < -neg_tol_ha
#
# Returns list(r2_q75, r2_mean, n_pairs, rmse)
# -----------------------------------------------------------------------------

score_max_daily_area_variant <- function(pairs_all, max_gap_da, neg_tol_ha) {

  dat <- pairs_all %>%
    filter(
      gap_days      <= max_gap_da,
      delta_ha      >= -neg_tol_ha,
      delta_ha      >  0,
      is.finite(daily_area_ha),
      daily_area_ha > 0,
      is.finite(FWI),
      is.finite(EffectiveWind)
    )

  if (nrow(dat) < 5)
    return(list(r2_q75  = NA_real_,
                r2_mean = NA_real_,
                n_pairs = nrow(dat),
                rmse    = NA_real_))

  m_mean <- lm(daily_area_ha ~ FWI + EffectiveWind, data = dat)

  wts_q75 <- ifelse(
    dat$daily_area_ha >= quantile(dat$daily_area_ha, 0.75, na.rm = TRUE),
    0.75, 0.25
  ) + 1e-6
  m_q75 <- lm(daily_area_ha ~ FWI + EffectiveWind, data = dat, weights = wts_q75)

  ss_res <- sum((dat$daily_area_ha - predict(m_q75))^2)
  ss_tot <- sum((dat$daily_area_ha - mean(dat$daily_area_ha))^2)
  r2_q75 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

  list(
    r2_q75  = r2_q75,
    r2_mean = summary(m_mean)$r.squared,
    n_pairs = nrow(dat),
    rmse    = sqrt(mean((dat$daily_area_ha - predict(m_mean))^2))
  )
}


# -----------------------------------------------------------------------------
# search_threshold_grid
#
# Loops over all combinations of thresholds and sampling controls, scoring
# each against spread probability AUC and max daily area R² (Q75 fit).
#
# Arguments:
#   per_pair_samples  List from build_spread_samples()
#   pairs_all         Full pairs tibble from bind_climate_to_pairs()$pairs
#   gap_sp_vals       Vector of max_gap_days values to try for spread probability
#   gap_da_vals       Vector of max_gap_days values to try for daily area
#   neg_tol_vals      Vector of neg_tol_ha values to try
#   fail_ratio_vals   Vector of failure:success ratios to try
#   max_samples_pair  Cap on cells per pair during GLM fitting
#   w_auc             Weight on AUC component of combined score
#   w_r2              Weight on R² (Q75) component of combined score
#   progress_fn       Optional function(i, n) for progress reporting
#
# Returns tibble sorted by score descending, with rank column added.
# Score = w_auc * AUC + w_r2 * max(0, R²_Q75), treating NA as neutral values.
# -----------------------------------------------------------------------------

search_threshold_grid <- function(per_pair_samples, pairs_all,
                                   gap_sp_vals     = c(1, 2, 3, 5),
                                   gap_da_vals     = c(3, 7, 14),
                                   neg_tol_vals    = c(0.5, 1.0, 2.0),
                                   fail_ratio_vals = c(2, 3, 5),
                                   max_samples_pair = 25000,
                                   w_auc = 0.6,
                                   w_r2  = 0.4,
                                   progress_fn = NULL) {

  grid <- expand.grid(
    gap_sp     = gap_sp_vals,
    gap_da     = gap_da_vals,
    neg_tol    = neg_tol_vals,
    fail_ratio = fail_ratio_vals,
    stringsAsFactors = FALSE
  )

  n_combos <- nrow(grid)
  message(sprintf("search_threshold_grid: evaluating %d combinations", n_combos))

  rows <- vector("list", n_combos)

  for (i in seq_len(n_combos)) {
    if (!is.null(progress_fn)) progress_fn(i, n_combos)

    g <- grid[i, ]

    sp_res <- fit_spread_prob_from_cache(
      per_pair_samples = per_pair_samples,
      max_gap_days     = g$gap_sp,
      neg_tol_ha       = g$neg_tol,
      failure_ratio    = g$fail_ratio,
      max_samples_pair = max_samples_pair
    )

    da_res <- score_max_daily_area_variant(
      pairs_all  = pairs_all,
      max_gap_da = g$gap_da,
      neg_tol_ha = g$neg_tol
    )

    auc     <- if (!is.null(sp_res)) sp_res$auc     else NA_real_
    n_sp    <- if (!is.null(sp_res)) sp_res$n_pairs  else 0L
    r2_q75  <- da_res$r2_q75
    n_da    <- da_res$n_pairs

    score <- w_auc * coalesce(auc, 0.5) + w_r2 * pmax(0, coalesce(r2_q75, 0))

    rows[[i]] <- tibble(
      gap_sp     = g$gap_sp,
      gap_da     = g$gap_da,
      neg_tol    = g$neg_tol,
      fail_ratio = g$fail_ratio,
      n_pairs_sp = n_sp,
      n_pairs_da = n_da,
      auc        = auc,
      r2_q75     = r2_q75,
      score      = score
    )
  }

  results <- bind_rows(rows) %>%
    arrange(desc(score)) %>%
    mutate(rank = seq_len(n()))

  best <- results[1, ]
  message(sprintf(
    "search_threshold_grid: best score=%.3f (gap_sp=%d, gap_da=%d, neg_tol=%.1f, fail_ratio=%d, AUC=%.3f, R2=%.3f)",
    best$score, best$gap_sp, best$gap_da, best$neg_tol, best$fail_ratio,
    coalesce(best$auc, NA_real_), coalesce(best$r2_q75, NA_real_)
  ))

  results
}
