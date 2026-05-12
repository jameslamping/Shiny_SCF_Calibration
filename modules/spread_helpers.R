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

# quantreg is the primary constrained fitter for MaxSpreadArea.
# nnls is the first fallback; plain OLS with clamping is the last resort.
for (.pkg in c("quantreg", "nnls")) {
  if (!requireNamespace(.pkg, quietly = TRUE))
    message(sprintf(paste0(
      "Optional package '%s' not found. Install it for better constrained ",
      "MaxSpreadArea regression: install.packages('%s')"), .pkg, .pkg))
}
rm(.pkg)


# =============================================================================
# Collinearity helpers (used by both fit_* functions)
# =============================================================================

# compute_vif: variance inflation factors for predictor columns (no intercept)
compute_vif <- function(X_pred) {
  if (ncol(X_pred) < 2)
    return(setNames(rep(NA_real_, ncol(X_pred)), colnames(X_pred)))
  tryCatch(
    diag(solve(cor(X_pred))),
    error = function(e) setNames(rep(NA_real_, ncol(X_pred)), colnames(X_pred))
  )
}

# condition_number: ratio of largest to smallest singular value of design matrix
condition_number <- function(X) {
  sv <- svd(X)$d
  max(sv) / min(sv[sv > .Machine$double.eps * max(sv)])
}

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
#   cal_vect        SpatVector of calibration boundary (any CRS)
#   working_crs     CRS string for output
#   spread_template SpatRaster defining the output grid (extent + resolution)
#   normalize       Logical (default TRUE). When TRUE, fuel loads are divided by
#                   max(FBFM40_FINE_FUELS) = 9.0 t/ac so the raster is 0-1
#                   scaled.  When FALSE, raw tons/acre values are returned
#                   (1-hr + 10-hr dead fuel, Scott & Burgan 2005).  Use FALSE
#                   to explore whether SCF should receive physical fuel loads
#                   rather than a dimensionless index.
#   progress_fn     optional function(msg) for status messages
#
# Returns a SpatRaster on spread_template grid.
#   normalize=TRUE : values 0-1 (NA = non-burnable)
#   normalize=FALSE: values 0-9 t/ac (NA = non-burnable)
# -----------------------------------------------------------------------------

fetch_landfire_fine_fuels <- function(cal_vect, working_crs, spread_template,
                                       normalize   = TRUE,
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

  # ---- Reclassify integer codes -> fine fuel load (tons/acre) ----------------
  if (!is.null(progress_fn)) progress_fn("Reclassifying FBFM40 to fine fuel load...")

  code_strs <- as.character(as.integer(round(raw_vals)))
  loads     <- as.numeric(FBFM40_FINE_FUELS[code_strs])   # tons/acre; NA for unknown codes

  max_load <- max(FBFM40_FINE_FUELS)                       # GR9 = 9.0 t/ac

  if (normalize) {
    loads <- loads / max_load                              # scale to 0-1
    if (!is.null(progress_fn)) progress_fn("Normalizing fine fuels to 0-1 index...")
    message(sprintf("FBFM40 normalized: divided by %.1f t/ac (GR9 max) -> 0-1 index", max_load))
  } else {
    message(sprintf("FBFM40 NOT normalized: values in tons/acre (range 0 – %.1f t/ac)", max_load))
  }

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
# Loads a user-provided raster, reprojects and resamples to the spread grid.
#
# Arguments:
#   path            Path to a GeoTIFF fine fuels raster (any numeric scale).
#   working_crs     CRS string for output.
#   spread_template SpatRaster defining the output grid (extent + resolution).
#   normalize       Logical (default TRUE). When TRUE, values are divided by
#                   the raster maximum so the output is 0-1 scaled.  When
#                   FALSE, the original values are preserved as-is.
# -----------------------------------------------------------------------------

load_local_fine_fuels <- function(path, working_crs, spread_template,
                                   normalize = TRUE) {

  if (!file.exists(path))
    stop("Fine fuels raster not found: ", path)

  r   <- rast(path)
  v   <- values(r, mat = FALSE)
  rng <- range(v, na.rm = TRUE)

  if (normalize) {
    mx <- max(v, na.rm = TRUE)
    if (is.finite(mx) && mx > 1) {
      message(sprintf("load_local_fine_fuels: normalizing by max (%.2f) -> 0-1", mx))
      v <- v / mx
    }
  } else {
    message(sprintf("load_local_fine_fuels: normalize=FALSE, values kept in original units (range %.2f – %.2f)",
                    rng[1], rng[2]))
  }

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


# =============================================================================
# load_fired_perimeters
#
# Loads FIRED (Fire Events Delineation) DAILY fire progression data as a
# drop-in replacement for load_geomac_perimeters().  Returns a terra SpatVector
# with the same schema (FIRE_ID, perimeter_date, AREA_HA, row_id) that
# bind_climate_to_pairs() expects.
#
# IMPORTANT — two FIRED products exist (Balch et al. 2020):
#   1. Events product  (fired_events_conus_*):  one row per fire, FINAL perimeter.
#                                                NOT usable here for spread prob.
#   2. Daily product   (fired_daily_conus_*):   one row per fire-day, daily patch.
#                                                THIS is what this function needs.
#
# Download FIRED CONUS DAILY GeoPackage from:
#   https://scholar.colorado.edu/concern/datasets/765372341
#   (Balch et al. 2020, Remote Sensing 12(21), 3498 — companion daily product)
#
# The events product (fired_events_conus_*) is useful for MaxSpreadArea
# calibration via load_fired_events_for_max_area() defined below.
#
# Daily polygon geometry represents NEWLY BURNED area on each specific day
# (incremental patches from MODIS MCD64A1 at 500 m).  This function builds
# CUMULATIVE perimeters by running union of daily patches within each event.
# The key advantage over GeoMAC: all polygons are from the same MODIS algorithm
# → consecutive perimeters are spatially consistent → high pct_proximal.
#
# Arguments:
#   fired_gpkg     path to the downloaded FIRED DAILY .gpkg file
#   cal_vect       SpatVector in working_crs (calibration boundary)
#   working_crs    target CRS string (from template_r)
#   cal_start      Date — earliest perimeter date to keep
#   cal_end        Date — latest perimeter date to keep
#   min_fire_ha    minimum FINAL fire size in ha; smaller events are dropped.
#                  Default 100 ha — removes tiny MODIS artefacts.
#   layer          GPKG layer name, or NULL to auto-detect.
# =============================================================================

load_fired_perimeters <- function(fired_gpkg, cal_vect, working_crs,
                                   cal_start, cal_end,
                                   min_fire_ha       = 100,
                                   min_detection_days = 3L,
                                   layer   = NULL,
                                   n_cores = NULL) {

  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' required for load_fired_perimeters(). install.packages('sf')")
  if (!file.exists(fired_gpkg))
    stop("FIRED GeoPackage not found: ", fired_gpkg)

  # ---- Detect layer -----------------------------------------------------------
  avail_layers <- tryCatch(sf::st_layers(fired_gpkg)$name, error = function(e) character(0))
  message("FIRED GeoPackage layers: ", paste(avail_layers, collapse = ", "))

  if (is.null(layer)) {
    # Prefer layers whose name contains "daily"; fall back to first available
    hit <- avail_layers[grepl("daily", tolower(avail_layers))]
    layer <- if (length(hit) > 0) hit[1] else avail_layers[1]
  }

  # ---- Guard: warn if this looks like the events layer, not daily -------------
  # The events layer name contains "events"; the daily layer contains "daily".
  # If the user points to the wrong file we tell them exactly where to download.
  if (grepl("event", tolower(layer)) && !grepl("daily", tolower(layer))) {
    stop(
      "The selected layer ('", layer, "') appears to be the FIRED EVENTS product\n",
      "(one row per fire event, final perimeter only).\n\n",
      "For spread probability calibration you need the FIRED DAILY product,\n",
      "which has one row per fire-day (daily burned area patches).\n\n",
      "Download the FIRED CONUS Daily GeoPackage from:\n",
      "  https://scholar.colorado.edu/concern/datasets/765372341\n\n",
      "Once downloaded, update FIRED_GPKG in Section 1 of the diagnostic script\n",
      "to point to the new file (fired_daily_conus_*.gpkg).\n\n",
      "To use the events file for MaxSpreadArea calibration instead, call:\n",
      "  load_fired_events_for_max_area(fired_gpkg, cal_vect, working_crs, ...)"
    )
  }

  message(sprintf("Reading FIRED daily layer: '%s' ...", layer))

  # ---- Disable S2 for all planar ops on MODIS Sinusoidal CRS -----------------
  # Must be set before any st_area() / st_union() / st_intersects() calls.
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE), add = TRUE)

  # ---- Read layer -------------------------------------------------------------
  fired_sf <- sf::st_read(fired_gpkg, layer = layer, quiet = TRUE)

  # Report raw columns so the user can verify below
  message("  Columns: ", paste(names(fired_sf), collapse = ", "))
  message(sprintf("  Raw rows: %d | Unique IDs: %d",
                  nrow(fired_sf), length(unique(fired_sf$id))))

  # ---- Standardise column names -----------------------------------------------
  # Actual FIRED CONUS daily columns (Balch et al. 2020, README confirmed):
  #   id              — fire event identifier
  #   date            — calendar date of this daily patch
  #   daily_area_km2  — new area burned this day in km² (incremental, NOT cumulative)
  #   cum_area_km2    — cumulative burned area through this day in km²
  #   total_area_km2  — total event area in km² (repeated for all rows of event)
  #   cum_pixels      — cumulative pixel count (useful for monotonicity check)
  #   event_day       — day number within event (1 = ignition day)
  #   ignition_date   — event start date (constant within event)
  #   last_date       — event end date (constant within event)
  #   duration        — event duration days (constant within event)
  .pick <- function(df, patterns) {
    nms <- tolower(names(df))
    for (p in patterns) {
      hit <- grep(p, nms, value = FALSE)[1]
      if (!is.na(hit)) return(names(df)[hit])
    }
    NA_character_
  }

  id_col      <- .pick(fired_sf, c("^id$", "^event_id$", "^fired_id$"))
  date_col    <- .pick(fired_sf, c("^date$", "^burn_date$", "^daily_date$"))
  # Incremental daily area (km² in the Balch et al. 2020 product)
  area_km2_col <- .pick(fired_sf, c("^daily_area_km2$", "^area_km2$",
                                     "^daily_area_ha$", "^new_area_ha$"))
  # Total event area for pre-filtering (km²)
  total_km2_col <- .pick(fired_sf, c("^total_area_km2$", "^total_km2$",
                                      "^total_area_ha$"))
  # Cumulative pixel count — used for monotonicity check (integers, no rounding noise)
  cum_pix_col  <- .pick(fired_sf, c("^cum_pixels$", "^cumulative_pixels$"))

  if (is.na(id_col))
    stop("Could not find event ID column in FIRED daily layer.\n",
         "Columns present: ", paste(names(fired_sf), collapse = ", "), "\n",
         "Expected: 'id' or 'event_id'")
  if (is.na(date_col))
    stop("Could not find daily date column in FIRED layer.\n",
         "Columns present: ", paste(names(fired_sf), collapse = ", "), "\n",
         "Expected: 'date' (the specific calendar date, varying per row).\n",
         "If only 'ignition_date' is present, you have the events product,\n",
         "not the daily product (layer='", layer, "').")

  message(sprintf(
    "  Mapped: id='%s'  date='%s'  daily_area='%s'  total_area='%s'",
    id_col, date_col,
    if (!is.na(area_km2_col))  area_km2_col  else "(none — will compute from geometry)",
    if (!is.na(total_km2_col)) total_km2_col else "(none)"
  ))

  fired_sf <- fired_sf %>%
    dplyr::rename(event_id = !!id_col, date = !!date_col)
  fired_sf$event_id <- as.character(fired_sf$event_id)
  fired_sf$date     <- as.Date(fired_sf$date)

  # ---- Events-layer guard: catch wrong file even if layer name check passes ---
  # Events product: one row per fire → rows == unique IDs
  rows_per_id <- max(table(fired_sf$event_id))
  if (rows_per_id == 1L) {
    stop(
      "Each event_id appears exactly once — this is the FIRED EVENTS product\n",
      "(one row per fire, final perimeter only), not the daily product.\n\n",
      "Download the FIRED CONUS Daily GeoPackage from:\n",
      "  https://scholar.colorado.edu/concern/datasets/765372341\n",
      "  Filename: fired_conus_daily_nov2001-jan2019.gpkg\n\n",
      "Update FIRED_GPKG in Section 1 of the diagnostic script to this path.\n",
      "The events file is still useful for MaxSpreadArea calibration via\n",
      "load_fired_events_for_max_area()."
    )
  }

  # ---- Date filter ------------------------------------------------------------
  fired_sf <- fired_sf %>%
    dplyr::filter(!is.na(date), date >= cal_start, date <= cal_end)
  message(sprintf("  After date filter [%s, %s]: %d rows", cal_start, cal_end, nrow(fired_sf)))

  if (nrow(fired_sf) == 0)
    stop("No FIRED records within the calibration date range (",
         cal_start, " to ", cal_end, ").")

  # ---- Pre-filter by event total size (attribute-based, fast) -----------------
  # Skip events whose total area is below min_fire_ha before the expensive
  # spatial filter and cumulative union loop.
  if (!is.na(total_km2_col)) {
    min_km2   <- min_fire_ha / 100          # ha → km²
    n_before  <- nrow(fired_sf)
    fired_sf  <- fired_sf %>%
      dplyr::filter(.data[[total_km2_col]] >= min_km2)
    message(sprintf(
      "  Pre-filter (total_area >= %.0f ha): kept %d of %d rows.",
      min_fire_ha, nrow(fired_sf), n_before
    ))
    if (nrow(fired_sf) == 0)
      stop("No events remain after size pre-filter. Reduce FIRED_MIN_FIRE_HA.")
  }

  # ---- Spatial filter: events intersecting calibration boundary ---------------
  message("Spatial filter: events intersecting calibration boundary...")

  fired_crs <- sf::st_crs(fired_sf)
  # CRS check: MODIS Sinusoidal reports input as "unnamed" — check wkt instead
  if (is.na(fired_crs) || is.null(fired_crs$wkt) || nchar(fired_crs$wkt) < 10)
    stop("FIRED layer CRS is missing or unreadable.")

  # Reproject cal_vect to MODIS Sinusoidal for the spatial join
  cal_sf <- sf::st_as_sf(cal_vect) %>%
    sf::st_transform(fired_crs) %>%
    sf::st_union()

  hits     <- sf::st_intersects(fired_sf, cal_sf, sparse = FALSE)[, 1]
  fired_sf <- fired_sf[hits, ]
  message(sprintf("  -> %d daily records for %d events intersect boundary.",
                  nrow(fired_sf), length(unique(fired_sf$event_id))))

  if (nrow(fired_sf) == 0)
    stop("No FIRED perimeters intersect the calibration boundary.\n",
         "  Check CAL_VECT_STATES / CAL_VECT_SHP covers the FIRED data extent.")

  # ---- Detect incremental vs cumulative geometry ------------------------------
  # FIRED daily geometry is INCREMENTAL (new pixels burned that day only).
  # Use cum_pixels attribute for monotonicity check — integer counts are more
  # reliable than area-from-geometry which can vary with polygon simplification.
  if (!is.na(cum_pix_col)) {
    sample_events <- head(unique(fired_sf$event_id), 200)
    check_df <- fired_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(event_id %in% sample_events) %>%
      dplyr::arrange(event_id, date) %>%
      dplyr::group_by(event_id) %>%
      dplyr::summarise(
        n_days   = dplyr::n(),
        monotone = all(diff(.data[[cum_pix_col]]) >= 0),
        .groups  = "drop"
      ) %>%
      dplyr::filter(n_days > 1)

    pct_monotone <- if (nrow(check_df) > 0) mean(check_df$monotone, na.rm = TRUE) else 0
    is_cumulative <- pct_monotone > 0.75
    message(sprintf(
      "  Monotonicity check via cum_pixels (%d multi-day events): %.0f%% non-decreasing -> %s",
      nrow(check_df), 100 * pct_monotone,
      if (is_cumulative) "CUMULATIVE perimeters (use directly)"
      else               "INCREMENTAL patches (building cumulative union)"
    ))
  } else {
    # Fallback: compute from geometry area
    fired_sf$geom_area_ha_tmp <- as.numeric(sf::st_area(fired_sf)) / 10000
    sample_events <- head(unique(fired_sf$event_id), 200)
    check_df <- fired_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(event_id %in% sample_events) %>%
      dplyr::arrange(event_id, date) %>%
      dplyr::group_by(event_id) %>%
      dplyr::summarise(
        n_days   = dplyr::n(),
        monotone = all(diff(geom_area_ha_tmp) >= -1),
        .groups  = "drop"
      ) %>%
      dplyr::filter(n_days > 1)
    pct_monotone  <- if (nrow(check_df) > 0) mean(check_df$monotone, na.rm = TRUE) else 0
    is_cumulative <- pct_monotone > 0.75
  }

  # ---- Build cumulative perimeters (incremental data only) --------------------
  if (!is_cumulative) {

    # ---- Resolve core count ----------------------------------------------------
    # n_cores = NULL  → auto (all physical cores minus 1, min 1)
    # n_cores = 1     → sequential (safe fallback)
    # n_cores = N     → use exactly N cores (capped at available - 1)
    avail_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
    if (is.null(n_cores)) n_cores <- avail_cores
    n_cores <- max(1L, min(as.integer(n_cores), avail_cores))

    all_events <- unique(fired_sf$event_id)
    n_ev       <- length(all_events)

    message(sprintf(
      "Building cumulative perimeters from incremental daily patches...\n  %d events | %d cores (PSOCK, cross-platform)",
      n_ev, n_cores
    ))

    # Pre-split sf into a named list — each worker receives only its own events,
    # not the full fired_sf, keeping per-worker memory overhead low.
    fired_split <- split(fired_sf, fired_sf$event_id)

    # Per-event worker function: builds the cumulative union across days and
    # returns a plain data.frame with a WKT geometry column (safe to serialize
    # across multisession workers).
    .process_event <- function(ev) {
      sf::sf_use_s2(FALSE)   # set per-worker — multisession workers are fresh R sessions
      ev       <- ev[order(ev$date), ]
      geoms    <- sf::st_geometry(ev)
      n_days   <- length(geoms)
      cum_rows <- vector("list", n_days)
      cum_geom <- geoms[[1]]
      eid      <- ev$event_id[1]

      for (d in seq_len(n_days)) {
        if (d > 1) cum_geom <- sf::st_union(cum_geom, geoms[[d]])
        cum_rows[[d]] <- data.frame(
          event_id      = eid,
          date          = ev$date[d],
          # MODIS Sinusoidal is equal-area (m²); divide by 10 000 for ha.
          cumul_area_ha = as.numeric(sf::st_area(cum_geom)) / 10000,
          wkt           = sf::st_as_text(cum_geom),
          stringsAsFactors = FALSE
        )
      }
      dplyr::bind_rows(cum_rows)
    }

    t0 <- proc.time()[["elapsed"]]

    if (n_cores == 1L) {

      # ---- Sequential fallback (n_cores=1 or for debugging) -------------------
      cum_list <- lapply(fired_split, .process_event)

    } else {

      # ---- future + furrr — cross-platform parallel (Mac / Windows / Linux) ---
      if (!requireNamespace("future", quietly = TRUE) ||
          !requireNamespace("furrr",  quietly = TRUE))
        stop("Packages 'future' and 'furrr' are required for parallel processing.\n",
             "Install with: install.packages(c('future', 'furrr'))")

      # Save and restore the caller's plan so we don't clobber it
      old_plan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(old_plan), add = TRUE)

      cum_list <- furrr::future_map(
        fired_split,
        .process_event,
        .options = furrr::furrr_options(seed = NULL)   # no RNG needed
      )
    }

    elapsed <- round(proc.time()[["elapsed"]] - t0)
    message(sprintf("  Cumulative union complete: %d events in %d s (~%.1f min).",
                    n_ev, elapsed, elapsed / 60))

    # Check for worker errors
    err_idx <- which(vapply(cum_list, inherits, logical(1), "try-error"))
    if (length(err_idx) > 0) {
      warning(sprintf(
        "%d events failed in parallel cumulative union (event IDs: %s). Dropping.",
        length(err_idx),
        paste(names(fired_split)[err_idx], collapse = ", ")
      ))
      cum_list <- cum_list[-err_idx]
    }

    cum_df       <- dplyr::bind_rows(cum_list)
    cum_geom_sfc <- sf::st_as_sfc(cum_df$wkt, crs = fired_crs)
    fired_sf <- sf::st_sf(
      event_id     = cum_df$event_id,
      date         = as.Date(cum_df$date),
      geom_area_ha = cum_df$cumul_area_ha,
      geometry     = cum_geom_sfc
    )

  } else {
    # Already cumulative — compute area from geometry for consistent schema
    fired_sf$geom_area_ha <- as.numeric(sf::st_area(fired_sf)) / 10000
    fired_sf <- fired_sf %>% dplyr::arrange(event_id, date)
  }

  # ---- Min fire size filter (post-union, on cumulative area) ------------------
  final_size <- fired_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(max_ha = max(geom_area_ha, na.rm = TRUE), .groups = "drop")

  keep_ids <- final_size$event_id[final_size$max_ha >= min_fire_ha]

  n_before <- length(unique(fired_sf$event_id))
  fired_sf <- fired_sf %>% dplyr::filter(event_id %in% keep_ids)
  message(sprintf(
    "Min fire size filter (>= %.0f ha): kept %d of %d events (%d daily records).",
    min_fire_ha, length(keep_ids), n_before, nrow(fired_sf)
  ))

  # ---- Minimum detection days filter ------------------------------------------
  # Require each fire event to have at least min_detection_days MODIS detections
  # (after cumulative union). A fire detected on only 1-2 days produces at most
  # 1 sequential pair, which is statistically fragile and more likely to be a
  # detection artifact than a real multi-day fire progression.
  if (!is.null(min_detection_days) && min_detection_days > 1L) {
    det_counts <- fired_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(event_id) %>%
      dplyr::summarise(n_det = dplyr::n(), .groups = "drop")

    keep_ids_det <- det_counts$event_id[det_counts$n_det >= as.integer(min_detection_days)]
    n_before_det <- length(unique(fired_sf$event_id))
    fired_sf     <- fired_sf %>% dplyr::filter(event_id %in% keep_ids_det)
    message(sprintf(
      "Min detection days filter (>= %d days): kept %d of %d events (%d daily records).",
      as.integer(min_detection_days),
      length(keep_ids_det), n_before_det, nrow(fired_sf)
    ))
    if (nrow(fired_sf) == 0)
      stop(sprintf(
        "No events remain after min_detection_days >= %d filter. ",
        "Reduce fired_min_detection_days.", as.integer(min_detection_days)
      ))
  }

  # ---- Standardise to GeoMAC output schema -----------------------------------
  fired_sf$FIRE_ID        <- paste0("FIRED-", fired_sf$event_id)
  fired_sf$perimeter_date <- as.Date(fired_sf$date)
  fired_sf$AREA_HA        <- fired_sf$geom_area_ha
  fired_sf$row_id         <- seq_len(nrow(fired_sf))

  # Reproject to working_crs
  result_v <- terra::vect(fired_sf[, c("FIRE_ID", "perimeter_date",
                                        "AREA_HA", "row_id")])
  result_v <- terra::project(result_v, working_crs)

  message(sprintf(
    "FIRED: %d cumulative perimeters for %d fire events.",
    nrow(result_v), length(unique(result_v$FIRE_ID))
  ))

  result_v
}


# =============================================================================
# load_fired_events_for_max_area
#
# Uses the FIRED EVENTS product (fired_events_conus_*.gpkg — one row per fire,
# final perimeter only) to build training data for MaxSpreadArea calibration.
#
# The MaxSpreadArea SCF parameter is:
#   MaxSpreadArea = B0 + B1*FWI + B2*EWS  [ha]
# Calibration target: 75th percentile of observed max daily growth per fire.
#
# FIRED events provide max_growth_ha (maximum single-day burned area) and
# max_growth_date (the date it occurred), which is exactly what is needed:
# pair each fire's max_growth_ha with FWI and EWS on max_growth_date.
#
# This function returns a tibble ready to pass to fit_max_daily_area().
#
# Arguments:
#   fired_events_gpkg  path to fired_events_conus_*.gpkg
#   cal_vect           SpatVector in working_crs
#   working_crs        target CRS (from template_r)
#   fwi_daily          tibble with columns date (Date), FWI (numeric)
#   wind_daily         tibble with columns date (Date), WindSpeed_kmh (numeric)
#   cal_start / end    Date range filter
#   min_fire_ha        minimum total_area_ha — default 1000 ha (large fires only)
# =============================================================================

load_fired_events_for_max_area <- function(fired_events_gpkg,
                                            cal_vect, working_crs,
                                            fwi_daily, wind_daily,
                                            cal_start, cal_end,
                                            min_fire_ha = 1000) {

  if (!file.exists(fired_events_gpkg))
    stop("FIRED events GeoPackage not found: ", fired_events_gpkg)

  message("Loading FIRED events for MaxSpreadArea calibration...")

  avail <- tryCatch(sf::st_layers(fired_events_gpkg)$name, error = function(e) character(0))
  lyr   <- avail[grepl("event", tolower(avail))][1]
  if (is.na(lyr)) lyr <- avail[1]

  ev_sf <- sf::st_read(fired_events_gpkg, layer = lyr, quiet = TRUE)
  message(sprintf("  Loaded %d events.", nrow(ev_sf)))

  # ---- Standardise columns ---------------------------------------------------
  ev_sf$ignition_date  <- as.Date(ev_sf$ignition_date)
  ev_sf$last_date      <- as.Date(ev_sf$last_date)
  ev_sf$max_growth_date <- as.Date(ev_sf$max_growth_date)

  # ---- Date filter -----------------------------------------------------------
  ev_sf <- ev_sf %>%
    dplyr::filter(
      !is.na(max_growth_date),
      max_growth_date >= cal_start,
      max_growth_date <= cal_end,
      total_area_ha   >= min_fire_ha,
      max_growth_ha   >  0
    )
  message(sprintf("  After date/size filter: %d events.", nrow(ev_sf)))

  # ---- Spatial filter --------------------------------------------------------
  fired_crs <- sf::st_crs(ev_sf)
  cal_sf    <- sf::st_as_sf(cal_vect) %>%
    sf::st_transform(fired_crs) %>%
    sf::st_union()

  hits  <- sf::st_intersects(ev_sf, cal_sf, sparse = FALSE)[, 1]
  ev_sf <- ev_sf[hits, ]
  message(sprintf("  After spatial filter (cal_vect): %d events.", nrow(ev_sf)))

  if (nrow(ev_sf) == 0)
    stop("No FIRED events intersect the calibration boundary.")

  # ---- Join FWI and EWS on max_growth_date -----------------------------------
  ev_df <- sf::st_drop_geometry(ev_sf) %>%
    dplyr::select(id, ignition_date, max_growth_date, total_area_ha,
                  max_growth_ha, fsr_ha_per_day) %>%
    dplyr::left_join(
      fwi_daily  %>%
        dplyr::rename_with(~"FWI", .cols = dplyr::any_of("value")) %>%
        dplyr::select(date, FWI),
      by = c("max_growth_date" = "date")
    ) %>%
    dplyr::left_join(
      wind_daily %>% dplyr::select(date, WindSpeed_kmh),
      by = c("max_growth_date" = "date")
    ) %>%
    dplyr::filter(is.finite(FWI), is.finite(WindSpeed_kmh)) %>%
    dplyr::rename(
      FIRE_ID       = id,
      date          = max_growth_date,
      daily_area_ha = max_growth_ha
    )

  message(sprintf(
    "  FWI join: %d / %d events matched.",
    sum(is.finite(ev_df$FWI)), nrow(ev_df)
  ))
  message(sprintf(
    "  daily_area_ha range: %.0f – %.0f ha  (max single-day growth)",
    min(ev_df$daily_area_ha, na.rm = TRUE),
    max(ev_df$daily_area_ha, na.rm = TRUE)
  ))

  ev_df
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
                                   max_growth_factor   = NULL,
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
      # Detection-gap artifact flag: daily increment > k × existing area is almost
      # certainly multiple days of MODIS-missed burning compressed into one record.
      # Only applied when max_growth_factor is set (NULL = no cap).
      drop_growth         = if (!is.null(max_growth_factor) &&
                                 is.finite(as.numeric(max_growth_factor)) &&
                                 as.numeric(max_growth_factor) > 0) {
                              (daily_area_ha / pmax(area_ha, 1)) >
                                as.numeric(max_growth_factor)
                            } else {
                              rep(FALSE, dplyr::n())
                            },
      use_for_spread_prob = !drop_neg & !drop_growth & gap_days <= max_gap_sp & delta_ha > 0,
      use_for_daily_area  = !drop_neg & !drop_growth & gap_days <= max_gap_da & delta_ha > 0
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

  n_growth_dropped <- sum(pairs$drop_growth, na.rm = TRUE)
  growth_msg <- if (!is.null(max_growth_factor) && n_growth_dropped > 0)
    sprintf(" | %d dropped (growth > %.1f× existing area)", n_growth_dropped, as.numeric(max_growth_factor))
  else ""

  message(sprintf(
    "Pairs: %d total | %d spread-prob (gap<=%d) | %d daily-area (gap<=%d)%s",
    nrow(pairs), sum(pairs$use_for_spread_prob, na.rm = TRUE), max_gap_sp,
    sum(pairs$use_for_daily_area, na.rm = TRUE), max_gap_da,
    growth_msg
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

fit_max_daily_area <- function(pairs      = NULL,
                               fired_dat  = NULL,
                               cal_months = NULL,
                               fwi_min    = NULL,
                               use_ews    = TRUE,
                               constrained = FALSE) {

  if (is.null(pairs) && is.null(fired_dat))
    stop("fit_max_daily_area: provide either 'pairs' or 'fired_dat'.")

  # ---- Assemble working dataset -----------------------------------------------
  if (!is.null(fired_dat)) {
    # fired_dat from load_fired_events_for_max_area(): FIRE_ID, date, daily_area_ha, FWI, ...
    # May include EffectiveWind if computed upstream; fall back to WindSpeed_kmh.
    dat <- fired_dat %>%
      filter(is.finite(daily_area_ha), daily_area_ha > 0, is.finite(FWI))

    if (!("EffectiveWind" %in% names(dat) && any(is.finite(dat$EffectiveWind)))) {
      if ("WindSpeed_kmh" %in% names(dat)) {
        dat <- dat %>% mutate(EffectiveWind = WindSpeed_kmh)
        message("MaxSpreadArea: no EffectiveWind in fired_dat — using WindSpeed_kmh as proxy.")
      } else {
        dat <- dat %>% mutate(EffectiveWind = NA_real_)
        message("MaxSpreadArea: no wind column found in fired_dat — EWS unavailable.")
      }
    }

    if (use_ews) dat <- dat %>% filter(is.finite(EffectiveWind))
  } else {
    dat <- pairs %>%
      filter(use_for_daily_area,
             is.finite(daily_area_ha), daily_area_ha > 0,
             is.finite(FWI))
    if (use_ews) dat <- dat %>% filter(is.finite(EffectiveWind))
  }

  # ---- Seasonal filter --------------------------------------------------------
  if (!is.null(cal_months) && length(cal_months) > 0 && length(cal_months) < 12) {
    dat <- dat %>% filter(lubridate::month(date) %in% as.integer(cal_months))
    message(sprintf(
      "Seasonal filter (max area): months [%s] → %d events remain.",
      paste(sort(as.integer(cal_months)), collapse = ","), nrow(dat)
    ))
  }

  # ---- FWI minimum filter -----------------------------------------------------
  if (!is.null(fwi_min) && is.finite(fwi_min) && fwi_min > 0) {
    dat <- dat %>% filter(FWI >= fwi_min)
    message(sprintf(
      "FWI minimum filter (max area): FWI >= %.1f → %d events remain.",
      fwi_min, nrow(dat)
    ))
  }

  if (nrow(dat) < 5)
    stop("Fewer than 5 usable pairs for max daily area fit.")

  X <- if (use_ews) model.matrix(~ FWI + EffectiveWind, data = dat) else
                    model.matrix(~ FWI, data = dat)
  y <- dat$daily_area_ha

  # ---- Collinearity diagnostics ----------------------------------------------
  X_pred    <- X[, -1, drop = FALSE]
  vif_vals  <- compute_vif(X_pred)
  kappa     <- condition_number(X)
  r_fwi_ews <- if (use_ews) cor(dat$FWI, dat$EffectiveWind, use = "complete.obs") else NA_real_

  if (kappa > 100 || any(vif_vals > 10, na.rm = TRUE)) {
    if (use_ews) {
      message(sprintf(paste0(
        "⚠ MaxSpreadArea design matrix: condition number = %.1f, VIF = [FWI:%.1f, EWS:%.1f], ",
        "FWI/EWS correlation = %.3f. Constrained quantile regression enforces ",
        "B0 ≥ 50 ha, B1 ≥ 0, B2 ≥ 0 to prevent collinearity artifacts."
      ), kappa, vif_vals[1], vif_vals[2], r_fwi_ews))
    } else {
      message(sprintf(paste0(
        "⚠ MaxSpreadArea design matrix: condition number = %.1f, VIF = [FWI:%.1f]. ",
        "Constrained quantile regression enforces B0 ≥ 50 ha, B1 ≥ 0 to prevent artifacts."
      ), kappa, vif_vals[1]))
    }
  }

  # ---- Quantile regression helper --------------------------------------------
  # When constrained = TRUE:
  #   Enforces B0 >= B0_MIN (50 ha), B1 >= 0, B2 >= 0 to prevent SCF warnings
  #   at low FWI/EWS and to guard against collinearity sign-flips.
  #   Priority: quantreg constrained → nnls reparameterization → OLS clamped.
  # When constrained = FALSE (default):
  #   Plain quantile regression (or OLS fallback) — coefficients are free.
  #   Negative signs are allowed; the user can inspect and re-run with constraints
  #   if sign violations occur.
  B0_MIN <- 50.0  # ha — only enforced when constrained = TRUE

  fit_quantile <- function(tau) {

    n_coef <- if (use_ews) 3L else 2L
    form   <- if (use_ews) daily_area_ha ~ FWI + EffectiveWind else
                           daily_area_ha ~ FWI

    # ---- Unconstrained path -------------------------------------------------
    if (!constrained) {
      if (requireNamespace("quantreg", quietly = TRUE)) {
        rq_fit <- tryCatch(
          quantreg::rq(form, tau = tau, data = dat, method = "br"),
          error = function(e) NULL
        )
        if (!is.null(rq_fit)) {
          raw_coefs <- coef(rq_fit)
          coefs <- if (use_ews) {
            setNames(raw_coefs, c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"))
          } else {
            c(setNames(raw_coefs, c("B0_Intercept", "B1_FWI")), B2_EffectiveWind = 0.0)
          }
          sum_rq    <- tryCatch(summary(rq_fit, se = "nid"), error = function(e) NULL)
          se_raw    <- if (!is.null(sum_rq)) sum_rq$coefficients[, 2] else rep(NA_real_, n_coef)
          pval_raw  <- if (!is.null(sum_rq)) sum_rq$coefficients[, 4] else rep(NA_real_, n_coef)
          se_vals   <- if (use_ews) se_raw else c(se_raw, NA_real_)
          pval_vals <- if (use_ews) pval_raw else c(pval_raw, NA_real_)
          return(list(coef   = coefs,
                      se     = setNames(se_vals,   names(coefs)),
                      pval   = setNames(pval_vals, names(coefs)),
                      method = "quantreg"))
        }
      }
      # OLS fallback (unconstrained)
      q_tgt <- quantile(y, tau, na.rm = TRUE)
      wts   <- ifelse(y >= q_tgt, tau, 1 - tau) + 1e-6
      m_ols <- if (use_ews) lm(daily_area_ha ~ FWI + EffectiveWind, data = dat, weights = wts) else
                            lm(daily_area_ha ~ FWI,               data = dat, weights = wts)
      raw_coefs <- coef(m_ols)
      coefs <- if (use_ews) {
        setNames(raw_coefs, c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"))
      } else {
        c(setNames(raw_coefs, c("B0_Intercept", "B1_FWI")), B2_EffectiveWind = 0.0)
      }
      ols_t     <- broom::tidy(m_ols) %>% rename(std_error = std.error, p_value = p.value)
      se_vals   <- if (use_ews) ols_t$std_error else c(ols_t$std_error, NA_real_)
      pval_vals <- if (use_ews) ols_t$p_value   else c(ols_t$p_value, NA_real_)
      return(list(coef  = coefs,
                  se    = setNames(se_vals,   names(coefs)),
                  pval  = setNames(pval_vals, names(coefs)),
                  method = "ols"))
    }

    # ---- Constrained path (B0 >= B0_MIN, B1 >= 0, B2 >= 0) -----------------
    # Path 1: quantreg with explicit inequality constraints
    if (requireNamespace("quantreg", quietly = TRUE)) {
      R_mat <- diag(n_coef)
      r_vec <- c(B0_MIN, rep(0.0, n_coef - 1L))

      rq_fit <- tryCatch(
        quantreg::rq(form, tau = tau, data = dat, method = "br",
                     R = R_mat, r = r_vec),
        error = function(e) NULL
      )
      if (!is.null(rq_fit)) {
        raw_coefs <- coef(rq_fit)
        coefs <- if (use_ews) {
          setNames(raw_coefs, c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"))
        } else {
          c(setNames(raw_coefs, c("B0_Intercept", "B1_FWI")), B2_EffectiveWind = 0.0)
        }
        sum_rq    <- tryCatch(summary(rq_fit, se = "nid"), error = function(e) NULL)
        se_raw    <- if (!is.null(sum_rq)) sum_rq$coefficients[, 2] else rep(NA_real_, n_coef)
        pval_raw  <- if (!is.null(sum_rq)) sum_rq$coefficients[, 4] else rep(NA_real_, n_coef)
        se_vals   <- if (use_ews) se_raw else c(se_raw, NA_real_)
        pval_vals <- if (use_ews) pval_raw else c(pval_raw, NA_real_)
        return(list(coef   = coefs,
                    se     = setNames(se_vals,   names(coefs)),
                    pval   = setNames(pval_vals, names(coefs)),
                    method = "constrained_qr"))
      }
    }

    # Path 2: nnls with reparameterization
    if (requireNamespace("nnls", quietly = TRUE)) {
      y_adj  <- y - B0_MIN
      q_tgt  <- quantile(y, tau, na.rm = TRUE)
      wts    <- sqrt(ifelse(y >= q_tgt, tau, 1 - tau) + 1e-6)
      fit_nn <- tryCatch(nnls::nnls(X * wts, y_adj * wts), error = function(e) NULL)
      if (!is.null(fit_nn)) {
        alpha <- fit_nn$x
        coefs <- if (use_ews) {
          setNames(alpha + c(B0_MIN, 0, 0), c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"))
        } else {
          c(setNames(alpha + c(B0_MIN, 0), c("B0_Intercept", "B1_FWI")), B2_EffectiveWind = 0.0)
        }
        return(list(coef   = coefs,
                    se     = rep(NA_real_, 3),
                    pval   = rep(NA_real_, 3),
                    method = "constrained_nnls"))
      }
    }

    # Path 3: weighted OLS with post-hoc clamping (last resort)
    q_tgt <- quantile(y, tau, na.rm = TRUE)
    wts   <- ifelse(y >= q_tgt, tau, 1 - tau) + 1e-6
    m_ols <- if (use_ews) lm(daily_area_ha ~ FWI + EffectiveWind, data = dat, weights = wts) else
                          lm(daily_area_ha ~ FWI,               data = dat, weights = wts)
    raw_coefs <- coef(m_ols)
    coefs <- if (use_ews) {
      setNames(raw_coefs, c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"))
    } else {
      c(setNames(raw_coefs, c("B0_Intercept", "B1_FWI")), B2_EffectiveWind = 0.0)
    }
    any_clamped <- coefs["B1_FWI"] < 0 || (use_ews && coefs["B2_EffectiveWind"] < 0) ||
                   coefs["B0_Intercept"] < B0_MIN
    if (any_clamped) {
      warning(paste0(
        "MaxSpreadArea: OLS produced constraint violations. Clamping to constraints. ",
        "Install 'quantreg' for proper constrained fitting."
      ))
      coefs["B0_Intercept"]     <- max(B0_MIN, coefs["B0_Intercept"])
      coefs["B1_FWI"]           <- max(0.0,    coefs["B1_FWI"])
      coefs["B2_EffectiveWind"] <- max(0.0,    coefs["B2_EffectiveWind"])
    }
    ols_t     <- broom::tidy(m_ols) %>% rename(std_error = std.error, p_value = p.value)
    se_vals   <- if (use_ews) ols_t$std_error else c(ols_t$std_error, NA_real_)
    pval_vals <- if (use_ews) ols_t$p_value   else c(ols_t$p_value, NA_real_)
    return(list(coef  = coefs,
                se    = setNames(se_vals,   names(coefs)),
                pval  = setNames(pval_vals, names(coefs)),
                method = if (any_clamped) "ols_clamped" else "ols"))
  }

  # ---- Three quantile fits --------------------------------------------------
  fit_med <- fit_quantile(0.50)
  fit_q75 <- fit_quantile(0.75)
  fit_q90 <- fit_quantile(0.90)

  method_used <- fit_q75$method
  message(sprintf(
    "MaxSpreadArea fit: n=%d events | range [%.1f, %.1f] ha | mean=%.0f ha | method=%s | cond=%.1f | use_ews=%s",
    nrow(dat), min(y), max(y), mean(y),
    method_used, kappa, use_ews
  ))
  message("  Recommended target: 75th percentile fit (MaxSpreadArea is a ceiling, not a mean)")

  # ---- Build tidy coefficient table -----------------------------------------
  method_note <- function(m) switch(m,
    quantreg         = "Unconstrained quantile regression (quantreg). Coefficients are free — check signs.",
    ols              = "Unconstrained weighted OLS (install quantreg for QR). Coefficients are free — check signs.",
    constrained_qr   = "Constrained quantile regression (quantreg, B0≥50 ha, B1≥0, B2≥0).",
    constrained_nnls = "Constrained nnls (B0≥50 ha, B1≥0, B2≥0); SE not available.",
    ols_clamped      = "OLS with clamping (install quantreg for proper constrained QR); SE from OLS before clamping.",
    "OLS."
  )

  build_coef_row <- function(fit_result, label) {
    tibble(
      term          = c("B0_Intercept", "B1_FWI", "B2_EffectiveWind"),
      estimate      = as.numeric(fit_result$coef),
      std_error     = as.numeric(fit_result$se),
      statistic     = as.numeric(fit_result$coef /
                        pmax(abs(fit_result$se), .Machine$double.eps)),
      p_value       = as.numeric(fit_result$pval),
      fit_target    = label,
      constrained   = fit_result$method %in% c("constrained_qr", "constrained_nnls",
                                                "ols_clamped"),
      method        = fit_result$method,
      scf_parameter = c("MaximumSpreadAreaB0", "MaximumSpreadAreaB1",
                        "MaximumSpreadAreaB2"),
      note          = paste0(method_note(fit_result$method),
                             " Raw linear (ha/day) — enter directly in SCF parameter file.")
    )
  }

  coef_df <- bind_rows(
    build_coef_row(fit_med, "Median (QR)"),
    build_coef_row(fit_q75, "75th percentile (recommended)"),
    build_coef_row(fit_q90, "90th percentile")
  )

  list(coef       = coef_df,
       model_mean = NULL,
       model_q75  = NULL,
       model_q90  = NULL,
       data       = dat,
       cond_num   = kappa,
       vif        = vif_vals,
       fwi_ews_cor = r_fwi_ews)
}


# -----------------------------------------------------------------------------
# fit_spread_probability
# Logistic model from rasterized perimeter pairs.
# -----------------------------------------------------------------------------

fit_spread_probability <- function(pairs2, geomac_v, template_r, park_mask,
                                    fine_fuels_r     = NULL,
                                    failure_ratio    = 3, max_samples = 25000,
                                    dil_kernel       = NULL,
                                    min_pct_proximal = 0,
                                    n_cores          = NULL,
                                    progress_fn      = NULL,
                                    cal_months       = NULL,
                                    fwi_min          = NULL,
                                    predictors       = c("FWI", "FineFuels", "EWS")) {

  # Default dilation kernel: 3x3 all-ones (rook-equivalent for max focal)
  if (is.null(dil_kernel)) dil_kernel <- matrix(1, 3, 3)

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

  # ---- Seasonal filter --------------------------------------------------------
  if (!is.null(cal_months) && length(cal_months) > 0 && length(cal_months) < 12) {
    eligible <- eligible %>%
      filter(lubridate::month(date) %in% as.integer(cal_months))
    message(sprintf(
      "Seasonal filter (spread prob): months [%s] → %d pairs remain.",
      paste(sort(as.integer(cal_months)), collapse = ","), nrow(eligible)
    ))
  }

  # ---- FWI minimum filter -----------------------------------------------------
  if (!is.null(fwi_min) && is.finite(fwi_min) && fwi_min > 0) {
    eligible <- eligible %>% filter(FWI >= fwi_min)
    message(sprintf(
      "FWI minimum filter (spread prob): FWI >= %.1f → %d pairs remain.",
      fwi_min, nrow(eligible)
    ))
  }

  if (nrow(eligible) == 0)
    stop("No eligible perimeter pairs after seasonal/FWI filters for spread probability.")

  # ---- Write terra objects to temp files for parallel worker access -----------
  .tmp_geomac <- tempfile(fileext = ".gpkg")
  .tmp_ff     <- if (!is.null(fine_fuels_r)) tempfile(fileext = ".tif") else NULL
  terra::writeVector(geomac_v, .tmp_geomac, filetype = "GPKG", overwrite = TRUE)
  if (!is.null(.tmp_ff)) terra::writeRaster(fine_fuels_r, .tmp_ff, overwrite = TRUE)

  # Use spread_template (the actual rasterization grid, passed in as template_r)
  spread_template <- template_r
  .template_meta <- list(
    nrows = terra::nrow(spread_template), ncols = terra::ncol(spread_template),
    xmin  = terra::xmin(spread_template),  xmax  = terra::xmax(spread_template),
    ymin  = terra::ymin(spread_template),  ymax  = terra::ymax(spread_template),
    crs   = terra::crs(spread_template)
  )

  # ---- Resolve worker count ---------------------------------------------------
  avail     <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  n_workers <- if (is.null(n_cores)) avail else max(1L, min(as.integer(n_cores), avail))
  n_workers <- min(n_workers, nrow(eligible))

  # ---- Split eligible rows into chunks ----------------------------------------
  eligible$.orig_idx <- seq_len(nrow(eligible))
  chunk_ids <- cut(seq_len(nrow(eligible)), breaks = n_workers,
                   labels = FALSE, include.lowest = TRUE)
  chunks <- split(eligible, chunk_ids)

  # ---- Worker function ---------------------------------------------------------
  .process_chunk <- function(chunk_rows,
                              geomac_path, template_meta, ff_path,
                              dil_kernel, wind_col,
                              MAX_SAMPLES, FAILURE_RATIO) {
    gv   <- terra::vect(geomac_path)
    tmpl <- terra::rast(
      nrows = template_meta$nrows, ncols = template_meta$ncols,
      xmin  = template_meta$xmin,  xmax  = template_meta$xmax,
      ymin  = template_meta$ymin,  ymax  = template_meta$ymax,
      crs   = template_meta$crs
    )
    ff <- if (!is.null(ff_path)) terra::rast(ff_path) else NULL

    .rasterize_perim <- function(v, row_id_val, template) {
      vv <- v[v$row_id == row_id_val, ]
      if (nrow(vv) == 0) return(NULL)
      r <- terra::rasterize(vv, template, field = 1, touches = TRUE)
      r[!is.na(r)] <- 1
      r[is.na(r)]  <- 0
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

      r_t  <- .rasterize_perim(gv, row$row_id,      tmpl)
      r_t1 <- .rasterize_perim(gv, row$row_id_next, tmpl)

      if (is.null(r_t) || is.null(r_t1)) {
        pair_stats_l[[j]] <- dplyr::tibble(
          pair_idx = i, FIRE_ID = row$FIRE_ID, date = row$date,
          FWI = fwi, EWS = ews, delta_ha = row$delta_ha,
          n_new_burned = 0L, n_fail = 0L, pct_proximal = NA_real_,
          status = "null_raster"
        )
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

      if (n_succ == 0) next

      succ_cap   <- as.integer(MAX_SAMPLES / (1L + FAILURE_RATIO))
      succ_cells <- if (succ_cap < n_succ) sample(succ_cells, succ_cap) else succ_cells
      n_succ     <- length(succ_cells)

      n_fail_keep <- min(n_fail, n_succ * FAILURE_RATIO)
      fail_keep   <- if (n_fail_keep > 0) sample(fail_cells, n_fail_keep) else integer(0)
      n_f         <- length(fail_keep)
      if (n_succ + n_f == 0) next

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
    }

    list(
      pair_stats = dplyr::bind_rows(pair_stats_l),
      samples    = dplyr::bind_rows(samples_l)
    )
  }

  # ---- Dispatch parallel or sequential ----------------------------------------
  if (n_workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("furrr",  quietly = TRUE))
      stop("Packages 'future' and 'furrr' required for parallel processing.\n",
           "Install with: install.packages(c('future', 'furrr'))")
    future::plan(future::multisession, workers = n_workers)
  } else {
    future::plan(future::sequential)
  }

  # ---- Parallel threshold -------------------------------------------------------
  # Require at least 4 pairs per worker to justify the process-spawn overhead.
  # GeoMAC datasets typically have few pairs — forcing sequential avoids both
  # startup cost and the externalptr serialisation issues that arise when
  # terra/GDAL handles are captured in worker closures.
  MIN_PAIRS_PER_WORKER <- 4L
  if (nrow(eligible) < n_workers * MIN_PAIRS_PER_WORKER) {
    n_workers <- 1L
    message(sprintf(
      "fit_spread_probability: only %d eligible pairs — running sequentially (threshold: %d pairs/worker).",
      nrow(eligible), MIN_PAIRS_PER_WORKER
    ))
    # Re-split into a single chunk
    chunks <- list(eligible)
    future::plan(future::sequential)
  }

  # ---- Isolate worker function from terra externalptr objects ------------------
  # furrr's global scanner walks up the full parent-environment chain.
  # .process_chunk and the dispatch lambda are both defined inside
  # fit_spread_probability, so their default parent env is this frame, which
  # contains geomac_v (a multi-GB SpatVector with GDAL externalptr handles that
  # cannot cross process boundaries).
  #
  # Fix: set parent = baseenv() on BOTH function environments.  baseenv() has
  # no variables, so the scanner stops there and never reaches fit_spread_probability.
  # Workers still work correctly: every value arrives as an explicit argument.
  environment(.process_chunk) <- new.env(parent = baseenv())

  # Build the dispatch environment manually (not via local(), which uses
  # parent.frame() and would reconnect the chain).
  .dispatch_env <- new.env(parent = baseenv())
  .dispatch_env$f  <- .process_chunk
  .dispatch_env$gp <- .tmp_geomac
  .dispatch_env$tm <- .template_meta
  .dispatch_env$fp <- .tmp_ff
  .dispatch_env$dk <- dil_kernel
  .dispatch_env$wc <- wind_col
  .dispatch_env$ms <- max_samples
  .dispatch_env$fr <- failure_ratio

  .dispatch_fn <- function(chunk) f(
    chunk_rows    = chunk,
    geomac_path   = gp,
    template_meta = tm,
    ff_path       = fp,
    dil_kernel    = dk,
    wind_col      = wc,
    MAX_SAMPLES   = ms,
    FAILURE_RATIO = fr
  )
  environment(.dispatch_fn) <- .dispatch_env

  message(sprintf("fit_spread_probability: %d pairs across %d worker(s)",
                  nrow(eligible), n_workers))

  results <- if (n_workers > 1L) {
    furrr::future_map(
      chunks,
      .dispatch_fn,
      .options = furrr::furrr_options(seed = NULL)
    )
  } else {
    lapply(chunks, .dispatch_fn)
  }

  future::plan(future::sequential)
  unlink(c(.tmp_geomac, .tmp_ff))

  # ---- Reassemble in original pair order --------------------------------------
  pair_diag   <- dplyr::bind_rows(lapply(results, `[[`, "pair_stats")) %>%
    dplyr::arrange(pair_idx)
  all_samples <- dplyr::bind_rows(lapply(results, `[[`, "samples")) %>%
    dplyr::arrange(pair_idx)

  # ---- Apply min_pct_proximal filter ------------------------------------------
  if (min_pct_proximal > 0 && nrow(all_samples) > 0) {
    keep_idx    <- pair_diag %>%
      dplyr::filter(status == "ok", !is.na(pct_proximal),
                    pct_proximal >= min_pct_proximal) %>%
      dplyr::pull(pair_idx)
    n_before    <- nrow(all_samples)
    p_before    <- sum(pair_diag$status == "ok")
    all_samples <- all_samples %>% dplyr::filter(pair_idx %in% keep_idx)
    message(sprintf(
      "min_pct_proximal=%d%%: kept %d / %d ok pairs (%d / %d samples)",
      min_pct_proximal, length(keep_idx), p_before, nrow(all_samples), n_before
    ))
  }

  if (nrow(all_samples) == 0)
    stop("No raster samples generated. Check that perimeters fall within template extent.")

  # ---- Fine fuels variance check --------------------------------------------
  # B2_FineFuels is only identifiable when FF varies across the sample.
  # When fine_fuels_r is NULL, FF=1.0 everywhere → zero variance → B2 unidentified.
  # When FF has real spatial variation, check whether that variation is meaningful.
  ff_var  <- var(all_samples$FineFuels, na.rm = TRUE)
  drop_b2 <- is.null(fine_fuels_r) || ff_var < 0.001 || !"FineFuels" %in% predictors
  drop_b3 <- !"EWS" %in% predictors

  if (drop_b2) {
    reason <- if (!"FineFuels" %in% predictors)
      "FineFuels excluded by user (predictor zeroed out)"
    else if (is.null(fine_fuels_r))
      "FineFuels placeholder (1.0 everywhere) — FF has zero variance"
    else
      sprintf("FineFuels variance = %.6f across samples (< 0.001 threshold)", ff_var)
    message(sprintf(paste0(
      "⚠ SpreadProbability: %s. ",
      "B2_FineFuels is unidentifiable or excluded — dropping from model and setting B2=0. ",
      "B0 from the reduced model (no FF) already absorbs the mean FF effect, ",
      "so it is the correct value to enter as SpreadProbabilityB0."
    ), reason))
  }

  if (drop_b3) {
    message("⚠ SpreadProbability: EWS excluded by user — B3_EffectiveWind set to 0.")
  }

  # ---- VIF for spread probability predictors --------------------------------
  sp_pred_cols <- c("FWI", if (!drop_b2) "FineFuels", if (!drop_b3) "EffectiveWind")
  sp_pred_cols <- sp_pred_cols[sp_pred_cols %in% names(all_samples)]
  X_sp  <- if (length(sp_pred_cols) >= 2) as.matrix(all_samples[, sp_pred_cols]) else
            matrix(NA_real_, nrow = nrow(all_samples), ncol = 1,
                   dimnames = list(NULL, sp_pred_cols[1]))
  vif_sp <- if (length(sp_pred_cols) >= 2) compute_vif(X_sp) else
            setNames(NA_real_, if (length(sp_pred_cols) == 1) sp_pred_cols else "none")

  # ---- GLM fit --------------------------------------------------------------
  # Build formula from active predictors (user-selected and auto-dropped).
  # Logit-link coefficients go DIRECTLY into the SCF parameter file.
  # SCF CanSpread(): Pspread = 1/(1+exp(-(B0 + B1*FWI + B2*FF + B3*EWS)))
  active_terms <- c(
    "FWI",
    if (!drop_b2) "FineFuels",
    if (!drop_b3) "EffectiveWind"
  )
  if (length(active_terms) == 0) active_terms <- "1"
  m_formula <- as.formula(paste("y ~", paste(active_terms, collapse = " + ")))
  m <- glm(m_formula, data = all_samples, family = binomial(link = "logit"))

  wind_label <- if (drop_b3)    "B3_EffectiveWind (dropped — zeroed)" else
                if (has_eff_wind) "B3_EffectiveWind" else
                "B3_EffectiveWind (raw wind, no terrain)"

  coef_df_raw <- broom::tidy(m) %>%
    rename(std_error = std.error, p_value = p.value)

  # Helper: extract coefficient components from fitted model (zero-fill when absent)
  .get_b <- function(term_name) {
    row <- coef_df_raw[coef_df_raw$term == term_name, ]
    if (nrow(row) == 0)
      list(est = 0.0, se = NA_real_, z = NA_real_, p = NA_real_)
    else
      list(est = row$estimate, se = row$std_error, z = row$statistic, p = row$p_value)
  }
  b0 <- .get_b("(Intercept)")
  b1 <- .get_b("FWI")
  b2 <- .get_b("FineFuels")
  b3 <- .get_b("EffectiveWind")

  b2_note <- if (drop_b2) {
    if (!"FineFuels" %in% predictors)
      "FineFuels excluded by user. SpreadProbabilityB2=0 in SCF parameter file."
    else if (is.null(fine_fuels_r))
      "FineFuels placeholder (1.0 everywhere) — FF has zero variance. SpreadProbabilityB2=0."
    else
      sprintf("FineFuels variance=%.6f (< 0.001). SpreadProbabilityB2=0.", ff_var)
  } else
    "logit-link; enter directly in SCF parameter file (no back-transform needed)"

  b3_note <- if (drop_b3)
    "EWS excluded by user. SpreadProbabilityB3=0 in SCF parameter file."
  else
    "logit-link; enter directly in SCF parameter file (no back-transform needed)"

  coef_df <- tibble(
    term          = c("B0_Intercept", "B1_FWI",
                      if (drop_b2) "B2_FineFuels (dropped — zeroed)" else "B2_FineFuels",
                      wind_label),
    estimate      = c(b0$est, b1$est, b2$est, b3$est),
    std_error     = c(b0$se,  b1$se,  b2$se,  b3$se),
    statistic     = c(b0$z,   b1$z,   b2$z,   b3$z),
    p_value       = c(b0$p,   b1$p,   b2$p,   b3$p),
    scf_parameter = c("SpreadProbabilityB0", "SpreadProbabilityB1",
                      "SpreadProbabilityB2", "SpreadProbabilityB3"),
    note          = c(
      "logit-link; enter directly in SCF parameter file (no back-transform needed)",
      "logit-link; enter directly in SCF parameter file (no back-transform needed)",
      b2_note,
      b3_note
    )
  )

  message(sprintf(
    "Spread prob fit: %d samples (%d successes, %d failures) | drop_b2=%s | drop_b3=%s | FF var=%.6f | VIF: %s",
    nrow(all_samples), sum(all_samples$y == 1), sum(all_samples$y == 0),
    drop_b2, drop_b3, ff_var,
    paste(names(vif_sp), round(vif_sp, 2), sep = "=", collapse = ", ")
  ))

  # ---- Multi-model comparison (M1–M4) ----------------------------------------
  .auc_est <- function(fit) {
    pr <- predict(fit, type = "response")
    y  <- fit$y
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    mean(sample(pr[y == 1], min(n1, 5000), replace = n1 < 5000) >
         sample(pr[y == 0], min(n0, 5000), replace = n0 < 5000))
  }
  .mcf_r2 <- function(fit) 1 - fit$deviance / fit$null.deviance

  .safe_glm <- function(formula) {
    tryCatch(
      glm(formula, data = all_samples, family = binomial(link = "logit")),
      error = function(e) NULL,
      warning = function(w) {
        withCallingHandlers(
          glm(formula, data = all_samples, family = binomial(link = "logit")),
          warning = function(w) invokeRestart("muffleWarning"))
      }
    )
  }

  m1 <- .safe_glm(y ~ FWI)
  m2 <- if (!drop_b2) .safe_glm(y ~ FWI + FineFuels) else NULL
  m3 <- .safe_glm(y ~ FWI + EffectiveWind)

  .comp_row <- function(fit, label) {
    if (is.null(fit)) return(NULL)
    cf <- coef(fit)
    tibble(
      Model    = label,
      McF_R2   = round(.mcf_r2(fit), 4),
      AUC      = round(.auc_est(fit), 4),
      B1_FWI   = if ("FWI"           %in% names(cf)) round(cf["FWI"],           5) else NA_real_,
      B2_FF    = if ("FineFuels"     %in% names(cf)) round(cf["FineFuels"],     4) else NA_real_,
      B3_EWS   = if ("EffectiveWind" %in% names(cf)) round(cf["EffectiveWind"], 5) else NA_real_,
      FWI_sign = if ("FWI"           %in% names(cf)) ifelse(cf["FWI"] > 0, "+", "-") else NA_character_,
      FF_sign  = if ("FineFuels"     %in% names(cf)) ifelse(cf["FineFuels"] > 0, "+", "-") else NA_character_,
      EWS_sign = if ("EffectiveWind" %in% names(cf)) ifelse(cf["EffectiveWind"] > 0, "+", "-") else NA_character_
    )
  }

  model_comparison <- bind_rows(
    .comp_row(m1, "M1: FWI only"),
    .comp_row(m2, "M2: FWI + FineFuels"),
    .comp_row(m3, "M3: FWI + EWS"),
    .comp_row(m,  "M4: FWI + FineFuels + EWS (primary)")
  )

  # AUC for primary model
  auc_val <- .auc_est(m)
  r2_val  <- .mcf_r2(m)
  n_pairs_val <- sum(pair_diag$status == "ok", na.rm = TRUE)

  list(coef             = coef_df,
       model            = m,
       model_comparison = model_comparison,
       pair_diag        = pair_diag,
       samples          = all_samples,
       wind_predictor   = wind_col,
       ff_placeholder   = is.null(fine_fuels_r),
       ff_var           = ff_var,
       drop_b2          = drop_b2,
       drop_b3          = drop_b3,
       predictors       = predictors,
       vif              = vif_sp,
       auc              = auc_val,
       r2               = r2_val,
       n_pairs          = n_pairs_val,
       n_obs            = nrow(all_samples))
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

  # Drop B2 when FF has no meaningful variance (same logic as fit_spread_probability)
  ff_var  <- var(all_samples$FineFuels, na.rm = TRUE)
  drop_b2 <- ff_var < 0.001

  m <- tryCatch(
    if (drop_b2)
      glm(y ~ FWI + EffectiveWind,
          data   = all_samples,
          family = binomial(link = "logit"))
    else
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
    model     = m,
    drop_b2   = drop_b2
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

# =============================================================================
# render_fire_animation
# Builds per-pair ggplot frames (PNG) and assembles a GIF for a selected fire.
# Uses annotation_raster for hillshade + geom_raster (alpha) for categories,
# so no ggnewscale dependency is needed.
# =============================================================================

render_fire_animation <- function(fire_id,
                                   pairs,
                                   geomac_v,
                                   spread_template,
                                   dil_cells   = 5L,
                                   dem_r       = NULL,
                                   fps         = 1,
                                   progress_fn = NULL) {

  pacman::p_load(ggplot2, terra, dplyr, gifski)

  # ---- Fire pairs in date order -----------------------------------------------
  fire_pairs <- pairs %>%
    dplyr::filter(FIRE_ID == !!fire_id) %>%
    dplyr::arrange(date)

  if (nrow(fire_pairs) == 0)
    stop(sprintf("No pairs found for fire ID: %s", fire_id))

  # ---- Fire perimeter geometries ---------------------------------------------
  fire_v <- geomac_v[geomac_v$FIRE_ID == fire_id, ]
  if (nrow(fire_v) == 0)
    stop(sprintf("No perimeters found for fire ID: %s", fire_id))

  # ---- Fire extent with 10% buffer -------------------------------------------
  fire_ext <- terra::ext(fire_v)
  dx       <- (fire_ext[2] - fire_ext[1]) * 0.10
  dy       <- (fire_ext[4] - fire_ext[3]) * 0.10
  fire_ext_buf <- terra::ext(fire_ext[1] - dx, fire_ext[2] + dx,
                              fire_ext[3] - dy, fire_ext[4] + dy)

  # Crop spread_template to fire extent
  fire_template <- tryCatch(
    terra::crop(spread_template, fire_ext_buf),
    error = function(e) {
      stop(sprintf(
        "Failed to crop spread_template to fire extent.\n  Fire ext: %s\n  Template ext: %s\n  Error: %s",
        paste(round(as.vector(fire_ext_buf), 1), collapse=", "),
        paste(round(as.vector(terra::ext(spread_template)), 1), collapse=", "),
        conditionMessage(e)
      ))
    }
  )
  nr <- terra::nrow(fire_template)
  nc <- terra::ncol(fire_template)

  if (nr == 0L || nc == 0L) {
    stop(sprintf(
      "Fire template has 0 cells after crop (%d rows x %d cols).\n  Fire ext: %s\n  Template ext: %s\n  CRS match: %s",
      nr, nc,
      paste(round(as.vector(fire_ext_buf), 1), collapse=", "),
      paste(round(as.vector(terra::ext(spread_template)), 1), collapse=", "),
      isTRUE(terra::same.crs(fire_v, spread_template))
    ))
  }

  xmin_t <- terra::xmin(fire_template)
  xmax_t <- terra::xmax(fire_template)
  ymin_t <- terra::ymin(fire_template)
  ymax_t <- terra::ymax(fire_template)

  message(sprintf(
    "render_fire_animation: fire=%s  pairs=%d  template=%dx%d cells  CRS_match=%s",
    fire_id, nrow(fire_pairs), nr, nc,
    isTRUE(terra::same.crs(fire_v, spread_template))
  ))

  # ---- Cell coordinate grid (must come BEFORE hillshade — shade_df uses xy_df) --
  # na.rm=FALSE is critical: spread_template is a geometry-only raster (no data
  # values), so all cells are NA.  The default na.rm=TRUE would drop every row,
  # returning a 0-row data frame and causing "replacement has N rows, data has 0".
  xy_df  <- as.data.frame(fire_template, xy = TRUE, na.rm = FALSE)[, c("x", "y")]
  cell_w <- terra::xres(fire_template)   # used by geom_tile to avoid uneven-spacing warning
  cell_h <- terra::yres(fire_template)

  # ---- Hillshade --------------------------------------------------------------
  # Build shade_df: xy_df + shade_col (gray hex).  Using terra coordinates
  # directly avoids all matrix-orientation/nativeRaster issues that caused
  # horizontal banding with annotation_raster.  terra::values() returns values
  # in the same row-major order as as.data.frame(raster, xy=TRUE), so
  # shade_df[k, "shade_col"] is always the correct color for xy_df[k, x/y].
  shade_df <- NULL
  if (!is.null(dem_r)) {
    shade_df <- tryCatch({

      dem_crop <- tryCatch({
        dem_proj <- terra::project(dem_r, terra::crs(fire_template))
        terra::crop(dem_proj, fire_ext_buf)
      }, error = function(e) NULL)
      if (is.null(dem_crop)) stop("DEM crop failed")

      # Step 1: clamp ocean/water fill values (elevatr encodes water as ≤0)
      dem_w <- terra::clamp(dem_crop, lower = -400, values = NA)

      # Step 2: iterative focal fill — grows window until NAs are gone
      for (w_sz in c(5L, 11L, 21L, 41L)) {
        n_na <- sum(is.na(terra::values(dem_w)))
        if (n_na == 0L) break
        dem_w <- tryCatch(
          terra::focal(dem_w, w = w_sz, fun = "mean", na.rm = TRUE,
                       na.policy = "only"),
          error = function(e)
            terra::focal(dem_w, w = w_sz, fun = "mean", na.rm = TRUE)
        )
      }

      # Step 3: smooth DEM (w=5 ≈ 1 km at 215 m) to dampen river-bank spikes
      dem_sm <- terra::focal(dem_w, w = 5L, fun = "mean", na.rm = TRUE)

      # Step 4: hillshade at native DEM resolution — do NOT resample to fire_template.
      # Keeping native resolution preserves terrain detail; the category layer is
      # plotted on top at fire_template resolution, so the two geom_tile layers can
      # have different cell sizes without any conflict.
      slope_r   <- terra::terrain(dem_sm, "slope",  unit = "radians")
      aspect_r  <- terra::terrain(dem_sm, "aspect", unit = "radians")
      shade_raw <- terra::shade(slope_r, aspect_r, angle = 60, direction = 315)
      # Light blur at native DEM resolution (w=3) to soften acquisition noise
      shade_res <- terra::focal(shade_raw, w = 3L, fun = "mean", na.rm = TRUE)
      sv        <- terra::values(shade_res, mat = FALSE)

      # Step 5: quantile normalization (2nd–98th percentile)
      sv_finite <- sv[is.finite(sv) & !is.na(sv)]
      if (length(sv_finite) < 10L) stop("too few finite shade values")
      sv_lo  <- quantile(sv_finite, 0.02)
      sv_hi  <- quantile(sv_finite, 0.98)
      sv_med <- median(sv_finite)
      sv[!is.finite(sv) | is.na(sv)] <- sv_med
      sv_norm <- pmax(0, pmin(1, (sv - sv_lo) / max(sv_hi - sv_lo, 1e-9)))

      # Step 6: build shade_df from native DEM coordinates (NOT fire_template xy_df).
      # terra::values() row-major order matches as.data.frame(xy=TRUE, na.rm=FALSE).
      shade_xy  <- as.data.frame(shade_res, xy = TRUE, na.rm = FALSE)[, c("x", "y")]
      shade_cw  <- terra::xres(shade_res)   # DEM cell width  (finer than fire_template)
      shade_ch  <- terra::yres(shade_res)   # DEM cell height
      sdf <- shade_xy
      sdf$shade_col  <- grDevices::gray(sv_norm)
      sdf$shade_cw   <- shade_cw
      sdf$shade_ch   <- shade_ch
      sdf

    }, error = function(e) {
      message("render_fire_animation: hillshade failed (", conditionMessage(e),
              ") — using plain gray background")
      NULL
    })
  }

  # ---- Dilation kernel --------------------------------------------------------
  dc         <- max(1L, as.integer(dil_cells))
  r_sz       <- 2L * dc + 1L
  dil_kernel <- outer(seq_len(r_sz), seq_len(r_sz),
                      function(i, j) as.integer(
                        sqrt((i - dc - 1L)^2 + (j - dc - 1L)^2) <= dc))

  # ---- Category colors (matching diagnostic script) --------------------------
  cat_hex <- c(
    burned         = "#B0B0B0",
    frontier_fail  = "#E8601C",
    new_distant    = "#F0C000",
    new_proximal   = "#4575B4"
  )
  cat_labels <- c(
    burned         = "Already burned (t)",
    frontier_fail  = "Frontier unburned (fail)",
    new_distant    = paste0("New area — distant"),
    new_proximal   = paste0("New area — proximal")
  )

  # ---- Generate frames -------------------------------------------------------
  n_frames   <- nrow(fire_pairs)
  frame_pngs <- character(n_frames)
  frame_meta <- vector("list", n_frames)
  tmp_dir    <- tempfile(pattern = "scf_anim_")
  dir.create(tmp_dir, recursive = TRUE)

  # Cumulative burned accumulator — guarantees the burned footprint only ever
  # grows across frames.  For GeoMAC (truly cumulative perimeters) this is
  # identical to using vt directly; for FIRED (where detection gaps or minor
  # perimeter redraws can cause vt to shrink slightly) this fills those gaps
  # and makes the animation show a monotonically growing fire.
  n_cells <- terra::ncell(fire_template)
  vcum    <- rep(0L, n_cells)   # nothing burned yet before frame 1

  for (i in seq_len(n_frames)) {
    if (!is.null(progress_fn)) progress_fn(i, n_frames)

    row <- fire_pairs[i, ]

    # ---- Row-id lookup (outside tryCatch so `next` works reliably) ------------
    # Use which() + integer index for robust terra SpatVector subsetting.
    # Terra stores integer columns as doubles internally; direct logical
    # subsetting (fire_v[fire_v$row_id == row$row_id, ]) can silently return
    # 0 rows due to type mismatch, causing every frame to be skipped.
    fv_ids <- as.numeric(fire_v$row_id)
    idx_t  <- which(fv_ids == as.numeric(row$row_id))[1L]
    idx_t1 <- which(fv_ids == as.numeric(row$row_id_next))[1L]

    if (is.na(idx_t) || is.na(idx_t1)) {
      message(sprintf(
        "  [frame %d] missing geometry: row_id=%s row_id_next=%s (first 10 available: %s)",
        i, row$row_id, row$row_id_next,
        paste(head(fv_ids, 10L), collapse = ", ")
      ))
      frame_meta[[i]] <- list(pair_idx = i, status = "missing_geometry",
                               date_t = row$date, date_t1 = row$date_next)
      next
    }
    v_t  <- fire_v[idx_t,  ]
    v_t1 <- fire_v[idx_t1, ]

    # ---- Rasterize + plot (tryCatch: errors log and skip frame, not crash) ----
    frame_result <- tryCatch({

    # Rasterize → pull values into R vectors immediately.
    # terra::rasterize() returns a file-backed SpatRaster; ANY subsequent terra
    # [<-/ifel on a file-backed raster throws "replacement has N rows, data has 0".
    # Extract to plain R integer vectors and work there; only push back to an
    # in-memory SpatRaster (via setValues) for the focal() dilation step.
    tmp_t  <- terra::rasterize(v_t,  fire_template, field = 1L,
                                touches = TRUE, background = 0L)
    tmp_t1 <- terra::rasterize(v_t1, fire_template, field = 1L,
                                touches = TRUE, background = 0L)

    vt_raw  <- terra::values(tmp_t,  mat = FALSE)
    vt1_raw <- terra::values(tmp_t1, mat = FALSE)

    # NA → 0 in plain R (no terra [<- needed)
    vt_raw[is.na(vt_raw)]   <- 0
    vt1_raw[is.na(vt1_raw)] <- 0

    vt  <- as.integer(vt_raw)
    vt1 <- as.integer(vt1_raw)

    # ---- Apply cumulative accumulator -----------------------------------------
    # vt_cum = union of all prior frames' vt1 with the current frame's vt.
    # This fills detection gaps and ensures burned area never visually shrinks.
    vt_cum <- pmax(vcum, vt)
    vt1_cum <- pmax(vt_cum, vt1)

    # ---- Classify cells -------------------------------------------------------
    r_new     <- (vt1_cum == 1L) & (vt_cum == 0L)

    # Dilation: build an in-memory raster (setValues → never file-backed) for focal
    r_t_mem   <- terra::setValues(fire_template, vt_cum)
    r_dil_raw <- terra::focal(r_t_mem, w = dil_kernel, fun = "max", na.rm = TRUE)
    vdil_raw  <- terra::values(r_dil_raw, mat = FALSE)
    vdil_raw[is.na(vdil_raw)] <- 0
    vdil      <- as.integer(vdil_raw)
    r_frontier     <- (vdil == 1L) & (vt_cum == 0L)

    # Priority: burned > proximal > distant > frontier_fail > NA
    cat_vec <- rep(NA_character_, length(vt_cum))
    cat_vec[r_frontier]                         <- "frontier_fail"
    cat_vec[r_new & !r_frontier]                <- "new_distant"
    cat_vec[r_new &  r_frontier]                <- "new_proximal"
    cat_vec[vt_cum == 1L]                       <- "burned"

    # ---- Advance cumulative accumulator ---------------------------------------
    # Must happen AFTER classification so this frame's new cells show as "new"
    # (not "burned") in the plot, but ARE included in the next frame's base.
    vcum <<- vt1_cum

    # ---- Statistics -----------------------------------------------------------
    n_new      <- sum(r_new)
    n_prox     <- sum(r_new & r_frontier)
    pct_prox   <- if (n_new > 0L) round(100 * n_prox / n_new, 1) else NA_real_
    n_ff       <- sum(r_frontier & !r_new)

    fwi_val    <- if ("FWI"      %in% names(row)) row$FWI      else NA_real_
    delta_val  <- if ("delta_ha" %in% names(row)) row$delta_ha else NA_real_

    frame_meta[[i]] <- list(
      pair_idx        = i,
      date_t          = row$date,
      date_t1         = row$date_next,
      fwi             = fwi_val,
      delta_ha        = delta_val,
      pct_proximal    = pct_prox,
      n_new           = n_new,
      n_proximal      = n_prox,
      n_frontier_fail = n_ff,
      status          = "ok"
    )

    # ---- Plot data frame ------------------------------------------------------
    plot_df <- xy_df
    plot_df$category <- factor(cat_vec,
                                levels = names(cat_hex),
                                labels = cat_labels)

    cat_df_nona <- plot_df[!is.na(plot_df$category), , drop = FALSE]

    # ---- Title / subtitle -----------------------------------------------------
    cum_ha    <- sum(vt1_cum) * terra::xres(fire_template) *
                   terra::yres(fire_template) / 10000  # cells → ha
    new_ha    <- sum(r_new)   * terra::xres(fire_template) *
                   terra::yres(fire_template) / 10000

    parts <- character(0)
    if (!is.na(fwi_val)) parts <- c(parts, sprintf("FWI: %.1f", fwi_val))
    parts <- c(parts, sprintf("Total: %.0f ha", cum_ha))
    parts <- c(parts, sprintf("+%.0f ha this step", new_ha))
    if (!is.na(pct_prox)) parts <- c(parts, sprintf("%.1f%% proximal", pct_prox))
    subtitle_str <- paste(parts, collapse = "  |  ")

    title_str <- sprintf(
      "Fire: %s   |   %s -> %s   (Pair %d / %d)",
      fire_id,
      format(row$date,      "%Y-%m-%d"),
      format(row$date_next, "%Y-%m-%d"),
      i, n_frames
    )

    # ---- Build ggplot ---------------------------------------------------------
    p <- ggplot2::ggplot() +
      ggplot2::coord_equal(
        xlim   = c(xmin_t, xmax_t),
        ylim   = c(ymin_t, ymax_t),
        expand = FALSE
      )

    # Hillshade basemap via geom_tile + I() (identity fill, no scale mapping).
    # shade_df uses native DEM coordinates and cell sizes (finer than fire_template),
    # so terrain detail is preserved.  The category layer on top uses fire_template
    # cell sizes — the two layers co-exist without conflict.
    if (!is.null(shade_df)) {
      p <- p + ggplot2::geom_tile(
        data    = shade_df,
        mapping = ggplot2::aes(x = x, y = y, fill = I(shade_col)),
        width   = shade_df$shade_cw[1L],
        height  = shade_df$shade_ch[1L],
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = xmin_t, xmax = xmax_t, ymin = ymin_t, ymax = ymax_t,
        fill = "gray90"
      )
    }

    # Category layer (semi-transparent so hillshade shows through).
    # geom_tile with explicit width/height avoids the "uneven horizontal intervals"
    # warning that geom_raster emits for projected rasters (slightly non-uniform spacing).
    if (nrow(cat_df_nona) > 0) {
      p <- p + ggplot2::geom_tile(
        data    = cat_df_nona,
        mapping = ggplot2::aes(x = x, y = y, fill = category),
        width   = cell_w,
        height  = cell_h,
        alpha   = 0.82
      )
    }

    p <- p +
      ggplot2::scale_fill_manual(
        values   = setNames(unname(cat_hex), cat_labels),
        name     = NULL,
        drop     = FALSE,
        na.value = NA
      ) +
      ggplot2::labs(
        title    = title_str,
        subtitle = subtitle_str,
        x = NULL, y = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position   = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = ggplot2::element_rect(fill  = "white",
                                                   colour = "gray50",
                                                   linewidth = 0.4),
        legend.key        = ggplot2::element_rect(fill = "white"),
        legend.margin     = ggplot2::margin(4, 6, 4, 6),
        legend.text       = ggplot2::element_text(size = 9),
        plot.title        = ggplot2::element_text(size = 9.5, face = "bold"),
        plot.subtitle     = ggplot2::element_text(size = 8.5),
        axis.text         = ggplot2::element_text(size = 7),
        axis.ticks.length = ggplot2::unit(2, "pt"),
        panel.grid        = ggplot2::element_blank()
      )

    # ---- Save PNG -------------------------------------------------------------
    png_path <- file.path(tmp_dir, sprintf("frame_%04d.png", i))
    ggplot2::ggsave(png_path, p,
                    width = 7, height = 6.5, dpi = 120, bg = "white")
    frame_pngs[i] <- png_path
    "ok"

    }, error = function(e) {
      message(sprintf("  [frame %d] ERROR: %s", i, conditionMessage(e)))
      frame_meta[[i]] <<- list(pair_idx = i, status = paste0("error: ", conditionMessage(e)),
                                date_t = row$date, date_t1 = row$date_next)
      "error"
    })
    if (identical(frame_result, "error")) next
  }

  # ---- Assemble GIF -----------------------------------------------------------
  valid_pngs <- frame_pngs[nchar(frame_pngs) > 0 & file.exists(frame_pngs)]
  gif_path   <- NULL

  if (length(valid_pngs) > 0) {
    gif_file <- file.path(tmp_dir,
                           paste0("fire_", gsub("[^A-Za-z0-9]", "_", fire_id), ".gif"))
    if (requireNamespace("gifski", quietly = TRUE)) {
      gifski::gifski(valid_pngs, gif_file,
                     width = 840L, height = 780L,
                     delay = 1 / max(fps, 0.1))
      gif_path <- gif_file
      message(sprintf("Animation GIF written: %s (%d frames, %.1f fps)",
                       basename(gif_file), length(valid_pngs), fps))
    } else {
      warning("Package 'gifski' not available — GIF not generated. ",
              "Install with: install.packages('gifski')")
    }
  }

  # ---- Return -----------------------------------------------------------------
  list(
    frame_pngs = valid_pngs,
    gif_path   = gif_path,
    n_frames   = length(valid_pngs),
    frame_meta = dplyr::bind_rows(lapply(frame_meta, function(x) {
      as.data.frame(lapply(x, function(v) if (is.null(v)) NA else v))
    })),
    fire_id    = fire_id,
    dil_cells  = dil_cells
  )
}
