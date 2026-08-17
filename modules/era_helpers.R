# =============================================================================
# era_helpers.R
# Functions to extract park-mean daily time series from ERA-Interim / ERA5
# NetCDF files and compute startup values / climatologies.
#
# ERA-Interim FWI data (ffmc, dmc, dc, fwi):
#   Vitolo et al. (2019) -- https://doi.org/10.5281/zenodo.1065400
#   Longitudes stored 0->360 in raw files; terra handles this transparently.
#
# ERA5 wind (u10, v10):
#   Copernicus CDS -- https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
#   Can be stored as two separate files or one combined file.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(climateR)
})

# -----------------------------------------------------------------------------
# extract_era_timeseries
#
# Crops an ERA NetCDF to the park boundary, computes daily spatial mean.
# Returns data.frame(date, value).
#
# Handles:
#   - ERA-Interim 0->360 longitude convention (terra reprojects transparently)
#   - Missing time dimension (falls back to layer index)
#   - Date filtering
# -----------------------------------------------------------------------------

extract_era_timeseries <- function(nc_path, park_vect, cal_start = NULL, cal_end = NULL) {

  if (!file.exists(nc_path)) stop("ERA NetCDF not found: ", nc_path)

  r  <- rast(nc_path)
  tt <- time(r)

  if (is.null(tt) || all(is.na(tt))) {
    stop("No time dimension detected in: ", nc_path)
  }

  dates <- as.Date(tt)

  # Filter to calibration window
  if (!is.null(cal_start)) {
    keep <- which(dates >= cal_start & (!is.null(cal_end) & dates <= cal_end |
                                        is.null(cal_end)))
    if (length(keep) == 0) stop("No ERA layers fall within calibration window.")
    r     <- r[[keep]]
    dates <- dates[keep]
  }

  # Project park boundary to match raster CRS
  park_proj <- project(park_vect, crs(r))

  # Crop + mask
  r2 <- crop(r, park_proj, snap = "out")
  r2 <- mask(r2, park_proj)

  # Spatial mean per layer
  g <- global(r2, "mean", na.rm = TRUE)

  tibble(
    date  = dates,
    value = as.numeric(g[, 1])
  )
}


# -----------------------------------------------------------------------------
# extract_era5_wind
#
# Extracts daily park-mean wind speed (m/s) and direction (degrees) from
# ERA5 u10/v10 component files.
#
# UNITS: wind speed is returned in m/s — the native ERA5 and LANDIS Climate
# Library unit.  The LANDIS climate CSV stores windspeed as m/s (confirmed
# from LTB test data); calibration must use the same units as simulation.
#
# Supports:
#   - Two separate files (u10_path and v10_path both specified)
#   - One combined file containing both u10 and v10 variables (only u10_path
#     specified; the function auto-detects both variables)
#
# Returns data.frame(date, WindSpeed_kmh, WindDir_deg)
# -----------------------------------------------------------------------------

extract_era5_wind <- function(u10_path, v10_path = NULL, park_vect,
                              cal_start = NULL, cal_end = NULL) {

  if (!file.exists(u10_path)) stop("u10 NetCDF not found: ", u10_path)

  # ---- Load u10 -------------------------------------------------------
  r_u10 <- rast(u10_path)

  # Auto-detect if this is a combined file containing both u10 and v10
  combined <- FALSE
  if (is.null(v10_path) || !nchar(v10_path)) {
    var_names <- names(r_u10)
    u_idx <- grep("u10|u_10|u-component|10u", tolower(var_names))
    v_idx <- grep("v10|v_10|v-component|10v", tolower(var_names))
    if (length(u_idx) > 0 && length(v_idx) > 0) {
      combined <- TRUE
      # Separate into u and v stacks (assumes same number of time steps)
      # terra stores multi-variable NetCDF as separate SpatRasters per variable
      r_v10 <- r_u10[[v_idx]]
      r_u10 <- r_u10[[u_idx]]
    } else {
      stop("v10_path not provided and could not detect v10 variable in: ", u10_path)
    }
  } else {
    if (!file.exists(v10_path)) stop("v10 NetCDF not found: ", v10_path)
    r_v10 <- rast(v10_path)
  }

  # ---- Time dimension -------------------------------------------------
  tt_u <- time(r_u10)
  tt_v <- time(r_v10)

  if (is.null(tt_u) || all(is.na(tt_u))) stop("No time dimension in u10 file.")
  if (is.null(tt_v) || all(is.na(tt_v))) stop("No time dimension in v10 file.")

  dates_u <- as.Date(tt_u)
  dates_v <- as.Date(tt_v)

  # Align on common dates
  common_dates <- intersect(as.character(dates_u), as.character(dates_v))
  if (length(common_dates) == 0) stop("u10 and v10 files share no common dates.")

  idx_u <- which(as.character(dates_u) %in% common_dates)
  idx_v <- which(as.character(dates_v) %in% common_dates)

  r_u10 <- r_u10[[idx_u]]
  r_v10 <- r_v10[[idx_v]]
  dates <- as.Date(common_dates)

  # ---- Filter to calibration window -----------------------------------
  if (!is.null(cal_start)) {
    keep  <- which(dates >= cal_start & (is.null(cal_end) | dates <= cal_end))
    if (length(keep) == 0) stop("No ERA5 wind layers fall within calibration window.")
    r_u10 <- r_u10[[keep]]
    r_v10 <- r_v10[[keep]]
    dates <- dates[keep]
  }

  # ---- Project park boundary ------------------------------------------
  park_proj <- project(park_vect, crs(r_u10))

  # ---- Crop + mask + spatial mean ------------------------------------
  u_crop <- crop(r_u10, park_proj, snap = "out")
  u_crop <- mask(u_crop, park_proj)
  u_mean <- as.numeric(global(u_crop, "mean", na.rm = TRUE)[, 1])

  v_crop <- crop(r_v10, park_proj, snap = "out")
  v_crop <- mask(v_crop, park_proj)
  v_mean <- as.numeric(global(v_crop, "mean", na.rm = TRUE)[, 1])

  # ---- Derive speed and direction -------------------------------------
  # Speed: sqrt(u^2 + v^2) converted m/s -> km/h (* 3.6).
  # SCF User Guide (p.5, sections 2.43/2.47/2.51):
  #   "The climate library converts all wind speed units into kilometers / hour.
  #    Be sure to convert your wind data into the correct units when inputting
  #    into the climate library."
  # ERA5 native units are m/s; multiply by 3.6 so calibrated EWS coefficients
  # match the km/h values that SCF reads at runtime from DailyWindSpeed.
  wind_speed_kmh <- sqrt(u_mean^2 + v_mean^2) * 3.6

  # Meteorological direction: direction FROM which wind blows (degrees, 0=N)
  # atan2(u, v) gives the direction the wind is going TO; subtract 180 to get FROM
  wind_dir_deg <- (atan2(u_mean, v_mean) * 180 / pi + 180) %% 360

  tibble(
    date          = dates,
    WindSpeed_kmh = wind_speed_kmh,
    WindDir_deg   = wind_dir_deg
  )
}


# -----------------------------------------------------------------------------
# extract_gridmet_wind
#
# Fetches daily wind speed and direction from GRIDMET via OPeNDAP (no local
# files required). Data are fetched year-by-year to support progress reporting
# and graceful partial failure.
#
# GRIDMET variables used:
#   vs  -- daily mean 10m wind speed (m/s) — kept in m/s to match LANDIS
#   th  -- daily mean 10m wind direction (degrees clockwise from north)
#
# Coverage: contiguous US (CONUS) only, 1979-present, ~4 km resolution.
# Source:   https://www.climatologylab.org/gridmet.html
#
# Arguments:
#   cal_vect    SpatVector of calibration boundary (any CRS; reprojected internally)
#   cal_start   start Date
#   cal_end     end Date
#   progress_fn optional function(i, n_total, year) for progress callbacks
#
# Returns data.frame(date, WindSpeed_kmh, WindDir_deg)
# -----------------------------------------------------------------------------

extract_gridmet_wind <- function(cal_vect, cal_start, cal_end,
                                  progress_fn = NULL) {

  if (!requireNamespace("climateR", quietly = TRUE))
    stop("climateR package required. Install with: install.packages('climateR')")

  # GRIDMET uses WGS84 geographic coordinates
  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  cal_sf    <- sf::st_as_sf(cal_wgs84)

  years   <- seq(year(cal_start), year(cal_end))
  results <- vector("list", length(years))

  for (i in seq_along(years)) {
    yr <- years[i]
    if (!is.null(progress_fn)) progress_fn(i, length(years), yr)

    yr_start <- max(cal_start, as.Date(sprintf("%d-01-01", yr)))
    yr_end   <- min(cal_end,   as.Date(sprintf("%d-12-31", yr)))

    gm <- tryCatch(
      climateR::getGridMET(
        AOI       = cal_sf,
        varname   = c("vs", "th"),
        startDate = as.character(yr_start),
        endDate   = as.character(yr_end)
      ),
      error = function(e) {
        warning(sprintf("GRIDMET fetch failed for %d: %s", yr, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(gm)) next

    # climateR variable names vary by version — match flexibly by pattern
    nm <- names(gm)
    message(sprintf("GRIDMET year %d: returned variables: %s", yr, paste(nm, collapse = ", ")))

    vs_idx <- grep("^vs$|wind.?speed|speed", nm, ignore.case = TRUE)[1]
    th_idx <- grep("^th$|wind.?dir|direction|bearing", nm, ignore.case = TRUE)[1]

    if (is.na(vs_idx))
      stop("Cannot identify wind speed layer in GRIDMET response. ",
           "Available names: ", paste(nm, collapse = ", "))
    if (is.na(th_idx))
      stop("Cannot identify wind direction layer in GRIDMET response. ",
           "Available names: ", paste(nm, collapse = ", "))

    vs_r <- gm[[vs_idx]]
    th_r <- gm[[th_idx]]

    # Mask to actual boundary polygon (climateR returns bounding-box crop).
    # Reproject the boundary to exactly match the raster CRS to avoid
    # minor WGS84 string mismatches between terra and climateR.
    cal_v  <- project(vect(cal_sf), crs(vs_r))
    vs_r   <- mask(vs_r, cal_v)
    th_r   <- mask(th_r, cal_v)

    dates <- tryCatch(as.Date(time(vs_r)),
                      error = function(e) seq(yr_start, yr_end, by = "day"))

    vs_mean <- as.numeric(global(vs_r, "mean", na.rm = TRUE)[, 1])
    th_mean <- as.numeric(global(th_r, "mean", na.rm = TRUE)[, 1])

    results[[i]] <- tibble(
      date          = dates,
      WindSpeed_kmh = vs_mean * 3.6,   # GRIDMET vs is m/s; convert to km/h per SCF User Guide
      WindDir_deg   = th_mean
    )
  }

  out <- bind_rows(results)

  if (nrow(out) == 0)
    stop(paste0(
      "No GRIDMET wind data returned for the specified region and date range.\n",
      "  Check that the calibration boundary is within the contiguous US (CONUS)\n",
      "  and that the GRIDMET THREDDS server is reachable."
    ))

  out %>%
    filter(date >= cal_start, date <= cal_end) %>%
    arrange(date)
}


# =============================================================================
# fetch_gridmet_met_year   (internal helper for compute_gridmet_fwi)
#
# Fetches the four GRIDMET variables needed by the Canadian FWI System for a
# single year, aggregates to a coarser grid, and returns a long data.frame of
# per-cell daily values.  Results are cached to disk as .rds so an interrupted
# multi-year fetch can resume without re-downloading.
#
# GRIDMET -> FWI input mapping (US daily gridded fire-danger convention):
#   tmmx (daily max temperature, K)      -> temp (deg C)
#   rmin (daily min relative humidity %) -> rh   (%)
#   vs   (daily mean 10 m wind, m/s)     -> ws   (km/h)   [x 3.6]
#   pr   (daily precipitation, mm)       -> prec (mm)
#
# Max-temp / min-RH approximates mid-afternoon peak burning conditions rather
# than true local noon.  This runs slightly hotter than a noon-referenced FWI
# but better represents the peak-of-day behaviour SCF's spread logic models.
# -----------------------------------------------------------------------------

fetch_gridmet_met_year <- function(cal_sf, yr, yr_start, yr_end,
                                    agg_fact  = 6L,
                                    cache_dir = NULL,
                                    aoi_key   = "aoi") {

  cache_file <- if (!is.null(cache_dir)) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    file.path(cache_dir, sprintf("gridmet_met_%s_%d.rds", aoi_key, yr))
  } else NULL

  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    # Reject caches that are unusable: all-NA dates (pre time-axis fix), or the
    # older 4-variable schema that stored derived temp/rh instead of the raw
    # tmax/tmin/rmin/rmax needed to switch FWI convention.
    need_cols <- c("tmax", "tmin", "rmin", "rmax", "ws", "prec")
    if (!is.null(cached) && nrow(cached) > 0 && !all(is.na(cached$date)) &&
        all(need_cols %in% names(cached))) {
      message(sprintf("  GRIDMET %d: loaded from cache (%s)", yr, basename(cache_file)))
      return(cached)
    }
    reason <- if (is.null(cached) || nrow(cached) == 0) "unreadable"
              else if (all(is.na(cached$date))) "bad dates"
              else "old 4-variable schema"
    message(sprintf("  GRIDMET %d: cache invalid (%s) — re-fetching", yr, reason))
    unlink(cache_file)
  }

  # Fetch each variable in its own request with retry.  The GRIDMET THREDDS
  # server intermittently returns HTTP 500 on individual aggregate files
  # (e.g. "readAgg agg_met_pr_1979_CurrentYear_CONUS.nc"); requesting all four
  # at once means one flaky variable kills the whole year, which is costly on a
  # multi-hour multi-year pull.  climateR returns long descriptive layer names
  # rather than the short variable codes, so match by pattern.
  fetch_var <- function(vn, pat, label, n_try = 3L) {
    for (k in seq_len(n_try)) {
      r <- tryCatch(
        climateR::getGridMET(AOI = cal_sf, varname = vn,
                             startDate = as.character(yr_start),
                             endDate   = as.character(yr_end)),
        error = function(e) {
          message(sprintf("    %s attempt %d/%d failed: %s",
                          vn, k, n_try, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(r)) {
        i <- grep(pat, names(r), ignore.case = TRUE)[1]
        if (is.na(i)) i <- 1L
        return(r[[i]])
      }
      if (k < n_try) Sys.sleep(5 * k)   # linear backoff
    }
    stop(sprintf("GRIDMET %s (%s) unavailable for %d after %d attempts.",
                 label, vn, yr, n_try))
  }

  # Fetch BOTH temperature and BOTH humidity extremes so the cache can serve
  # either FWI convention (see compute_gridmet_fwi) without re-downloading:
  #   "landis" — mean temp / mean RH, matching the LANDIS Climate Library
  #   "peak"   — max temp / min RH, the US gridded fire-danger convention
  r_tmax <- fetch_var("tmmx", "^tmmx$|maximum.*temp",        "max temperature")
  r_tmin <- fetch_var("tmmn", "^tmmn$|minimum.*temp",        "min temperature")
  r_rmin <- fetch_var("rmin", "^rmin$|minimum.*humid",       "min relative humidity")
  r_rmax <- fetch_var("rmax", "^rmax$|maximum.*humid",       "max relative humidity")
  r_ws   <- fetch_var("vs",   "^vs$|wind.?speed|mean.*wind", "wind speed")
  r_pr   <- fetch_var("pr",   "^pr$|precip",                 "precipitation")
  r_temp <- r_tmax   # reference layer for the time axis / geometry

  # Capture the time axis BEFORE aggregating — terra::aggregate() does not carry
  # the time dimension through, which silently yields all-NA dates downstream.
  n_lyr <- terra::nlyr(r_temp)
  dates <- tryCatch(as.Date(terra::time(r_temp)), error = function(e) NULL)
  if (is.null(dates) || all(is.na(dates)) || length(dates) != n_lyr) {
    dates <- seq(yr_start, by = "day", length.out = n_lyr)
    message(sprintf("  GRIDMET %d: time axis unusable — assuming daily from %s",
                    yr, yr_start))
  }

  # Aggregate to a coarser grid before computing FWI.  FWI is nonlinear and the
  # moisture codes carry day-to-day memory, so we preserve per-location code
  # evolution (unlike averaging the met fields first) without paying the full
  # 4 km cost.  agg_fact = 6 gives ~24 km cells.
  .agg <- function(r) if (agg_fact > 1L)
    terra::aggregate(r, fact = agg_fact, fun = "mean", na.rm = TRUE) else r
  r_tmax <- .agg(r_tmax); r_tmin <- .agg(r_tmin)
  r_rmin <- .agg(r_rmin); r_rmax <- .agg(r_rmax)
  r_ws   <- .agg(r_ws);   r_pr   <- .agg(r_pr)

  cal_v  <- project(vect(cal_sf), crs(r_tmax))
  .msk <- function(r) mask(r, cal_v)
  r_tmax <- .msk(r_tmax); r_tmin <- .msk(r_tmin)
  r_rmin <- .msk(r_rmin); r_rmax <- .msk(r_rmax)
  r_ws   <- .msk(r_ws);   r_pr   <- .msk(r_pr)

  xy <- terra::xyFromCell(r_tmax, seq_len(terra::ncell(r_tmax)))
  m_tmax <- terra::values(r_tmax); m_tmin <- terra::values(r_tmin)
  m_rmin <- terra::values(r_rmin); m_rmax <- terra::values(r_rmax)
  m_ws   <- terra::values(r_ws);   m_pr   <- terra::values(r_pr)

  # Keep only cells with complete data across all six variables
  keep <- which(
    rowSums(is.na(m_tmax)) == 0 & rowSums(is.na(m_tmin)) == 0 &
    rowSums(is.na(m_rmin)) == 0 & rowSums(is.na(m_rmax)) == 0 &
    rowSums(is.na(m_ws))   == 0 & rowSums(is.na(m_pr))   == 0
  )
  if (length(keep) == 0)
    stop(sprintf("GRIDMET %d: no cells with complete data inside boundary.", yr))

  nd <- length(dates)
  # Cache the RAW variables.  Deriving temp/rh here would bake one convention
  # into the cache and force a full re-download to switch.
  out <- tibble(
    id    = rep(keep, times = nd),
    long  = rep(xy[keep, 1], times = nd),
    lat   = rep(xy[keep, 2], times = nd),
    date  = rep(dates, each = length(keep)),
    tmax  = as.numeric(m_tmax[keep, ]) - 273.15,   # K -> deg C
    tmin  = as.numeric(m_tmin[keep, ]) - 273.15,   # K -> deg C
    rmin  = as.numeric(m_rmin[keep, ]),
    rmax  = as.numeric(m_rmax[keep, ]),
    ws    = as.numeric(m_ws[keep, ]) * 3.6,        # m/s -> km/h
    prec  = as.numeric(m_pr[keep, ])               # mm (Van Wagner expects mm)
  )

  if (!is.null(cache_file)) saveRDS(out, cache_file)
  out
}


# =============================================================================
# compute_gridmet_fwi
#
# Derives daily FFMC / DMC / DC / FWI time series from GRIDMET meteorology
# using the Canadian Forest Fire Weather Index System (cffdrs package).
#
# Motivation: the ERA-Interim FWI product (Vitolo et al. 2019) ends in 2018,
# but NIROPS perimeters run 2020-2024 and GRIDMET runs 1979-present.  Deriving
# FWI from GRIDMET gives a single internally consistent FWI definition that
# spans both the ignition record (FPA FOD, 1992-2020) and the modern perimeter
# record — which matters because SCF applies one FWI stream at runtime.
#
# Moisture codes are re-initialised each year on `reset_mmdd` (default 1 April).
# Running continuously through winter would let DC accumulate unrealistically.
#
# Returns list(FWI, FFMC, DMC, DC), each a tibble(date, value) matching the
# schema produced by extract_era_timeseries() so downstream code is unchanged.
# -----------------------------------------------------------------------------

compute_gridmet_fwi <- function(cal_vect, cal_start, cal_end,
                                 agg_fact    = 6L,
                                 reset_mmdd  = "04-01",
                                 init_ffmc   = 85,
                                 init_dmc    = 6,
                                 init_dc     = 15,
                                 convention  = c("landis", "peak"),
                                 precip_scale = 1.0,
                                 cache_dir   = NULL,
                                 progress_fn = NULL) {

  # FWI convention.  This matters more than any coefficient:
  #
  #   "landis" — daily MEAN temperature and MEAN relative humidity.  This is
  #     what the LANDIS-II Climate Library computes (ClimateDataParser.cs:314
  #     sets Temp = 0.5*(MinTemp+MaxTemp); :329 sets RH = 0.5*(MinRH+MaxRH)),
  #     and therefore what SCF actually consumes at runtime.
  #
  #   "peak" — daily MAX temperature and MIN relative humidity, the US gridded
  #     fire-danger convention approximating mid-afternoon burning conditions.
  #
  # Calibrating on "peak" and running SCF on Climate-Library FWI applies the
  # fitted coefficients to a systematically different (much higher) predictor.
  # Default to "landis" so calibration and simulation share one definition.
  convention <- match.arg(convention)

  if (!requireNamespace("climateR", quietly = TRUE))
    stop("climateR package required. Install with: install.packages('climateR')")
  if (!requireNamespace("cffdrs", quietly = TRUE))
    stop("cffdrs package required for GRIDMET-derived FWI. ",
         "Install with: install.packages('cffdrs')")

  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  cal_sf    <- sf::st_as_sf(cal_wgs84)

  # Stable short key so cached files are tied to this specific boundary
  bb      <- as.vector(terra::ext(cal_wgs84))
  aoi_key <- substr(
    paste0(format(round(bb, 3), trim = TRUE), collapse = "_"), 1, 40
  )
  aoi_key <- gsub("[^0-9A-Za-z_-]", "", aoi_key)

  years   <- seq(year(cal_start), year(cal_end))
  results <- vector("list", length(years))

  for (i in seq_along(years)) {
    yr <- years[i]
    if (!is.null(progress_fn)) progress_fn(i, length(years), yr)

    yr_start <- max(cal_start, as.Date(sprintf("%d-01-01", yr)))
    yr_end   <- min(cal_end,   as.Date(sprintf("%d-12-31", yr)))

    met <- tryCatch(
      fetch_gridmet_met_year(cal_sf, yr, yr_start, yr_end,
                             agg_fact = agg_fact, cache_dir = cache_dir,
                             aoi_key  = aoi_key),
      error = function(e) {
        warning(sprintf("GRIDMET fetch/compute failed for %d: %s",
                        yr, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(met)) next

    # Restrict to the fire season starting at reset_mmdd so the moisture codes
    # begin from the standard spring startup state each year.
    reset_date <- as.Date(sprintf("%d-%s", yr, reset_mmdd))
    met <- met %>% filter(date >= reset_date) %>% arrange(date, id)
    if (nrow(met) == 0) next

    # Derive temp / RH per the selected convention.  The cache stores raw
    # tmax/tmin/rmin/rmax so switching convention costs nothing.
    if (convention == "landis") {
      met_temp <- 0.5 * (met$tmax + met$tmin)   # LANDIS: 0.5*(MinTemp+MaxTemp)
      met_rh   <- 0.5 * (met$rmin + met$rmax)   # LANDIS: 0.5*(MinRH+MaxRH)
    } else {
      met_temp <- met$tmax                      # peak burning conditions
      met_rh   <- met$rmin
    }
    # FWI equations require RH strictly within (0, 100]
    met_rh <- pmax(pmin(met_rh, 100), 0.01)

    fwi_in <- data.frame(
      id   = met$id,   lat  = met$lat,  long = met$long,
      yr   = as.integer(format(met$date, "%Y")),
      mon  = as.integer(format(met$date, "%m")),
      day  = as.integer(format(met$date, "%d")),
      temp = met_temp, rh = met_rh, ws = met$ws,
      # precip_scale replicates a LANDIS climate file whose precipitation column
      # is in CENTIMETRES.  The Climate Library feeds Precip straight into the
      # Van Wagner rain equations, which are calibrated in MILLIMETRES and use
      # thresholds of 1.5 / 2.8 mm with no unit conversion (FireClimate.cs:95,
      # :122).  A cm file is therefore treated as 1/10 the true rainfall, which
      # inflates runtime FWI.  Set precip_scale = 0.1 to reproduce that, so the
      # calibration FWI matches what SCF actually consumes.
      prec = met$prec * precip_scale
    )

    n_cells <- length(unique(met$id))
    fwi_out <- tryCatch(
      cffdrs::fwi(
        fwi_in,
        init  = data.frame(
          ffmc = rep(init_ffmc, n_cells),
          dmc  = rep(init_dmc,  n_cells),
          dc   = rep(init_dc,   n_cells),
          lat  = unique(met[, c("id", "lat")])$lat
        ),
        batch = TRUE, out = "all", lat.adjust = TRUE, uppercase = TRUE
      ),
      error = function(e) {
        warning(sprintf("cffdrs::fwi failed for %d: %s", yr, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fwi_out)) next

    # Spatial mean per date — matches how extract_era_timeseries() reduces the
    # gridded ERA FWI product to a single regional series.
    results[[i]] <- tibble(
      date = met$date,
      FWI  = as.numeric(fwi_out$FWI),
      FFMC = as.numeric(fwi_out$FFMC),
      DMC  = as.numeric(fwi_out$DMC),
      DC   = as.numeric(fwi_out$DC)
    ) %>%
      group_by(date) %>%
      summarise(across(c(FWI, FFMC, DMC, DC), ~ mean(.x, na.rm = TRUE)),
                .groups = "drop")
  }

  all_df <- bind_rows(results)
  if (nrow(all_df) == 0)
    stop(paste0(
      "No GRIDMET-derived FWI produced for the specified region and period.\n",
      "  Check that the calibration boundary is within CONUS and that the\n",
      "  GRIDMET THREDDS server is reachable."
    ))

  all_df <- all_df %>%
    filter(date >= cal_start, date <= cal_end) %>%
    arrange(date)

  message(sprintf(
    "GRIDMET FWI: %d daily values | FWI range %.1f - %.1f | %d year(s)",
    nrow(all_df), min(all_df$FWI, na.rm = TRUE), max(all_df$FWI, na.rm = TRUE),
    length(unique(format(all_df$date, "%Y")))
  ))

  list(
    FWI  = all_df %>% transmute(date, value = FWI),
    FFMC = all_df %>% transmute(date, value = FFMC),
    DMC  = all_df %>% transmute(date, value = DMC),
    DC   = all_df %>% transmute(date, value = DC)
  )
}


# =============================================================================
# compute_landis_climate_fwi
#
# Computes FWI directly from a LANDIS-II daily climate input file, replicating
# the Climate Library's own conventions.  This is the FWI SCF actually consumes
# at runtime, so it is the reference against which any calibration FWI should
# be checked.
#
# WHY THIS EXISTS: spread and ignition coefficients are fit against whatever FWI
# series the calibration uses.  SCF then applies them to the Climate Library's
# FWI.  If the two are computed from different variables the coefficients are
# applied to a systematically different predictor and no amount of tuning fixes
# it.  This function makes that comparison possible.
#
# Expected file format (LANDIS daily climate CSV):
#   Year,Month,Day,Variable,eco1,eco2,...
# with Variable rows including tmin/tmax, ppt, and either RH, or dewpoint, or
# specific humidity.  Wind (windSpeed) is optional but improves ISI.
#
# Climate Library conventions replicated (commit a10d6c7):
#   Temp = 0.5*(MinTemp + MaxTemp)                  ClimateDataParser.cs:314
#   RH   = 0.5*(MinRH + MaxRH), else from dewpoint  ClimateDataParser.cs:329,333
#   RH from dewpoint uses the mean temp             ClimateDataParser.cs:453-462
#   Wind [m/s] * 3.6 -> [km/h]                      ClimateDataParser.cs:12
#   Codes reset each year from the config startup   FireClimate.cs:19-21
#   Only days SpringStart..WinterStart are computed FireClimate.cs:26-30
#
# PRECIPITATION UNITS — read carefully, the sign of this is easy to get wrong.
#   The Van Wagner rain equations expect MILLIMETRES.  LANDIS climate files
#   commonly supply CENTIMETRES (ClimateRecord.Precip is documented [cm]).
#   An UNPATCHED Climate Library passes the raw column straight through, so a
#   cm file is silently read as 1/10 the true rainfall.
#
#   precip_units    "cm" or "mm" — what the FILE contains.
#   library_patched TRUE  = compute physically correct FWI (converts cm -> mm).
#                   FALSE = replicate the unpatched library, i.e. feed the raw
#                           column in regardless of units.
#
#   Resulting multiplier applied to the file's values:
#       cm file + patched   -> x10   (correct)
#       cm file + unpatched -> x1    (reproduces the bug)
#       mm file, either     -> x1    (no bug to reproduce)
#
# NOTE: cffdrs implements the FFMC rainfall routine, which the unpatched
# Climate Library omits entirely.  Against an unpatched library this function
# will therefore read somewhat LOWER than the library's true output.  Against a
# patched library the two should agree.
#
# Returns list(FWI, FFMC, DMC, DC), each tibble(date, value) — same schema as
# extract_era_timeseries() so it drops into the rest of the app unchanged.
# =============================================================================

compute_landis_climate_fwi <- function(csv_path,
                                        ecoregions   = NULL,
                                        init_ffmc    = 85,
                                        init_dmc     = 6,
                                        init_dc      = 15,
                                        first_day    = 1,
                                        last_day     = 365,
                                        precip_units    = c("cm", "mm"),
                                        library_patched = TRUE,
                                        lat          = 46.85,
                                        progress_fn  = NULL) {

  precip_units <- match.arg(precip_units)
  precip_scale <- if (isTRUE(library_patched) && precip_units == "cm") 10.0 else 1.0

  if (!file.exists(csv_path)) stop("LANDIS climate file not found: ", csv_path)
  if (!requireNamespace("cffdrs", quietly = TRUE))
    stop("cffdrs package required. install.packages('cffdrs')")

  rd <- if (requireNamespace("data.table", quietly = TRUE))
    data.table::fread(csv_path, showProgress = FALSE) else
    utils::read.csv(csv_path, stringsAsFactors = FALSE)
  rd <- as.data.frame(rd)

  nm <- names(rd)
  .find <- function(pat) nm[grepl(pat, nm, ignore.case = TRUE)][1]
  c_year <- .find("^year$"); c_mon <- .find("^month$")
  c_day  <- .find("^day$");  c_var <- .find("^variable$")
  if (any(is.na(c(c_year, c_mon, c_day, c_var))))
    stop("Expected Year, Month, Day and Variable columns. Found: ",
         paste(head(nm, 12), collapse = ", "))

  eco_cols <- nm[grepl("^eco", nm, ignore.case = TRUE)]
  if (length(eco_cols) == 0) stop("No ecoregion columns (eco1, eco2, ...) found.")
  if (!is.null(ecoregions) && length(ecoregions) > 0) {
    sel <- eco_cols[tolower(eco_cols) %in% tolower(ecoregions)]
    if (length(sel) > 0) eco_cols <- sel
  }
  if (!is.null(progress_fn))
    progress_fn(sprintf("Averaging %d ecoregion column(s)", length(eco_cols)))

  # Landscape mean across the selected ecoregions, then one row per day/variable
  rd$.val <- rowMeans(
    suppressWarnings(sapply(rd[eco_cols], as.numeric)), na.rm = TRUE)

  long <- data.frame(
    year = suppressWarnings(as.integer(rd[[c_year]])),
    mon  = suppressWarnings(as.integer(rd[[c_mon]])),
    day  = suppressWarnings(as.integer(rd[[c_day]])),
    var  = tolower(trimws(as.character(rd[[c_var]]))),
    val  = rd$.val,
    stringsAsFactors = FALSE
  )
  long <- long[stats::complete.cases(long[, c("year", "mon", "day")]), ]

  w <- tidyr::pivot_wider(long, names_from = var, values_from = val,
                          values_fn = function(x) mean(x, na.rm = TRUE))
  w <- as.data.frame(w)
  w$date <- as.Date(sprintf("%04d-%02d-%02d", w$year, w$mon, w$day))
  w <- w[!is.na(w$date), ]
  w <- w[order(w$date), ]

  g <- function(pat) { h <- names(w)[grepl(pat, names(w), ignore.case = TRUE)][1]
                       if (is.na(h)) NULL else suppressWarnings(as.numeric(w[[h]])) }
  tmin <- g("^tmin$|min.*temp"); tmax <- g("^tmax$|max.*temp")
  tmean_direct <- g("^temp$|^tmean$")
  rmin <- g("^rmin$|min.*(rh|humid)"); rmax <- g("^rmax$|max.*(rh|humid)")
  rh_direct <- g("^rh$|^relativehumidity$")
  dewp <- g("^dewpoint$|^tdew$|dew")
  ppt  <- g("^ppt$|^precip"); wsp <- g("^windspeed$|^ws$|^wind$")

  if (is.null(ppt)) stop("No precipitation column (ppt / precip) found.")

  # Temp: prefer an explicit column, else the Climate Library's 0.5*(min+max)
  temp <- if (!is.null(tmean_direct)) tmean_direct
          else if (!is.null(tmin) && !is.null(tmax)) 0.5 * (tmin + tmax)
          else stop("Need either a temp column or both tmin and tmax.")

  # RH: explicit, else mean of min/max, else derived from dewpoint + mean temp
  rh <- if (!is.null(rh_direct)) rh_direct
        else if (!is.null(rmin) && !is.null(rmax)) 0.5 * (rmin + rmax)
        else if (!is.null(dewp))
          100 * pmin(1, exp(17.625 * (dewp / (243.04 + dewp) -
                                      temp / (243.04 + temp))))
        else stop("Need RH, or MinRH+MaxRH, or dewpoint.")
  rh <- pmax(pmin(rh, 100), 0.01)

  # Wind: LANDIS input is m/s and the library multiplies by 3.6
  ws <- if (!is.null(wsp)) wsp * 3.6 else rep(0, nrow(w))
  if (is.null(wsp))
    message("compute_landis_climate_fwi: no wind column — ISI computed with ws = 0.")

  prec <- ppt * precip_scale
  doy  <- as.integer(format(w$date, "%j"))

  # Codes reset each year; only SpringStart..WinterStart are computed
  out <- list(); i <- 1
  yrs <- sort(unique(w$year))
  for (k in seq_along(yrs)) {
    y <- yrs[k]
    if (!is.null(progress_fn))
      progress_fn(sprintf("FWI year %d (%d of %d)", y, k, length(yrs)))
    sel <- which(w$year == y & doy >= first_day & doy <= last_day)
    if (length(sel) < 2) next
    inp <- data.frame(
      yr = w$year[sel], mon = w$mon[sel], day = w$day[sel],
      lat = lat, long = -120,
      temp = temp[sel], rh = rh[sel], ws = ws[sel], prec = prec[sel]
    )
    inp <- inp[stats::complete.cases(inp), ]
    if (nrow(inp) < 2) next
    r <- tryCatch(
      suppressWarnings(cffdrs::fwi(
        inp, init = data.frame(ffmc = init_ffmc, dmc = init_dmc,
                               dc = init_dc, lat = lat),
        batch = TRUE, out = "all", uppercase = TRUE)),
      error = function(e) { warning(sprintf("FWI failed for %d: %s", y,
                                            conditionMessage(e))); NULL })
    if (is.null(r)) next
    out[[i]] <- tibble(
      date = as.Date(sprintf("%04d-%02d-%02d", inp$yr, inp$mon, inp$day)),
      FWI = as.numeric(r$FWI), FFMC = as.numeric(r$FFMC),
      DMC = as.numeric(r$DMC), DC = as.numeric(r$DC)
    )
    i <- i + 1
  }

  all_df <- bind_rows(out)
  if (nrow(all_df) == 0)
    stop("No FWI produced from the LANDIS climate file. Check the day window ",
         "(first_day/last_day) and that required variables are present.")
  all_df <- all_df %>% arrange(date)

  message(sprintf("  precip: file in %s, library_patched=%s -> multiplier x%g",
                  precip_units, library_patched, precip_scale))
  message(sprintf(
    "LANDIS climate FWI: %d days | %d years | FWI median %.1f, p90 %.1f, max %.1f",
    nrow(all_df), length(unique(format(all_df$date, "%Y"))),
    median(all_df$FWI, na.rm = TRUE),
    quantile(all_df$FWI, 0.9, na.rm = TRUE), max(all_df$FWI, na.rm = TRUE)))

  list(
    FWI  = all_df %>% transmute(date, value = FWI),
    FFMC = all_df %>% transmute(date, value = FFMC),
    DMC  = all_df %>% transmute(date, value = DMC),
    DC   = all_df %>% transmute(date, value = DC)
  )
}


# -----------------------------------------------------------------------------
# plot_fwi_comparison
#
# Side-by-side comparison of the CALIBRATION FWI (whatever series the models are
# fit against) and the RUNTIME FWI (computed from the LANDIS climate file with
# the Climate Library's conventions).
#
# This is the check that catches the single most damaging calibration error:
# fitting coefficients against one FWI definition and running SCF on another.
# The two series should overlap.  Where they do not, every fitted coefficient is
# being applied to a systematically different predictor, and no amount of
# parameter tuning will reconcile them.
#
# Panel A — distribution over the fire season (the days that matter).
# Panel B — monthly median cycle, which exposes seasonal-timing differences that
#           a pooled distribution can hide.
#
# Arguments:
#   cal_fwi, runtime_fwi  tibble(date, value) each
#   season_months         months treated as fire season for Panel A
# -----------------------------------------------------------------------------

plot_fwi_comparison <- function(cal_fwi, runtime_fwi = NULL,
                                 cal_label     = "Calibration FWI",
                                 runtime_label = "Runtime FWI (LANDIS)",
                                 season_months = 6:9) {

  if (is.null(cal_fwi) || !all(c("date", "value") %in% names(cal_fwi))) return(NULL)

  mk <- function(df, lab) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    tibble(date = as.Date(df$date), value = suppressWarnings(as.numeric(df$value)),
           series = lab) %>% filter(is.finite(value))
  }
  dd <- bind_rows(mk(cal_fwi, cal_label), mk(runtime_fwi, runtime_label))
  if (nrow(dd) == 0) return(NULL)
  dd$series <- factor(dd$series, levels = c(cal_label, runtime_label))

  cols <- setNames(c("#2c3e50", "#a0522d"), c(cal_label, runtime_label))
  seas <- dd %>% filter(lubridate::month(date) %in% season_months)

  # Summary line: the numbers that decide whether these are the same variable
  s <- seas %>% group_by(series) %>%
    summarise(n = n(),
              p50 = median(value, na.rm = TRUE),
              p90 = quantile(value, 0.90, na.rm = TRUE),
              mx  = max(value, na.rm = TRUE), .groups = "drop")
  sub <- paste(sprintf("%s: p50 %.1f, p90 %.1f, max %.1f",
                       s$series, s$p50, s$p90, s$mx), collapse = "   |   ")
  if (nrow(s) == 2) {
    r50 <- s$p50[2] / max(s$p50[1], 1e-9)
    r90 <- s$p90[2] / max(s$p90[1], 1e-9)
    sub <- paste0(sub, sprintf("\nRuntime / calibration: p50 %.2f\u00d7, p90 %.2f\u00d7",
                               r50, r90))
  } else {
    sub <- paste0(sub, "\nLoad the LANDIS climate file to add the runtime series.")
  }

  p_a <- ggplot(seas, aes(value, fill = series, colour = series)) +
    geom_density(alpha = 0.32, linewidth = 0.7) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(title = sprintf("A. Fire-season distribution (months %s)",
                         paste(range(season_months), collapse = "\u2013")),
         x = "FWI", y = "Density", fill = NULL, colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  mo <- dd %>%
    mutate(mon = lubridate::month(date)) %>%
    group_by(series, mon) %>%
    summarise(med = median(value, na.rm = TRUE), .groups = "drop")

  p_b <- ggplot(mo, aes(mon, med, colour = series)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.1) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "B. Monthly median cycle",
         x = NULL, y = "Median FWI", colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    (p_a / p_b) + patchwork::plot_annotation(
      title    = "Calibration vs Runtime FWI",
      subtitle = sub,
      caption  = paste0(
        "Calibration FWI = the series the ignition and spread models are fit ",
        "against.  Runtime FWI = computed from the LANDIS climate file using the\n",
        "Climate Library's conventions, i.e. what SCF actually consumes.  ",
        "These should overlap; if they do not, fitted coefficients are being\n",
        "applied to a different predictor than they were estimated from.")
    )
  } else {
    p_a + labs(subtitle = sub)
  }
}


# -----------------------------------------------------------------------------
# compute_startup_values
#
# Derives recommended SCF startup values (FFMC, DMC, DC) from ERA climatologies.
# Also returns overall means and FWI for reference.
#
# Returns a data.frame with one row per index.
# -----------------------------------------------------------------------------

compute_startup_values <- function(ffmc_ts, dmc_ts, dc_ts, fwi_ts, ref_date) {

  make_clim <- function(ts_df) {
    ts_df %>%
      mutate(
        doy  = yday(date),
        mmdd = format(date, "%m-%d")
      ) %>%
      group_by(doy, mmdd) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(doy)
  }

  pick_value <- function(clim_df, ref) {
    target_mmdd <- format(as.Date(ref), "%m-%d")
    hit <- clim_df %>% filter(mmdd == target_mmdd)
    if (nrow(hit) > 0) return(hit$mean_value[1])
    target_doy <- yday(as.Date(ref))
    clim_df %>% slice(which.min(abs(doy - target_doy))) %>% pull(mean_value)
  }

  codes <- list(FFMC = ffmc_ts, DMC = dmc_ts, DC = dc_ts, FWI = fwi_ts)
  result <- lapply(names(codes), function(nm) {
    ts   <- codes[[nm]]
    clim <- make_clim(ts)
    tibble(
      index         = nm,
      overall_mean  = mean(ts$value, na.rm = TRUE),
      startup_value = pick_value(clim, ref_date),
      ref_date      = as.character(ref_date)
    )
  })
  bind_rows(result)
}


# -----------------------------------------------------------------------------
# plot_era_climatologies
#
# Produces a 4-panel seasonal climatology plot (FFMC, DMC, DC, FWI)
# with a vertical line at ref_date.
# -----------------------------------------------------------------------------

plot_era_climatologies <- function(ffmc_ts, dmc_ts, dc_ts, fwi_ts, ref_date) {

  make_clim_df <- function(ts_df, code_name) {
    ts_df %>%
      mutate(
        doy       = yday(date),
        mmdd      = format(date, "%m-%d"),
        plot_date = as.Date(paste0("2001-", format(date, "%m-%d")))
      ) %>%
      group_by(doy, mmdd, plot_date) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      mutate(code = code_name)
  }

  all_clim <- bind_rows(
    make_clim_df(ffmc_ts, "FFMC"),
    make_clim_df(dmc_ts,  "DMC"),
    make_clim_df(dc_ts,   "DC"),
    make_clim_df(fwi_ts,  "FWI")
  )

  ref_plot_date <- as.Date(paste0("2001-", format(as.Date(ref_date), "%m-%d")))

  ggplot(all_clim, aes(plot_date, mean_value)) +
    geom_line(linewidth = 0.8, color = "#2c5f2e") +
    geom_vline(xintercept = as.numeric(ref_plot_date),
               linetype = "dashed", color = "#a0522d", linewidth = 0.8) +
    facet_wrap(~code, scales = "free_y", ncol = 2) +
    scale_x_date(date_labels = "%b", date_breaks = "2 months") +
    labs(
      title    = "ERA Seasonal Climatologies",
      subtitle = paste0("Dashed line = startup reference date (", ref_date, ")"),
      x        = NULL,
      y        = "Mean Value"
    ) +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "#2c5f2e"),
          strip.text = element_text(color = "white", face = "bold"))
}
