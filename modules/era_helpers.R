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
# Extracts daily park-mean wind speed (km/h) and direction (degrees) from
# ERA5 u10/v10 component files.
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
  # Speed: sqrt(u^2 + v^2), convert m/s -> km/h (*3.6)
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
#   vs  -- daily mean 10m wind speed (m/s)  -> converted to km/h
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

    vs_r <- gm$vs
    th_r <- gm$th

    # Mask to actual boundary polygon (climateR returns bounding-box crop)
    cal_v  <- vect(cal_sf)
    vs_r   <- mask(vs_r, cal_v)
    th_r   <- mask(th_r, cal_v)

    dates <- tryCatch(as.Date(time(vs_r)),
                      error = function(e) seq(yr_start, yr_end, by = "day"))

    vs_mean <- as.numeric(global(vs_r, "mean", na.rm = TRUE)[, 1])
    th_mean <- as.numeric(global(th_r, "mean", na.rm = TRUE)[, 1])

    results[[i]] <- tibble(
      date          = dates,
      WindSpeed_kmh = vs_mean * 3.6,   # m/s -> km/h
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
