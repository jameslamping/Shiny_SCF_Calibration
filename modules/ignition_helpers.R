# =============================================================================
# ignition_helpers.R
# Functions for SCF ignition parameter calibration.
#
# Key change from v1:
#   load_fpa_fod() now accepts cal_vect (calibration boundary, already in
#   working_crs) and working_crs separately. FPA FOD points are projected
#   to working_crs BEFORE the spatial filter, eliminating the CRS mismatch
#   that caused "replacement has N rows, data has 0" errors.
#
#   make_ignition_surfaces() accepts working_crs explicitly and reprojects
#   ignition points to match the LANDIS template grid.
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(terra, dplyr, lubridate, tidyr, ggplot2, glue)

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

pick_first_col <- function(nm, patterns) {
  for (p in patterns) {
    hit <- nm[grepl(p, tolower(nm))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

parse_any_date <- function(x) {
  if (inherits(x, "Date"))                 return(x)
  if (inherits(x, c("POSIXct","POSIXt"))) return(as.Date(x))
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  d <- suppressWarnings(as.Date(x, "%m/%d/%Y"))  # FPA FOD format
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y/%m/%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, "%Y%m%d"))
  d
}

# -----------------------------------------------------------------------------
# classify_fpa_ignition_type
# -----------------------------------------------------------------------------

classify_fpa_ignition_type <- function(general_cause, cause_classification) {
  gc  <- toupper(trimws(as.character(general_cause)))
  cc  <- toupper(trimws(as.character(cause_classification)))
  out <- rep(NA_character_, length(gc))

  out[gc == "LIGHTNING"]                            <- "Lightning"
  idx_na <- is.na(out)
  out[idx_na & cc == "NATURAL"]                     <- "Lightning"
  out[idx_na & !is.na(cc) & cc != "NATURAL"]       <- "Accidental"
  out[is.na(out)]                                   <- "Accidental"
  out
}

# -----------------------------------------------------------------------------
# load_fpa_fod
#
# Loads FPA FOD geodatabase and filters to the calibration boundary.
#
# Arguments:
#   gdb_path    path to FPA_FOD_*.gdb
#   cal_vect    calibration boundary SpatVector, ALREADY in working_crs
#   working_crs CRS string (from template raster) — FPA points projected here
#   year_min/max calibration year window
#
# Returns data.frame with columns:
#   date, ignition_type, FIRE_SIZE_HA, x, y  (coords in working_crs)
# -----------------------------------------------------------------------------

load_fpa_fod <- function(gdb_path, cal_vect, working_crs,
                          year_min = 1992, year_max = 2018) {

  if (!file.exists(gdb_path)) stop("FPA FOD GDB not found: ", gdb_path)

  # ---- Detect layer name --------------------------------------------------
  lyrs <- tryCatch(terra::vector_layers(gdb_path)$layer, error = function(e) NULL)
  lyr  <- if (!is.null(lyrs)) {
    hit <- lyrs[grepl("fires|fod|fpa", tolower(lyrs))]
    if (length(hit) > 0) hit[1] else lyrs[1]
  } else "Fires"

  message("Loading FPA FOD layer: ", lyr)
  fpa_v <- vect(gdb_path, layer = lyr)

  # ---- Validate source CRS ------------------------------------------------
  src_crs <- crs(fpa_v)
  if (is.na(src_crs) || nchar(src_crs) == 0) {
    stop("FPA FOD has no CRS. Cannot reproject to working CRS.")
  }

  # ---- Reproject FPA to working CRS BEFORE spatial filter ----------------
  # This is the critical step — always project the data to a common CRS
  # before any relate() / intersects() call.
  message("Reprojecting FPA FOD to working CRS...")
  fpa_proj <- project(fpa_v, working_crs)

  # ---- Spatial filter: FPA points within calibration boundary -------------
  # Use a bounding box pre-filter first (fast) then exact within test.
  # terra::relate() with "within" is reliable for points-in-polygon when
  # both objects share the same CRS.
  message("Filtering FPA FOD to calibration boundary...")

  # Step A: coarse bounding box filter (avoids relate() on millions of points)
  bb      <- ext(cal_vect)
  bb_poly <- as.polygons(bb, crs = working_crs)
  in_bb   <- relate(fpa_proj, bb_poly, "within")
  fpa_bb  <- fpa_proj[in_bb, ]

  if (nrow(fpa_bb) == 0) {
    stop(paste0(
      "No FPA FOD points fall within the calibration boundary bounding box.\n",
      "  FPA extent (working CRS): ", as.character(ext(fpa_proj)), "\n",
      "  Calibration boundary extent: ", as.character(ext(cal_vect)), "\n",
      "  Check that the calibration shapefile covers the correct region."
    ))
  }

  # Step B: exact within test on the bbox-filtered subset
  in_boundary <- relate(fpa_bb, cal_vect, "within")
  fpa_park    <- fpa_bb[in_boundary, ]

  if (nrow(fpa_park) == 0) {
    stop(paste0(
      nrow(fpa_bb), " FPA points were in the bounding box but 0 passed the ",
      "exact boundary test. The shapefile may have complex geometry — try ",
      "simplifying it or using a convex hull."
    ))
  }

  message(sprintf("FPA FOD: %d points within calibration boundary.", nrow(fpa_park)))

  # ---- Column detection ---------------------------------------------------
  nm       <- names(fpa_park)
  date_col <- pick_first_col(nm, c("discovery_date", "disc_date", "date", "fire_date"))
  gc_col   <- pick_first_col(nm, c("nwcg_general_cause", "general_cause", "cause"))
  cc_col   <- pick_first_col(nm, c("nwcg_cause_classification", "cause_classification",
                                    "stat_cause_descr"))
  size_col <- pick_first_col(nm, c("fire_size", "gis_acres", "acres", "size"))

  if (is.na(date_col)) stop("Could not find a date column in FPA FOD.")

  # ---- Build data.frame ---------------------------------------------------
  xy <- crds(fpa_park)

  df <- as.data.frame(fpa_park) %>%
    mutate(
      x    = xy[, 1],
      y    = xy[, 2],
      
      # Use FIRE_YEAR directly from the GDB — it's a clean integer field
      FIRE_YEAR = as.integer(FIRE_YEAR),
      
      # Parse date separately; NAs here won't kill the year filter
      date = parse_any_date(.data[[date_col]]),
      
      ignition_type = classify_fpa_ignition_type(
        if (!is.na(gc_col)) .data[[gc_col]] else rep(NA_character_, n()),
        if (!is.na(cc_col)) .data[[cc_col]] else rep(NA_character_, n())
      ),
      
      FIRE_SIZE_ACRES = if (!is.na(size_col))
        suppressWarnings(as.numeric(.data[[size_col]])) else NA_real_,
      FIRE_SIZE_HA = FIRE_SIZE_ACRES * 0.40468564224
    ) %>%
    # Filter on FIRE_YEAR first — doesn't depend on date parsing
    filter(FIRE_YEAR >= year_min, FIRE_YEAR <= year_max) %>%
    # Only drop records where date is NA after year filter passes
    filter(!is.na(date))

  if (nrow(df) == 0) {
    stop(sprintf("FPA FOD: 0 records remain after filtering to years %d-%d.",
                 year_min, year_max))
  }

  message(sprintf("FPA FOD: %d records retained (%d Lightning, %d Accidental).",
                  nrow(df),
                  sum(df$ignition_type == "Lightning",  na.rm = TRUE),
                  sum(df$ignition_type == "Accidental", na.rm = TRUE)))

  df
}


# -----------------------------------------------------------------------------
# build_ignition_model_data
#
# Aggregates FPA FOD to daily counts per ignition type, joined to daily FWI.
# Zero-fills days with no ignitions.
# Returns data.frame(date, FWI, ignition_type, n)
# -----------------------------------------------------------------------------

build_ignition_model_data <- function(ign_df, fwi_daily) {

  counts <- ign_df %>%
    group_by(date, ignition_type) %>%
    summarise(n = n(), .groups = "drop")

  grid <- expand.grid(
    date          = fwi_daily$date,
    ignition_type = c("Lightning", "Accidental"),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))

  grid %>%
    left_join(counts, by = c("date", "ignition_type")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    left_join(fwi_daily %>% rename(FWI = value), by = "date") %>%
    filter(!is.na(FWI))
}


# -----------------------------------------------------------------------------
# fit_ignition_models
#
# Fits Poisson or ZIP per ignition type.
# Returns data.frame(ignition_type, distribution, b0, b1, bz0, bz1)
# -----------------------------------------------------------------------------

fit_ignition_models <- function(model_data, distribution = c("ZIP", "Poisson")) {

  distribution <- match.arg(distribution)
  out <- list()

  for (typ in c("Lightning", "Accidental")) {

    d <- model_data %>% filter(ignition_type == typ, !is.na(FWI))

    if (nrow(d) < 10) {
      message("Skipping ", typ, ": fewer than 10 usable days.")
      next
    }

    if (distribution == "Poisson") {

      m <- glm(n ~ FWI, data = d, family = poisson(link = "log"))
      out[[typ]] <- tibble(
        ignition_type = typ, distribution = "Poisson",
        b0 = unname(coef(m)[1]), b1 = unname(coef(m)[2]),
        bz0 = NA_real_, bz1 = NA_real_
      )

    } else {

      if (!requireNamespace("pscl", quietly = TRUE))
        stop("ZIP requires package 'pscl'. Install with install.packages('pscl').")

      m  <- pscl::zeroinfl(n ~ FWI | FWI, data = d, dist = "poisson", link = "logit")
      co <- coef(m)

      # coef names vary by pscl version; handle both "count_(Intercept)" and "count.(Intercept)"
      get_coef <- function(pattern) {
        idx <- grep(pattern, names(co), perl = TRUE)
        if (length(idx) == 0) NA_real_ else unname(co[idx[1]])
      }

      out[[typ]] <- tibble(
        ignition_type = typ, distribution = "ZeroInflatedPoisson",
        b0  = get_coef("count.*Intercept"),
        b1  = get_coef("count.*FWI"),
        bz0 = get_coef("zero.*Intercept"),
        bz1 = get_coef("zero.*FWI")
      )
    }
  }

  bind_rows(out)
}


# -----------------------------------------------------------------------------
# make_ignition_surfaces
#
# Generates Lightning, Accidental, and Rx ignition allocation rasters
# on the LANDIS template grid using a Gaussian distance-decay kernel.
#
# Arguments:
#   ign_df      data.frame from load_fpa_fod() (x,y in working_crs)
#   template    LANDIS SpatRaster (defines output grid)
#   working_crs CRS string — both ign_df coords and template are in this CRS
#   bw_m        kernel bandwidth (metres)
#   maxdist_m   max influence distance (metres)
# -----------------------------------------------------------------------------

make_ignition_surfaces <- function(ign_df, template, working_crs,
                                    bw_m = 7000, maxdist_m = 25000) {
  surfaces <- list()

  for (typ in c("Lightning", "Accidental")) {
    pts <- ign_df %>% filter(ignition_type == typ) %>% select(x, y)
    surfaces[[typ]] <- make_kernel_surface(pts, template, bw_m, maxdist_m)
  }

  # Rx: uniform surface (no empirical Rx data in FPA FOD)
  rx_r <- rast(template)
  tv   <- values(template, mat = FALSE)
  values(rx_r) <- ifelse(is.na(tv), NA_real_, 1.0)
  surfaces[["Rx"]] <- rx_r

  surfaces
}


# -----------------------------------------------------------------------------
# make_kernel_surface  (internal)
# Gaussian distance-decay kernel on LANDIS template grid.
# -----------------------------------------------------------------------------

make_kernel_surface <- function(points_df, template, bandwidth_m, maxdist_m) {
  
  if (nrow(points_df) == 0) {
    r <- rast(template); values(r) <- NA_real_; return(r)
  }
  
  # ---- Step 1: rasterize point counts onto template grid ------------------
  # Convert points to SpatVector, project to template CRS
  pts_v <- vect(points_df, geom = c("x", "y"), crs = crs(template))
  count_r <- rasterize(pts_v, template, fun = "count", background = 0)
  
  # Mask to active cells only
  active_mask <- !is.na(values(template, mat = FALSE))
  count_vals  <- values(count_r, mat = FALSE)
  count_vals[!active_mask] <- NA
  values(count_r) <- count_vals
  
  # ---- Step 2: build Gaussian kernel weight matrix ------------------------
  # Cell size in metres (assumes projected CRS)
  cell_m <- mean(res(template))
  
  # Kernel radius in cells
  radius_cells <- ceiling(maxdist_m / cell_m)
  
  # Build kernel weight matrix
  # Odd-sized matrix centred on focal cell
  k_size <- 2 * radius_cells + 1
  kernel  <- matrix(0, nrow = k_size, ncol = k_size)
  
  for (row in seq_len(k_size)) {
    for (col in seq_len(k_size)) {
      dx <- (col - radius_cells - 1) * cell_m
      dy <- (row - radius_cells - 1) * cell_m
      d2 <- dx^2 + dy^2
      if (sqrt(d2) <= maxdist_m) {
        kernel[row, col] <- exp(-d2 / (2 * bandwidth_m^2))
      }
    }
  }
  
  # ---- Step 3: apply focal filter (convolution) ---------------------------
  smooth_r <- focal(count_r, w = kernel, fun = "sum",
                    na.policy = "omit", na.rm = TRUE)
  
  # ---- Step 4: normalize to 0-1 -------------------------------------------
  sv  <- values(smooth_r, mat = FALSE)
  mx  <- max(sv, na.rm = TRUE)
  if (is.finite(mx) && mx > 0) sv <- sv / mx
  sv[!active_mask] <- NA
  values(smooth_r) <- sv
  
  smooth_r
}


# -----------------------------------------------------------------------------
# plot_ignition_diagnostics
# Observed vs. predicted ignition counts by FWI bin.
# -----------------------------------------------------------------------------

plot_ignition_diagnostics <- function(model_data, coef_df) {

  predict_ign <- function(fwi_vec, row) {
    if (row$distribution == "Poisson") {
      exp(row$b0 + row$b1 * fwi_vec)
    } else {
      lambda <- exp(row$b0 + row$b1 * fwi_vec)
      p_zero <- 1 / (1 + exp(row$bz0 + row$bz1 * fwi_vec))
      (1 - p_zero) * lambda
    }
  }

  fwi_seq <- seq(0, max(model_data$FWI, na.rm = TRUE), length.out = 100)

  pred_df <- map_dfr(seq_len(nrow(coef_df)), function(i) {
    row <- coef_df[i, ]
    tibble(FWI = fwi_seq,
           pred = predict_ign(fwi_seq, row),
           ignition_type = row$ignition_type)
  })

  obs_binned <- model_data %>%
    mutate(FWI_bin = round(FWI)) %>%
    group_by(FWI_bin, ignition_type) %>%
    summarise(mean_n = mean(n), .groups = "drop")

  ggplot() +
    geom_col(data = obs_binned, aes(FWI_bin, mean_n),
             fill = "#c8d8c8", width = 1) +
    geom_line(data = pred_df, aes(FWI, pred),
              color = "#a0522d", linewidth = 1) +
    facet_wrap(~ignition_type, scales = "free_y") +
    labs(
      title    = "Ignition Model Diagnostics",
      subtitle = "Bars = observed mean daily ignitions per FWI bin; line = fitted model",
      x = "Fire Weather Index (FWI)", y = "Mean daily ignitions"
    ) +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "#2c5f2e"),
          strip.text = element_text(color = "white", face = "bold"))
}
