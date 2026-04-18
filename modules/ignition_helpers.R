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
pacman::p_load(terra, dplyr, lubridate, tidyr, ggplot2, glue, scales)

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

  # Rx: uniform surface (no empirical Rx data in FPA FOD).
  # Set all active cells to 100 so every cell carries equal weight after
  # SCF's integer cast — consistent with the 0-100 scale used for the
  # empirical Lightning and Accidental surfaces.
  rx_r <- rast(template)
  tv   <- values(template, mat = FALSE)
  values(rx_r) <- ifelse(is.na(tv), NA_real_, 100)
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
  
  # ---- Step 4: scale to 0-100 integers ---------------------------------------
  # SCF's MapUtility.ReadMap() casts pixel values to (int) at load time, and
  # PreShuffleEther() casts again when building the weighted selector.  Any
  # value in [0, 1) truncates to 0 and the cell is never selected as an ignition
  # site.  Scaling to the range 0-100 (rounded integers) preserves the full
  # relative variation in ignition probability while surviving both int casts.
  sv  <- values(smooth_r, mat = FALSE)
  mx  <- max(sv, na.rm = TRUE)
  if (is.finite(mx) && mx > 0) sv <- round(sv / mx * 100)
  sv[!active_mask] <- NA
  values(smooth_r) <- sv

  smooth_r
}


# Shared theme element for all ignition diagnostic plots
.ign_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "#2c5f2e"),
      strip.text       = element_text(color = "white", face = "bold"),
      legend.position  = "bottom"
    )
}

# Internal: compute full E[n] = (1 - pi) * lambda for one coef row
.predict_ign_row <- function(fwi_vec, row) {
  lambda <- exp(row$b0 + row$b1 * fwi_vec)
  if (row$distribution == "Poisson") {
    lambda
  } else {
    # plogis() is the standard inverse-logit: 1 / (1 + exp(-x))
    # pi is the structural-zero probability; (1 - pi) * lambda is E[n]
    pi <- plogis(row$bz0 + row$bz1 * fwi_vec)
    (1 - pi) * lambda
  }
}

# Internal: observed mean count and zero-fraction per integer FWI bin
.obs_binned <- function(model_data) {
  model_data %>%
    mutate(FWI_bin = round(FWI)) %>%
    group_by(FWI_bin, ignition_type) %>%
    summarise(mean_n    = mean(n),
              zero_frac = mean(n == 0),
              n_days    = n(),
              .groups   = "drop")
}


# -----------------------------------------------------------------------------
# plot_ignition_diagnostics
# Observed mean daily ignitions (bars) vs. full model E[n] (line).
# -----------------------------------------------------------------------------

plot_ignition_diagnostics <- function(model_data, coef_df) {

  fwi_seq <- seq(0, max(model_data$FWI, na.rm = TRUE), length.out = 200)

  pred_df <- map_dfr(seq_len(nrow(coef_df)), function(i) {
    row <- coef_df[i, ]
    tibble(FWI           = fwi_seq,
           pred          = .predict_ign_row(fwi_seq, row),
           ignition_type = row$ignition_type,
           distribution  = row$distribution)
  })

  obs <- .obs_binned(model_data)

  dist_label <- unique(pred_df$distribution)
  subtitle <- if (length(dist_label) == 1 && dist_label == "Poisson") {
    "Bars = observed mean daily ignitions per FWI bin  |  Line = Poisson E[n]"
  } else {
    "Bars = observed mean daily ignitions per FWI bin  |  Line = ZIP E[n] = (1\u2212\u03c0)\u00d7\u03bb"
  }

  ggplot() +
    geom_col(data = obs, aes(FWI_bin, mean_n),
             fill = "#c8d8c8", width = 1) +
    geom_line(data = pred_df, aes(FWI, pred),
              color = "#a0522d", linewidth = 1.1) +
    facet_wrap(~ignition_type, scales = "free_y") +
    labs(title    = "Observed vs. Fitted Ignition Rate",
         subtitle = subtitle,
         x = "Fire Weather Index (FWI)", y = "Mean daily ignitions") +
    .ign_theme()
}


# -----------------------------------------------------------------------------
# plot_zero_inflation
# Shows how the structural-zero probability (pi) varies with FWI alongside
# the observed fraction of zero-count days.  ZIP only.
# -----------------------------------------------------------------------------

plot_zero_inflation <- function(model_data, coef_df) {

  zip_rows <- coef_df %>% filter(distribution != "Poisson")
  if (nrow(zip_rows) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Zero-inflation panel only available for ZIP model.",
                      size = 5, hjust = 0.5) +
             theme_void())
  }

  fwi_seq <- seq(0, max(model_data$FWI, na.rm = TRUE), length.out = 200)

  pi_df <- map_dfr(seq_len(nrow(zip_rows)), function(i) {
    row <- zip_rows[i, ]
    tibble(FWI           = fwi_seq,
           pi            = plogis(row$bz0 + row$bz1 * fwi_seq),
           ignition_type = row$ignition_type)
  })

  obs <- .obs_binned(model_data) %>% filter(n_days >= 5)

  ggplot() +
    geom_col(data = obs, aes(FWI_bin, zero_frac),
             fill = "#e8d5b0", width = 1) +
    geom_line(data = pi_df, aes(FWI, pi),
              color = "#6a3d9a", linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
    facet_wrap(~ignition_type) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      title    = "Zero-Inflation Component: Structural-Zero Probability (\u03c0) vs FWI",
      subtitle = paste0(
        "Bars = observed fraction of zero-count days (bins with \u22655 days)\n",
        "Line = \u03c0 = plogis(bz0 + bz1 \u00d7 FWI)  |  Values above 50% line: ",
        "most days are structurally zero"
      ),
      x = "Fire Weather Index (FWI)",
      y = "Probability of structural zero (\u03c0)"
    ) +
    .ign_theme()
}


# -----------------------------------------------------------------------------
# plot_count_component
# Isolates lambda = exp(b0 + b1 * FWI) for Poisson and ZIP side-by-side.
# Helps explain how b0/b1 differ when zeros are modelled separately.
# -----------------------------------------------------------------------------

plot_count_component <- function(model_data, coef_df) {

  fwi_seq <- seq(0, max(model_data$FWI, na.rm = TRUE), length.out = 200)

  lambda_df <- map_dfr(seq_len(nrow(coef_df)), function(i) {
    row   <- coef_df[i, ]
    label <- if (row$distribution == "Poisson") "Poisson \u03bb" else "ZIP \u03bb (count component)"
    tibble(FWI           = fwi_seq,
           lambda        = exp(row$b0 + row$b1 * fwi_seq),
           model         = label,
           ignition_type = row$ignition_type)
  })

  obs <- .obs_binned(model_data)

  ggplot() +
    geom_col(data = obs, aes(FWI_bin, mean_n),
             fill = "#c8d8c8", width = 1) +
    geom_line(data = lambda_df,
              aes(FWI, lambda, color = model, linetype = model),
              linewidth = 1.1) +
    scale_color_manual(
      values = c("Poisson \u03bb" = "#a0522d",
                 "ZIP \u03bb (count component)" = "#2980b9")
    ) +
    scale_linetype_manual(
      values = c("Poisson \u03bb" = "solid",
                 "ZIP \u03bb (count component)" = "dashed")
    ) +
    facet_wrap(~ignition_type, scales = "free_y") +
    labs(
      title    = "Count Component (\u03bb) by Model",
      subtitle = paste0(
        "\u03bb = exp(b0 + b1 \u00d7 FWI)  |  SCF uses B0/B1 from this component\n",
        "ZIP \u03bb > Poisson \u03bb at low FWI because ZIP attributes zeros to \u03c0, ",
        "not to a low count rate"
      ),
      x = "Fire Weather Index (FWI)", y = "\u03bb (expected count | non-structural-zero day)",
      color = NULL, linetype = NULL
    ) +
    .ign_theme()
}


# -----------------------------------------------------------------------------
# plot_ignition_residuals
# Mean (observed - predicted) per 2-unit FWI bin.  Flags systematic bias.
# -----------------------------------------------------------------------------

plot_ignition_residuals <- function(model_data, coef_df) {

  fwi_seq_full <- model_data$FWI

  resid_df <- map_dfr(seq_len(nrow(coef_df)), function(i) {
    row   <- coef_df[i, ]
    d     <- model_data %>% filter(ignition_type == row$ignition_type, !is.na(FWI))
    mu    <- .predict_ign_row(d$FWI, row)
    label <- if (row$distribution == "Poisson") "Poisson" else "ZIP"
    tibble(FWI           = d$FWI,
           residual      = d$n - mu,
           model         = label,
           ignition_type = row$ignition_type)
  }) %>%
    mutate(FWI_bin = round(FWI / 2) * 2) %>%
    group_by(FWI_bin, ignition_type, model) %>%
    summarise(mean_r = mean(residual),
              sd_r   = sd(residual),
              .groups = "drop")

  ggplot(resid_df, aes(FWI_bin, mean_r, color = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r, fill = model),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c("Poisson" = "#a0522d", "ZIP" = "#2980b9")) +
    scale_fill_manual(values  = c("Poisson" = "#a0522d", "ZIP" = "#2980b9")) +
    facet_wrap(~ignition_type, scales = "free_y") +
    labs(
      title    = "Model Residuals by FWI Bin",
      subtitle = "Mean \u00b1 SD of (observed \u2212 predicted) per 2-unit FWI bin  |  Zero line = perfect fit",
      x = "Fire Weather Index (FWI)", y = "Mean residual (observed \u2212 predicted)",
      color = NULL, fill = NULL
    ) +
    .ign_theme()
}


# =============================================================================
# SCALE ADJUSTMENT HELPERS
# Functions for adjusting ignition model coefficients from calibration-region
# scale to LANDIS landscape scale, and for summarising expected annual
# ignitions at the landscape scale.
#
# Two approaches:
#   Option A — Area offset: adjust b0 by log(template_area / cal_area).
#              Assumes uniform ignition density across the calibration region.
#   Option B — Park-specific intercept: re-estimate b0 from FPA FOD points
#              within the template extent, holding b1 fixed from regional fit.
#              Requires enough fire history within the template (~15+ years).
# =============================================================================

# -----------------------------------------------------------------------------
# compute_landscape_areas
# Returns calibration boundary area, template landscape area, and the log
# area ratio used as the b0 offset in Option A.
# Both inputs must be in the same projected CRS (working_crs).
# -----------------------------------------------------------------------------

compute_landscape_areas <- function(cal_vect_proj, template_r) {
  cal_area_ha  <- sum(expanse(cal_vect_proj, unit = "ha"))

  cell_ha      <- prod(res(template_r)) / 1e4        # m² -> ha
  n_active     <- sum(!is.na(values(template_r, mat = FALSE)))
  tmpl_area_ha <- n_active * cell_ha

  list(
    cal_area_ha  = cal_area_ha,
    tmpl_area_ha = tmpl_area_ha,
    ratio        = tmpl_area_ha / cal_area_ha,
    log_ratio    = log(tmpl_area_ha / cal_area_ha)
  )
}


# -----------------------------------------------------------------------------
# apply_area_offset  (Option A)
# Adjusts b0 for each ignition type by log(template_area / cal_area).
# b1, bz0, bz1 are unchanged — only the intercept scales with area.
# -----------------------------------------------------------------------------

apply_area_offset <- function(coef_df, cal_area_ha, tmpl_area_ha) {
  log_offset <- log(tmpl_area_ha / cal_area_ha)
  coef_df %>%
    mutate(
      b0     = b0 + log_offset,
      method = sprintf("Area offset: log(%.0f / %.0f) = %.4f added to b0",
                       tmpl_area_ha, cal_area_ha, log_offset)
    )
}


# -----------------------------------------------------------------------------
# filter_fpa_to_template  (Option B helper)
# Keeps only FPA FOD records whose locations fall within the active (non-NA)
# cells of the LANDIS template raster.  Uses terra::extract() for exactness.
# ign_df must have columns x, y in the same CRS as template_r.
# -----------------------------------------------------------------------------

filter_fpa_to_template <- function(ign_df, template_r) {
  if (nrow(ign_df) == 0) return(ign_df)

  pts_v     <- vect(ign_df, geom = c("x", "y"), crs = crs(template_r))
  tmpl_vals <- extract(template_r, pts_v)[, 2]   # NA where outside landscape
  ign_df[!is.na(tmpl_vals), ]
}


# -----------------------------------------------------------------------------
# fit_park_intercept  (Option B)
# Re-estimates b0 for each ignition type using FPA FOD fires within the
# template extent, holding b1 fixed at the regional value.  Always returns
# Poisson (ZIP is not feasible at park scale with sparse data).
# Falls back to area-offset-adjusted coef if <3 fires exist for a type.
# -----------------------------------------------------------------------------

fit_park_intercept <- function(park_model_data, coef_regional,
                                cal_area_ha, tmpl_area_ha) {
  out <- list()

  for (typ in c("Lightning", "Accidental")) {
    d   <- park_model_data %>% filter(ignition_type == typ, !is.na(FWI))
    reg <- coef_regional   %>% filter(ignition_type == typ)
    n_fires <- sum(d$n)

    if (n_fires < 3) {
      message(sprintf(
        "%s: only %d fires in template — too few to fit. Using area-offset b0.",
        typ, n_fires))
      fallback <- apply_area_offset(reg, cal_area_ha, tmpl_area_ha)
      out[[typ]] <- fallback %>%
        mutate(method = paste0("Fallback (area offset): ", method,
                               sprintf(" [only %d fires in template]", n_fires)))
      next
    }

    b1_fixed <- reg$b1
    m        <- glm(n ~ 1, offset = b1_fixed * FWI,
                    data = d, family = poisson(link = "log"))
    b0_park  <- unname(coef(m)[1])

    out[[typ]] <- tibble(
      ignition_type = typ,
      distribution  = "Poisson",
      b0  = b0_park,
      b1  = b1_fixed,
      bz0 = NA_real_,
      bz1 = NA_real_,
      method = sprintf(
        "Park intercept: b0 from %d fires in template; b1=%.4f fixed from regional fit",
        n_fires, b1_fixed)
    )
  }

  bind_rows(out)
}


# -----------------------------------------------------------------------------
# predict_annual_ignitions
# Applies adjusted coefficients to the ERA FWI daily time series to produce
# expected annual ignition counts at the landscape scale.
# Returns tibble(year, ignition_type, expected_annual, scenario).
# -----------------------------------------------------------------------------

predict_annual_ignitions <- function(coef_df, fwi_daily,
                                      scenario_label = "Landscape-adjusted") {
  d <- fwi_daily %>% mutate(year = year(date))

  map_dfr(seq_len(nrow(coef_df)), function(i) {
    row <- coef_df[i, ]
    d %>%
      mutate(E_n = .predict_ign_row(value, row)) %>%
      group_by(year) %>%
      summarise(expected_annual = sum(E_n, na.rm = TRUE), .groups = "drop") %>%
      mutate(ignition_type = row$ignition_type,
             scenario      = scenario_label)
  })
}


# -----------------------------------------------------------------------------
# summarise_expected_ignitions
# Summary statistics table from predict_annual_ignitions() output.
# -----------------------------------------------------------------------------

summarise_expected_ignitions <- function(annual_df) {
  type_rows <- annual_df %>%
    group_by(scenario, ignition_type) %>%
    summarise(
      mean_annual   = round(mean(expected_annual),              1),
      median_annual = round(median(expected_annual),            1),
      sd_annual     = round(sd(expected_annual),                1),
      p10_annual    = round(quantile(expected_annual, 0.10),    1),
      p90_annual    = round(quantile(expected_annual, 0.90),    1),
      .groups = "drop"
    )

  # Total row per scenario
  total_rows <- annual_df %>%
    group_by(scenario, year) %>%
    summarise(expected_annual = sum(expected_annual), .groups = "drop") %>%
    group_by(scenario) %>%
    summarise(
      ignition_type = "TOTAL",
      mean_annual   = round(mean(expected_annual),           1),
      median_annual = round(median(expected_annual),         1),
      sd_annual     = round(sd(expected_annual),             1),
      p10_annual    = round(quantile(expected_annual, 0.10), 1),
      p90_annual    = round(quantile(expected_annual, 0.90), 1),
      .groups = "drop"
    )

  bind_rows(type_rows, total_rows) %>%
    arrange(scenario, ignition_type)
}


# -----------------------------------------------------------------------------
# plot_expected_ignitions
# Distribution of expected annual ignitions across the FWI climatology period.
# Shows both regional (unadjusted) and landscape-adjusted scenarios together.
# -----------------------------------------------------------------------------

plot_expected_ignitions <- function(annual_df) {
  # Add total column per year/scenario
  totals <- annual_df %>%
    group_by(scenario, year) %>%
    summarise(expected_annual = sum(expected_annual), .groups = "drop") %>%
    mutate(ignition_type = "Total")

  plot_df <- bind_rows(annual_df, totals) %>%
    mutate(ignition_type = factor(ignition_type,
                                   levels = c("Lightning", "Accidental", "Total")))

  means_df <- plot_df %>%
    group_by(scenario, ignition_type) %>%
    summarise(mean_val = mean(expected_annual), .groups = "drop")

  scenario_colors <- c(
    "Regional (unadjusted)" = "#a0522d",
    "Landscape-adjusted"    = "#2980b9",
    "Park intercept"        = "#27ae60"
  )

  ggplot(plot_df, aes(expected_annual, fill = scenario, color = scenario)) +
    geom_histogram(bins = 12, alpha = 0.55, position = "identity") +
    geom_vline(data = means_df,
               aes(xintercept = mean_val, color = scenario),
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values  = scenario_colors, na.value = "grey70") +
    scale_color_manual(values = scenario_colors, na.value = "grey70") +
    facet_wrap(~ignition_type, scales = "free", nrow = 1) +
    labs(
      title    = "Expected Annual Ignitions — Calibration FWI Climatology",
      subtitle = "Histograms across calibration years  |  Dashed line = mean  |  Landscape-adjusted = what SCF will simulate",
      x = "Expected ignitions per year",
      y = "Number of years",
      fill = NULL, color = NULL
    ) +
    .ign_theme(base_size = 12)
}
