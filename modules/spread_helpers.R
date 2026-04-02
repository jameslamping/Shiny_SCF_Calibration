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
               lubridate, glue, broom, scales)

ROOK_W <- matrix(c(0,1,0, 1,1,1, 0,1,0), nrow = 3, byrow = TRUE)

# Reuse utilities from ignition_helpers (sourced before this file in app.R)
# pick_first_col(), parse_any_date() are defined there.


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
# bind_climate_to_pairs
# -----------------------------------------------------------------------------

bind_climate_to_pairs <- function(geomac_v, fwi_daily, wind_daily,
                                   max_gap_sp = 1, max_gap_da = 3,
                                   neg_tol_ha = 1.0) {

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

  # Join climate
  climate <- fwi_daily %>%
    rename(FWI = value) %>%
    left_join(wind_daily, by = "date")

  pairs <- pairs %>%
    left_join(climate, by = "date") %>%
    mutate(
      FWI           = suppressWarnings(as.numeric(FWI)),
      WindSpeed_kmh = suppressWarnings(as.numeric(WindSpeed_kmh)),
      WindDir_deg   = suppressWarnings(as.numeric(WindDir_deg)),
      EffectiveWind = WindSpeed_kmh  # placeholder: swap for topographic formula later
    )

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
# log1p(daily_area_ha) ~ FWI + EffectiveWind
# -----------------------------------------------------------------------------

fit_max_daily_area <- function(pairs) {

  dat <- pairs %>%
    filter(use_for_daily_area,
           is.finite(daily_area_ha), daily_area_ha > 0,
           is.finite(FWI), is.finite(EffectiveWind))

  if (nrow(dat) < 5)
    stop("Fewer than 5 usable pairs for max daily area fit.")

  m <- lm(log1p(daily_area_ha) ~ FWI + EffectiveWind, data = dat)

  coef_df <- broom::tidy(m) %>%
    rename(std_error = std.error, p_value = p.value) %>%
    mutate(
      term = case_when(
        term == "(Intercept)"  ~ "B0_Intercept",
        term == "FWI"          ~ "B1_FWI",
        term == "EffectiveWind" ~ "B2_EffectiveWind",
        TRUE ~ term
      ),
      note = "Fit on log1p(ha/day) scale"
    )

  list(coef = coef_df, model = m, data = dat)
}


# -----------------------------------------------------------------------------
# fit_spread_probability
# Logistic model from rasterized perimeter pairs.
# -----------------------------------------------------------------------------

fit_spread_probability <- function(pairs2, geomac_v, template_r, park_mask,
                                    failure_ratio = 3, max_samples = 25000,
                                    progress_fn = NULL) {

  eligible <- pairs2 %>%
    filter(use_for_spread_prob, is.finite(WindSpeed_kmh))

  if (nrow(eligible) == 0)
    stop("No eligible perimeter pairs for spread probability fit.")


  samples_list <- vector("list", nrow(eligible))

  for (i in seq_len(nrow(eligible))) {
    if (!is.null(progress_fn)) progress_fn(i, nrow(eligible))

    row  <- eligible[i, ]
    fwi  <- row$FWI
    wspd <- row$WindSpeed_kmh
    if (!is.finite(wspd)) next

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

    samples_list[[i]] <- tibble(
      FIRE_ID       = rep(row$FIRE_ID, n_succ + n_f),
      date          = rep(row$date,    n_succ + n_f),
      y             = c(rep(1L, n_succ), rep(0L, n_f)),
      FWI           = rep(fwi,  n_succ + n_f),
      WindSpeed_kmh = rep(wspd, n_succ + n_f),
      FineFuels     = rep(1.0,  n_succ + n_f)
    )
  }

  all_samples <- bind_rows(samples_list)
  if (nrow(all_samples) == 0)
    stop("No raster samples generated. Check that perimeters fall within template extent.")

  m <- glm(y ~ FWI + FineFuels + WindSpeed_kmh,
           data   = all_samples,
           family = binomial(link = "logit"))

  coef_df <- broom::tidy(m) %>%
    rename(std_error = std.error, p_value = p.value) %>%
    mutate(
      term = case_when(
        term == "(Intercept)"   ~ "B0_Intercept",
        term == "FWI"           ~ "B1_FWI",
        term == "FineFuels"     ~ "B2_FineFuels (placeholder)",
        term == "WindSpeed_kmh" ~ "B3_EffectiveWind",
        TRUE ~ term
      ),
      scf_parameter = case_when(
        grepl("B0", term) ~ "SpreadProbabilityB0",
        grepl("B1", term) ~ "SpreadProbabilityB1",
        grepl("B2", term) ~ "SpreadProbabilityB2",
        grepl("B3", term) ~ "SpreadProbabilityB3",
        TRUE ~ NA_character_
      )
    )

  message(sprintf("Spread prob fit: %d samples (%d successes, %d failures)",
                  nrow(all_samples),
                  sum(all_samples$y == 1),
                  sum(all_samples$y == 0)))

  list(coef = coef_df, model = m, samples = all_samples)
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

plot_spread_prob_diagnostics <- function(pairs, coef_df) {

  b  <- setNames(coef_df$estimate, coef_df$term)
  b0 <- b["B0_Intercept"]
  b1 <- b["B1_FWI"]
  b3 <- b["B3_EffectiveWind"]

  if (any(is.na(c(b0, b1))))
    return(ggplot() + labs(title = "Spread probability not yet fitted.") + theme_bw())

  fwi_seq   <- seq(0, max(pairs$FWI, na.rm = TRUE), length.out = 100)
  wspd_mean <- mean(pairs$WindSpeed_kmh, na.rm = TRUE)

  pred_df <- tibble(
    FWI  = fwi_seq,
    prob = plogis(b0 + b1 * FWI + (if (!is.na(b3)) b3 * wspd_mean else 0) + 1.0)
  )

  ggplot(pred_df, aes(FWI, prob)) +
    geom_line(color = "#2c5f2e", linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    labs(title   = "Fitted Spread Probability vs FWI",
         subtitle = sprintf("At mean wind speed (%.1f km/h), fine fuels = 1.0", wspd_mean),
         x = "FWI", y = "P(spread to neighbor cell)") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw(base_size = 12)
}


plot_max_area_diagnostics <- function(pairs, coef_df) {

  dat <- pairs %>%
    filter(use_for_daily_area,
           is.finite(daily_area_ha), daily_area_ha > 0,
           is.finite(FWI), is.finite(EffectiveWind))

  if (nrow(dat) == 0 || is.null(coef_df))
    return(ggplot() + labs(title = "No max daily area data.") + theme_bw())

  b <- setNames(coef_df$estimate, coef_df$term)
  dat <- dat %>%
    mutate(pred_ha = expm1(b["B0_Intercept"] +
                           b["B1_FWI"] * FWI +
                           b["B2_EffectiveWind"] * EffectiveWind))

  ggplot(dat, aes(pred_ha, daily_area_ha)) +
    geom_point(alpha = 0.4, color = "#2c5f2e", size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#a0522d") +
    scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
    labs(title    = "Max Daily Area: Observed vs Predicted",
         subtitle = "Log-log scale; dashed = 1:1 line",
         x = "Predicted (ha/day)", y = "Observed (ha/day)") +
    theme_bw(base_size = 12)
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

write_scf_snippet <- function(row) {
  paste(c(
    sprintf(">> Candidate ID: %d", row$candidate_id),
    "",
    ">> Cell-to-cell spread probability (SCF Eq. 5)",
    sprintf("SpreadProbabilityB0    %.6f", row$B0),
    sprintf("SpreadProbabilityB1    %.6f", row$B1_FWI),
    sprintf("SpreadProbabilityB2    %.6f", row$B2_FineFuels),
    sprintf("SpreadProbabilityB3    %.6f", row$B3_WindSpeed),
    "",
    ">> Maximum daily spread area (SCF Eq. 6)",
    sprintf("MaximumSpreadAreaB0    %.6f", row$MaxAreaB0),
    sprintf("MaximumSpreadAreaB1    %.6f", row$MaxAreaB1_FWI),
    sprintf("MaximumSpreadAreaB2    %.6f", row$MaxAreaB2_Wind)
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
