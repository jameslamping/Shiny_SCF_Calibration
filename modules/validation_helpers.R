# =============================================================================
# validation_helpers.R
# Publication-ready model validation figures for the SCF calibration app.
#
# Compares SCF socialclimatefire-events-log.csv output against observed
# datasets (FPA FOD fire occurrence/size, ERA FWI climate) to assess how
# well the calibrated model reproduces key fire regime characteristics.
#
# Panels produced:
#   1. Annual ignition counts   — simulated vs expected from ERA+model
#   2. Fire size distribution   — ECDF, simulated vs FPA FOD historical
#   3. Annual area burned        — distribution comparison, simulated vs FPA FOD
#   4. Fire-climate relationship — log(fire_size) vs FWI
#   5. Severity distribution    — MeanDNBR from SCF events log
#   6. Fire duration             — NumberOfDays distribution
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, tidyr, purrr, scales, lubridate, stringr)

# -----------------------------------------------------------------------------
# Publication theme  (white bg, clean axes, bold strip labels, legend bottom)
# -----------------------------------------------------------------------------

.pub_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      axis.title         = element_text(face = "bold", size = base_size),
      axis.text          = element_text(size = base_size - 1),
      strip.background   = element_rect(fill = "grey95", colour = "grey70"),
      strip.text         = element_text(face = "bold", size = base_size),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = base_size - 1),
      plot.title         = element_text(face = "bold", size = base_size + 1,
                                        margin = margin(b = 4)),
      plot.subtitle      = element_text(colour = "grey40", size = base_size - 1,
                                        margin = margin(b = 8)),
      plot.caption       = element_text(colour = "grey50", size = base_size - 2,
                                        hjust = 0, margin = margin(t = 6))
    )
}

# Shared colour constants
.COL_OBS <- "#2c3e50"   # observed / FPA FOD
.COL_SIM <- "#c0392b"   # simulated / SCF
.COL_EXP <- "#2c5f2e"   # expected / model range


# -----------------------------------------------------------------------------
# load_scf_events_full
#
# Reads the full socialclimatefire-events-log.csv and returns a tidy
# event-level tibble.  cell_area_ha converts TotalSitesBurned to ha.
# -----------------------------------------------------------------------------

load_scf_events_full <- function(path, cell_area_ha = 1.0) {
  if (!file.exists(path)) stop("Events log not found: ", path)

  df <- read.csv(path, stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))

  # Standardise ignition type label
  if ("IgnitionType" %in% names(df)) {
    df$ignition_type <- case_when(
      grepl("lightning",  df$IgnitionType, ignore.case = TRUE) ~ "Lightning",
      grepl("accidental", df$IgnitionType, ignore.case = TRUE) ~ "Accidental",
      grepl("rx|prescribed", df$IgnitionType, ignore.case = TRUE) ~ "Rx",
      TRUE ~ trimws(df$IgnitionType)
    )
  } else {
    df$ignition_type <- NA_character_
  }

  # Convert burned cells to ha
  if ("TotalSitesBurned" %in% names(df)) {
    df$fire_size_ha <- df$TotalSitesBurned * cell_area_ha
  } else {
    df$fire_size_ha <- NA_real_
  }

  # Normalise FWI column names
  fwi_col <- intersect(
    c("InitialFireWeatherIndex", "MeanFWI", "FWI", "InitFWI"),
    names(df)
  )[1]
  if (!is.na(fwi_col)) df$FWI_event <- df[[fwi_col]]

  # Normalise DNBR column
  dnbr_col <- intersect(c("MeanDNBR", "DNBR", "dNBR", "MeanDnbr"), names(df))[1]
  if (!is.na(dnbr_col)) df$MeanDNBR <- df[[dnbr_col]]

  # Normalise duration column
  dur_col <- intersect(c("NumberOfDays", "FireDuration", "Duration",
                          "FireDays", "TotalDays"), names(df))[1]
  if (!is.na(dur_col)) df$duration_days <- df[[dur_col]]

  as_tibble(df)
}


# -----------------------------------------------------------------------------
# plot_val_ignitions
#
# Panel 1: Annual ignition counts.
# Boxplot of expected annual counts (ERA FWI × model) overlaid with simulated
# SCF years as jittered points.  Faceted by ignition type.
# -----------------------------------------------------------------------------

plot_val_ignitions <- function(expected_df, simulated_events_df) {

  if (is.null(simulated_events_df)) return(NULL)

  # Aggregate simulated events to annual counts
  sim_annual <- simulated_events_df %>%
    filter(!is.na(ignition_type), ignition_type %in% c("Lightning", "Accidental")) %>%
    group_by(SimulationYear, ignition_type) %>%
    summarise(n_fires = n(), .groups = "drop")

  # Fill any missing type × year combos with 0
  full_grid <- expand.grid(
    SimulationYear = unique(sim_annual$SimulationYear),
    ignition_type  = c("Lightning", "Accidental"),
    stringsAsFactors = FALSE
  )
  sim_annual <- full_grid %>%
    left_join(sim_annual, by = c("SimulationYear", "ignition_type")) %>%
    mutate(n_fires = replace_na(n_fires, 0L))

  # Totals across types
  sim_total <- sim_annual %>%
    group_by(SimulationYear) %>%
    summarise(n_fires = sum(n_fires), .groups = "drop") %>%
    mutate(ignition_type = "Total")

  sim_plot <- bind_rows(sim_annual, sim_total) %>%
    mutate(ignition_type = factor(ignition_type,
                                   levels = c("Lightning", "Accidental", "Total")))

  p <- ggplot()

  # Add expected range when available (requires ignition calibration)
  has_expected <- !is.null(expected_df) && nrow(expected_df) > 0
  if (has_expected) {
    exp_total <- expected_df %>%
      group_by(year) %>%
      summarise(expected_annual = sum(expected_annual), .groups = "drop") %>%
      mutate(ignition_type = "Total")

    exp_plot <- bind_rows(expected_df, exp_total) %>%
      mutate(ignition_type = factor(ignition_type,
                                     levels = c("Lightning", "Accidental", "Total")))
    p <- p +
      geom_boxplot(
        data  = exp_plot,
        aes(x = ignition_type, y = expected_annual,
            fill = "Expected (ERA FWI climatology)"),
        width = 0.45, outlier.shape = NA, alpha = 0.6
      ) +
      scale_fill_manual(values = c("Expected (ERA FWI climatology)" = .COL_EXP))
  }

  exp_note <- if (!has_expected)
    "\nRun Ignition Calibration + Scale Adjustment to add the expected ERA range." else ""

  p +
    geom_jitter(
      data  = sim_plot,
      aes(x = ignition_type, y = n_fires,
          colour = "SCF simulated years"),
      width = 0.12, size = 2.5, alpha = 0.85
    ) +
    stat_summary(
      data  = sim_plot,
      aes(x = ignition_type, y = n_fires, colour = "SCF simulated years"),
      fun   = mean, geom = "crossbar",
      width = 0.35, linewidth = 0.7, fatten = 1
    ) +
    scale_colour_manual(values = c("SCF simulated years" = .COL_SIM)) +
    facet_wrap(~ignition_type, scales = "free_y", nrow = 1) +
    labs(
      title    = "Annual Ignition Count Validation",
      subtitle = paste0(
        if (has_expected)
          "Green box = expected range from ERA FWI climatology (IQR); "
        else "",
        "red = individual SCF simulation years", exp_note
      ),
      x = NULL, y = "Annual ignitions",
      caption = "Data: FPA FOD (Short 2022); ERA-Interim FWI (Vitolo et al. 2019)"
    ) +
    .pub_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}


# -----------------------------------------------------------------------------
# plot_val_size_ecdf
#
# Panel 2: Empirical cumulative distribution function (ECDF) of fire sizes.
# Observed = FPA FOD fires within calibration boundary.
# Simulated = TotalSitesBurned × cell_area_ha per SCF event.
# Log x-axis; both distributions exclude fires < 1 ha.
# -----------------------------------------------------------------------------

plot_val_size_ecdf <- function(fpa_sizes_ha, sim_events_df,
                                min_ha = 0.4) {

  if (is.null(sim_events_df)) return(NULL)
  if (!"fire_size_ha" %in% names(sim_events_df)) return(NULL)

  sim_fires <- sim_events_df$fire_size_ha[
    !is.na(sim_events_df$fire_size_ha) & sim_events_df$fire_size_ha >= min_ha]

  sim_df <- tibble(
    fire_size_ha = sim_fires,
    source       = sprintf("Simulated \u2014 SCF  (n = %s)",
                           scales::comma(length(sim_fires)))
  )

  colour_vals <- c(setNames(.COL_SIM, unique(sim_df$source)))
  plot_df     <- sim_df

  # Add observed if available
  if (!is.null(fpa_sizes_ha) && length(fpa_sizes_ha) > 0) {
    obs_fires <- fpa_sizes_ha[fpa_sizes_ha >= min_ha & is.finite(fpa_sizes_ha)]
    obs_df    <- tibble(
      fire_size_ha = obs_fires,
      source       = sprintf("Observed \u2014 FPA FOD  (n = %s)",
                             scales::comma(length(obs_fires)))
    )
    plot_df     <- bind_rows(obs_df, sim_df)
    colour_vals <- c(colour_vals, setNames(.COL_OBS, unique(obs_df$source)))
  }

  obs_note <- if (is.null(fpa_sizes_ha))
    "\nRun Ignition Calibration to add observed FPA FOD comparison." else ""

  ggplot(plot_df, aes(fire_size_ha, colour = source)) +
    stat_ecdf(linewidth = 1.0, pad = FALSE) +
    scale_x_log10(
      labels = scales::label_comma(),
      breaks = 10^(0:6)
    ) +
    scale_colour_manual(values = colour_vals) +
    labs(
      title    = "Fire Size Distribution \u2014 ECDF",
      subtitle = paste0(
        "Empirical cumulative distribution of individual fire sizes (\u2265",
        min_ha, " ha).  Curves should overlap if the model reproduces",
        " the observed size distribution.", obs_note
      ),
      x = "Fire size (ha, log\u2081\u2080 scale)",
      y = "Cumulative proportion",
      caption = "Data: FPA FOD (Short 2022); SCF events log"
    ) +
    .pub_theme() +
    theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35))
}


# -----------------------------------------------------------------------------
# plot_val_annual_area
#
# Panel 3: Distribution of annual area burned.
# Boxplot side-by-side: observed FPA FOD vs simulated SCF.
# FPA FOD: sum of FIRE_SIZE_HA per year; SCF: sum of fire_size_ha per SimulationYear.
# -----------------------------------------------------------------------------

plot_val_annual_area <- function(fpa_df, sim_events_df,
                                  year_min = 1992, year_max = 2018) {

  if (is.null(sim_events_df)) return(NULL)
  if (!"fire_size_ha" %in% names(sim_events_df)) return(NULL)

  # SCF annual totals (always available)
  sim_annual <- sim_events_df %>%
    filter(is.finite(fire_size_ha)) %>%
    group_by(year = SimulationYear) %>%
    summarise(area_ha = sum(fire_size_ha, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = "Simulated \u2014 SCF")

  sim_med   <- median(sim_annual$area_ha, na.rm = TRUE)
  plot_df   <- sim_annual
  fill_vals <- c(setNames(.COL_SIM, "Simulated \u2014 SCF"))
  subtitle  <- sprintf("Simulated median: %s ha/yr", scales::comma(round(sim_med)))
  obs_note  <- if (is.null(fpa_df))
    "\nRun Ignition Calibration to add observed FPA FOD comparison." else ""

  # Add observed if available
  if (!is.null(fpa_df) && "FIRE_YEAR" %in% names(fpa_df) &&
      "FIRE_SIZE_HA" %in% names(fpa_df)) {
    obs_annual <- fpa_df %>%
      filter(FIRE_YEAR >= year_min, FIRE_YEAR <= year_max,
             is.finite(FIRE_SIZE_HA)) %>%
      group_by(year = FIRE_YEAR) %>%
      summarise(area_ha = sum(FIRE_SIZE_HA, na.rm = TRUE), .groups = "drop") %>%
      mutate(source = "Observed \u2014 FPA FOD")
    obs_med   <- median(obs_annual$area_ha, na.rm = TRUE)
    plot_df   <- bind_rows(obs_annual, sim_annual)
    fill_vals <- c(fill_vals, setNames(.COL_OBS, "Observed \u2014 FPA FOD"))
    subtitle  <- sprintf(
      "Observed median: %s ha  |  Simulated median: %s ha  |  Ratio: %.2f\u00d7",
      scales::comma(round(obs_med)), scales::comma(round(sim_med)),
      sim_med / obs_med
    )
  }

  ggplot(plot_df, aes(x = source, y = area_ha, fill = source)) +
    geom_boxplot(width = 0.45, outlier.shape = 21, outlier.size = 2,
                 alpha = 0.65) +
    geom_jitter(width = 0.1, size = 1.8, alpha = 0.6, colour = "grey30") +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      title    = "Annual Area Burned",
      subtitle = paste0(subtitle, obs_note),
      x = NULL, y = "Area burned (ha)",
      caption = sprintf(
        "Observed: FPA FOD %d\u2013%d within calibration boundary.  Simulated: all SCF events.",
        year_min, year_max
      )
    ) +
    .pub_theme() +
    theme(legend.position = "none")
}


# -----------------------------------------------------------------------------
# plot_val_fire_climate
#
# Panel 4: Fire size vs FWI relationship — binned median ± IQR approach.
#
# Two sub-panels stacked vertically:
#   Top   : Observed — FPA FOD fires joined to ERA FWI by discovery date.
#            Fire size expressed as % of calibration boundary area so the
#            result is landscape-scale-neutral.
#   Bottom: Simulated — SCF events using InitialFireWeatherIndex column.
#            Fire size expressed as % of LANDIS template area.
#
# Using % of respective landscape removes the spatial scale mismatch
# (calibration region = 14 M ha; LANDIS landscape = 95–200 k ha) and lets
# the FWI slope be compared directly.  The median + IQR per FWI bin is far
# more robust than LOESS on individual scattered fires.
#
# Arguments:
#   fpa_df         data.frame from load_fpa_fod()
#   sim_events_df  from load_scf_events_full()
#   era_fwi_daily  tibble(date, value) from ERA extraction
#   cal_area_ha    calibration boundary area (ha) — for normalisation
#   tmpl_area_ha   LANDIS template area (ha) — for normalisation
#   min_ha         minimum fire size to include
#   fwi_breaks     bin boundaries for FWI; default roughly 0,5,10,...,60+
# -----------------------------------------------------------------------------

plot_val_fire_climate <- function(fpa_df, sim_events_df,
                                   era_fwi_daily,
                                   cal_area_ha  = NULL,
                                   tmpl_area_ha = NULL,
                                   min_ha       = 4,
                                   fwi_breaks   = c(0, 5, 10, 15, 20, 25, 30, 40, 55, Inf)) {

  if (is.null(sim_events_df)) return(NULL)

  # Decide whether to normalise
  normalise <- !is.null(cal_area_ha) && !is.null(tmpl_area_ha) &&
               is.finite(cal_area_ha) && is.finite(tmpl_area_ha) &&
               cal_area_ha > 0 && tmpl_area_ha > 0

  y_label <- if (normalise) "Median fire size (% of landscape, IQR bar)"  else
                             "Median fire size (ha, IQR bar)"
  size_fn_obs <- if (normalise) function(ha) ha / cal_area_ha * 100 else identity
  size_fn_sim <- if (normalise) function(ha) ha / tmpl_area_ha * 100 else identity

  # FWI bin labels
  fwi_labs <- paste0(
    fwi_breaks[-length(fwi_breaks)], "\u2013",
    ifelse(is.finite(fwi_breaks[-1]),
           as.character(fwi_breaks[-1]),
           "+")
  )

  # Helper: bin a data.frame into FWI bins and compute median + IQR
  .bin_fwi <- function(df, fwi_col, size_col, size_fn, source_label) {
    df %>%
      mutate(
        size_norm = size_fn(.data[[size_col]]),
        fwi_bin   = cut(.data[[fwi_col]], breaks = fwi_breaks,
                        labels = fwi_labs, right = FALSE, include.lowest = TRUE)
      ) %>%
      filter(!is.na(fwi_bin), is.finite(size_norm)) %>%
      group_by(fwi_bin) %>%
      summarise(
        median_size = median(size_norm, na.rm = TRUE),
        q25         = quantile(size_norm, 0.25, na.rm = TRUE),
        q75         = quantile(size_norm, 0.75, na.rm = TRUE),
        n           = n(),
        .groups     = "drop"
      ) %>%
      filter(n >= 3) %>%
      mutate(source = source_label)
  }

  # ---- Observed ---------------------------------------------------------------
  obs_binned <- NULL
  if (!is.null(fpa_df) && !is.null(era_fwi_daily) &&
      "FIRE_SIZE_HA" %in% names(fpa_df) && "date" %in% names(fpa_df)) {
    obs_raw <- fpa_df %>%
      filter(is.finite(FIRE_SIZE_HA), FIRE_SIZE_HA >= min_ha, !is.na(date)) %>%
      left_join(era_fwi_daily %>% rename(FWI_obs = value), by = "date") %>%
      filter(is.finite(FWI_obs))

    if (nrow(obs_raw) > 0) {
      obs_binned <- .bin_fwi(obs_raw, "FWI_obs", "FIRE_SIZE_HA",
                              size_fn_obs, "Observed \u2014 FPA FOD")
    }
  }

  # ---- Simulated --------------------------------------------------------------
  sim_binned <- NULL
  if ("fire_size_ha" %in% names(sim_events_df) &&
      "FWI_event" %in% names(sim_events_df)) {
    sim_raw <- sim_events_df %>%
      filter(is.finite(fire_size_ha), fire_size_ha >= min_ha,
             is.finite(FWI_event))

    if (nrow(sim_raw) > 0) {
      sim_binned <- .bin_fwi(sim_raw, "FWI_event", "fire_size_ha",
                              size_fn_sim, "Simulated \u2014 SCF")
    }
  }

  if (is.null(obs_binned) && is.null(sim_binned)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 5, hjust = 0.5,
                      label = paste0(
                        "No FWI data available.\n",
                        "Ensure ERA FWI is loaded and events log has FWI_event column.")) +
             theme_void())
  }

  plot_df <- bind_rows(obs_binned, sim_binned) %>%
    mutate(source = factor(source,
                           levels = c("Observed \u2014 FPA FOD",
                                      "Simulated \u2014 SCF")))

  # Y-axis transform (log10 for ha scale; linear for % scale)
  if (normalise) {
    y_trans  <- "log10"
    y_labels <- function(x) paste0(scales::comma(x, accuracy = 0.001), "%")
  } else {
    y_trans  <- "log10"
    y_labels <- function(x) scales::comma(x)
  }

  scale_label <- if (normalise)
    sprintf("Sizes normalised as %% of landscape (%s = %.0f ha, %s = %.0f ha)",
            "FPA FOD boundary", cal_area_ha, "LANDIS template", tmpl_area_ha)
  else
    paste0("Fire sizes in ha. Note: observed fires from calibration boundary (",
           "~14 M ha); simulated from LANDIS landscape. Use area normalisation for fair comparison.")

  ggplot(plot_df,
         aes(x = fwi_bin, y = median_size, colour = source, group = source)) +
    geom_linerange(aes(ymin = q25, ymax = q75, colour = source),
                   linewidth = 4.5, alpha = 0.25,
                   position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4), linewidth = 0.9) +
    geom_point(size = 3, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c(
      "Observed \u2014 FPA FOD" = .COL_OBS,
      "Simulated \u2014 SCF"    = .COL_SIM
    )) +
    scale_y_continuous(trans = y_trans, labels = y_labels) +
    labs(
      title    = "Fire Size \u2013 Climate Relationship",
      subtitle = paste0(
        "Median fire size per FWI bin (IQR = thick bar).  Only bins with \u22653 fires shown.\n",
        scale_label
      ),
      x       = "Fire Weather Index (FWI) at ignition",
      y       = y_label,
      caption = paste0(
        "Observed FWI: ERA-Interim spatially averaged, matched to FPA FOD discovery date.  ",
        "Simulated FWI: InitialFireWeatherIndex from SCF events log."
      )
    ) +
    .pub_theme() +
    theme(
      axis.text.x    = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}


# -----------------------------------------------------------------------------
# fetch_mtbs_fires
#
# Loads the MTBS burn severity fire occurrence point shapefile — either from
# a user-supplied local path or downloaded from USGS EROS (~30 MB, cached).
#
# The MTBS BSP FOD shapefile (bsp_FODpoints_DD.shp or mtbs_fod_pts_data.shp)
# includes per-fire severity metrics in ALL LOWERCASE field names:
#   dnbr_val   — mean dNBR of the burned area per fire (primary severity metric,
#                directly comparable to SCF's MeanDNBR output)
#   high_t     — dNBR threshold above which pixels are "High" severity (0 for
#                CBI-assessed fires; ignore if 0)
#   mod_t      — moderate severity threshold
#   low_t      — low severity threshold
#   dnbr_offst — offset correction applied to dNBR
#   dnbr_stddv — within-fire dNBR standard deviation
#   burnbndlat / burnbndlon — centroid coordinates (used for spatial filter)
#   ig_date    — ignition date (ISO format "YYYY-MM-DD")
#   incid_type — "Wildfire", "Prescribed Fire", etc.
#   burnbndac  — burned area in acres
#
# There is no state field; spatial filtering is done via lat/lon.
#
# Arguments:
#   year_min / year_max  year filter applied to ig_date
#   cal_vect   SpatVector in any CRS — fires are filtered to within this polygon
#   local_path path to any sidecar file of the shapefile (any of .shp/.dbf/.prj
#              etc. accepted — the function auto-corrects to .shp)
#   cache_dir  directory for download cache
#   progress_fn function(msg) for status messages
# -----------------------------------------------------------------------------

fetch_mtbs_fires <- function(year_min    = 1984,
                              year_max    = 2018,
                              cal_vect    = NULL,
                              local_path  = NULL,
                              cache_dir   = tempdir(),
                              progress_fn = message) {

  # ---- 1. Locate or download the shapefile ----------------------------------
  if (!is.null(local_path) && nchar(trimws(local_path)) > 0) {
    # Auto-correct any shapefile sidecar extension to .shp
    shp_path <- sub("\\.(prj|dbf|shx|cpg|qpj|sbn|sbx)$", ".shp",
                    trimws(local_path), ignore.case = TRUE)
    if (!file.exists(shp_path))
      stop("MTBS file not found: ", shp_path)
    progress_fn(paste("Reading local MTBS file:", basename(shp_path)))
    v <- tryCatch(terra::vect(shp_path),
                  error = function(e) stop("Could not read MTBS file: ", e$message))
  } else {
    zip_url  <- paste0(
      "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/",
      "data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip"
    )
    zip_path <- file.path(cache_dir, "mtbs_fod_pts_data.zip")

    if (!file.exists(zip_path)) {
      progress_fn("Downloading MTBS fire occurrence data (~30 MB)...")
      tryCatch(
        download.file(zip_url, zip_path, mode = "wb", quiet = TRUE),
        error = function(e) stop("MTBS download failed: ", e$message,
                                 "\nCheck internet connection or provide local path.")
      )
    } else {
      progress_fn("Using cached MTBS download...")
    }

    out_dir <- file.path(cache_dir, "mtbs_pts")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    unzip(zip_path, exdir = out_dir, overwrite = FALSE)

    shp  <- list.files(out_dir, pattern = "\\.shp$", recursive = TRUE,
                       full.names = TRUE)[1]
    gpkg <- list.files(out_dir, pattern = "\\.gpkg$", recursive = TRUE,
                       full.names = TRUE)[1]
    src  <- if (!is.na(shp)) shp else if (!is.na(gpkg)) gpkg else
      stop("No shapefile or GeoPackage found in MTBS download.")

    progress_fn(paste("Reading MTBS data:", basename(src)))
    v <- tryCatch(terra::vect(src),
                  error = function(e) stop("Could not read MTBS file: ", e$message))
  }

  nm <- tolower(names(v))          # normalise to lowercase for lookup
  names(v) <- nm
  df <- as.data.frame(v)

  # ---- 2. Standardise / extract key columns --------------------------------

  # Ignition date — already ISO "YYYY-MM-DD" in this file
  date_col <- intersect(c("ig_date", "ignition_date", "start_date"), nm)[1]
  if (!is.na(date_col)) {
    df$Ig_Date <- suppressWarnings(as.Date(as.character(df[[date_col]])))
    df$year    <- lubridate::year(df$Ig_Date)
  }

  # Burned area
  ac_col <- intersect(c("burnbndac", "acres"), nm)[1]
  if (!is.na(ac_col)) {
    df$burn_area_ha <- suppressWarnings(as.numeric(df[[ac_col]])) * 0.404686
  }

  # Primary severity: mean dNBR of the burned area (directly comparable to SCF MeanDNBR)
  if ("dnbr_val" %in% nm) {
    df$MeanDNBR_obs <- suppressWarnings(as.numeric(df[["dnbr_val"]]))
  }

  # Severity thresholds (present in standard MTBS; may be 0 for CBI-only fires)
  for (fld in c("low_t", "mod_t", "high_t", "dnbr_offst", "dnbr_stddv")) {
    if (fld %in% nm)
      df[[fld]] <- suppressWarnings(as.numeric(df[[fld]]))
  }
  # Expose threshold fields with expected capitalised names for downstream code
  if ("high_t"     %in% nm) df$High_T     <- df[["high_t"]]
  if ("mod_t"      %in% nm) df$Mod_T      <- df[["mod_t"]]
  if ("low_t"      %in% nm) df$Low_T      <- df[["low_t"]]
  if ("dnbr_offst" %in% nm) df$dNBR_offst <- df[["dnbr_offst"]]
  if ("dnbr_stddv" %in% nm) df$dNBR_stdDv <- df[["dnbr_stddv"]]

  # Incidence type (filter to wildfires)
  if ("incid_type" %in% nm) {
    df$incid_type_clean <- trimws(df[["incid_type"]])
  }

  # ---- 3. Filter to wildfires only -----------------------------------------
  if ("incid_type_clean" %in% names(df)) {
    wf <- df %>%
      filter(grepl("wildfire", incid_type_clean, ignore.case = TRUE) |
             incid_type_clean == "Unknown")
    if (nrow(wf) > 0) {
      df <- wf
      progress_fn(sprintf("After wildfire filter: %d fires", nrow(df)))
    } else {
      progress_fn("Note: no records matched 'Wildfire' incid_type — keeping all records.")
    }
  }

  # ---- 4. Filter by year ---------------------------------------------------
  if ("year" %in% names(df)) {
    df <- df %>% filter(!is.na(year), year >= year_min, year <= year_max)
    progress_fn(sprintf("After year filter (%d-%d): %d fires",
                        year_min, year_max, nrow(df)))
  }

  # ---- 5. Spatial filter to calibration boundary ---------------------------
  if (!is.null(cal_vect) && nrow(df) > 0) {
    lat_col <- intersect(c("burnbndlat", "latitude", "lat"), nm)[1]
    lon_col <- intersect(c("burnbndlon", "longitude", "lon"), nm)[1]

    if (!is.na(lat_col) && !is.na(lon_col)) {
      lats  <- suppressWarnings(as.numeric(df[[lat_col]]))
      lons  <- suppressWarnings(as.numeric(df[[lon_col]]))
      valid <- is.finite(lats) & is.finite(lons)

      if (any(valid)) {
        pts_v   <- terra::vect(data.frame(lon = lons[valid], lat = lats[valid]),
                               geom = c("lon", "lat"), crs = "EPSG:4326")
        cal_wgs <- terra::project(cal_vect, "EPSG:4326")
        in_bnd  <- terra::relate(pts_v, cal_wgs, "within")
        keep    <- which(valid)[in_bnd]
        df      <- df[keep, ]
        progress_fn(sprintf("After spatial filter: %d fires within calibration boundary",
                            nrow(df)))
      } else {
        progress_fn("Warning: lat/lon columns are all NA — cannot apply spatial filter.")
      }
    } else {
      progress_fn("Warning: lat/lon columns not found — skipping spatial filter.")
    }
  }

  as_tibble(df)
}


# -----------------------------------------------------------------------------
# plot_val_severity
#
# Panel 5: Burn severity distribution — two-panel comparison.
#
# Left panel (always shown):
#   Histogram of SCF MeanDNBR per fire event, coloured by severity class
#   (Miller & Thode 2007 thresholds).  Median annotated.
#
# Right panel (when mtbs_df provided):
#   Distribution of MTBS fire-level High_T (dNBR threshold above which pixels
#   are classified as "High" severity).  A LOWER High_T means severe burning
#   conditions required less extreme dNBR — i.e., more of the landscape reached
#   high severity.  The shading matches the same severity class zones.
#   Where SCF MeanDNBR sits relative to observed MTBS High_T values tells you
#   whether simulated fires are burning at realistic severity levels.
#
# Arguments:
#   sim_events_df  from load_scf_events_full()
#   mtbs_df        optional tibble from fetch_mtbs_fires()
# -----------------------------------------------------------------------------

plot_val_severity <- function(sim_events_df, mtbs_df = NULL) {

  # Miller & Thode (2007) severity class shading (shared by both panels)
  thresholds <- tibble(
    xmin  = c(-Inf, 100, 270, 440, 660),
    xmax  = c(100,  270, 440, 660, Inf),
    label = c("Unburned/\nVery Low", "Low", "Moderate", "High", "Very High"),
    fill  = c("#ffffcc", "#c2e699", "#78c679", "#31a354", "#006837")
  )
  sev_vlines <- c(100, 270, 440, 660)
  sev_labs_x <- c(50, 185, 355, 550, 730)

  .sev_shading <- function() {
    list(
      geom_rect(data = thresholds,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                alpha = 0.15, inherit.aes = FALSE),
      scale_fill_identity(),
      geom_vline(xintercept = sev_vlines,
                 linetype = "dashed", colour = "grey55", linewidth = 0.5)
    )
  }

  # ---- Panel A: SCF MeanDNBR -----------------------------------------------
  has_scf <- !is.null(sim_events_df) && "MeanDNBR" %in% names(sim_events_df)

  if (!has_scf) {
    p_scf <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, size = 5, hjust = 0.5,
               label = "MeanDNBR not found in events log.") +
      theme_void()
  } else {
    sev_df   <- sim_events_df %>% filter(is.finite(MeanDNBR))
    n_events <- nrow(sev_df)
    med_dnbr <- median(sev_df$MeanDNBR, na.rm = TRUE)

    sev_pct <- sev_df %>%
      mutate(class = cut(MeanDNBR,
                         breaks = c(-Inf, 100, 270, 440, 660, Inf),
                         labels = c("Unburned/Very Low","Low","Moderate",
                                    "High","Very High"),
                         right = FALSE)) %>%
      count(class) %>%
      mutate(pct = round(n / sum(n) * 100, 1))

    # Build severity-class % label for subtitle
    pct_str <- sev_pct %>%
      filter(!is.na(class)) %>%
      mutate(lab = paste0(class, " ", pct, "%")) %>%
      pull(lab) %>%
      paste(collapse = "  |  ")

    p_scf <- ggplot(sev_df, aes(MeanDNBR)) +
      .sev_shading() +
      geom_histogram(bins = 40, fill = .COL_SIM, colour = "white",
                     alpha = 0.80, linewidth = 0.3) +
      annotate("text", x = sev_labs_x, y = Inf, vjust = 1.6,
               label = thresholds$label, size = 2.8, colour = "grey30") +
      geom_vline(xintercept = med_dnbr, colour = .COL_SIM,
                 linewidth = 1.0, linetype = "solid") +
      annotate("text", x = med_dnbr + 12, y = Inf, vjust = 2.2,
               label = sprintf("Median\n%.0f", med_dnbr),
               size = 3.0, colour = .COL_SIM, hjust = 0) +
      scale_x_continuous(limits = c(NA, 1300),
                         labels = scales::comma) +
      labs(
        title    = "A. SCF Simulated Fire Severity",
        subtitle = sprintf("%s events  |  %s", scales::comma(n_events), pct_str),
        x = "Mean \u0394NBR per SCF fire event",
        y = "Number of SCF fire events"
      ) +
      .pub_theme()
  }

  # ---- Panel B: MTBS High_T distribution -----------------------------------
  # ---- Panel B: MTBS observed mean dNBR (dnbr_val field) --------------------
  # dnbr_val is the mean dNBR of the burned area per fire — directly comparable
  # to SCF's MeanDNBR.  This is much better than High_T for the comparison.
  has_mtbs <- !is.null(mtbs_df) && "MeanDNBR_obs" %in% names(mtbs_df) &&
              any(is.finite(mtbs_df$MeanDNBR_obs))

  if (!has_mtbs) {
    p_mtbs <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, size = 4.5, hjust = 0.5,
               colour = "grey40",
               label = paste0(
                 "Load MTBS data to enable observed severity comparison.\n",
                 "Provide local path or click 'Fetch MTBS Data'."
               )) +
      theme_void()
  } else {
    mtbs_sev <- mtbs_df %>%
      filter(is.finite(MeanDNBR_obs), MeanDNBR_obs > 0)

    n_mtbs   <- nrow(mtbs_sev)
    med_obs  <- median(mtbs_sev$MeanDNBR_obs, na.rm = TRUE)

    obs_pct <- mtbs_sev %>%
      mutate(class = cut(MeanDNBR_obs,
                         breaks = c(-Inf, 100, 270, 440, 660, Inf),
                         labels = c("Unburned/Very Low","Low","Moderate",
                                    "High","Very High"),
                         right = FALSE)) %>%
      count(class) %>%
      mutate(pct = round(n / sum(n) * 100, 1))

    obs_pct_str <- obs_pct %>%
      filter(!is.na(class)) %>%
      mutate(lab = paste0(class, " ", pct, "%")) %>%
      pull(lab) %>%
      paste(collapse = "  |  ")

    p_mtbs <- ggplot(mtbs_sev, aes(MeanDNBR_obs)) +
      .sev_shading() +
      geom_histogram(bins = 40, fill = .COL_OBS, colour = "white",
                     alpha = 0.80, linewidth = 0.3) +
      annotate("text", x = sev_labs_x, y = Inf, vjust = 1.6,
               label = thresholds$label, size = 2.8, colour = "grey30") +
      geom_vline(xintercept = med_obs, colour = .COL_OBS,
                 linewidth = 1.0, linetype = "solid") +
      annotate("text", x = med_obs + 12, y = Inf, vjust = 2.2,
               label = sprintf("Median\n%.0f", med_obs),
               size = 3.0, colour = .COL_OBS, hjust = 0) +
      scale_x_continuous(limits = c(NA, 1300),
                         labels = scales::comma) +
      labs(
        title    = "B. Observed Burn Severity (MTBS mean dNBR per fire)",
        subtitle = sprintf("%s MTBS fires  |  %s", scales::comma(n_mtbs), obs_pct_str),
        x = "MTBS mean \u0394NBR per fire (dnbr_val)",
        y = "Number of MTBS fires"
      ) +
      .pub_theme()
  }

  # ---- Panel C: Combined violin + boxplot (distribution comparison) ----------
  # Build a combined data frame with both sources for side-by-side comparison.
  # This makes distributions directly comparable regardless of sample-size
  # differences between many short SCF simulations and the MTBS record.
  p_box <- NULL
  if (has_scf && has_mtbs) {
    # n_events / n_mtbs already computed above in Panel A / B respectively
    sev_df_box <- sim_events_df %>%
      filter(is.finite(MeanDNBR)) %>%
      transmute(dNBR   = MeanDNBR,
                source = sprintf("SCF Simulated\n(n\u2009=\u2009%s)",
                                 scales::comma(n_events)))
    mtbs_box <- mtbs_sev %>%
      transmute(dNBR   = MeanDNBR_obs,
                source = sprintf("MTBS Observed\n(n\u2009=\u2009%s)",
                                 scales::comma(n_mtbs)))

    box_df <- bind_rows(mtbs_box, sev_df_box) %>%
      mutate(source = factor(source, levels = c(unique(mtbs_box$source),
                                                 unique(sev_df_box$source))))

    # Severity threshold lines (horizontal since coord_flip)
    sev_hlines <- tibble(
      yintercept = sev_vlines,
      label      = c("Low", "Moderate", "High", "Very High")
    )

    p_box <- ggplot(box_df, aes(x = source, y = dNBR, fill = source)) +
      geom_violin(alpha = 0.35, width = 0.75, colour = NA, trim = TRUE) +
      geom_boxplot(width = 0.18, outlier.shape = NA,
                   colour = "grey20", alpha = 0.75) +
      geom_hline(data = sev_hlines,
                 aes(yintercept = yintercept),
                 linetype = "dashed", colour = "grey55", linewidth = 0.45,
                 inherit.aes = FALSE) +
      geom_text(data = sev_hlines,
                aes(x = Inf, y = yintercept, label = label),
                hjust = 1.08, vjust = -0.4, size = 2.8, colour = "grey40",
                inherit.aes = FALSE) +
      coord_flip() +
      scale_fill_manual(values = c(
        setNames(.COL_OBS, unique(mtbs_box$source)),
        setNames(.COL_SIM, unique(sev_df_box$source))
      )) +
      scale_y_continuous(limits = c(NA, 1400), labels = scales::comma) +
      labs(
        title    = "C. Severity Distribution Comparison (Violin + Boxplot)",
        subtitle = paste0(
          "Box shows IQR and median; violin shows full density.  ",
          "Dashed lines = severity class boundaries (Miller & Thode 2007)."
        ),
        x = NULL,
        y = "Mean \u0394NBR per fire"
      ) +
      .pub_theme() +
      theme(legend.position = "none",
            axis.text.y     = element_text(size = 10, face = "bold"))
  }

  # ---- Combine panels -------------------------------------------------------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    top_row <- patchwork::wrap_plots(p_scf, p_mtbs, ncol = 2)

    if (!is.null(p_box)) {
      combined <- patchwork::wrap_plots(top_row, p_box, ncol = 1,
                                         heights = c(1.4, 1))
    } else {
      combined <- top_row
    }

    combined <- combined +
      patchwork::plot_annotation(
        title    = "Burn Severity Distribution: Simulated vs Observed",
        subtitle = paste0(
          "A & B: Histograms of mean \u0394NBR per fire on the same axis and severity scale.  ",
          "C: Direct distribution comparison regardless of sample-size differences.  ",
          "If distributions align, the model is simulating realistic severity levels."
        ),
        caption = paste0(
          "Severity thresholds: Miller & Thode (2007) Ecol. Indic.  ",
          "MTBS data: Eidenshink et al. (2007); USGS EROS Center."
        ),
        theme = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey40", size = 10),
          plot.caption  = element_text(colour = "grey50", size = 9, hjust = 0)
        )
      )
    return(combined)
  }

  # Fallback: return SCF panel only (patchwork not installed)
  p_scf + labs(
    title   = "Burn Severity Distribution (SCF Simulated)",
    caption = paste0(
      "Install 'patchwork' for the MTBS comparison.  ",
      "dNBR thresholds: <100 Unburned/VL, 100\u2013270 Low, 270\u2013440 Mod, 440\u2013660 High, >660 VH."
    )
  )
}


# -----------------------------------------------------------------------------
# plot_val_duration
#
# Panel 6: Fire duration distribution.
# Simulated: histogram + ECDF of NumberOfDays per SCF event.
# Faceted by ignition type.
# -----------------------------------------------------------------------------

plot_val_duration <- function(sim_events_df) {

  if (is.null(sim_events_df)) return(NULL)
  if (!"duration_days" %in% names(sim_events_df)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 5, hjust = 0.5,
                      label = "NumberOfDays column not found in events log.") +
             theme_void())
  }

  dur_df <- sim_events_df %>%
    filter(is.finite(duration_days), duration_days > 0)

  if (nrow(dur_df) == 0) return(NULL)

  # Summary stats by ignition type
  sum_df <- dur_df %>%
    group_by(ignition_type) %>%
    summarise(
      med   = median(duration_days),
      p90   = quantile(duration_days, 0.90),
      n     = n(),
      .groups = "drop"
    )

  ggplot(dur_df, aes(duration_days)) +
    geom_histogram(bins = 30, fill = .COL_SIM, colour = "white",
                   alpha = 0.80, linewidth = 0.3) +
    geom_vline(data = sum_df, aes(xintercept = med),
               colour = .COL_EXP, linewidth = 1.1, linetype = "solid") +
    geom_vline(data = sum_df, aes(xintercept = p90),
               colour = "#e67e22", linewidth = 0.9, linetype = "dashed") +
    geom_text(data = sum_df,
              aes(x = med, y = Inf, label = sprintf("Median: %.1f d", med)),
              vjust = 2, hjust = -0.1, size = 3.2, colour = .COL_EXP) +
    geom_text(data = sum_df,
              aes(x = p90, y = Inf, label = sprintf("P90: %.1f d", p90)),
              vjust = 3.5, hjust = -0.1, size = 3.2, colour = "#e67e22") +
    facet_wrap(~ignition_type, scales = "free_y") +
    labs(
      title    = "Fire Duration Distribution",
      subtitle = "SCF simulated fires \u2014 green line = median, orange dashed = 90th percentile",
      x = "Fire duration (days)",
      y = "Number of fire events",
      caption  = "Source: SCF socialclimatefire-events-log.csv"
    ) +
    .pub_theme()
}


# -----------------------------------------------------------------------------
# plot_val_fire_count_vs_fwi
#
# Supplemental: Annual fire count vs mean summer FWI.
# Shows whether the climate-fire relationship is preserved.
# Observed: FPA FOD annual counts vs mean annual FWI.
# Simulated: SCF annual counts vs mean annual FWI (by simulation year index).
#
# When cal_area_ha and tmpl_area_ha are supplied, BOTH series are expressed as
# fires per million ha of their respective landscapes — this removes the spatial
# scale mismatch and lets the FWI slope be compared directly on a single axis.
# Without normalisation the two series are on completely different numeric
# scales (e.g. ~1,500 fires/yr in a 14 M ha calibration region vs ~25/yr in
# a 95 k ha LANDIS landscape) and the simulated line will appear flat at 0.
# -----------------------------------------------------------------------------

plot_val_fire_count_vs_fwi <- function(fpa_df, sim_events_df,
                                        era_fwi_daily,
                                        year_min     = 1992, year_max = 2018,
                                        cal_area_ha  = NULL,
                                        tmpl_area_ha = NULL) {

  if (is.null(era_fwi_daily)) return(NULL)

  normalise <- !is.null(cal_area_ha) && !is.null(tmpl_area_ha) &&
               is.finite(cal_area_ha) && is.finite(tmpl_area_ha) &&
               cal_area_ha > 0 && tmpl_area_ha > 0

  # Mean summer (Jun-Sep) FWI per year from ERA
  summer_fwi <- era_fwi_daily %>%
    mutate(year = lubridate::year(date),
           month = lubridate::month(date)) %>%
    filter(month %in% 6:9, year >= year_min, year <= year_max) %>%
    group_by(year) %>%
    summarise(mean_summer_fwi = mean(value, na.rm = TRUE), .groups = "drop")

  # Observed annual fire counts joined to FWI
  obs_df <- NULL
  if (!is.null(fpa_df) && "FIRE_YEAR" %in% names(fpa_df)) {
    obs_df <- fpa_df %>%
      filter(FIRE_YEAR >= year_min, FIRE_YEAR <= year_max) %>%
      count(year = FIRE_YEAR) %>%
      left_join(summer_fwi, by = "year") %>%
      filter(is.finite(mean_summer_fwi)) %>%
      mutate(
        count_raw = n,
        n         = if (normalise) n / (cal_area_ha / 1e6) else n,
        source    = "Observed \u2014 FPA FOD"
      )
  }

  # Simulated: assign ERA years cyclically to simulation years
  if (!is.null(sim_events_df) && "SimulationYear" %in% names(sim_events_df)) {
    sim_yrs    <- sort(unique(sim_events_df$SimulationYear))
    era_yrs    <- summer_fwi$year
    era_cycled <- era_yrs[(seq_along(sim_yrs) - 1) %% length(era_yrs) + 1]

    sim_yr_map <- tibble(SimulationYear = sim_yrs,
                         era_year       = era_cycled)

    sim_annual <- sim_events_df %>%
      count(SimulationYear) %>%
      left_join(sim_yr_map, by = "SimulationYear") %>%
      left_join(summer_fwi, by = c("era_year" = "year")) %>%
      filter(is.finite(mean_summer_fwi)) %>%
      transmute(
        year             = SimulationYear,
        count_raw        = n,
        n                = if (normalise) n / (tmpl_area_ha / 1e6) else n,
        mean_summer_fwi  = mean_summer_fwi,
        source           = "Simulated \u2014 SCF"
      )

    plot_df <- bind_rows(obs_df, sim_annual)
  } else {
    plot_df <- obs_df
  }

  if (is.null(plot_df) || nrow(plot_df) == 0) return(NULL)

  y_label   <- if (normalise) "Annual fires per million ha" else "Annual fire count"
  scale_note <- if (normalise)
    sprintf(
      "Both series normalised to fires / M ha of respective landscape\n(FPA FOD boundary = %.1f M ha; LANDIS template = %.0f ha)",
      cal_area_ha / 1e6, tmpl_area_ha
    )
  else
    paste0(
      "Note: FPA FOD fires cover the full calibration boundary (~14 M ha); ",
      "SCF fires cover only the LANDIS landscape.\n",
      "Run Ignition Calibration \u2192 Scale Adjustment to enable area normalisation."
    )

  ggplot(plot_df, aes(mean_summer_fwi, n, colour = source)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.15) +
    scale_colour_manual(values = c(
      "Observed \u2014 FPA FOD" = .COL_OBS,
      "Simulated \u2014 SCF"    = .COL_SIM
    )) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      title    = "Annual Fire Count vs Mean Summer FWI",
      subtitle = paste0(
        "Observed (", year_min, "\u2013", year_max, ") and simulated years.  ",
        "Regression lines show climate\u2013fire sensitivity.\n",
        scale_note
      ),
      x = "Mean summer (Jun\u2013Sep) FWI",
      y = y_label,
      caption = "Summer FWI: ERA-Interim spatial mean over calibration boundary."
    ) +
    .pub_theme()
}


# -----------------------------------------------------------------------------
# save_validation_figures
# Saves all validation figures to a temporary directory at 300 dpi.
# Returns a named list of file paths.
# -----------------------------------------------------------------------------

save_validation_figures <- function(
    fig_list,          # named list of ggplot objects (NULL elements are skipped)
    width_in  = 7,
    height_in = 5,
    dpi       = 300
  ) {

  out_dir <- tempfile("val_figs_")
  dir.create(out_dir, recursive = TRUE)

  paths <- list()
  for (nm in names(fig_list)) {
    g <- fig_list[[nm]]
    if (is.null(g)) next
    fname <- file.path(out_dir, paste0(nm, ".png"))
    tryCatch(
      ggplot2::ggsave(fname, g,
                      width = width_in, height = height_in,
                      dpi   = dpi, bg = "white"),
      error = function(e) message("Could not save ", nm, ": ", e$message)
    )
    if (file.exists(fname)) paths[[nm]] <- fname
  }
  paths
}
