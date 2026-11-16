# =============================================================================
# tuning_helpers.R
# Interactive parameter tuning plots for the SCF calibration app.
#
# All functions take raw numeric coefficients (no reactive objects) and return
# ggplot objects or plain R values, making them easy to test independently.
#
# Functions:
#   simulate_fire_simple()   One SCF fire event, FireEvent.cs spread rules
#   plot_msa_surface()       MaxSpreadArea vs FWI at multiple EWS levels
#   plot_sp_surface()        SpreadProbability vs FWI
#   plot_sp_ff_sensitivity() P vs FWI at FF=0 vs FF=1 (shows mismatch risk)
#   plot_ff_scaling()        FF_scaled vs raw g C/m² and resulting P(spread)
#                            — shows the saturation problem when MaxFF is too low
#   plot_dnbr_surface()      dNBR vs EWS at multiple LadderFuel levels
#   dnbr_sensitivity_table() Denominator / dNBR tabulation across LF values
#   plot_cm_curve()          P(cohort mort) vs BarkThickness at multiple dNBR
#   plot_fire_map()          Fire grid from simulate_fire_simple() as ggplot
#   plot_fire_timeline()     Area burned per day vs MaxSpreadArea
#   run_fire_ensemble_weather() Fires on real daily weather, SCF-weighted ignition days
#   plot_fire_ensemble_weather() Size CDF vs FPA FOD + size vs ignition FWI
#   sweep_b0_weather()       SpreadProbabilityB0 by simulation (size + annual area)
#   plot_b0_sweep()          Scores and fire sizes across the B0 sweep
#   ign_tune_coef_df()       Hand-typed ignition coefficients -> coef tibble
#   build_scf_tuning_snippet() SCF parameter file text from tuning inputs
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, tidyr, scales)


# =============================================================================
# Fire spread simulation
# =============================================================================

# -----------------------------------------------------------------------------
# simulate_fire_simple
#
# Single SCF fire event on a square grid of active cells, replicating the
# spread loop in FireEvent.cs (Extension-Social-Climate-Fire master, 2026-04):
#
#   * FIFO queue.  The ignition cell is queued first and marked Disturbed
#     (FireEvent.cs:105).  Each popped cell burns, then tries its 4 cardinal
#     neighbours (N, E, S, W) and then its 4 diagonals, skipping any neighbour
#     already Disturbed (:679, :707, :735).
#   * ONE TRIAL PER CELL.  CanSpread() marks the neighbour Disturbed whether or
#     not it ignites (:456), so a cell that fails is never tried again.  This
#     is single-trial (site) percolation: fire only propagates when P clears
#     roughly 0.4 for this 8-neighbour lattice.
#   * Diagonal tries use P * 0.71 (:537).  Which neighbour tries a cell, and so
#     which weight applies, is decided by queue order exactly as in SCF.
#   * MaxSpreadArea is PACING, not a ceiling (:204-216).  After each burned
#     cell the day's area grows by one cell; once it exceeds
#     (int)(B0 + B1*FWI + B2*EWS) the day advances and the fire keeps burning.
#     So each day burns floor(MSA / cell_area) + 1 cells.  Nothing is discarded.
#   * The loop runs while day <= max_days and a cell only spreads while
#     day < max_days (:166, :226; SCF uses DaysPerYear = 364).
#   * FWI and EWS may be single values or DAILY VECTORS (element d is day d;
#     the last value is reused past the end), as SCF reads daily climate.  The
#     weather of the day a cell burns sets both its spread trials and MSA.
#
# Not modelled: suppression, BurningSitesThreshold (absent from SCF master),
# Cb escalation of EWS with intensity, spatially varying fuels/terrain, and
# cells already Disturbed by earlier fires in the same year.
#
# Implementation: the queue is processed one generation at a time (all cells
# queued by the previous generation, in queue order).  Each untried neighbour
# is tried once, by the earliest cell in queue order that reaches it, which is
# exactly the first CanSpread() call SCF would make for it.  Same result as a
# per-cell loop, at vector speed.
#
# Returns a list:
#   total_ha, days, msa_ha (day-1 MSA), p_spread (day-1 cardinal P), p_diag,
#   cells_per_day (day 1), n_burned, n_tried, hit_edge (fire reached the grid
#   edge, so the size is truncated by the grid), timeline
#   data.frame(day, new_cells, total_cells), grid (0 = untouched,
#   1 = burned on the final day, 2 = burned, 3 = tried but did not ignite),
#   fwi, ews, ff (day-1 values), grid_dim, cell_area_ha.
# -----------------------------------------------------------------------------

simulate_fire_simple <- function(grid_dim     = 100,
                                  fwi          = 20,
                                  ews          = 15,
                                  ff           = 0.5,
                                  msa_b0, msa_b1, msa_b2,
                                  sp_b0,  sp_b1,  sp_b2,  sp_b3,
                                  cell_area_ha = 0.09,
                                  max_days     = 60,
                                  seed         = 42) {

  if (!is.null(seed)) set.seed(seed)
  grid_dim <- as.integer(grid_dim)
  n_cells  <- grid_dim * grid_dim
  max_days <- max(1L, as.integer(max_days))

  # Daily weather lookups (recycle the last value past the end)
  wx <- function(v, d) v[pmin(d, length(v))]
  p_of_day <- function(d) {
    1 / (1 + exp(-(sp_b0 + sp_b1 * wx(fwi, d) + sp_b2 * ff + sp_b3 * wx(ews, d))))
  }
  msa_of_day <- function(d) {
    trunc(msa_b0 + msa_b1 * wx(fwi, d) + msa_b2 * wx(ews, d))   # C# (int) cast
  }
  # Cells burned per day: smallest k with k * cell_area > MSA  (k >= 1)
  k_of_day <- function(d) {
    m <- msa_of_day(d)
    ifelse(m < 0, 1L, as.integer(floor(m / cell_area_ha)) + 1L)
  }

  # Neighbour offsets in SCF's order: N, E, S, W, then the four diagonals
  # (Get4CardinalActiveNeighbors / Get4DiagonalNeighbors).
  d_row  <- c(-1L, 0L, 1L,  0L, -1L, 1L,  1L, -1L)
  d_col  <- c( 0L, 1L, 0L, -1L,  1L, 1L, -1L, -1L)
  weight <- c(1, 1, 1, 1, 0.71, 0.71, 0.71, 0.71)

  disturbed <- logical(n_cells)
  burn_day  <- integer(n_cells)        # 0 = not burned
  ign       <- (grid_dim %/% 2 - 1L) * grid_dim + grid_dim %/% 2   # centre
  disturbed[ign] <- TRUE

  gen      <- ign         # current generation, in queue order
  day      <- 1L          # SCF's `day`
  day_area <- 0L          # cells burned so far today
  n_tried  <- 0L
  hit_edge <- FALSE

  while (length(gen) > 0 && day <= max_days) {

    # ---- Pop the generation: assign burn days, advance `day` ---------------
    m         <- length(gen)
    pre_day   <- integer(m)      # day the cell burns (weather for its trials)
    post_day  <- integer(m)      # day after any increment (gates spreading)
    i <- 1L
    while (i <= m && day <= max_days) {
      k    <- k_of_day(day)
      take <- min(m - i + 1L, k - day_area)   # cells that fit in today
      idx  <- i:(i + take - 1L)
      pre_day[idx]  <- day
      post_day[idx] <- day
      day_area <- day_area + take
      i <- i + take
      if (day_area >= k) {                    # area > MSA after the k-th cell
        post_day[idx[take]] <- day + 1L
        day      <- day + 1L
        day_area <- 0L
      }
    }
    popped <- seq_len(i - 1L)                 # cells actually processed
    if (length(popped) == 0) break
    gen      <- gen[popped]
    pre_day  <- pre_day[popped]
    post_day <- post_day[popped]
    burn_day[gen] <- pre_day

    # ---- Spread: each untried neighbour gets one try, from its first reacher
    spreaders <- which(post_day < max_days)
    if (length(spreaders) == 0) break
    src  <- gen[spreaders]
    srow <- (src - 1L) %% grid_dim + 1L
    scol <- (src - 1L) %/% grid_dim + 1L

    # All (source, direction) pairs in SCF's processing order:
    # source queue order first, then direction rank 1..8.
    nr  <- rep(srow, each = 8) + rep(d_row, times = length(src))
    nc  <- rep(scol, each = 8) + rep(d_col, times = length(src))
    dir <- rep(1:8, times = length(src))
    who <- rep(spreaders, each = 8)
    inside <- nr >= 1L & nr <= grid_dim & nc >= 1L & nc <= grid_dim
    if (!all(inside)) hit_edge <- TRUE
    nb  <- (nc[inside] - 1L) * grid_dim + nr[inside]
    dir <- dir[inside]; who <- who[inside]

    keep <- !disturbed[nb] & !duplicated(nb)      # first reacher wins
    nb   <- nb[keep]; dir <- dir[keep]; who <- who[keep]
    if (length(nb) == 0) { gen <- integer(0); break }

    disturbed[nb] <- TRUE                         # tried, burn or not
    n_tried <- n_tried + length(nb)
    p   <- p_of_day(pre_day[who]) * weight[dir]
    ok  <- runif(length(nb)) < p
    gen <- nb[ok]                                 # already in queue order
  }

  burned  <- which(burn_day > 0L)
  n_burn  <- length(burned)
  last_d  <- if (n_burn) max(burn_day) else 0L
  tl <- if (n_burn) {
    nd <- tabulate(burn_day[burned], nbins = last_d)
    data.frame(day = seq_len(last_d), new_cells = nd, total_cells = cumsum(nd))
  } else data.frame(day = integer(0), new_cells = integer(0), total_cells = integer(0))

  grid <- matrix(0L, grid_dim, grid_dim)
  grid[disturbed]  <- 3L
  grid[burned]     <- 2L
  grid[burn_day == last_d & burn_day > 0L] <- 1L

  list(
    total_ha      = n_burn * cell_area_ha,
    days          = last_d,
    msa_ha        = msa_of_day(1L),
    p_spread      = p_of_day(1L),
    p_diag        = 0.71 * p_of_day(1L),
    cells_per_day = k_of_day(1L),
    n_burned      = n_burn,
    n_tried       = n_tried,
    hit_edge      = hit_edge,
    timeline      = tl,
    grid          = grid,
    fwi           = fwi[1],
    ews           = ews[1],
    ff            = ff,
    grid_dim      = grid_dim,
    cell_area_ha  = cell_area_ha
  )
}


# =============================================================================
# Plot: MaxSpreadArea surface
# =============================================================================

plot_msa_surface <- function(msa_b0, msa_b1, msa_b2,
                               fwi_max = 50,
                               ews_vals = c(5, 15, 25, 35)) {

  df <- expand.grid(FWI = seq(0, fwi_max, 1), EWS = ews_vals) |>
    dplyr::mutate(
      MSA_raw     = msa_b0 + msa_b1 * FWI + msa_b2 * EWS,
      MSA         = pmax(MSA_raw, 0),
      is_negative = MSA_raw < 0,
      EWS_label   = factor(paste0("EWS = ", EWS, " km/h"))
    )

  b1_sign <- if (msa_b1 < 0) "⚠ B1 NEGATIVE — higher FWI shrinks ceiling (wrong)" else
             if (msa_b1 == 0) "⚠ B1 = 0" else
             "✓ B1 positive"
  b2_sign <- if (msa_b2 < 0) "  |  ⚠ B2 NEGATIVE" else "  |  ✓ B2 positive"
  sub_txt <- paste0(
    sprintf("B0=%.1f  B1=%.3f  B2=%.3f   ", msa_b0, msa_b1, msa_b2),
    b1_sign, b2_sign
  )

  p <- ggplot(df, aes(FWI, MSA, colour = EWS_label)) +
    # plausible range shading
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 50, ymax = 800,
             fill = "#2c5f2e", alpha = 0.07) +
    annotate("text", x = fwi_max * 0.97, y = 425,
             label = "Plausible range\n50–800 ha/day",
             hjust = 1, size = 3, colour = "#2c5f2e") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 1.1) +
    # mark where MSA was clamped from negative
    geom_point(data = dplyr::filter(df, is_negative),
               aes(y = 0), shape = 4, size = 3, colour = "#c0392b") +
    scale_colour_manual(values = c("#9ecae1", "#4292c6", "#2171b5", "#084594")) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
    labs(
      title    = "Maximum Spread Area vs FWI",
      subtitle = sub_txt,
      x        = "FWI (fire weather index)",
      y        = "Max daily spread area (ha)",
      colour   = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9,
            colour = if (msa_b1 < 0) "#c0392b" else "grey40"))
  p
}


# =============================================================================
# Plot: SpreadProbability vs FWI
# =============================================================================

plot_sp_surface <- function(sp_b0, sp_b1, sp_b2, sp_b3,
                              ff_val  = 0.5,
                              fwi_max = 50,
                              ews_vals = c(5, 15, 25, 35)) {

  df <- expand.grid(FWI = seq(0, fwi_max, 0.5), EWS = ews_vals) |>
    dplyr::mutate(
      logit    = sp_b0 + sp_b1 * FWI + sp_b2 * ff_val + sp_b3 * EWS,
      Pspread  = 1 / (1 + exp(-logit)),
      EWS_label = factor(paste0("EWS = ", EWS, " km/h"))
    )

  # FWI at which P crosses 0.5 for mean EWS
  mid_ews <- ews_vals[ceiling(length(ews_vals) / 2)]
  fwi_50  <- tryCatch({
    b_eff <- sp_b0 + sp_b2 * ff_val + sp_b3 * mid_ews
    if (abs(sp_b1) < 1e-10) NA_real_ else -b_eff / sp_b1
  }, error = function(e) NA_real_)

  sub_txt <- sprintf(
    "B0=%.3f  B1=%.4f  B2=%.4f  B3=%.4f   FF=%.2f   P=0.5 at FWI≈%.1f (EWS=%d)",
    sp_b0, sp_b1, sp_b2, sp_b3, ff_val,
    if (is.finite(fwi_50)) fwi_50 else 999, mid_ews
  )

  ggplot(df, aes(FWI, Pspread, colour = EWS_label)) +
    # Plausible Cascades range
    annotate("rect", xmin = 10, xmax = 40, ymin = 0.05, ymax = 0.35,
             fill = "#2c5f2e", alpha = 0.07) +
    annotate("text", x = 25, y = 0.36,
             label = "Plausible range (Cascades)", size = 3, colour = "#2c5f2e") +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = c("#fdae6b", "#f16913", "#d94801", "#7f2704")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title    = "Spread Probability vs FWI",
      subtitle = sub_txt,
      x        = "FWI",
      y        = "P(spread to neighbour cell)",
      colour   = NULL,
      caption  = sprintf("FineFuels fixed at %.2f", ff_val)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 9, colour = "grey40"))
}


# =============================================================================
# Plot: SpreadProbability FF sensitivity (FF=0 vs FF=1)
# =============================================================================

plot_sp_ff_sensitivity <- function(sp_b0, sp_b1, sp_b2, sp_b3,
                                    fwi_max  = 50,
                                    ews_fixed = 15) {

  ff_levels <- c(0, 0.25, 0.5, 0.75, 1.0)
  df <- expand.grid(FWI = seq(0, fwi_max, 0.5), FF = ff_levels) |>
    dplyr::mutate(
      Pspread  = 1 / (1 + exp(-(sp_b0 + sp_b1 * FWI + sp_b2 * FF + sp_b3 * ews_fixed))),
      FF_label = factor(sprintf("FF = %.2f", FF))
    )

  # Compute the shift caused by FF going 0 → 1
  p_ff0 <- 1 / (1 + exp(-(sp_b0 + sp_b1 * 20 + sp_b2 * 0   + sp_b3 * ews_fixed)))
  p_ff1 <- 1 / (1 + exp(-(sp_b0 + sp_b1 * 20 + sp_b2 * 1.0 + sp_b3 * ews_fixed)))
  delta_p <- p_ff1 - p_ff0

  sub_txt <- sprintf(
    "EWS fixed at %d km/h | B2=%.4f | At FWI=20: P(FF=0)=%.3f, P(FF=1)=%.3f, Δ=%.3f",
    ews_fixed, sp_b2, p_ff0, p_ff1, delta_p
  )
  warn_col <- if (abs(delta_p) > 0.15 || sp_b2 < 0) "#c0392b" else "grey40"

  ggplot(df, aes(FWI, Pspread, colour = FF_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 1.1) +
    scale_colour_brewer(palette = "YlOrRd", direction = 1) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title    = "Spread Probability — FineFuels Sensitivity",
      subtitle = sub_txt,
      x        = "FWI",
      y        = "P(spread to neighbour cell)",
      colour   = "Fine Fuels",
      caption  = paste0(
        if (sp_b2 < 0) "⚠ B2 is negative — higher fuels REDUCE spread (likely unidentified from FF=1 placeholder). " else "",
        "Wide spread between FF=0 and FF=1 lines means FF matters a lot at runtime. ",
        "If calibrated with FF=1 placeholder, set B2=0 and adjust B0 = B0 + B2×1.0."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, colour = warn_col),
          plot.caption  = element_text(size = 8, colour = warn_col))
}


# =============================================================================
# Plot: MaximumFineFuels scaling — raw g C/m² → FF_scaled (0-1)
# =============================================================================

# -----------------------------------------------------------------------------
# plot_ff_scaling
#
# Two-panel figure showing:
#   Top: FF_scaled = min(FF_raw / MaxFF, 1.0) as a function of raw g C/m².
#        Vertical line at MaxFF (saturation point). Optional line at the
#        user-reported maximum observed LANDIS output.
#        Annotation shows what % of the landscape is already at FF_scaled=1.0
#        if nearly all cells exceed MaxFF.
#
#   Bottom: P(spread) at fixed FWI and EWS across the raw g C/m² range.
#           Shows how much spread probability actually varies given the
#           current MaxFF setting vs. a higher setting.
#
# Key diagnostic:
#   If max_ff << observed_landis_max, most cells will have FF_scaled ≈ 1.0.
#   B2 then acts as a constant offset (same as the FF=1 calibration placeholder
#   problem) and provides no real fuel loading signal at runtime.
#   Recommendation: set MaxFF to the 95th-99th percentile of your LANDIS
#   FineFuels output so the full observed range maps to 0-1.
#
# Arguments:
#   max_ff           Current MaximumFineFuels setting (g C/m²)
#   sp_b0..sp_b3     SpreadProbability coefficients
#   observed_max     User-reported max raw LANDIS FineFuels (g C/m²); optional
#   fwi_fixed        FWI for P(spread) bottom panel (default 20)
#   ews_fixed        EWS for P(spread) bottom panel (default 15 km/h)
# -----------------------------------------------------------------------------

plot_ff_scaling <- function(max_ff,
                             sp_b0, sp_b1, sp_b2, sp_b3,
                             observed_max = NULL,
                             fwi_fixed    = 20,
                             ews_fixed    = 15) {

  # Determine x-axis range: at least 2× MaxFF, or up to 1.1× observed max
  x_max <- max(max_ff * 2.5,
               if (is.numeric(observed_max) && observed_max > 0) observed_max * 1.05
               else 0,
               100)
  x_seq <- seq(0, x_max, length.out = 600)

  df <- data.frame(
    FF_raw    = x_seq,
    FF_scaled = pmin(x_seq / max_ff, 1.0),
    P_spread  = 1 / (1 + exp(-(sp_b0 + sp_b1 * fwi_fixed +
                                sp_b2 * pmin(x_seq / max_ff, 1.0) +
                                sp_b3 * ews_fixed)))
  )

  # Saturation diagnosis
  pct_saturated <- if (is.numeric(observed_max) && observed_max > 0 && max_ff < observed_max)
    round(100 * (1 - max_ff / observed_max), 1) else 0

  # P(spread) range across the raw-FF dimension
  p_at_ff0   <- 1 / (1 + exp(-(sp_b0 + sp_b1 * fwi_fixed + sp_b2 * 0   + sp_b3 * ews_fixed)))
  p_at_maxff <- 1 / (1 + exp(-(sp_b0 + sp_b1 * fwi_fixed + sp_b2 * 1.0 + sp_b3 * ews_fixed)))

  warn_col <- if (pct_saturated > 50) "#c0392b" else
              if (pct_saturated > 20) "#e67e22" else "grey40"

  sub_top <- sprintf(
    "MaximumFineFuels = %.0f g C/m²  |  Saturation at FF_scaled = 1.0 (FF_raw ≥ MaxFF)%s",
    max_ff,
    if (pct_saturated > 0)
      sprintf("  |  ⚠ %.1f%% of LANDIS range [0–%.0f] is already saturated",
              pct_saturated, observed_max)
    else ""
  )

  sub_bot <- sprintf(
    "FWI=%.0f  EWS=%.0f km/h  |  P spans %.3f (FF_raw=0) → %.3f (FF_raw=MaxFF)  |  ΔP = %.3f",
    fwi_fixed, ews_fixed, p_at_ff0, p_at_maxff, abs(p_at_maxff - p_at_ff0)
  )

  # ---- Top panel: scaling curve -------------------------------------------
  p_top <- ggplot(df, aes(FF_raw, FF_scaled)) +
    geom_line(colour = "#2171b5", linewidth = 1.2) +
    # Saturation cap line
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    # MaxFF vertical line
    geom_vline(xintercept = max_ff, linetype = "longdash",
               colour = "#2171b5", linewidth = 0.7) +
    annotate("text", x = max_ff, y = 0.5,
             label = sprintf("MaxFF = %.0f", max_ff),
             hjust = -0.08, size = 3.2, colour = "#2171b5") +
    # Observed max line (if supplied)
    {if (is.numeric(observed_max) && observed_max > 0)
      list(
        geom_vline(xintercept = observed_max, linetype = "dotted",
                   colour = "#c0392b", linewidth = 0.7),
        annotate("text", x = observed_max, y = 0.25,
                 label = sprintf("LANDIS max ≈ %.0f", observed_max),
                 hjust = 1.05, size = 3.2, colour = "#c0392b")
      )
    } +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title    = "Fine Fuels Scaling: Raw g C/m² → FF_scaled (0–1)",
      subtitle = sub_top,
      x        = "Raw LANDIS FineFuels (g C/m²)",
      y        = "FF_scaled = min(FF / MaxFF, 1)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, colour = warn_col))

  # ---- Bottom panel: P(spread) across the raw-FF range --------------------
  p_bot <- ggplot(df, aes(FF_raw, P_spread)) +
    geom_line(colour = "#d94801", linewidth = 1.2) +
    geom_vline(xintercept = max_ff, linetype = "longdash",
               colour = "#2171b5", linewidth = 0.7) +
    {if (is.numeric(observed_max) && observed_max > 0)
      geom_vline(xintercept = observed_max, linetype = "dotted",
                 colour = "#c0392b", linewidth = 0.7)
    } +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      subtitle = sub_bot,
      x        = "Raw LANDIS FineFuels (g C/m²)",
      y        = "P(spread)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey40"))

  # Stack with patchwork if available, otherwise just return top panel
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p_top, p_bot, ncol = 1, heights = c(1.1, 1))
  } else {
    p_top
  }
}


# =============================================================================
# Table + Plot: Site Mortality / dNBR sensitivity to LadderFuels
# =============================================================================

# Returns a data.frame showing denominator and dNBR at each LF value
dnbr_sensitivity_table <- function(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                    clay, et_mm, ews_kmh, cwd_mm, ff,
                                    lf_vals = c(0, 1, 5, 20, 100, 300, 500, 1000)) {

  base <- sm_b0 + sm_b1 * clay + sm_b2 * et_mm + sm_b3 * ews_kmh +
          sm_b4 * cwd_mm + sm_b5 * ff

  data.frame(
    LadderFuels = lf_vals,
    Denominator = base + sm_b6 * lf_vals,
    dNBR        = pmin(ifelse(base + sm_b6 * lf_vals <= 0, 2000,
                              1 / (base + sm_b6 * lf_vals)), 2000),
    Capped      = (base + sm_b6 * lf_vals) <= 0
  ) |>
    dplyr::mutate(
      Status = dplyr::case_when(
        Capped      ~ "⚠ dNBR capped at 2000 (denominator <= 0)",
        dNBR > 1000 ~ "⚠ Very high severity",
        dNBR > 600  ~ "High severity",
        dNBR > 355  ~ "Moderate severity",
        dNBR > 185  ~ "Low-Moderate severity",
        TRUE        ~ "Low severity"
      )
    )
}


plot_dnbr_surface <- function(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                               clay, et_mm, cwd_mm, ff,
                               ews_max = 40,
                               lf_scenarios = c(0, 5, 50, 200, 500)) {

  df <- expand.grid(EWS = seq(0, ews_max, 0.5), LF = lf_scenarios) |>
    dplyr::mutate(
      eta  = sm_b0 + sm_b1 * clay + sm_b2 * et_mm +
             sm_b3 * EWS + sm_b4 * cwd_mm + sm_b5 * ff + sm_b6 * LF,
      dNBR = pmin(ifelse(eta <= 0, 2000, 1 / eta), 2000),
      LF_label = factor(paste0("LF = ", LF, " g/m²"))
    )

  # Check if any scenario is always capped
  always_capped <- df |>
    dplyr::group_by(LF_label) |>
    dplyr::summarise(all_cap = all(dNBR >= 1999), .groups = "drop")

  n_capped <- sum(always_capped$all_cap)
  sub_txt  <- sprintf(
    "Clay=%.3f  ET=%.0f mm  CWD=%.0f mm  FF=%.2f   %s",
    clay, et_mm, cwd_mm, ff,
    if (n_capped > 0)
      sprintf("⚠ %d LF scenario(s) give dNBR=2000 at all EWS (denominator <= 0)", n_capped)
    else "✓ All LF scenarios give sub-capped dNBR at low EWS"
  )

  ggplot(df, aes(EWS, dNBR, colour = LF_label)) +
    # Severity reference lines
    geom_hline(yintercept = c(185, 355, 600), linetype = "dashed",
               colour = "grey65", linewidth = 0.4) +
    annotate("text", x = ews_max * 0.98, y = c(185, 355, 600),
             label = c("Low", "Moderate", "High"),
             hjust = 1, vjust = -0.5, size = 3, colour = "grey45") +
    geom_hline(yintercept = 2000, linetype = "solid", colour = "#c0392b", linewidth = 0.6) +
    annotate("text", x = ews_max * 0.98, y = 2000,
             label = "Cap = 2000", hjust = 1, vjust = -0.5, size = 3, colour = "#c0392b") +
    geom_line(linewidth = 1.1) +
    scale_colour_brewer(palette = "YlOrBr", direction = 1) +
    scale_y_continuous(limits = c(0, 2100), labels = scales::comma) +
    labs(
      title    = "Predicted dNBR vs Effective Wind Speed",
      subtitle = sub_txt,
      x        = "Effective Wind Speed (km/h)",
      y        = "Predicted dNBR",
      colour   = "LadderFuels"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9,
            colour = if (n_capped > 0) "#c0392b" else "grey40"))
}


# =============================================================================
# Site Mortality — Eta Decomposition Bar Chart
#
# At one set of predictor values, shows the additive contribution of each
# term (Bn * predictor) to eta.  Positive contributions (blue) raise eta and
# reduce severity; negative contributions (red) lower eta and raise dNBR.
# The total eta and resulting dNBR are shown in the subtitle.
# =============================================================================

plot_sm_eta_decomposition <- function(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                       clay, et_mm, ews_kmh, cwd_mm, ff, lf) {

  terms <- data.frame(
    term  = c("B0 (Intercept)",
              "B1 x Clay",
              "B2 x ET",
              "B3 x EWS",
              "B4 x CWD",
              "B5 x FineFuels",
              "B6 x LadderFuels"),
    coef  = c(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6),
    pred  = c(1, clay, et_mm, ews_kmh, cwd_mm, ff, lf)
  ) |>
    dplyr::mutate(
      contribution = coef * pred,
      direction    = ifelse(contribution >= 0,
                            "Increases eta (reduces severity)",
                            "Decreases eta (increases severity)"),
      term         = factor(term, levels = rev(term))  # bottom-to-top display
    )

  eta_total <- sum(terms$contribution)
  dnbr_val  <- if (eta_total <= 0) Inf else 1 / eta_total
  dnbr_show <- min(dnbr_val, 2000)

  sev_label <- dplyr::case_when(
    is.infinite(dnbr_val) | dnbr_val > 1999 ~ "CAPPED at 2000 (eta <= 0)",
    dnbr_val > 1000 ~ sprintf("%.0f  — Very High severity",  dnbr_show),
    dnbr_val > 600  ~ sprintf("%.0f  — High severity",       dnbr_show),
    dnbr_val > 355  ~ sprintf("%.0f  — Moderate severity",   dnbr_show),
    dnbr_val > 185  ~ sprintf("%.0f  — Low-Moderate",        dnbr_show),
    TRUE            ~ sprintf("%.0f  — Low severity",         dnbr_show)
  )

  eta_ok <- is.finite(dnbr_val) && dnbr_val < 2000
  sub_col <- if (eta_ok) "grey40" else "#c0392b"

  ggplot(terms, aes(x = contribution, y = term, fill = direction)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2e", contribution),
                  hjust = ifelse(contribution >= 0, -0.1, 1.1)),
              size = 3.2, colour = "grey20") +
    scale_fill_manual(
      values = c("Increases eta (reduces severity)"  = "#2980b9",
                 "Decreases eta (increases severity)" = "#c0392b"),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = 0.3)) +
    labs(
      title    = "Eta Decomposition — term-by-term contributions",
      subtitle = sprintf("eta = %.5f   =>   dNBR = %s", eta_total, sev_label),
      x        = "Contribution to eta  (eta = sum of all terms = 1/dNBR)",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.subtitle   = element_text(size = 9.5, colour = sub_col, face = "bold")
    )
}


# =============================================================================
# Site Mortality — Marginal Effects Grid
#
# Holds all predictors at their current slider values and varies one at a
# time, showing the isolated dNBR response.  A red diamond marks the
# current slider position on each panel.
# =============================================================================

plot_sm_marginal_effects <- function(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                      clay, et_mm, ews_kmh, cwd_mm, ff, lf) {

  # Baseline eta at current slider values
  eta_base <- sm_b0 + sm_b1*clay + sm_b2*et_mm + sm_b3*ews_kmh +
              sm_b4*cwd_mm + sm_b5*ff + sm_b6*lf

  # For each predictor, replace that term with a range of values
  make_curve <- function(b_coef, pred_range, current_val, panel_label, x_label) {
    eta_vec <- eta_base - b_coef * current_val + b_coef * pred_range
    data.frame(
      panel   = panel_label,
      x_label = x_label,
      x       = pred_range,
      eta     = eta_vec,
      dNBR    = pmin(ifelse(eta_vec <= 0, 2000, 1 / eta_vec), 2000),
      is_current = abs(pred_range - current_val) == min(abs(pred_range - current_val))
    )
  }

  df <- dplyr::bind_rows(
    make_curve(sm_b1, seq(0,    0.7, 0.005), clay,     "Clay",         "Clay (fraction 0-1)"),
    make_curve(sm_b2, seq(0,  1500, 5),       et_mm,   "ET",           "ET (mm/yr)"),
    make_curve(sm_b3, seq(0,    50, 0.25),    ews_kmh, "EWS",          "Effective Wind Speed (km/h)"),
    make_curve(sm_b4, seq(0,   700, 2),       cwd_mm,  "CWD",          "Climatic Water Deficit (mm)"),
    make_curve(sm_b5, seq(0,     1, 0.01),    ff,      "FineFuels",    "Fine Fuels (0-1 scaled)"),
    make_curve(sm_b6, seq(0,  1000, 5),       lf,      "LadderFuels",  "Ladder Fuels (g C/m2)")
  ) |>
    dplyr::mutate(panel = factor(panel,
                                  levels = c("Clay","ET","EWS","CWD","FineFuels","LadderFuels")))

  current_pts <- df |> dplyr::filter(is_current)

  ggplot(df, aes(x, dNBR)) +
    geom_hline(yintercept = c(185, 355, 600, 1000), linetype = "dashed",
               colour = "grey72", linewidth = 0.35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 185,
             fill = "#ffffcc", alpha = 0.35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 185, ymax = 355,
             fill = "#fed976", alpha = 0.35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 355, ymax = 600,
             fill = "#fd8d3c", alpha = 0.30) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 600, ymax = 2100,
             fill = "#e31a1c", alpha = 0.12) +
    geom_line(colour = "#2c3e50", linewidth = 1.0) +
    geom_point(data = current_pts, aes(x, dNBR),
               colour = "#c0392b", size = 4, shape = 18) +
    facet_wrap(~ panel, scales = "free_x", nrow = 2,
               labeller = labeller(panel = c(
                 Clay        = "Clay  (B1)",
                 ET          = "ET  (B2)",
                 EWS         = "EWS  (B3)",
                 CWD         = "CWD  (B4)",
                 FineFuels   = "FineFuels  (B5)",
                 LadderFuels = "LadderFuels  (B6)"
               ))) +
    scale_y_continuous(limits = c(0, 2100),
                       labels = scales::comma,
                       breaks = c(0, 185, 355, 600, 1000, 2000)) +
    labs(
      title    = "Marginal Effects — isolated dNBR response per predictor",
      subtitle = "All other predictors held at current slider values.  Red diamond = current value.  Shading = severity class.",
      x        = NULL,
      y        = "Predicted dNBR"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#e8eaf6"),
      strip.text       = element_text(face = "bold", size = 9),
      plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
      axis.text.x      = element_text(size = 8)
    )
}


# =============================================================================
# Site Mortality — EWS x LadderFuels Severity Heatmap
#
# 2D color map of predicted severity class over the full EWS x LF space
# at fixed Clay, ET, CWD, FF values from sliders.  Replaces guesswork
# about which parameter combinations produce each severity class.
# =============================================================================

plot_sm_ews_lf_heatmap <- function(sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                    clay, et_mm, cwd_mm, ff,
                                    ews_max = 50, lf_max = 800) {

  sev_breaks <- c(0, 185, 355, 600, 1000, 2001)
  sev_labels <- c("Low (<185)", "Low-Mod (185-355)",
                  "Moderate (355-600)", "High (600-1000)", "Very High (>1000)")
  sev_colors <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c", "#800026")

  df <- expand.grid(EWS = seq(0, ews_max, 0.5),
                    LF  = seq(0, lf_max,  5)) |>
    dplyr::mutate(
      eta      = sm_b0 + sm_b1*clay + sm_b2*et_mm + sm_b3*EWS +
                 sm_b4*cwd_mm + sm_b5*ff + sm_b6*LF,
      dNBR     = pmin(ifelse(eta <= 0, 2000, 1/eta), 2000),
      severity = cut(dNBR, breaks = sev_breaks, labels = sev_labels,
                     right = FALSE, include.lowest = TRUE)
    )

  # Contour lines at key dNBR thresholds
  df_cont <- expand.grid(EWS = seq(0, ews_max, 0.1),
                          LF  = seq(0, lf_max,  1)) |>
    dplyr::mutate(
      eta  = sm_b0 + sm_b1*clay + sm_b2*et_mm + sm_b3*EWS + sm_b4*cwd_mm + sm_b5*ff + sm_b6*LF,
      dNBR = pmin(ifelse(eta <= 0, 2000, 1/eta), 2000)
    )

  ggplot(df, aes(EWS, LF, fill = severity)) +
    geom_raster(interpolate = FALSE) +
    geom_contour(data        = df_cont,
                 mapping     = aes(x = EWS, y = LF, z = dNBR),
                 inherit.aes = FALSE,
                 breaks      = c(185, 355, 600, 1000),
                 colour      = "white", linewidth = 0.5, linetype = "solid") +
    scale_fill_manual(
      values = setNames(sev_colors, sev_labels),
      name   = "Severity class",
      drop   = FALSE
    ) +
    labs(
      title    = "EWS x LadderFuels Severity Map",
      subtitle = sprintf("Clay=%.3f  ET=%.0f mm  CWD=%.0f mm  FF=%.2f  |  White contours at dNBR 185/355/600/1000",
                         clay, et_mm, cwd_mm, ff),
      x        = "Effective Wind Speed (km/h)",
      y        = "LadderFuels (g C/m2)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "right",
      plot.subtitle    = element_text(size = 8.5, colour = "grey40")
    )
}


# =============================================================================
# Plot: Cohort Mortality
# =============================================================================

plot_cm_curve <- function(cm_b0, cm_b1, cm_b2,
                           bt_max   = 6,
                           dnbr_vals = c(185, 355, 600, 1000)) {

  df <- expand.grid(BT = seq(0, bt_max, 0.05), dNBR = dnbr_vals) |>
    dplyr::mutate(
      Pmort    = 1 / (1 + exp(-(cm_b0 + cm_b1 * BT + cm_b2 * dNBR))),
      sev      = factor(dNBR,
                        levels = dnbr_vals,
                        labels = sprintf("dNBR=%d", dnbr_vals))
    )

  b1_ok <- cm_b1 < 0
  b2_ok <- cm_b2 > 0
  sign_note <- paste0(
    if (!b1_ok) "⚠ B1 > 0 (thicker bark should PROTECT, not kill) " else "✓ B1 negative ",
    if (!b2_ok) "| ⚠ B2 < 0 (higher dNBR should KILL, not protect) " else "| ✓ B2 positive"
  )
  warn_col <- if (!b1_ok || !b2_ok) "#c0392b" else "grey40"

  ggplot(df, aes(BT, Pmort, colour = sev)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = c("#fee08b", "#fc8d59", "#e34a33", "#b30000")) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title    = "Cohort Mortality vs Bark Thickness",
      subtitle = paste0(
        sprintf("B0=%.3f  B1=%.4f  B2=%.5f   ", cm_b0, cm_b1, cm_b2),
        sign_note
      ),
      x        = "Bark Thickness (cm)",
      y        = "P(cohort mortality)",
      colour   = "Severity",
      caption  = "Correct behaviour: P decreases with BT; P increases from Low → High severity."
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, colour = warn_col))
}


# =============================================================================
# Plots: Fire simulation output
# =============================================================================

plot_fire_timeline <- function(sim, cell_area_ha = sim$cell_area_ha) {

  tl <- sim$timeline
  if (nrow(tl) <= 1 && sim$n_burned <= 1) {
    return(ggplot() +
      labs(title = sprintf(
        "Fire went out at the ignition cell (P(spread) = %.3f; percolation needs ~0.45+)",
        sim$p_spread)) +
      theme_bw())
  }

  tl$area_ha <- tl$new_cells * cell_area_ha
  msa_line   <- max(sim$msa_ha, 0)

  ggplot(tl, aes(day, area_ha)) +
    geom_col(fill = "#fc8d59", alpha = 0.75, width = 1) +
    geom_hline(yintercept = msa_line, linetype = "dashed",
               colour = "#c0392b", linewidth = 0.7) +
    annotate("text", x = max(tl$day), y = msa_line,
             label = sprintf("MaxSpreadArea = %.0f ha/day", sim$msa_ha),
             hjust = 1.05, vjust = -0.5, colour = "#c0392b", size = 3) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = sprintf("Area Burned per Day  |  FWI=%.0f  EWS=%.0f  FF=%.2f",
                         sim$fwi, sim$ews, sim$ff),
      subtitle = sprintf(
        "Total %s ha over %d days  |  P(spread) %.3f cardinal, %.3f diagonal  |  Grid %d\u00d7%d%s",
        scales::comma(round(sim$total_ha, 1)), sim$days, sim$p_spread, sim$p_diag,
        sim$grid_dim, sim$grid_dim, if (isTRUE(sim$hit_edge)) "  |  HIT GRID EDGE" else ""),
      x        = "Fire day",
      y        = "Area burned that day (ha)",
      caption  = "MaxSpreadArea paces the fire (cells per day); it does not cap total size."
    ) +
    theme_bw(base_size = 12)
}


plot_fire_map <- function(sim) {

  grid <- sim$grid
  nr   <- nrow(grid); nc <- ncol(grid)

  df <- data.frame(
    row   = rep(seq_len(nr), times = nc),
    col   = rep(seq_len(nc), each  = nr),
    state = as.integer(as.vector(grid))
  ) |>
    dplyr::mutate(
      State = factor(state, levels = c(0L, 3L, 2L, 1L),
                     labels = c("Untouched", "Tried, did not ignite",
                                "Burned", "Burned on final day"))
    )

  ggplot(df, aes(col, row, fill = State)) +
    geom_raster() +
    scale_fill_manual(values = c("Untouched"             = "#f0f0f0",
                                  "Tried, did not ignite" = "#9ecae1",
                                  "Burned"                = "#67001f",
                                  "Burned on final day"   = "#ff6b00"),
                      drop = FALSE) +
    coord_equal() +
    labs(
      title = "Final Fire Footprint",
      x = NULL, y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom")
}


# =============================================================================
# SpreadProbability response-surface heatmap (FWI × EWS)
# =============================================================================

# -----------------------------------------------------------------------------
# plot_sp_heatmap
#
# 2-D heatmap of P(spread) across the full FWI × EWS space. White contour
# lines mark the 10 / 25 / 50 / 75 % probability levels. The shaded region
# marks the "plausible Cascades operating window" (FWI 10-40, EWS 5-35 km/h).
# More informative than line plots because it shows the joint response and
# makes it immediately visible if P is already near 1 across the whole space.
# -----------------------------------------------------------------------------

plot_sp_heatmap <- function(sp_b0, sp_b1, sp_b2, sp_b3, ff_val = 0.5) {

  df <- expand.grid(FWI = seq(0, 50, 0.5), EWS = seq(0, 50, 0.5)) |>
    dplyr::mutate(
      Pspread = 1 / (1 + exp(-(sp_b0 + sp_b1 * FWI + sp_b2 * ff_val + sp_b3 * EWS)))
    )

  # P(spread) at the center of the plausible operating window
  p_center <- 1 / (1 + exp(-(sp_b0 + sp_b1 * 25 + sp_b2 * ff_val + sp_b3 * 20)))
  sub_txt  <- sprintf(
    "B0=%.3f  B1=%.4f  B2=%.4f  B3=%.4f  FF=%.2f  |  P at FWI=25,EWS=20: %.2f",
    sp_b0, sp_b1, sp_b2, sp_b3, ff_val, p_center
  )

  ggplot(df, aes(FWI, EWS, fill = Pspread)) +
    geom_tile() +
    # Operating window outline
    annotate("rect", xmin = 10, xmax = 40, ymin = 5, ymax = 35,
             fill = NA, colour = "white", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = 10.5, y = 34.5,
             label = "Typical\noperating\nwindow",
             hjust = 0, vjust = 1, size = 3, colour = "white") +
    # Contour lines at key probabilities
    geom_contour(aes(z = Pspread), breaks = c(0.10, 0.25, 0.50, 0.75),
                 colour = "white", linewidth = 0.5, linetype = "solid") +
    geom_contour(aes(z = Pspread), breaks = 0.50,
                 colour = "white", linewidth = 1.2) +
    scale_fill_gradientn(
      name   = "P(spread)",
      colours = c("#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
      limits  = c(0, 1),
      labels  = scales::percent,
      guide   = guide_colorbar(barwidth = 12, barheight = 0.8,
                               title.position = "top", title.hjust = 0.5)
    ) +
    labs(
      title    = "Spread Probability — FWI × EWS Response Surface",
      subtitle = sub_txt,
      x        = "FWI (fire weather index)",
      y        = "Effective Wind Speed (km/h)",
      caption  = sprintf(
        "White contours at P = 10%%, 25%%, 50%% (bold), 75%%.  FineFuels fixed at %.2f.",
        ff_val
      )
    ) +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}


# =============================================================================
# Fire size distribution — simulation ensemble
# =============================================================================

# -----------------------------------------------------------------------------
# run_fire_ensemble
#
# Runs the cellular-automaton fire simulation across a grid of FWI × EWS
# conditions, with multiple random seeds per cell, to build a fire size
# distribution under the current parameter set.
#
# Returns a data.frame with columns:
#   FWI, EWS, rep, fire_ha, days, p_spread, msa_ha, hit_edge, grid_dim
# -----------------------------------------------------------------------------

run_fire_ensemble <- function(fwi_vals     = c(5, 15, 25, 35),
                               ews_vals     = c(5, 15, 25),
                               ff_val       = 0.5,
                               n_reps       = 12,
                               grid_dim     = 200L,
                               msa_b0, msa_b1, msa_b2,
                               sp_b0, sp_b1, sp_b2, sp_b3,
                               cell_area_ha = 0.09,
                               max_days     = 150L,
                               progress_fn  = NULL) {

  grid_dim  <- max(50L, as.integer(grid_dim))
  scenarios <- expand.grid(FWI = fwi_vals, EWS = ews_vals,
                            stringsAsFactors = FALSE)
  n_total   <- nrow(scenarios) * n_reps
  results   <- vector("list", n_total)
  k <- 0L

  for (i in seq_len(nrow(scenarios))) {
    fwi <- scenarios$FWI[i]
    ews <- scenarios$EWS[i]
    for (rep in seq_len(n_reps)) {
      k <- k + 1L
      if (!is.null(progress_fn)) progress_fn(k, n_total)
      sim <- simulate_fire_simple(
        grid_dim     = grid_dim,
        fwi          = fwi,
        ews          = ews,
        ff           = ff_val,
        seed         = rep * 997L + i,   # deterministic but varied seeds
        msa_b0       = msa_b0, msa_b1 = msa_b1, msa_b2 = msa_b2,
        sp_b0        = sp_b0,  sp_b1  = sp_b1,
        sp_b2        = sp_b2,  sp_b3  = sp_b3,
        cell_area_ha = cell_area_ha,
        max_days     = max_days
      )
      results[[k]] <- data.frame(
        FWI      = fwi,
        EWS      = ews,
        rep      = rep,
        fire_ha  = sim$total_ha,
        days     = sim$days,
        p_spread = sim$p_spread,
        msa_ha   = sim$msa_ha,
        hit_edge = sim$hit_edge,
        grid_dim = grid_dim
      )
    }
  }
  do.call(rbind, results)
}


# -----------------------------------------------------------------------------
# ignition_day_weights
#
# Expected SCF ignitions on each day, summed over ignition types, from the
# ignition coefficients (PlugIn.cs): lambda = exp(b0 + b1*FWI), and for ZIP
# the Poisson draw only happens when runif() >= plogis(bz0 + bz1*FWI), so
# E[n] = (1 - plogis(bz0 + bz1*FWI)) * lambda.  Used to pick ignition days in
# proportion to how often SCF would ignite fires on them.
# ign_coef: data.frame(ignition_type, distribution, b0, b1, bz0, bz1), or NULL
# for equal weights.
# -----------------------------------------------------------------------------

ignition_day_weights <- function(fwi, ign_coef = NULL) {
  if (is.null(ign_coef) || nrow(ign_coef) == 0) return(rep(1, length(fwi)))
  w <- numeric(length(fwi))
  for (i in seq_len(nrow(ign_coef))) {
    r   <- ign_coef[i, ]
    lam <- exp(r$b0 + r$b1 * fwi)
    zip <- !identical(as.character(r$distribution), "Poisson") &&
           isTRUE(is.finite(r$bz0) && is.finite(r$bz1))
    if (zip) lam <- lam * (1 - stats::plogis(r$bz0 + r$bz1 * fwi))
    w <- w + lam
  }
  w[!is.finite(w) | w < 0] <- 0
  w
}


# -----------------------------------------------------------------------------
# run_fire_ensemble_weather
#
# Fire-size distribution under REAL daily weather, the way SCF experiences it:
#   1. Ignition days are drawn from the weather record in proportion to the
#      expected number of SCF ignitions that day (ignition_day_weights), so
#      fire weather at ignition has SCF's own distribution.
#   2. Each fire then burns on that year's actual daily FWI and wind from the
#      ignition day to day `days_per_year` (SCF DaysPerYear = 364), via
#      simulate_fire_simple().  Fires go out when cooler days drop P(spread)
#      below the percolation threshold, as they do in SCF.
#
# weather: data.frame(date, FWI, EWS).  EWS in km/h.  Days missing from the
#   record get FWI = 0 (the Climate Library only computes FWI inside the fire
#   season) and the median EWS.
#
# Returns one row per fire: fire, start_date, doy, ign_FWI, ign_EWS, fire_ha,
# n_cells, days, hit_edge, p_spread (ignition day, cardinal), msa_ha.
# attr(, "mode") = "weather"; attr(, "info") = list of run settings.
# -----------------------------------------------------------------------------

run_fire_ensemble_weather <- function(weather,
                                      ign_coef      = NULL,
                                      n_fires       = 100L,
                                      grid_dim      = 200L,
                                      ff_val        = 0.5,
                                      msa_b0, msa_b1, msa_b2,
                                      sp_b0,  sp_b1,  sp_b2,  sp_b3,
                                      cell_area_ha  = 0.09,
                                      days_per_year = 364L,
                                      seed          = 1L,
                                      progress_fn   = NULL) {

  wx <- weather %>%
    transmute(date = as.Date(date), FWI = as.numeric(FWI), EWS = as.numeric(EWS)) %>%
    filter(!is.na(date)) %>%
    distinct(date, .keep_all = TRUE) %>%
    arrange(date)
  ews_fill <- stats::median(wx$EWS, na.rm = TRUE)
  if (!is.finite(ews_fill)) ews_fill <- 0
  doy <- as.integer(format(wx$date, "%j"))

  # Candidate ignition days: FWI known, and at least one day left to spread
  cand <- which(is.finite(wx$FWI) & doy < days_per_year)
  if (length(cand) == 0) stop("No usable weather days for ignition.")
  w <- ignition_day_weights(wx$FWI[cand], ign_coef)
  weighting <- if (is.null(ign_coef) || nrow(ign_coef) == 0) "equal (no ignition coefficients)"
               else "expected SCF ignitions per day"
  if (sum(w) <= 0) { w <- rep(1, length(cand)); weighting <- "equal (weights were all zero)" }

  set.seed(seed)
  picks <- cand[sample.int(length(cand), n_fires, replace = TRUE, prob = w)]
  date_key <- as.integer(wx$date)

  rows <- vector("list", n_fires)
  for (k in seq_len(n_fires)) {
    if (!is.null(progress_fn)) progress_fn(k, n_fires)
    d0  <- wx$date[picks[k]]
    yr  <- as.integer(format(d0, "%Y"))
    d1  <- as.Date(sprintf("%d-01-01", yr)) + (days_per_year - 1L)
    dts <- seq(d0, d1, by = "day")
    idx <- match(as.integer(dts), date_key)
    fwi_v <- wx$FWI[idx]; fwi_v[!is.finite(fwi_v)] <- 0
    ews_v <- wx$EWS[idx]; ews_v[!is.finite(ews_v)] <- ews_fill

    sim <- simulate_fire_simple(
      grid_dim = grid_dim, fwi = fwi_v, ews = ews_v, ff = ff_val,
      msa_b0 = msa_b0, msa_b1 = msa_b1, msa_b2 = msa_b2,
      sp_b0  = sp_b0,  sp_b1  = sp_b1,  sp_b2  = sp_b2, sp_b3 = sp_b3,
      cell_area_ha = cell_area_ha, max_days = length(dts),
      seed = seed * 100003L + k)

    rows[[k]] <- data.frame(
      fire = k, start_date = d0, doy = doy[picks[k]],
      ign_FWI = fwi_v[1], ign_EWS = ews_v[1],
      fire_ha = sim$total_ha, n_cells = sim$n_burned, days = sim$days,
      hit_edge = sim$hit_edge, p_spread = sim$p_spread, msa_ha = sim$msa_ha)
  }
  out <- do.call(rbind, rows)
  attr(out, "mode") <- "weather"
  attr(out, "info") <- list(
    weighting = weighting, n_weather_days = nrow(wx),
    years = range(as.integer(format(wx$date, "%Y"))),
    grid_dim = grid_dim, cell_area_ha = cell_area_ha,
    days_per_year = days_per_year)
  out
}


# -----------------------------------------------------------------------------
# plot_fire_ensemble_weather
#
# A: size CDF of simulated fires vs observed FPA FOD, both restricted to fires
#    larger than one cell and no larger than the simulation grid (the same
#    truncation as candidate scoring, with the grid standing in for the
#    landscape).  Subtitle gives the size score (mean |log quantile ratio|).
# B: fire size against ignition-day FWI; open points hit the grid edge.
# -----------------------------------------------------------------------------

plot_fire_ensemble_weather <- function(ens, obs_sizes = NULL) {
  info    <- attr(ens, "info")
  cell    <- info$cell_area_ha
  grid_ha <- info$grid_dim^2 * cell

  sim_sz <- ens$fire_ha[ens$fire_ha > cell]
  obs_sz <- if (!is.null(obs_sizes)) {
    x <- obs_sizes[is.finite(obs_sizes)]
    x[x > cell & x <= grid_ha]
  } else numeric(0)

  .cdf <- function(x, lab) if (length(x))
    data.frame(size_ha = sort(x), cdf = seq_along(x) / length(x), source = lab)
  cdf_df <- rbind(.cdf(sim_sz, "Simulated (SCF rules, real weather)"),
                  .cdf(obs_sz, "Observed (FPA FOD)"))

  score_txt <- if (length(sim_sz) && length(obs_sz))
    sprintf("  |  size score %.2f (0 = match, 0.69 = 2x off)",
            score_size_quantiles(obs_sz, sim_sz)) else ""
  pct_out <- 100 * mean(ens$n_cells <= 1)

  pA <- if (is.null(cdf_df)) {
    ggplot() + labs(title = "No fire spread beyond its ignition cell") + theme_bw()
  } else {
    ggplot(cdf_df, aes(size_ha, cdf, colour = source)) +
      geom_step(linewidth = 1) +
      scale_x_log10(labels = scales::comma) +
      scale_colour_manual(values = c("Simulated (SCF rules, real weather)" = "#d94801",
                                     "Observed (FPA FOD)" = "#2166ac")) +
      labs(title = "Fire size distribution",
           subtitle = sprintf("%.0f%% of ignitions never left the ignition cell%s",
                              pct_out, score_txt),
           x = "Fire size (ha, log scale)", y = "Cumulative probability",
           colour = NULL,
           caption = sprintf("Fires > 1 cell (%.2f ha) and <= grid area (%s ha)",
                             cell, scales::comma(round(grid_ha)))) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  }

  pB <- ggplot(ens, aes(ign_FWI, pmax(fire_ha, cell), shape = hit_edge)) +
    geom_point(colour = "#d94801", alpha = 0.6, size = 2) +
    scale_y_log10(labels = scales::comma) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c(`FALSE` = "Contained in grid",
                                  `TRUE` = "Reached grid edge")) +
    labs(title = "Size vs ignition-day FWI",
         subtitle = sprintf("%d fires  |  ignition days weighted by %s",
                            nrow(ens), info$weighting),
         x = "FWI on ignition day", y = "Fire size (ha, log scale)", shape = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(pA, pB, ncol = 2)
  else pA
}


# -----------------------------------------------------------------------------
# sweep_b0_weather
#
# Sets SpreadProbabilityB0 the way CLAUDE.md says it has to be set: by
# simulation, because the regression intercept is not identifiable (single-
# trial percolation + case-control sampling).  For each candidate B0 the
# real-weather ensemble (run_fire_ensemble_weather) is run with the SAME seed,
# so every B0 sees the same ignition days and the same random draws (common
# random numbers): score differences come from B0, not from noise.
#
# Scores (same scale as candidate scoring: mean |log ratio|, 0 = perfect,
# 0.69 = off by a factor of two):
#   score_size  fire-size distribution vs observed, both restricted to
#               fires > 1 cell and <= the simulation grid area.
#   score_aab   expected annual area burned as a fraction of the landscape,
#               (SCF ignitions per year on the landscape) x (mean area per
#               ignition), vs the observed regional fraction per year
#               (obs_rate).  Ignitions per year use the ignition coefficients
#               with the log(landscape / boundary) area offset applied
#               (area_ratio = tmpl_area_ha / cal_area_ha), summed over days with
#               known FWI and divided by the years in the record.
#   score_total weighted mean of the available scores.
#
# Coarse pass over b0_values, stopping after `patience` consecutive values
# worse than the best so far; then (refine = TRUE) a fine pass at step/4 within
# one coarse step of the coarse best.  progress_fn(fraction, detail) is called
# per fire within each B0 run.
#
# Returns list(table, ensembles (named by B0), best_b0, info).
# -----------------------------------------------------------------------------

sweep_b0_weather <- function(b0_values, weather,
                             obs_sizes,
                             obs_rate      = NULL,
                             ign_coef      = NULL,
                             area_ratio    = NULL,
                             tmpl_area_ha  = NULL,
                             n_fires       = 150L,
                             grid_dim      = 200L,
                             ff_val        = 0.5,
                             msa_b0, msa_b1, msa_b2,
                             sp_b1, sp_b2, sp_b3,
                             cell_area_ha  = 0.09,
                             w_size        = 0.5,
                             w_aab         = 0.5,
                             patience      = 3L,
                             refine        = TRUE,
                             seed          = 1L,
                             progress_fn   = NULL) {

  b0_values <- sort(unique(b0_values))
  grid_ha   <- grid_dim^2 * cell_area_ha
  obs_sz    <- obs_sizes[is.finite(obs_sizes) & obs_sizes > cell_area_ha &
                         obs_sizes <= grid_ha]

  # Expected SCF ignitions per year on the landscape
  ign_per_yr <- NA_real_
  if (!is.null(ign_coef) && nrow(ign_coef) > 0 && isTRUE(is.finite(area_ratio))) {
    ign_land <- ign_coef
    ign_land$b0 <- ign_land$b0 + log(area_ratio)
    fwi_ok <- weather$FWI[is.finite(weather$FWI)]
    n_yrs  <- length(unique(format(as.Date(weather$date[is.finite(weather$FWI)]), "%Y")))
    ign_per_yr <- sum(ignition_day_weights(fwi_ok, ign_land)) / n_yrs
  }

  run_one <- function(b0, label) {
    prog <- if (is.null(progress_fn)) NULL else
      function(k, n) progress_fn(k / n, sprintf("%s: B0 = %.3f, fire %d of %d", label, b0, k, n))
    ens <- run_fire_ensemble_weather(
      weather = weather, ign_coef = ign_coef, n_fires = n_fires,
      grid_dim = grid_dim, ff_val = ff_val,
      msa_b0 = msa_b0, msa_b1 = msa_b1, msa_b2 = msa_b2,
      sp_b0 = b0, sp_b1 = sp_b1, sp_b2 = sp_b2, sp_b3 = sp_b3,
      cell_area_ha = cell_area_ha, seed = seed, progress_fn = prog)
    ens_list[[sprintf("%.4f", b0)]] <<- ens

    sim_sz   <- ens$fire_ha[ens$fire_ha > cell_area_ha]
    s_size   <- score_size_quantiles(obs_sz, sim_sz)
    sim_rate <- if (is.finite(ign_per_yr) && isTRUE(tmpl_area_ha > 0))
                  ign_per_yr * mean(ens$fire_ha) / tmpl_area_ha
                else NA_real_
    s_aab    <- if (!is.null(obs_rate) && is.finite(sim_rate)) score_annual_area(obs_rate, sim_rate)
                else NA_real_
    s <- c(s_size, s_aab); w <- c(w_size, w_aab); keep <- !is.na(s) & w > 0
    data.frame(
      B0 = b0, pass = label,
      pct_no_spread   = 100 * mean(ens$n_cells <= 1),
      median_ha       = if (length(sim_sz)) stats::median(sim_sz) else NA_real_,
      p90_ha          = if (length(sim_sz)) stats::quantile(sim_sz, 0.9, names = FALSE) else NA_real_,
      mean_ha_per_ign = mean(ens$fire_ha),
      median_days     = if (length(sim_sz)) stats::median(ens$days[ens$fire_ha > cell_area_ha]) else NA_real_,
      pct_edge        = 100 * mean(ens$hit_edge),
      sim_annual_pct  = 100 * sim_rate,
      score_size      = s_size,
      score_aab       = s_aab,
      score_total     = if (any(keep)) sum(w[keep] * s[keep]) / sum(w[keep]) else NA_real_)
  }

  # Coarse pass.  Stop after `patience` consecutive values worse than the best
  # so far: beyond the percolation cliff fires only get bigger (or saturate at
  # MaxSpreadArea pace x days left), so the score cannot recover.
  rows <- list(); ens_list <- list(); stopped_at <- NA_real_
  best_sc <- Inf; worse <- 0L
  for (j in seq_along(b0_values)) {
    r <- run_one(b0_values[j], sprintf("Coarse %d/%d", j, length(b0_values)))
    rows[[length(rows) + 1]] <- r
    if (is.finite(r$score_total) && r$score_total < best_sc) {
      best_sc <- r$score_total; worse <- 0L
    } else if (is.finite(best_sc)) {
      worse <- worse + 1L
      if (worse >= patience && j < length(b0_values)) { stopped_at <- b0_values[j]; break }
    }
  }

  # Fine pass around the coarse best (the cliff is steep, so a coarse grid
  # can step right over the best value).
  tab0 <- do.call(rbind, rows)
  if (refine && length(b0_values) >= 3 && any(is.finite(tab0$score_total))) {
    step  <- stats::median(diff(b0_values))
    b_c   <- tab0$B0[which.min(tab0$score_total)]
    fine  <- seq(b_c - step, b_c + step, by = step / 4)
    fine  <- fine[!sapply(fine, function(x) any(abs(x - tab0$B0) < 1e-9))]
    for (j in seq_along(fine))
      rows[[length(rows) + 1]] <- run_one(fine[j], sprintf("Fine %d/%d", j, length(fine)))
  }

  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$B0), ]
  rownames(tab) <- NULL
  best <- if (any(is.finite(tab$score_total))) tab$B0[which.min(tab$score_total)] else NA_real_

  list(table = tab, ensembles = ens_list, best_b0 = best,
       info = list(ign_per_yr = ign_per_yr, obs_rate = obs_rate,
                   obs_median_ha = if (length(obs_sz)) stats::median(obs_sz) else NA_real_,
                   obs_p90_ha = if (length(obs_sz)) stats::quantile(obs_sz, 0.9, names = FALSE) else NA_real_,
                   n_obs = length(obs_sz), grid_ha = grid_ha, n_fires = n_fires,
                   stopped_at = stopped_at, w_size = w_size, w_aab = w_aab,
                   best_b0_size = if (any(is.finite(tab$score_size)))
                     tab$B0[which.min(tab$score_size)] else NA_real_,
                   best_b0_aab  = if (any(is.finite(tab$score_aab)))
                     tab$B0[which.min(tab$score_aab)] else NA_real_))
}


# -----------------------------------------------------------------------------
# plot_b0_sweep
#
# A: the three scores against B0, best B0 marked.
# B: simulated median / P90 fire size (fires > 1 cell) against B0, with the
#    observed values as horizontal references, and % of ignitions that never
#    left their cell on a secondary scale (as labels).
# -----------------------------------------------------------------------------

plot_b0_sweep <- function(sweep) {
  tab  <- sweep$table
  info <- sweep$info
  sc <- rbind(
    data.frame(B0 = tab$B0, score = tab$score_size,  what = "Fire size"),
    data.frame(B0 = tab$B0, score = tab$score_aab,   what = "Annual area"),
    data.frame(B0 = tab$B0, score = tab$score_total, what = "Total"))
  sc <- sc[is.finite(sc$score), ]

  pA <- ggplot(sc, aes(B0, score, colour = what, linewidth = what == "Total")) +
    geom_line() + geom_point(size = 1.8) +
    scale_linewidth_manual(values = c(`TRUE` = 1.4, `FALSE` = 0.7), guide = "none") +
    scale_colour_manual(values = c("Fire size" = "#d94801", "Annual area" = "#2166ac",
                                   "Total" = "black")) +
    { if (is.finite(sweep$best_b0))
        geom_vline(xintercept = sweep$best_b0, linetype = "dashed") } +
    labs(title = sprintf("B0 sweep: best B0 = %.2f", sweep$best_b0),
         subtitle = sprintf(
           "Size alone favours %.2f, annual area alone %.2f%s",
           info$best_b0_size, info$best_b0_aab,
           if (is.finite(info$stopped_at))
             sprintf("  |  stopped after %.2f (no further improvement)", info$stopped_at) else ""),
         x = "SpreadProbabilityB0", y = "Score (mean |log ratio|)", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  sz <- rbind(data.frame(B0 = tab$B0, ha = tab$median_ha, what = "Median"),
              data.frame(B0 = tab$B0, ha = tab$p90_ha,    what = "P90"))
  sz <- sz[is.finite(sz$ha) & sz$ha > 0, ]
  ref <- data.frame(what = c("Median", "P90"),
                    ha = c(info$obs_median_ha, info$obs_p90_ha))
  pB <- ggplot(sz, aes(B0, ha, colour = what)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    geom_hline(data = ref[is.finite(ref$ha), ], aes(yintercept = ha, colour = what),
               linetype = "dashed") +
    geom_text(data = tab[!grepl("^Fine", tab$pass), ],
              aes(B0, Inf, label = sprintf("%.0f%%", pct_no_spread)),
              inherit.aes = FALSE, vjust = 1.4, size = 2.8, colour = "grey40") +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(values = c(Median = "#d94801", P90 = "#7f2704")) +
    labs(title = "Simulated fire size vs B0",
         subtitle = "Dashed = observed FPA FOD.  Top labels: % of ignitions that never spread",
         x = "SpreadProbabilityB0", y = "Fire size, fires > 1 cell (ha, log)", colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pA, pB, ncol = 2)
  else pA
}


# -----------------------------------------------------------------------------
# plot_fire_ensemble
#
# Violin + box + jitter plot of fire sizes from run_fire_ensemble(), faceted
# by EWS level. Plausible size reference lines are shown in green.
# Optional reference_p50 / reference_p90 draw observed FPA FOD quantiles.
# -----------------------------------------------------------------------------

plot_fire_ensemble <- function(ensemble_df,
                                reference_p50 = NULL,
                                reference_p90 = NULL,
                                obs_fires_df  = NULL,
                                min_obs_ha    = 4) {

  if (is.null(ensemble_df) || nrow(ensemble_df) == 0)
    return(ggplot() + labs(title = "Click 'Run Ensemble' to generate fire size distribution.") + theme_bw())

  fwi_vals <- sort(unique(ensemble_df$FWI))
  ews_vals <- sort(unique(ensemble_df$EWS))
  fwi_lvls <- paste0("FWI = ", fwi_vals)
  ews_lvls <- paste0("EWS ", ews_vals, " km/h")

  # Grid size (stored per row; grab from first row)
  gd <- if ("grid_dim" %in% names(ensemble_df)) ensemble_df$grid_dim[1] else 100L

  sim_df <- ensemble_df |>
    dplyr::mutate(
      FWI_label    = factor(paste0("FWI = ", FWI),        levels = fwi_lvls),
      EWS_label    = factor(paste0("EWS ", EWS, " km/h"), levels = ews_lvls),
      fire_ha_plot = pmax(fire_ha, 0.01),
      source       = "Simulated"
    )

  med_ha <- round(stats::median(sim_df$fire_ha), 1)
  p90_ha <- round(stats::quantile(sim_df$fire_ha, 0.90), 1)

  # ---- Prepare observed (FPA FOD) — filter to landscape-scale fires ----------
  # obs_fires_df columns: fire_ha, FWI_bin, EWS_bin (numeric, binned to ensemble values)
  # Filter to >= min_obs_ha to exclude tiny escaped ignitions (86% of FOD < 4 ha)
  obs_plot <- NULL
  if (!is.null(obs_fires_df) && nrow(obs_fires_df) > 0 &&
      all(c("fire_ha", "FWI_bin", "EWS_bin") %in% names(obs_fires_df))) {
    obs_plot <- obs_fires_df |>
      dplyr::filter(is.finite(fire_ha), fire_ha >= min_obs_ha) |>
      dplyr::mutate(
        fire_ha_plot = pmax(fire_ha, 0.01),
        FWI_label    = factor(paste0("FWI = ",  FWI_bin),       levels = fwi_lvls),
        EWS_label    = factor(paste0("EWS ", EWS_bin, " km/h"), levels = ews_lvls),
        source       = "Observed (FPA FOD)"
      ) |>
      dplyr::filter(!is.na(FWI_label), !is.na(EWS_label))
    if (nrow(obs_plot) == 0) obs_plot <- NULL
  }

  n_obs <- if (!is.null(obs_plot)) nrow(obs_plot) else 0L

  # Shared layers used in both modes
  y_scale <- scale_y_log10(
    labels = scales::comma,
    breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000),
    limits = c(0.01, 200000)
  )
  ref_band <- list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 5000,
             fill = "#2c5f2e", alpha = 0.07),
    annotate("text", x = 0.6, y = 50, label = "Plausible\n10-5000 ha",
             size = 2.8, colour = "#2c5f2e", hjust = 0)
  )
  ref_lines <- list(
    if (!is.null(reference_p50))
      geom_hline(yintercept = reference_p50, linetype = "dashed",
                 colour = "#2171b5", linewidth = 0.7) else NULL,
    if (!is.null(reference_p90))
      geom_hline(yintercept = reference_p90, linetype = "dotted",
                 colour = "#2171b5", linewidth = 0.7) else NULL
  )
  base_theme <- list(
    theme_bw(base_size = 12),
    theme(
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 9),
      strip.background = element_rect(fill = "#f0f0f0"),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )
  )

  # ============================================================
  # MODE A: observed data present — dodged boxplot comparison
  # ============================================================
  if (!is.null(obs_plot) && nrow(obs_plot) >= 3) {

    # Colours: navy = observed, orange = simulated
    col_obs <- "#2c4f8f"
    col_sim <- "#e07b54"

    all_df <- dplyr::bind_rows(
      sim_df |> dplyr::select(FWI_label, EWS_label, fire_ha_plot, source),
      obs_plot |> dplyr::select(FWI_label, EWS_label, fire_ha_plot, source)
    ) |>
      dplyr::mutate(source = factor(source,
                                     levels = c("Observed (FPA FOD)", "Simulated")))

    obs_med <- round(stats::median(obs_plot$fire_ha), 1)
    obs_p90 <- round(stats::quantile(obs_plot$fire_ha, 0.90), 1)

    p <- ggplot(all_df, aes(x = FWI_label, y = fire_ha_plot, fill = source)) +
      geom_boxplot(
        position      = position_dodge(0.7),
        outlier.shape = 21, outlier.size = 1.0, outlier.alpha = 0.45,
        alpha = 0.75, width = 0.55
      ) +
      facet_wrap(~ EWS_label, nrow = 1) +
      y_scale +
      scale_fill_manual(
        values = c("Observed (FPA FOD)" = col_obs, "Simulated" = col_sim),
        name   = NULL
      ) +
      ref_band +
      ref_lines +
      labs(
        title   = "Fire Size Distribution — Simulation Ensemble vs Observed (FPA FOD)",
        subtitle = sprintf(
          paste0("%dx%d grid, %d runs/scenario.  ",
                 "Obs (>=%.0f ha): median=%.0f ha, P90=%.0f ha  |  ",
                 "Sim: median=%.1f ha, P90=%.1f ha  |  n_obs=%d"),
          gd, gd, dplyr::n_distinct(sim_df$rep),
          min_obs_ha, obs_med, obs_p90,
          med_ha, p90_ha, n_obs
        ),
        x       = NULL,
        y       = "Fire size (ha, log scale)",
        caption = paste0(
          "FPA FOD fires binned to nearest FWI/EWS scenario (wind speed on discovery date as EWS proxy). ",
          sprintf("Min size: %.0f ha.  ", min_obs_ha),
          if (!is.null(reference_p50)) sprintf("Blue dashed = FPA P50 (%.0f ha, >=%.0f ha fires). ", reference_p50, min_obs_ha) else "",
          "Green band = plausible range."
        )
      ) +
      base_theme +
      theme(legend.position = "bottom")

  # ============================================================
  # MODE B: no observed data — violin + boxplot + jitter (simulated only)
  # ============================================================
  } else {

    p <- ggplot(sim_df, aes(FWI_label, fire_ha_plot, fill = FWI_label)) +
      geom_violin(scale = "width", alpha = 0.55, colour = NA) +
      geom_boxplot(width = 0.18, outlier.shape = NA,
                   colour = "grey20", fill = "white", alpha = 0.8) +
      geom_jitter(width = 0.12, alpha = 0.4, size = 1.4, colour = "grey50") +
      facet_wrap(~ EWS_label, nrow = 1) +
      y_scale +
      scale_fill_brewer(palette = "YlOrRd", guide = "none") +
      ref_band +
      ref_lines +
      labs(
        title    = "Fire Size Distribution — Simulation Ensemble",
        subtitle = sprintf(
          "Each scenario: %d runs, %dx%d cell grid.  Sim median = %.1f ha  |  Sim P90 = %.1f ha",
          dplyr::n_distinct(sim_df$rep), gd, gd, med_ha, p90_ha
        ),
        x       = NULL,
        y       = "Fire size (ha, log scale)",
        caption = paste0(
          "Load FPA FOD ignition data to overlay observed fire size distribution. ",
          if (!is.null(reference_p50)) sprintf("Blue dashed = FPA P50 (%.0f ha). ", reference_p50) else "",
          "Green band = plausible range."
        )
      ) +
      base_theme +
      theme(legend.position = "none")
  }

  p
}


# =============================================================================
# Ignition tuning
#
# SCF draws a daily ignition count per type from
#     lambda = exp(b0 + b1 * FWI_landscape_mean)
# and, under the ZIP option, skips the draw entirely with probability
#     pi = plogis(bz0 + bz1 * FWI)        [PlugIn.cs:816-818]
# so the expected daily count is  E[n] = (1 - pi) * lambda.
#
# The tuning tab reuses the Ignition tab's own diagnostic plots rather than
# defining new ones, so the only thing needed here is a way to feed hand-typed
# coefficients into those functions in the shape they already expect.
# =============================================================================

# -----------------------------------------------------------------------------
# ign_tune_coef_df
# Assembles a coefficient tibble in the same shape fit_ignition_models()
# returns, from the flat per-type values held in the tuning inputs.
# -----------------------------------------------------------------------------

ign_tune_coef_df <- function(l_b0, l_b1, l_bz0, l_bz1,
                              a_b0, a_b1, a_bz0, a_bz1,
                              distribution = c("ZeroInflatedPoisson", "Poisson")) {
  distribution <- match.arg(distribution)
  dplyr::tibble(
    ignition_type = c("Lightning", "Accidental"),
    distribution  = distribution,
    b0            = c(l_b0,  a_b0),
    b1            = c(l_b1,  a_b1),
    bz0           = c(l_bz0, a_bz0),
    bz1           = c(l_bz1, a_bz1)
  )
}


# =============================================================================
# SCF snippet builder
# =============================================================================

build_scf_tuning_snippet <- function(msa_b0, msa_b1, msa_b2,
                                      sp_b0, sp_b1, sp_b2, sp_b3,
                                      sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                      cm_b0, cm_b1, cm_b2,
                                      max_ff = 500,
                                      ign_coef = NULL) {

  # ---- Ignition block (omitted entirely when no ignition coefficients set) ----
  ign_lines <- character(0)
  if (!is.null(ign_coef) && nrow(ign_coef) > 0) {
    zip <- !identical(unique(ign_coef$distribution)[1], "Poisson")
    ign_lines <- c(
      ">> ---- Ignitions (lambda = exp(b0 + b1*FWI); landscape-mean FWI) -----",
      unlist(lapply(seq_len(nrow(ign_coef)), function(i) {
        r  <- ign_coef[i, ]
        kw <- if (identical(r$ignition_type, "Lightning")) "Lightning" else "Accidental"
        c(
          sprintf("%-30s %12.6f", paste0(kw, "IgnitionsB0"), r$b0),
          sprintf("%-30s %12.6f", paste0(kw, "IgnitionsB1"), r$b1),
          if (zip && is.finite(r$bz0) && is.finite(r$bz1)) c(
            sprintf("%-30s %12.6f", paste0(kw, "IgnitionsBinomialB0"), r$bz0),
            sprintf("%-30s %12.6f", paste0(kw, "IgnitionsBinomialB1"), r$bz1)
          )
        )
      })),
      if (zip)
        ">> Binomial keywords are the ZIP structural-zero terms: SCF draws from"
      else
        ">> Poisson only — no Binomial keywords written.",
      if (zip)
        ">> the Poisson only when runif() >= plogis(bz0 + bz1*FWI).",
      "",
      ">> NOTE: these are the REGIONAL (calibration boundary) coefficients as",
      ">> tuned.  Before use in SCF, b0 must carry the log(landscape / boundary)",
      ">> area offset from the Scale Adjustment tab.  b1 transfers unchanged.",
      ""
    )
  }

  paste(c(
    ">> =====================================================================",
    ">> SCF Parameter Tuning Export",
    sprintf(">> Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    ">> =====================================================================",
    "",
    ign_lines,
    ">> ---- Max Spread Area (B0 + B1*FWI + B2*EWS, result in hectares) ----",
    sprintf("MaximumSpreadAreaB0    %10.4f",  msa_b0),
    sprintf("MaximumSpreadAreaB1    %10.4f",  msa_b1),
    sprintf("MaximumSpreadAreaB2    %10.4f",  msa_b2),
    "",
    ">> ---- Spread Probability (logistic; EWS and FF in km/h and 0-1) -----",
    sprintf("SpreadProbabilityB0    %12.6f",  sp_b0),
    sprintf("SpreadProbabilityB1    %12.6f",  sp_b1),
    sprintf("SpreadProbabilityB2    %12.6f",  sp_b2),
    sprintf("SpreadProbabilityB3    %12.6f",  sp_b3),
    "",
    sprintf("MaximumFineFuels       %10.1f   >> g C m^-2; divides LANDIS fine fuels to 0-1", max_ff),
    "",
    ">> ---- Site Mortality (Gamma GLM inverse link; dNBR = 1/eta) ---------",
    sprintf("SiteMortalityB0        %14.8f",  sm_b0),
    sprintf("SiteMortalityB1        %14.8f",  sm_b1),
    sprintf("SiteMortalityB2        %14.8f",  sm_b2),
    sprintf("SiteMortalityB3        %14.8f",  sm_b3),
    sprintf("SiteMortalityB4        %14.8f",  sm_b4),
    sprintf("SiteMortalityB5        %14.8f",  sm_b5),
    sprintf("SiteMortalityB6        %14.8f",  sm_b6),
    "",
    ">> ---- Cohort Mortality (logistic; BarkThickness in cm) --------------",
    sprintf("CohortMortalityB0      %10.6f",  cm_b0),
    sprintf("CohortMortalityB1      %10.6f",  cm_b1),
    sprintf("CohortMortalityB2      %10.6f",  cm_b2),
    ""
  ), collapse = "\n")
}
