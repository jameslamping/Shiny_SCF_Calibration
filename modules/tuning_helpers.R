# =============================================================================
# tuning_helpers.R
# Interactive parameter tuning plots for the SCF calibration app.
#
# All functions take raw numeric coefficients (no reactive objects) and return
# ggplot objects or plain R values, making them easy to test independently.
#
# Functions:
#   simulate_fire_simple()   Vectorized cellular automaton (no terra deps)
#   plot_msa_surface()       MaxSpreadArea vs FWI at multiple EWS levels
#   plot_sp_surface()        SpreadProbability vs FWI
#   plot_sp_ff_sensitivity() P vs FWI at FF=0 vs FF=1 (shows mismatch risk)
#   plot_ff_scaling()        FF_scaled vs raw g C/m² and resulting P(spread)
#                            — shows the saturation problem when MaxFF is too low
#   plot_dnbr_surface()      dNBR vs EWS at multiple LadderFuel levels
#   dnbr_sensitivity_table() Denominator / dNBR tabulation across LF values
#   plot_cm_curve()          P(cohort mort) vs BarkThickness at multiple dNBR
#   plot_fire_map()          Fire grid from simulate_fire_simple() as ggplot
#   plot_fire_timeline()     Cumulative area burned per day
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
# Vectorized cellular automaton that mimics SCF's daily fire spread loop.
#
# At each day:
#   1. Find cells adjacent (8-connectivity) to the current burn front.
#   2. Each candidate cell ignites with probability P(spread).
#   3. If total new cells > MaxSpreadArea/cell_area_ha, randomly subsample.
#   4. Burning cells advance to burned; new ignitions become the new front.
#
# Returns a list with:
#   total_ha   : total area burned (ha)
#   days       : number of simulation days
#   msa_ha     : MaxSpreadArea ceiling used
#   p_spread   : spread probability used
#   timeline   : data.frame (day, new_cells, total_cells)
#   grid       : integer matrix — 0=unburned, 1=still burning (last front), 2=burned
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

  set.seed(seed)

  msa_ha    <- max(0, msa_b0 + msa_b1 * fwi + msa_b2 * ews)
  msa_cells <- max(1L, ceiling(msa_ha / cell_area_ha))
  p_spread  <- 1 / (1 + exp(-(sp_b0 + sp_b1 * fwi + sp_b2 * ff + sp_b3 * ews)))
  p_spread  <- max(0, min(1, p_spread))   # clamp to [0,1]

  grid <- matrix(0L, grid_dim, grid_dim)
  cr <- grid_dim %/% 2
  cc <- grid_dim %/% 2
  grid[cr, cc] <- 1L

  n_total  <- 1L
  timeline <- data.frame(day = integer(max_days),
                          new_cells  = integer(max_days),
                          total_cells = integer(max_days))
  n_days <- 0L

  # Vectorized 8-neighbour dilation: shift burning mask in all 8 directions
  dilate8 <- function(burning) {
    nr <- nrow(burning); nc <- ncol(burning)
    Reduce(`|`, list(
      # cardinal
      rbind(burning[-1, ],  matrix(FALSE, 1, nc)),  # shift up
      rbind(matrix(FALSE, 1, nc), burning[-nr, ]),  # shift down
      cbind(burning[, -1],  matrix(FALSE, nr, 1)),  # shift left
      cbind(matrix(FALSE, nr, 1), burning[, -nc]),  # shift right
      # diagonal
      rbind(cbind(burning[-1, -1],  matrix(FALSE, nr-1, 1)), matrix(FALSE, 1, nc)),
      rbind(cbind(matrix(FALSE, nr-1, 1), burning[-1, -nc]), matrix(FALSE, 1, nc)),
      rbind(matrix(FALSE, 1, nc), cbind(burning[-nr, -1],  matrix(FALSE, nr-1, 1))),
      rbind(matrix(FALSE, 1, nc), cbind(matrix(FALSE, nr-1, 1), burning[-nr, -nc]))
    ))
  }

  for (day in seq_len(max_days)) {
    burning <- grid == 1L

    if (!any(burning)) break

    # Candidates: adjacent to burning, currently unburned
    cand <- dilate8(burning) & (grid == 0L)
    n_cand <- sum(cand)
    if (n_cand == 0L) break

    # Spread probabilistically (fully vectorized)
    spread_mask <- cand & (matrix(runif(grid_dim^2), grid_dim, grid_dim) < p_spread)

    # Enforce MaxSpreadArea ceiling
    n_new <- sum(spread_mask)
    if (n_new > msa_cells) {
      new_idx  <- which(spread_mask)
      kept_idx <- sample(new_idx, msa_cells)
      spread_mask[] <- FALSE
      spread_mask[kept_idx] <- TRUE
      n_new <- msa_cells
    }

    # Advance state
    grid[burning]     <- 2L   # burning -> burned
    grid[spread_mask] <- 1L   # new front

    n_total <- n_total + n_new
    n_days  <- n_days + 1L
    timeline$day[n_days]         <- day
    timeline$new_cells[n_days]   <- n_new
    timeline$total_cells[n_days] <- n_total

    if (n_new == 0L || n_total >= grid_dim^2) break
  }

  timeline <- timeline[seq_len(n_days), , drop = FALSE]

  list(
    total_ha  = n_total * cell_area_ha,
    days      = n_days,
    msa_ha    = msa_ha,
    p_spread  = p_spread,
    timeline  = timeline,
    grid      = grid,
    fwi       = fwi,
    ews       = ews,
    ff        = ff
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
    geom_raster(interpolate = TRUE) +
    geom_contour(data = df_cont, aes(EWS, LF, z = dNBR),
                 breaks = c(185, 355, 600, 1000),
                 colour = "white", linewidth = 0.5, linetype = "solid") +
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

plot_fire_timeline <- function(sim, cell_area_ha) {

  tl <- sim$timeline
  if (nrow(tl) == 0) {
    return(ggplot() +
      labs(title = "Fire went out immediately (P(spread) ≈ 0 or MSA = 0)") +
      theme_bw())
  }

  tl$area_ha <- tl$total_cells * cell_area_ha

  ggplot(tl, aes(day, area_ha)) +
    geom_area(fill = "#fc8d59", alpha = 0.55) +
    geom_line(colour = "#d94801", linewidth = 0.8) +
    geom_hline(yintercept = sim$msa_ha, linetype = "dashed",
               colour = "#c0392b", linewidth = 0.7) +
    annotate("text", x = max(tl$day), y = sim$msa_ha,
             label = sprintf("MaxSpreadArea = %.0f ha", sim$msa_ha),
             hjust = 1.05, vjust = -0.5, colour = "#c0392b", size = 3) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = sprintf("Simulated Fire Growth  |  FWI=%.0f  EWS=%.0f  FF=%.2f",
                         sim$fwi, sim$ews, sim$ff),
      subtitle = sprintf("Total: %.1f ha over %d days  |  P(spread)=%.3f  |  Grid: 100×100 cells",
                         sim$total_ha, sim$days, sim$p_spread),
      x        = "Simulation day",
      y        = "Cumulative area burned (ha)"
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
      State = factor(state, levels = 0:2,
                     labels = c("Unburned", "Final front", "Burned"))
    )

  ggplot(df, aes(col, row, fill = State)) +
    geom_raster() +
    scale_fill_manual(values = c("Unburned"    = "#f0f0f0",
                                  "Final front" = "#ff6b00",
                                  "Burned"      = "#67001f")) +
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
#   FWI, EWS, rep, fire_ha, days, p_spread, msa_ha
# -----------------------------------------------------------------------------

run_fire_ensemble <- function(fwi_vals     = c(5, 15, 25, 35),
                               ews_vals     = c(5, 15, 25),
                               ff_val       = 0.5,
                               n_reps       = 12,
                               grid_dim     = 200L,
                               msa_b0, msa_b1, msa_b2,
                               sp_b0, sp_b1, sp_b2, sp_b3,
                               cell_area_ha = 0.09,
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
        max_days     = 90
      )
      results[[k]] <- data.frame(
        FWI      = fwi,
        EWS      = ews,
        rep      = rep,
        fire_ha  = sim$total_ha,
        days     = sim$days,
        p_spread = sim$p_spread,
        msa_ha   = sim$msa_ha,
        grid_dim = grid_dim
      )
    }
  }
  do.call(rbind, results)
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
# SCF snippet builder
# =============================================================================

build_scf_tuning_snippet <- function(msa_b0, msa_b1, msa_b2,
                                      sp_b0, sp_b1, sp_b2, sp_b3,
                                      sm_b0, sm_b1, sm_b2, sm_b3, sm_b4, sm_b5, sm_b6,
                                      cm_b0, cm_b1, cm_b2,
                                      max_ff = 500) {
  paste(c(
    ">> =====================================================================",
    ">> SCF Parameter Tuning Export",
    sprintf(">> Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    ">> =====================================================================",
    "",
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
