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
