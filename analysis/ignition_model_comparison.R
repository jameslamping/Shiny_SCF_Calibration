# =============================================================================
# ignition_model_comparison.R
# Standalone script to compare Poisson vs ZIP ignition model fits using
# the same data pipeline as the Shiny app.
#
# Outputs (written to analysis/plots/):
#   01_zip_sign_bug.png        — shows the bug vs corrected ZIP prediction
#   02_poisson_vs_zip.png      — side-by-side Poisson vs corrected ZIP
#   03_zero_inflation_curve.png — how pi (structural zero prob) varies with FWI
#   04_count_components.png    — lambda (count) component comparison
#   05_residuals.png           — deviance residuals, both models, both types
#   06_model_summary.txt       — coefficient table + AIC comparison
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(pscl)
  library(purrr)
  library(scales)
})

source("modules/ignition_helpers.R")
source("modules/era_helpers.R")

dir.create("analysis/plots", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA (same paths as app sidebar defaults)
# =============================================================================

cat("--- Loading calibration boundary ---\n")
shp_path <- "/Users/jlamping/Documents/Anthropic_Projects/scf_calibration_app/data_raw/landscape/S_USA.ECOSYS_ECOMAPPROVINCES_2025/ECOSYS_CoastalRange.shp"
cal_vect_raw <- vect(shp_path)

cat("--- Loading LANDIS template raster ---\n")
tif_path <- "/Users/jlamping/Documents/Anthropic_Projects/scf_calibration_app/data_raw/landscape/NOCA_ClimateRegions_09122025.tif"
template_r  <- rast(tif_path)
working_crs <- crs(template_r)

# Reproject calibration boundary to working CRS
cal_vect_proj <- project(cal_vect_raw, working_crs)
cat("  Template CRS:", substr(working_crs, 1, 60), "\n")

cat("--- Extracting ERA-Interim FWI time series ---\n")
fwi_nc <- "/Users/jlamping/University of Oregon Dropbox/James Lamping/Lamping/NPS_postdoc/Spatial/Climate/ERA-Interim/data/fwi.nc"
fwi_daily <- extract_era_timeseries(
  nc_path   = fwi_nc,
  park_vect = cal_vect_proj,
  cal_start = as.Date("1992-01-01"),
  cal_end   = as.Date("2018-12-31")
)
cat(sprintf("  FWI: %d daily values (%s to %s)\n",
            nrow(fwi_daily), min(fwi_daily$date), max(fwi_daily$date)))

cat("--- Loading FPA FOD ignitions ---\n")
fpa_gdb <- "/Users/jlamping/Documents/Anthropic_Projects/scf_calibration_app/data_raw/FPA_FOD/FPA_FOD_20221014.gdb"
ign_df <- load_fpa_fod(
  gdb_path    = fpa_gdb,
  cal_vect    = cal_vect_proj,
  working_crs = working_crs,
  year_min    = 1992,
  year_max    = 2018
)

cat("--- Building model data ---\n")
model_data <- build_ignition_model_data(ign_df, fwi_daily)
cat(sprintf("  %d type × day rows; zero days: %.1f%%\n",
            nrow(model_data),
            100 * mean(model_data$n == 0)))

# =============================================================================
# 2. FIT ALL MODELS
# =============================================================================

fit_all <- function(d_type) {
  d <- model_data %>% filter(ignition_type == d_type, !is.na(FWI))

  pois <- glm(n ~ FWI, data = d, family = poisson(link = "log"))
  zip  <- zeroinfl(n ~ FWI | FWI, data = d, dist = "poisson", link = "logit")

  list(poisson = pois, zip = zip, data = d, type = d_type)
}

cat("--- Fitting models ---\n")
fits <- map(c("Accidental", "Lightning"), fit_all)
names(fits) <- c("Accidental", "Lightning")

# =============================================================================
# 3. EXTRACT COEFFICIENTS & PRINT SUMMARY
# =============================================================================

coef_summary <- map_dfr(fits, function(f) {
  p <- f$poisson
  z <- f$zip
  bind_rows(
    tibble(
      type         = f$type,
      model        = "Poisson",
      b0           = coef(p)[1],
      b1           = coef(p)[2],
      bz0          = NA_real_,
      bz1          = NA_real_,
      AIC          = AIC(p),
      log_lik      = as.numeric(logLik(p))
    ),
    tibble(
      type         = f$type,
      model        = "ZIP",
      b0           = coef(z)["count_(Intercept)"],
      b1           = coef(z)["count_FWI"],
      bz0          = coef(z)["zero_(Intercept)"],
      bz1          = coef(z)["zero_FWI"],
      AIC          = AIC(z),
      log_lik      = as.numeric(logLik(z))
    )
  )
})

# Fix NA coef names (pscl version differences: "count_(Intercept)" vs "count.(Intercept)")
fix_coefs <- function(f) {
  z   <- f$zip
  co  <- coef(z)
  nms <- names(co)
  get_c <- function(pat) { idx <- grep(pat, nms, perl = TRUE); if (!length(idx)) NA_real_ else unname(co[idx[1]]) }
  list(b0 = get_c("count.*Intercept"), b1 = get_c("count.*FWI"),
       bz0 = get_c("zero.*Intercept"), bz1 = get_c("zero.*FWI"))
}

coef_zip <- map(fits, fix_coefs)
coef_pois <- map(fits, function(f) list(b0 = unname(coef(f$poisson)[1]), b1 = unname(coef(f$poisson)[2])))

cat("\n========== COEFFICIENT SUMMARY ==========\n")
print(coef_summary %>% mutate(across(where(is.numeric), ~round(., 4))), n = Inf)

sink("analysis/plots/06_model_summary.txt")
cat("SCF Ignition Model Comparison — Poisson vs ZIP\n")
cat("Data: FPA FOD 1992–2018, ERA-Interim FWI, WA North Cascades calibration region\n\n")
cat("ZIP model formula: n ~ FWI | FWI  (pscl::zeroinfl, Poisson count, logit zero)\n\n")
print(coef_summary %>% mutate(across(where(is.numeric), ~round(., 4))), n = Inf)
cat("\nNOTE: Lower AIC = better fit. AIC penalizes ZIP for 2 extra parameters.\n")
cat("\nZIP SIGN CONVENTION:\n")
cat("  pi (structural zero prob) = plogis(bz0 + bz1 * FWI)\n")
cat("  i.e. bz1 > 0 means MORE structural zeros at HIGH FWI (physically odd)\n")
cat("       bz1 < 0 means FEWER structural zeros at HIGH FWI (more plausible)\n")
sink()

# =============================================================================
# 4. PREDICTION HELPER — CORRECT ZIP FORMULA
# =============================================================================

fwi_seq <- seq(0, max(model_data$FWI, na.rm = TRUE), length.out = 200)

pred_poisson <- function(b0, b1, fwi) exp(b0 + b1 * fwi)

# CORRECT ZIP prediction (note: plogis is the standard logistic / inverse-logit)
pred_zip_correct <- function(b0, b1, bz0, bz1, fwi) {
  lambda <- exp(b0 + b1 * fwi)
  pi     <- plogis(bz0 + bz1 * fwi)        # structural zero probability
  (1 - pi) * lambda
}

# BUGGY ZIP prediction as currently coded in the app (sign error)
pred_zip_bug <- function(b0, b1, bz0, bz1, fwi) {
  lambda <- exp(b0 + b1 * fwi)
  p_zero <- 1 / (1 + exp(bz0 + bz1 * fwi)) # missing negation — computes 1-pi
  (1 - p_zero) * lambda
}

# =============================================================================
# 5. PLOT 01: THE SIGN BUG — what the app shows vs what it should show
# =============================================================================

cat("--- Plot 01: ZIP sign bug ---\n")

obs_binned <- model_data %>%
  mutate(FWI_bin = round(FWI)) %>%
  group_by(FWI_bin, ignition_type) %>%
  summarise(mean_n = mean(n), .groups = "drop")

bug_df <- map_dfr(c("Accidental", "Lightning"), function(typ) {
  cz <- coef_zip[[typ]]
  tibble(
    FWI           = fwi_seq,
    ignition_type = typ,
    `App (buggy)`     = pred_zip_bug(cz$b0, cz$b1, cz$bz0, cz$bz1, fwi_seq),
    `Corrected ZIP`   = pred_zip_correct(cz$b0, cz$b1, cz$bz0, cz$bz1, fwi_seq)
  ) %>%
    pivot_longer(c(`App (buggy)`, `Corrected ZIP`),
                 names_to = "version", values_to = "pred")
})

p01 <- ggplot() +
  geom_col(data = obs_binned, aes(FWI_bin, mean_n), fill = "#c8d8c8", width = 1) +
  geom_line(data = bug_df, aes(FWI, pred, color = version, linetype = version),
            linewidth = 1.1) +
  scale_color_manual(values = c("App (buggy)" = "#c0392b", "Corrected ZIP" = "#2980b9")) +
  scale_linetype_manual(values = c("App (buggy)" = "dashed", "Corrected ZIP" = "solid")) +
  facet_wrap(~ignition_type, scales = "free_y") +
  labs(title    = "ZIP Prediction: App Bug vs Corrected Formula",
       subtitle = "Bug: p_zero = 1/(1+exp(bz)) instead of 1/(1+exp(-bz))  |  Bars = observed mean daily ignitions",
       x = "Fire Weather Index (FWI)", y = "Mean daily ignitions",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "#2c5f2e"),
        strip.text       = element_text(color = "white", face = "bold"),
        legend.position  = "bottom")

ggsave("analysis/plots/01_zip_sign_bug.png", p01, width = 10, height = 5, dpi = 150)

# =============================================================================
# 6. PLOT 02: Poisson vs corrected ZIP side-by-side
# =============================================================================

cat("--- Plot 02: Poisson vs corrected ZIP ---\n")

model_pred_df <- map_dfr(c("Accidental", "Lightning"), function(typ) {
  cp <- coef_pois[[typ]]
  cz <- coef_zip[[typ]]
  bind_rows(
    tibble(FWI = fwi_seq, ignition_type = typ, model = "Poisson",
           pred = pred_poisson(cp$b0, cp$b1, fwi_seq)),
    tibble(FWI = fwi_seq, ignition_type = typ, model = "ZIP (corrected)",
           pred = pred_zip_correct(cz$b0, cz$b1, cz$bz0, cz$bz1, fwi_seq))
  )
})

aic_labels <- coef_summary %>%
  group_by(type) %>%
  mutate(label = sprintf("%s  AIC=%.0f", model, AIC)) %>%
  ungroup()

p02 <- ggplot() +
  geom_col(data = obs_binned, aes(FWI_bin, mean_n), fill = "#c8d8c8", width = 1) +
  geom_line(data = model_pred_df,
            aes(FWI, pred, color = model, linetype = model), linewidth = 1.1) +
  scale_color_manual(values = c("Poisson" = "#a0522d", "ZIP (corrected)" = "#2980b9")) +
  scale_linetype_manual(values = c("Poisson" = "solid", "ZIP (corrected)" = "dashed")) +
  facet_wrap(~ignition_type, scales = "free_y") +
  labs(title    = "Ignition Model Comparison: Poisson vs ZIP",
       subtitle = "Bars = observed mean daily ignitions per FWI bin",
       x = "Fire Weather Index (FWI)", y = "Mean daily ignitions",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "#2c5f2e"),
        strip.text       = element_text(color = "white", face = "bold"),
        legend.position  = "bottom") +
  geom_text(data = coef_summary %>% filter(model == "Poisson"),
            aes(x = 40, y = Inf,
                label = sprintf("Poisson AIC: %.0f", AIC)),
            vjust = 2, hjust = 1, size = 3.5, color = "#a0522d", inherit.aes = FALSE) +
  geom_text(data = coef_summary %>% filter(model == "ZIP"),
            aes(x = 40, y = Inf,
                label = sprintf("ZIP AIC: %.0f", AIC)),
            vjust = 3.8, hjust = 1, size = 3.5, color = "#2980b9", inherit.aes = FALSE)

ggsave("analysis/plots/02_poisson_vs_zip.png", p02, width = 10, height = 5, dpi = 150)

# =============================================================================
# 7. PLOT 03: Zero-inflation probability curve (pi vs FWI)
# =============================================================================

cat("--- Plot 03: Zero-inflation probability (pi) ---\n")

pi_df <- map_dfr(c("Accidental", "Lightning"), function(typ) {
  cz <- coef_zip[[typ]]
  tibble(
    FWI           = fwi_seq,
    ignition_type = typ,
    pi_structural = plogis(cz$bz0 + cz$bz1 * fwi_seq),
    lambda        = exp(cz$b0 + cz$b1 * fwi_seq)
  )
})

# How many observed zeros per FWI bin (for reference)
zero_frac <- model_data %>%
  mutate(FWI_bin = round(FWI)) %>%
  group_by(FWI_bin, ignition_type) %>%
  summarise(zero_frac = mean(n == 0), n_days = n(), .groups = "drop")

p03 <- ggplot() +
  geom_col(data = zero_frac %>% filter(n_days >= 5),
           aes(FWI_bin, zero_frac), fill = "#e8d5b0", width = 1) +
  geom_line(data = pi_df, aes(FWI, pi_structural),
            color = "#8e44ad", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  facet_wrap(~ignition_type) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(title    = "ZIP Zero-Inflation Component: π (structural zero probability) vs FWI",
       subtitle = "Bars = observed fraction of zero-count days per FWI bin (≥5 days)\nLine = model-fitted π = plogis(bz0 + bz1 × FWI)",
       x = "Fire Weather Index (FWI)", y = "Probability of structural zero") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "#2c5f2e"),
        strip.text       = element_text(color = "white", face = "bold"))

ggsave("analysis/plots/03_zero_inflation_curve.png", p03, width = 10, height = 5, dpi = 150)

# =============================================================================
# 8. PLOT 04: Count component (lambda) — the part ZIP passes to SCF
# =============================================================================

cat("--- Plot 04: Lambda (count component) comparison ---\n")

lambda_df <- map_dfr(c("Accidental", "Lightning"), function(typ) {
  cp <- coef_pois[[typ]]
  cz <- coef_zip[[typ]]
  bind_rows(
    tibble(FWI = fwi_seq, ignition_type = typ, model = "Poisson λ",
           value = exp(cp$b0 + cp$b1 * fwi_seq)),
    tibble(FWI = fwi_seq, ignition_type = typ, model = "ZIP λ (count component only)",
           value = exp(cz$b0 + cz$b1 * fwi_seq))
  )
})

p04 <- ggplot(lambda_df, aes(FWI, value, color = model, linetype = model)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Poisson λ" = "#a0522d",
                                 "ZIP λ (count component only)" = "#2980b9")) +
  scale_linetype_manual(values = c("Poisson λ" = "solid",
                                    "ZIP λ (count component only)" = "dashed")) +
  facet_wrap(~ignition_type, scales = "free_y") +
  labs(title    = "Count Component (λ) Comparison",
       subtitle = "λ = exp(b0 + b1 × FWI)  — this is what SCF's B0/B1 parameters control",
       x = "Fire Weather Index (FWI)", y = "λ (expected count | non-zero day)",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "#2c5f2e"),
        strip.text       = element_text(color = "white", face = "bold"),
        legend.position  = "bottom")

ggsave("analysis/plots/04_count_components.png", p04, width = 10, height = 5, dpi = 150)

# =============================================================================
# 9. PLOT 05: Deviance residuals — does ZIP actually improve the fit?
# =============================================================================

cat("--- Plot 05: Residuals ---\n")

resid_df <- map_dfr(c("Accidental", "Lightning"), function(typ) {
  f  <- fits[[typ]]
  d  <- f$data

  pois_mu  <- predict(f$poisson, type = "response")
  # For ZIP, compute fitted values manually (pscl predict can be inconsistent)
  cz      <- coef_zip[[typ]]
  pi_hat  <- plogis(cz$bz0 + cz$bz1 * d$FWI)
  lam_hat <- exp(cz$b0    + cz$b1    * d$FWI)
  zip_mu  <- (1 - pi_hat) * lam_hat

  bind_rows(
    tibble(FWI = d$FWI, ignition_type = typ, model = "Poisson",
           residual = d$n - pois_mu, n = d$n),
    tibble(FWI = d$FWI, ignition_type = typ, model = "ZIP (corrected)",
           residual = d$n - zip_mu, n = d$n)
  )
})

# Bin residuals by FWI for visibility
resid_binned <- resid_df %>%
  mutate(FWI_bin = round(FWI / 2) * 2) %>%
  group_by(FWI_bin, ignition_type, model) %>%
  summarise(mean_resid = mean(residual), sd_resid = sd(residual), .groups = "drop")

p05 <- ggplot(resid_binned, aes(FWI_bin, mean_resid, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_ribbon(aes(ymin = mean_resid - sd_resid, ymax = mean_resid + sd_resid,
                  fill = model), alpha = 0.15, color = NA) +
  scale_color_manual(values = c("Poisson" = "#a0522d", "ZIP (corrected)" = "#2980b9")) +
  scale_fill_manual(values  = c("Poisson" = "#a0522d", "ZIP (corrected)" = "#2980b9")) +
  facet_wrap(~ignition_type, scales = "free_y") +
  labs(title    = "Model Residuals by FWI Bin",
       subtitle = "Mean ± SD of (observed − predicted) per 2-unit FWI bin",
       x = "Fire Weather Index (FWI)", y = "Mean residual (observed − predicted)",
       color = NULL, fill = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "#2c5f2e"),
        strip.text       = element_text(color = "white", face = "bold"),
        legend.position  = "bottom")

ggsave("analysis/plots/05_residuals.png", p05, width = 10, height = 5, dpi = 150)

# =============================================================================
# 10. FINAL CONSOLE REPORT
# =============================================================================

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Plots written to analysis/plots/\n\n")

cat("--- AIC COMPARISON ---\n")
coef_summary %>%
  select(type, model, AIC, b0, b1, bz0, bz1) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print(n = Inf)

cat("\n--- ZIP ZERO-INFLATION COEFFICIENTS ---\n")
for (typ in c("Accidental", "Lightning")) {
  cz <- coef_zip[[typ]]
  pi_at_0   <- plogis(cz$bz0)
  pi_at_30  <- plogis(cz$bz0 + cz$bz1 * 30)
  pi_at_60  <- plogis(cz$bz0 + cz$bz1 * 60)
  cat(sprintf("\n%s:\n  bz0=%.4f  bz1=%.4f\n  pi(FWI=0)=%.3f  pi(FWI=30)=%.3f  pi(FWI=60)=%.3f\n",
              typ, cz$bz0, cz$bz1, pi_at_0, pi_at_30, pi_at_60))
}

cat("\n--- RECOMMENDATION ---\n")
for (typ in c("Accidental", "Lightning")) {
  aic_p <- coef_summary %>% filter(type == typ, model == "Poisson") %>% pull(AIC)
  aic_z <- coef_summary %>% filter(type == typ, model == "ZIP")     %>% pull(AIC)
  delta <- aic_p - aic_z
  cat(sprintf("%s: ΔAIC (Pois−ZIP) = %.1f  →  %s\n",
              typ, delta,
              ifelse(delta > 4,  "ZIP has substantially better fit",
              ifelse(delta < -4, "Poisson is better; ZIP not warranted",
                     "Models are essentially equivalent (|ΔAIC| < 4)"))))
}
