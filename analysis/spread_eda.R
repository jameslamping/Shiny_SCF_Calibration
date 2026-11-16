# =============================================================================
# spread_eda.R  —  Fire spread exploratory data analysis
#
# PURPOSE
#   Systematically examine the relationships between predictor variables
#   (FWI, Effective Wind Speed, Fine Fuels, slope, aspect) and observed fire
#   spread in the calibration dataset.  Outputs a multi-page PDF report and
#   a CSV of numerical summaries for discussion.
#
# HOW TO RUN
#   1. Set RDATA_PATH below to your most recent scf_cache_*.RData file.
#   2. Rscript analysis/spread_eda.R          (from the app root)
#      — or open in RStudio and source it.
#   3. Share  analysis/plots/spread_eda_report_<date>.pdf  and
#              analysis/plots/spread_eda_summary_<date>.csv  for review.
#
# WHAT IT EXPLORES
#   Part 1 : Data overview & quality
#   Part 2 : Spread probability — univariate marginal relationships
#   Part 3 : Spread probability — bivariate / interaction surfaces
#   Part 4 : Spread probability — linearity & alternative link diagnostics
#   Part 5 : Max daily area — distribution & predictor relationships
#   Part 6 : Max daily area — linear vs log-linear model comparison
#   Part 7 : Pair-level diagnostics (pct_proximal, gap days, seasonality)
#   Part 8 : Numerical summary table (saved to CSV)
# =============================================================================

# ---- 0. Configuration --------------------------------------------------------
RDATA_PATH  <- "scf_cache_20260508_211333.RData"   # ← update if needed
OUT_DIR     <- "analysis/plots"
STAMP       <- format(Sys.time(), "%Y%m%d_%H%M%S")
PDF_FILE    <- file.path(OUT_DIR, paste0("spread_eda_report_", STAMP, ".pdf"))
CSV_FILE    <- file.path(OUT_DIR, paste0("spread_eda_summary_", STAMP, ".csv"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Packages -------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales",
          "lubridate", "mgcv", "broom", "pROC", "grDevices", "grid", "gridExtra")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

theme_set(theme_bw(base_size = 11))
pal_spread <- c("0" = "#E8601C", "1" = "#4575B4")

# ---- 2. Load data ------------------------------------------------------------
cat("Loading", RDATA_PATH, "...\n")
suppressWarnings(load(RDATA_PATH))

# Data lives either at top level or inside a 'cache' list
.get <- function(nm) {
  if (exists(nm, envir = .GlobalEnv)) get(nm, envir = .GlobalEnv)
  else if (exists("cache") && !is.null(cache[[nm]])) cache[[nm]]
  else NULL
}

pairs   <- .get("pairs_clean")
samples <- .get("spread_prob_samples")
da_data <- .get("spread_da_data")
sp_coef <- .get("spread_prob_coef")
ma_coef <- .get("max_area_coef")
sp_fit  <- .get("spread_prob_fit")
ma_fit  <- .get("max_area_fit")

if (is.null(pairs))   stop("pairs_clean not found in ", RDATA_PATH)
if (is.null(samples)) stop("spread_prob_samples not found — run spread calibration first.")

cat(sprintf("  pairs_clean   : %d rows | %d fires\n",
            nrow(pairs), length(unique(pairs$FIRE_ID))))
cat(sprintf("  spread samples: %d rows | y=1: %.1f%%\n",
            nrow(samples), 100 * mean(samples$y == 1, na.rm = TRUE)))
if (!is.null(da_data))
  cat(sprintf("  daily area    : %d rows\n", nrow(da_data)))

# ---- Helper: bin-mean plot ---------------------------------------------------
bin_mean_plot <- function(df, x_col, y_col = "y", n_bins = 20,
                           x_lab = x_col, title = NULL, log_x = FALSE) {
  d <- df %>%
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
    mutate(.x = .data[[x_col]], .y = as.numeric(.data[[y_col]]))

  if (log_x) d <- d %>% filter(.x > 0)

  d <- d %>%
    mutate(.bin = cut(.x, breaks = n_bins, include.lowest = TRUE,
                      labels = FALSE)) %>%
    group_by(.bin) %>%
    summarise(x_mid  = median(.x),
              y_mean = mean(.y),
              y_se   = sd(.y) / sqrt(n()),
              n      = n(), .groups = "drop")

  p <- ggplot(d, aes(x = x_mid, y = y_mean)) +
    geom_ribbon(aes(ymin = pmax(0, y_mean - 1.96 * y_se),
                    ymax = pmin(1, y_mean + 1.96 * y_se)),
                fill = "steelblue", alpha = 0.25) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(aes(size = n), colour = "steelblue") +
    scale_size_area(max_size = 4, guide = "none") +
    labs(x = x_lab, y = "Mean P(spread)", title = title) +
    ylim(0, 1)

  if (log_x) p <- p + scale_x_log10()
  p
}

# ---- 3. Open PDF device ------------------------------------------------------
cat("Writing report to", PDF_FILE, "\n")
pdf(PDF_FILE, width = 11, height = 8.5)

# ============================================================================
# PART 1 — Data overview
# ============================================================================
cat("Part 1: overview\n")

# 1a. Pairs per fire histogram
p1a <- pairs %>%
  group_by(FIRE_ID) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  ggplot(aes(n_pairs)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
  labs(title = "Pairs per fire event",
       x = "Number of sequential pairs", y = "Fire count")

# 1b. FWI distribution in eligible pairs vs all
p1b <- pairs %>%
  filter(is.finite(FWI)) %>%
  mutate(used = ifelse(use_for_spread_prob, "Used (spread prob)", "Excluded")) %>%
  ggplot(aes(FWI, fill = used)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 40) +
  scale_fill_manual(values = c("Used (spread prob)" = "steelblue",
                                "Excluded" = "grey60")) +
  labs(title = "FWI distribution — all pairs", fill = NULL,
       x = "Fire Weather Index (FWI)", y = "Pair count")

# 1c. EWS distribution
p1c <- pairs %>%
  filter(is.finite(EffectiveWind)) %>%
  ggplot(aes(EffectiveWind)) +
  geom_histogram(bins = 40, fill = "tomato", colour = "white") +
  labs(title = "Effective Wind Speed distribution",
       x = "EWS (km/h)", y = "Pair count")

# 1d. Monthly pair counts
p1d <- pairs %>%
  filter(is.finite(FWI), use_for_spread_prob) %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  count(month) %>%
  ggplot(aes(month, n)) +
  geom_col(fill = "steelblue") +
  labs(title = "Spread probability pairs by month", x = NULL, y = "Pairs")

print((p1a | p1b) / (p1c | p1d) +
        plot_annotation(title = "Part 1 — Dataset overview",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))

# ============================================================================
# PART 2 — Spread probability: univariate marginal relationships
# ============================================================================
cat("Part 2: univariate spread probability\n")

sp <- samples %>%
  filter(is.finite(FWI), is.finite(y))

# 2a. FWI vs P(spread) — binned means + GLM line
glm_fwi <- glm(y ~ FWI, data = sp, family = binomial())
pred_fwi <- data.frame(FWI = seq(min(sp$FWI), max(sp$FWI), length.out = 200))
pred_fwi$y_glm <- predict(glm_fwi, pred_fwi, type = "response")

p2a <- bin_mean_plot(sp, "FWI", title = "FWI vs P(spread)") +
  geom_line(data = pred_fwi, aes(x = FWI, y = y_glm),
            colour = "firebrick", linewidth = 1, linetype = "dashed") +
  annotate("text", x = Inf, y = 0.05, hjust = 1.1,
           label = sprintf("GLM B1_FWI = %.4f", coef(glm_fwi)[2]),
           size = 3, colour = "firebrick")

# 2b. EWS vs P(spread)
p2b_data <- sp %>% filter(is.finite(EffectiveWind))
glm_ews <- glm(y ~ EffectiveWind, data = p2b_data, family = binomial())
pred_ews <- data.frame(EffectiveWind = seq(min(p2b_data$EffectiveWind),
                                            max(p2b_data$EffectiveWind), length.out = 200))
pred_ews$y_glm <- predict(glm_ews, pred_ews, type = "response")

p2b <- bin_mean_plot(p2b_data, "EffectiveWind",
                      x_lab = "Effective Wind Speed (km/h)",
                      title = "EWS vs P(spread)") +
  geom_line(data = pred_ews, aes(x = EffectiveWind, y = y_glm),
            colour = "firebrick", linewidth = 1, linetype = "dashed") +
  annotate("text", x = Inf, y = 0.05, hjust = 1.1,
           label = sprintf("GLM B3_EWS = %.4f", coef(glm_ews)[2]),
           size = 3, colour = "firebrick")

# 2c. Full model (FWI + EWS) vs each predictor (partial residual style)
if (all(c("FWI", "EffectiveWind") %in% names(sp)) && sum(is.finite(sp$EffectiveWind)) > 100) {
  glm_full <- glm(y ~ FWI + EffectiveWind, data = sp, family = binomial())
  sp$resid_full <- residuals(glm_full, type = "response")
  auc_full <- pROC::auc(pROC::roc(sp$y, fitted(glm_full), quiet = TRUE))

  p2c <- ggplot(sp %>% sample_n(min(5000, nrow(sp))),
                aes(x = FWI, y = resid_full)) +
    geom_point(alpha = 0.15, size = 0.6, colour = "steelblue") +
    geom_smooth(method = "loess", se = TRUE, colour = "firebrick",
                linewidth = 0.9, formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = sprintf("Full GLM residuals vs FWI  (AUC = %.3f)", auc_full),
         x = "FWI", y = "Residual (response scale)")

  p2d <- ggplot(sp %>% filter(is.finite(EffectiveWind)) %>%
                  sample_n(min(5000, nrow(.))),
                aes(x = EffectiveWind, y = resid_full)) +
    geom_point(alpha = 0.15, size = 0.6, colour = "tomato") +
    geom_smooth(method = "loess", se = TRUE, colour = "firebrick",
                linewidth = 0.9, formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Full GLM residuals vs EWS",
         x = "Effective Wind Speed (km/h)", y = "Residual")
} else {
  p2c <- ggplot() + labs(title = "Full GLM residuals: EWS not available")
  p2d <- p2c
  auc_full <- NA_real_
}

print((p2a | p2b) / (p2c | p2d) +
        plot_annotation(title = "Part 2 — Univariate spread probability relationships",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))

# ============================================================================
# PART 3 — Spread probability: bivariate interaction surface
# ============================================================================
cat("Part 3: bivariate surface\n")

sp_bi <- sp %>% filter(is.finite(FWI), is.finite(EffectiveWind))

# 3a. 2D heat map: FWI x EWS → mean P(spread)
n_fwi_bins <- 15; n_ews_bins <- 12
sp_bi2 <- sp_bi %>%
  mutate(fwi_cut = cut(FWI, n_fwi_bins, include.lowest = TRUE),
         ews_cut = cut(EffectiveWind, n_ews_bins, include.lowest = TRUE)) %>%
  group_by(fwi_cut, ews_cut) %>%
  summarise(mean_y = mean(y), n = n(), .groups = "drop") %>%
  filter(n >= 30)

p3a <- ggplot(sp_bi2, aes(x = fwi_cut, y = ews_cut, fill = mean_y)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1,
                        limits = c(0, 1), name = "P(spread)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7)) +
  labs(title = "Mean P(spread) by FWI × EWS bin  (n ≥ 30 per cell)",
       x = "FWI bin", y = "EWS bin (km/h)")

# 3b. Interaction: does the FWI slope change with EWS level?
sp_bi3 <- sp_bi %>%
  mutate(ews_tertile = cut(EffectiveWind,
                            quantile(EffectiveWind, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                            labels = c("Low EWS", "Mid EWS", "High EWS"),
                            include.lowest = TRUE))

p3b <- ggplot(sp_bi3, aes(x = FWI, y = as.numeric(y), colour = ews_tertile)) +
  geom_smooth(method = "glm", method.args = list(family = binomial()),
              se = TRUE, formula = y ~ x, linewidth = 1.1) +
  scale_colour_brewer(palette = "Set1", name = "EWS tertile") +
  labs(title = "FWI–spread slope by EWS tertile  (GLM fits)",
       x = "FWI", y = "P(spread)",
       caption = "Diverging slopes → interaction term may be warranted") +
  ylim(0, 1)

# 3c. Same but flipped: EWS slope by FWI tertile
sp_bi4 <- sp_bi %>%
  mutate(fwi_tertile = cut(FWI,
                            quantile(FWI, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                            labels = c("Low FWI", "Mid FWI", "High FWI"),
                            include.lowest = TRUE))

p3c <- ggplot(sp_bi4, aes(x = EffectiveWind, y = as.numeric(y), colour = fwi_tertile)) +
  geom_smooth(method = "glm", method.args = list(family = binomial()),
              se = TRUE, formula = y ~ x, linewidth = 1.1) +
  scale_colour_brewer(palette = "Set2", name = "FWI tertile") +
  labs(title = "EWS–spread slope by FWI tertile  (GLM fits)",
       x = "EWS (km/h)", y = "P(spread)") +
  ylim(0, 1)

print(p3a / (p3b | p3c) +
        plot_annotation(title = "Part 3 — Bivariate interaction surface",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))

# ============================================================================
# PART 4 — Linearity & alternative form diagnostics
# ============================================================================
cat("Part 4: linearity diagnostics\n")

# 4a. GAM vs GLM for FWI
gam_fwi <- mgcv::gam(y ~ s(FWI, k = 8), data = sp, family = binomial(), method = "REML")
pred_gam <- data.frame(FWI = seq(min(sp$FWI), max(sp$FWI), length.out = 300))
pred_gam$gam <- predict(gam_fwi, pred_gam, type = "response")
pred_gam$glm <- predict(glm_fwi, pred_gam, type = "response")

p4a <- ggplot(pred_gam, aes(x = FWI)) +
  geom_line(aes(y = glm, colour = "GLM (linear logit)"), linewidth = 1) +
  geom_line(aes(y = gam, colour = "GAM (spline)"),       linewidth = 1, linetype = "dashed") +
  scale_colour_manual(values = c("GLM (linear logit)" = "firebrick",
                                  "GAM (spline)"       = "steelblue"), name = NULL) +
  labs(title = "GLM vs GAM: FWI effect on P(spread)",
       subtitle = sprintf("GAM deviance explained: %.1f%%",
                           100 * summary(gam_fwi)$dev.expl),
       x = "FWI", y = "P(spread)") +
  ylim(0, 1)

# 4b. GAM partial effect on logit scale — is it linear?
plot_obj <- plot(gam_fwi, pages = 0, residuals = FALSE, rug = FALSE, select = 1)
gam_df <- data.frame(FWI = plot_obj[[1]]$x,
                      fit = plot_obj[[1]]$fit,
                      se  = plot_obj[[1]]$se)

p4b <- ggplot(gam_df, aes(x = FWI)) +
  geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = fit), colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "GAM smooth: FWI effect on logit(P(spread))",
       subtitle = "Linear if this is approximately a straight line",
       x = "FWI", y = "Partial effect (logit scale)")

# 4c. GAM for EWS
if (sum(is.finite(sp$EffectiveWind)) > 200) {
  gam_ews <- mgcv::gam(y ~ s(EffectiveWind, k = 8), data = sp %>% filter(is.finite(EffectiveWind)),
                         family = binomial(), method = "REML")
  plot_ews <- plot(gam_ews, pages = 0, residuals = FALSE, rug = FALSE, select = 1)
  gam_ews_df <- data.frame(EWS = plot_ews[[1]]$x,
                             fit  = plot_ews[[1]]$fit,
                             se   = plot_ews[[1]]$se)
  p4c <- ggplot(gam_ews_df, aes(x = EWS)) +
    geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se),
                fill = "tomato", alpha = 0.25) +
    geom_line(aes(y = fit), colour = "tomato", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "GAM smooth: EWS effect on logit(P(spread))",
         subtitle = "Linear if approximately straight",
         x = "EWS (km/h)", y = "Partial effect (logit scale)")
} else {
  p4c <- ggplot() + labs(title = "EWS GAM: insufficient data")
}

# 4d. AUC comparison table
models_sp <- list(
  "FWI only"     = glm(y ~ FWI,                          data = sp, family = binomial()),
  "EWS only"     = glm(y ~ EffectiveWind,                data = sp %>% filter(is.finite(EffectiveWind)), family = binomial()),
  "FWI + EWS"    = glm(y ~ FWI + EffectiveWind,          data = sp %>% filter(is.finite(EffectiveWind)), family = binomial()),
  "FWI × EWS"   = glm(y ~ FWI * EffectiveWind,          data = sp %>% filter(is.finite(EffectiveWind)), family = binomial()),
  "GAM(FWI)"     = gam_fwi
)

auc_df <- purrr::map_dfr(names(models_sp), function(nm) {
  m <- models_sp[[nm]]
  d <- if (inherits(m, "gam")) sp else model.frame(m)
  preds <- predict(m, type = "response")
  resp  <- if ("y" %in% names(d)) d$y else d[[1]]
  auc_v <- tryCatch(as.numeric(pROC::auc(pROC::roc(resp, preds, quiet = TRUE))), error = function(e) NA)
  dev_e <- if (inherits(m, "gam")) summary(m)$dev.expl else
    1 - m$deviance / m$null.deviance
  tibble(Model = nm, AUC = round(auc_v, 4),
         `Dev. explained` = scales::percent(dev_e, 0.1),
         AIC = round(AIC(m), 1))
})

p4d <- ggplot(auc_df, aes(x = reorder(Model, AUC), y = AUC, fill = AUC)) +
  geom_col() +
  geom_text(aes(label = round(AUC, 3)), hjust = -0.1, size = 3.5) +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") +
  coord_flip() +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(title = "Spread probability model comparison (AUC)",
       x = NULL, y = "AUC (higher = better discrimination)")

print((p4a | p4b) / (p4c | p4d) +
        plot_annotation(title = "Part 4 — Linearity & model form diagnostics",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))

# ============================================================================
# PART 5 — Max daily area: distribution & predictor relationships
# ============================================================================
cat("Part 5: max daily area\n")

da <- if (!is.null(da_data)) {
  da_data %>% filter(is.finite(daily_area_ha), daily_area_ha > 0, is.finite(FWI))
} else {
  pairs %>% filter(use_for_daily_area, is.finite(daily_area_ha), daily_area_ha > 0, is.finite(FWI))
}

if (nrow(da) < 10) {
  cat("  WARNING: insufficient daily area data for Part 5 plots\n")
  da <- NULL
}

if (!is.null(da)) {

  # 5a. Distribution of daily area
  p5a <- ggplot(da, aes(daily_area_ha)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
    labs(title = "Daily area increment distribution (raw scale)",
         x = "Daily area (ha)", y = "Count")

  p5b <- ggplot(da, aes(daily_area_ha)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
    scale_x_log10(labels = scales::comma) +
    labs(title = "Daily area increment (log scale)",
         x = "Daily area (ha, log)", y = "Count",
         caption = "Log-normal shape → log-linear model may fit better")

  # 5c. FWI vs daily area — scatter + quantile lines
  q_cols <- c("50th" = "steelblue", "75th" = "orange", "90th" = "firebrick")

  q_lines <- purrr::map_dfr(c(0.5, 0.75, 0.9), function(tau) {
    x_seq <- seq(min(da$FWI), max(da$FWI), length.out = 100)
    if (requireNamespace("quantreg", quietly = TRUE)) {
      m <- tryCatch(quantreg::rq(daily_area_ha ~ FWI, tau = tau, data = da),
                    error = function(e) NULL)
      if (!is.null(m))
        return(data.frame(FWI = x_seq, y = predict(m, data.frame(FWI = x_seq)),
                           quantile = paste0(tau*100, "th")))
    }
    data.frame(FWI = x_seq,
               y = coef(lm(daily_area_ha ~ FWI, data = da))[1] +
                   coef(lm(daily_area_ha ~ FWI, data = da))[2] * x_seq,
               quantile = paste0(tau*100, "th (OLS)"))
  })

  p5c <- ggplot(da %>% sample_n(min(3000, nrow(da))),
                aes(x = FWI, y = daily_area_ha)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey50") +
    geom_line(data = q_lines, aes(y = y, colour = quantile), linewidth = 1) +
    scale_colour_manual(values = c("50th"="steelblue","75th"="orange","90th"="firebrick"),
                         name = "Quantile") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "FWI vs daily area increment (ha)",
         x = "FWI", y = "Daily area (ha)")

  # 5d. Log(daily area) vs FWI — linearity check
  da <- da %>% mutate(log_area = log(daily_area_ha))
  p5d <- ggplot(da %>% sample_n(min(3000, nrow(da))),
                aes(x = FWI, y = log_area)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey50") +
    geom_smooth(method = "lm",  colour = "firebrick",  se = TRUE, formula = y ~ x,
                linewidth = 1, aes(linetype = "Linear")) +
    geom_smooth(method = "loess", colour = "steelblue", se = FALSE, formula = y ~ x,
                linewidth = 1, aes(linetype = "LOESS")) +
    scale_linetype_manual(values = c(Linear = "dashed", LOESS = "solid"), name = NULL) +
    labs(title = "FWI vs log(daily area ha)",
         subtitle = "LOESS ≈ linear → log-linear model appropriate",
         x = "FWI", y = "log(daily area ha)")

  print((p5a | p5b) / (p5c | p5d) +
          plot_annotation(title = "Part 5 — Max daily area distribution & FWI relationship",
                          theme = theme(plot.title = element_text(size = 14, face = "bold"))))

  # ============================================================================
  # PART 6 — Max daily area: linear vs log-linear model comparison
  # ============================================================================
  cat("Part 6: linear vs log-linear model comparison\n")

  da_ews <- da %>% filter(is.finite(EffectiveWind))

  fit_lin  <- lm(daily_area_ha ~ FWI + EffectiveWind, data = da_ews)
  fit_log  <- lm(log_area      ~ FWI + EffectiveWind, data = da_ews)

  # Residuals vs fitted
  da_ews$resid_lin <- residuals(fit_lin)
  da_ews$resid_log <- residuals(fit_log)
  da_ews$fit_lin   <- fitted(fit_lin)
  da_ews$fit_log   <- fitted(fit_log)

  p6a <- ggplot(da_ews, aes(x = fit_lin, y = resid_lin)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", formula = y ~ x) +
    labs(title = "Linear model: residuals vs fitted",
         subtitle = sprintf("R² = %.3f  RMSE = %.0f ha",
                             summary(fit_lin)$r.squared,
                             sqrt(mean(da_ews$resid_lin^2))),
         x = "Fitted (ha)", y = "Residual (ha)")

  p6b <- ggplot(da_ews, aes(x = fit_log, y = resid_log)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", formula = y ~ x) +
    labs(title = "Log-linear model: residuals vs fitted",
         subtitle = sprintf("R² = %.3f  (log scale)  Back-transform RMSE ≈ %.0f ha",
                             summary(fit_log)$r.squared,
                             sqrt(mean((exp(da_ews$fit_log) - da_ews$daily_area_ha)^2))),
         x = "Fitted (log ha)", y = "Residual (log ha)")

  # Predicted vs observed comparison
  p6c <- ggplot(da_ews %>% sample_n(min(2000, nrow(da_ews))),
                aes(x = daily_area_ha, y = fit_lin)) +
    geom_point(alpha = 0.3, size = 0.8, colour = "steelblue") +
    geom_abline(slope = 1, colour = "firebrick", linetype = "dashed") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Linear: observed vs predicted",
         x = "Observed ha", y = "Predicted ha")

  p6d <- ggplot(da_ews %>% sample_n(min(2000, nrow(da_ews))),
                aes(x = daily_area_ha, y = exp(fit_log))) +
    geom_point(alpha = 0.3, size = 0.8, colour = "tomato") +
    geom_abline(slope = 1, colour = "firebrick", linetype = "dashed") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Log-linear (back-transformed): observed vs predicted",
         x = "Observed ha", y = "Predicted ha")

  # Coefficient table comparison
  coef_comp <- bind_rows(
    broom::tidy(fit_lin) %>% mutate(model = "Linear (raw ha)"),
    broom::tidy(fit_log) %>% mutate(model = "Log-linear (log ha)")
  ) %>%
    select(model, term, estimate, std.error, p.value) %>%
    mutate(across(where(is.numeric), ~round(., 5)))

  p6e <- ggplot(coef_comp %>% filter(term != "(Intercept)"),
                aes(x = term, y = estimate, fill = model)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = estimate - 1.96*std.error,
                       ymax = estimate + 1.96*std.error),
                  position = position_dodge(0.9), width = 0.25) +
    scale_fill_brewer(palette = "Set1", name = NULL) +
    labs(title = "Coefficient comparison: linear vs log-linear",
         subtitle = "Log-linear B1, B2 are NOT directly usable in SCF (see Part 8 note)",
         x = NULL, y = "Coefficient estimate")

  print((p6a | p6b) / (p6c | p6d) +
          plot_annotation(title = "Part 6a — Max daily area: linear vs log-linear residuals",
                          theme = theme(plot.title = element_text(size = 14, face = "bold"))))
  print(p6e +
          plot_annotation(title = "Part 6b — Coefficient comparison",
                          theme = theme(plot.title = element_text(size = 14, face = "bold"))))
}

# ============================================================================
# PART 7 — Pair-level diagnostics: proximal spread, seasonality, gap days
# ============================================================================
cat("Part 7: pair-level diagnostics\n")

sp_pairs <- pairs %>%
  filter(use_for_spread_prob, is.finite(FWI))

# 7a. pct_proximal vs FWI — does fire "coherence" improve with danger?
if ("pct_proximal" %in% names(pairs)) {
  p7a <- sp_pairs %>%
    filter(is.finite(pct_proximal)) %>%
    ggplot(aes(x = FWI, y = pct_proximal)) +
    geom_point(alpha = 0.4, size = 0.8, colour = "steelblue") +
    geom_smooth(method = "loess", se = TRUE, colour = "firebrick",
                formula = y ~ x, linewidth = 1) +
    labs(title = "FWI vs % proximal spread",
         subtitle = "High pct_proximal = new area adjacent to existing perimeter (coherent spread)",
         x = "FWI", y = "% proximal")
} else {
  p7a <- ggplot() + labs(title = "pct_proximal not in pairs data")
}

# 7b. gap_days distribution
p7b <- sp_pairs %>%
  count(gap_days) %>%
  ggplot(aes(x = gap_days, y = n)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = 1:max(sp_pairs$gap_days, na.rm = TRUE)) +
  labs(title = "Distribution of gap days (spread probability pairs)",
       x = "Days between perimeters", y = "Pair count")

# 7c. Seasonal FWI vs mean P(spread)
sp_monthly <- samples %>%
  filter(is.finite(FWI), is.finite(y)) %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  group_by(month) %>%
  summarise(mean_y = mean(y), mean_fwi = mean(FWI),
            n = n(), .groups = "drop")

p7c <- ggplot(sp_monthly, aes(x = month)) +
  geom_col(aes(y = mean_y), fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = mean_fwi / max(mean_fwi, na.rm=TRUE), group = 1),
            colour = "firebrick", linewidth = 1) +
  geom_point(aes(y = mean_fwi / max(mean_fwi, na.rm=TRUE)), colour = "firebrick") +
  scale_y_continuous(
    name = "Mean P(spread)",
    sec.axis = sec_axis(~ . * max(sp_monthly$mean_fwi, na.rm=TRUE),
                         name = "Mean FWI (red line)")
  ) +
  labs(title = "Seasonal pattern: P(spread) and mean FWI by month", x = NULL)

# 7d. FWI distribution by data source (FIRED vs GeoMAC)
if ("FIRE_ID" %in% names(sp_pairs)) {
  sp_pairs <- sp_pairs %>%
    mutate(source = ifelse(grepl("^FIRED-", FIRE_ID), "FIRED", "GeoMAC"))
  p7d <- ggplot(sp_pairs, aes(x = FWI, fill = source)) +
    geom_density(alpha = 0.5) +
    scale_fill_brewer(palette = "Set1", name = "Source") +
    labs(title = "FWI distribution by perimeter data source",
         x = "FWI", y = "Density")
} else {
  p7d <- ggplot() + labs(title = "Source info not available")
}

print((p7a | p7b) / (p7c | p7d) +
        plot_annotation(title = "Part 7 — Pair-level diagnostics & seasonality",
                        theme = theme(plot.title = element_text(size = 14, face = "bold"))))

# ============================================================================
# PART 8 — Numerical summary table
# ============================================================================
cat("Part 8: numerical summary\n")

summary_rows <- list()

# Spread probability model AUCs
summary_rows[["sp_n_samples"]] <- tibble(
  Section = "Spread probability", Metric = "N samples (cell level)",
  Value = as.character(nrow(samples)), Note = ""
)
summary_rows[["sp_base_rate"]] <- tibble(
  Section = "Spread probability", Metric = "Base spread rate (y=1)",
  Value = sprintf("%.3f", mean(samples$y == 1, na.rm = TRUE)),
  Note  = "Imbalanced? Consider weighted GLM or oversampling"
)

for (nm in names(models_sp)) {
  row_auc <- auc_df %>% filter(Model == nm)
  if (nrow(row_auc) > 0)
    summary_rows[[paste0("auc_", nm)]] <- tibble(
      Section = "Spread probability — AUC", Metric = nm,
      Value = as.character(row_auc$AUC), Note = row_auc$`Dev. explained`
    )
}

# Daily area model R²
if (!is.null(da) && exists("fit_lin")) {
  summary_rows[["da_n"]] <- tibble(
    Section = "Max daily area", Metric = "N fire-days",
    Value = as.character(nrow(da)), Note = "")
  summary_rows[["da_lin_r2"]] <- tibble(
    Section = "Max daily area", Metric = "Linear OLS R²",
    Value = sprintf("%.3f", summary(fit_lin)$r.squared),
    Note  = "On raw ha scale — as used in SCF directly")
  summary_rows[["da_log_r2"]] <- tibble(
    Section = "Max daily area", Metric = "Log-linear OLS R²",
    Value = sprintf("%.3f", summary(fit_log)$r.squared),
    Note  = "On log(ha) scale — requires back-transform for SCF")
  summary_rows[["da_75pct"]] <- tibble(
    Section = "Max daily area", Metric = "Observed 75th percentile ha/day",
    Value   = sprintf("%.0f ha", quantile(da$daily_area_ha, 0.75, na.rm = TRUE)),
    Note    = "Recommended calibration target (MaxSpreadArea is a ceiling)")
  summary_rows[["da_90pct"]] <- tibble(
    Section = "Max daily area", Metric = "Observed 90th percentile ha/day",
    Value   = sprintf("%.0f ha", quantile(da$daily_area_ha, 0.90, na.rm = TRUE)),
    Note    = "Upper ceiling target")
}

summary_df <- bind_rows(summary_rows)

# Print to console
cat("\n=== EDA NUMERICAL SUMMARY ===\n")
print(as.data.frame(summary_df), row.names = FALSE)
write.csv(summary_df, CSV_FILE, row.names = FALSE)
cat("\nSaved:", CSV_FILE, "\n")

# Render summary as a ggplot table (last PDF page)
library(gridExtra)
p_table <- tableGrob(summary_df, rows = NULL,
                      theme = ttheme_minimal(base_size = 9))
grid.newpage()
grid.draw(p_table)
grid.text("Part 8 — Numerical summary (also saved to CSV)",
          x = 0.5, y = 0.98, just = "top",
          gp = grid::gpar(fontsize = 13, fontface = "bold"))

# ============================================================================
dev.off()
cat("\nDone. Report:", PDF_FILE, "\n")
cat("Share both the PDF and CSV for discussion.\n")
