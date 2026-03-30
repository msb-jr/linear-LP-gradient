# ==============================================================================
# Script: 8_substrate_lpr_fidelity.R
# Manuscript relevance: 3.5, Fig. 5
# ==============================================================================
# PURPOSE:
#   Elucidate the relationship between substrate delivery (glucose), the linear
#   gradient (m), and the LPR at within-segment and between-segment levels.
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#
# OUTPUT:
#   - 1_output/8_substrate_lpr_fidelity/__execution_time.csv: Runtime log
#   - 1_output/8_substrate_lpr_fidelity/_1_model_summaries.txt: Full lmer summaries
#   - 1_output/8_substrate_lpr_fidelity/_1_stats.csv: Model coefficients and statistics
#   - 1_output/8_substrate_lpr_fidelity/0_notes.txt: Figure legend
#   - 1_output/8_substrate_lpr_fidelity/1_figure.png: 4-panel figure
#
# DATA FILTERS:
#   - Retained segments only: p6e_m > 0 (positive gradient)
#   - Observations with non-NA glucose
#
# ANALYSES:
# ------------------------------------------------------------------------------
# Seven LME models. All use Satterthwaite approximation for t-tests.
# Models 1, 5 produce no figure panels.
#
#     Model 1: Pyruvate ~ Glucose (POINT-LEVEL, no panel)
#         - Establishes pyruvate–glucose association; reported in text.
#         Unit:       Point-level (individual datapoints)
#         Model:      pyruvate ~ glucose + (1|patient/p6e_seg_index)
#         Weighting:  None
#         Hypothesis: β_glucose > 0 (one-sided)
#         Test:       Satterthwaite approximation
#
#     Model 2: LPR ~ 1/Glucose in Type Pb segments (POINT-LEVEL → Panel A)
#         - Within-segment: LPR vs glucose when b > 0.
#           Modelled as lpr ~ inv_glucose to capture hyperbolic shape;
#           plotted on glucose x-axis.
#         Unit:       Point-level, Type Pb segments only
#         Model:      lpr ~ inv_glucose + (1|patient/p6e_seg_index)
#         Weighting:  None
#         Hypothesis: β > 0 (one-sided; theory: b > 0 → LPR ↓ as glucose ↑)
#         Test:       Satterthwaite approximation
#
#     Model 3: LPR ~ 1/Glucose in Type Nb segments (POINT-LEVEL → Panel B)
#         - Within-segment: LPR vs glucose when b < 0.
#         Unit:       Point-level, Type Nb segments only
#         Model:      lpr ~ inv_glucose + (1|patient/p6e_seg_index)
#         Weighting:  None
#         Hypothesis: β < 0 (one-sided; theory: b < 0 → LPR ↑ as glucose ↑)
#         Test:       Satterthwaite approximation
#
#     Model 4: LPR ~ m (SEGMENT-LEVEL → Panel C)
#         - How faithfully does LPR track the underlying gradient?
#         Unit:       Segment-level (median LPR, gradient per segment)
#         Model:      median_lpr ~ m + (1|patient)
#         Weighting:  n_points (segment size)
#         Hypothesis: β > 0 (one-sided); theoretical β = 1 if LPR = m
#         Test:       Satterthwaite approximation
#
#     Model 5: b ~ m (SEGMENT-LEVEL, no panel)
#         - Exploratory: quantifies the m–b coupling
#         Unit:       Segment-level
#         Model:      b ~ m + (1|patient)
#         Weighting:  n_points (segment size)
#         Hypothesis: Two-sided (exploratory)
#         Test:       Satterthwaite approximation
#
#     Model 6: m ~ Glucose (SEGMENT-LEVEL → Panel D)
#         - Is substrate availability associated with anaerobic activity?
#         Unit:       Segment-level (gradient, median glucose per segment)
#         Model:      m ~ glucose + (1|patient)
#         Weighting:  n_points (segment size)
#         Hypothesis: β < 0 (higher glucose → lower m; one-sided)
#         Test:       Satterthwaite approximation
#
#     Model 7: LPR ~ Glucose (SEGMENT-LEVEL → Panel D)
#         - Between-segment: does LPR track glucose at the segment level?
#         Unit:       Segment-level (median LPR, median glucose per segment)
#         Model:      median_lpr ~ glucose + (1|patient)
#         Weighting:  n_points (segment size)
#         Hypothesis: Two-sided (exploratory)
#         Test:       Satterthwaite approximation
#
# Figure panels:
#     Panel A: Within-segment LPR vs glucose, Type Pb (Model 2; hyperbolic fit)
#     Panel B: Within-segment LPR vs glucose, Type Nb (Model 3; hyperbolic fit)
#     Panel C: Between-segment median LPR vs m (Model 4)
#     Panel D: Between-segment m (Model 6) and LPR (Model 7) vs median glucose
#
# MULTIPLE COMPARISONS CORRECTION:
#   Applied only to Models 2 and 3 (Benjamini-Hochberg FDR, n=2). These two models jointly test 
#   the overarching hyperbolic hypothesis (LPR vs 1/Glucose) across two mutually 
#   exclusive subsets (Pb and Nb). The remaining models are uncorrected because they 
#   address distinct, independent physiological questions:
#     - Model 1: Validates the fundamental Pyruvate-Glucose link.
#     - Model 4: Independent assessment of LPR-m association.
#     - Model 5: Exploratory quantification of the m-b coupling.
#     - Model 6: Tests the primary hypothesis linking m to substrate availability.
#     - Model 7: Exploratory comparison of LPR tracking vs m tracking.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 8_substrate_lpr_fidelity.R\n")

start_time <- Sys.time()

save_plot_portable <- function(filename, plot, ...) {
  if (isTRUE(capabilities("cairo"))) {
    ggsave(filename = filename, plot = plot, type = "cairo", ...)
  } else {
    warning(sprintf("Cairo backend unavailable; saving %s with the default device.", basename(filename)))
    ggsave(filename = filename, plot = plot, ...)
  }
}

safe_divide <- function(numerator, denominator, tol = 1e-12, fill = NA_real_) {
  result <- numerator / denominator
  invalid <- !is.finite(result) | !is.finite(denominator) | abs(denominator) <= tol
  result[invalid] <- fill
  result
}

get_lmer_pred_se <- function(model, newdata) {
  fit <- predict(model, newdata = newdata, re.form = NA, allow.new.levels = TRUE)
  X <- model.matrix(delete.response(terms(model)), newdata)
  vcov_mat <- as.matrix(vcov(model))
  se <- sqrt(rowSums((X %*% vcov_mat) * X))
  list(fit = fit, se = se)
}

extract_coef <- function(model, predictor, one_sided = FALSE, expected_sign = 1) {
  coefs <- summary(model)$coefficients
  beta <- coefs[predictor, "Estimate"]
  se <- coefs[predictor, "Std. Error"]
  t_val <- coefs[predictor, "t value"]
  p_two <- coefs[predictor, "Pr(>|t|)"]
  if (one_sided) {
    p_val <- if ((expected_sign > 0 && beta > 0) || (expected_sign < 0 && beta < 0)) p_two / 2 else 1 - p_two / 2
  } else {
    p_val <- p_two
  }
  list(beta = beta, se = se, t_value = t_val, p_value = p_val)
}

format_p_text <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

# ==============================================================================
# Input / Output
# ==============================================================================
INPUT_FILE <- here::here("1_output", "2_linear_segmentation", "1_results.csv")

if (!file.exists(INPUT_FILE)) {
  cat("\n\nError: Results file not found. Run 2_linear_segmentation.py first.\n")
  stop()
}

OUTPUT_DIR <- here::here("1_output", "8_substrate_lpr_fidelity")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Data Loading and Preparation
# ==============================================================================

raw_df <- read_csv(INPUT_FILE, show_col_types = FALSE, progress = FALSE)

retained_df <- raw_df %>%
  filter(is.finite(p6e_m), p6e_m > 0)

# Point-level data
point_df <- retained_df %>%
  filter(!is.na(glucose)) %>%
  select(patient, p6e_seg_index, m = p6e_m, b = p6e_b, glucose, pyruvate, lpr) %>%
  mutate(
    type_b = case_when(b > 0 ~ "Pb", b < 0 ~ "Nb", TRUE ~ "Zb"),
    inv_glucose = safe_divide(1, glucose)
  )

point_pb <- point_df %>% filter(type_b == "Pb")
point_nb <- point_df %>% filter(type_b == "Nb")

# Segment-level data
segment_df <- point_df %>%
  group_by(patient, p6e_seg_index) %>%
  summarise(
    median_lpr = median(lpr, na.rm = TRUE),
    m = first(m),
    b = first(b),
    glucose = median(glucose, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  )

n_points_total <- nrow(point_df)
n_segments <- nrow(segment_df)
n_patients <- n_distinct(segment_df$patient)

n_points_pb <- nrow(point_pb)
n_seg_pb <- n_distinct(paste(point_pb$patient, point_pb$p6e_seg_index))
n_pat_pb <- n_distinct(point_pb$patient)

n_points_nb <- nrow(point_nb)
n_seg_nb <- n_distinct(paste(point_nb$patient, point_nb$p6e_seg_index))
n_pat_nb <- n_distinct(point_nb$patient)

# ==============================================================================
# Statistical Models
# ==============================================================================

# --- Model 1: Pyruvate (mM) ~ Glucose (POINT-LEVEL, no panel) ---
model_1 <- lmer(pyruvate ~ glucose + (1|patient/p6e_seg_index), data = point_df)
stats_1 <- extract_coef(model_1, "glucose", one_sided = TRUE, expected_sign = 1)
resid_sd_pyr <- sigma(model_1)

# --- Model 2: LPR ~ 1/Glucose in Type Pb (POINT-LEVEL → Panel A) ---
model_2 <- lmer(lpr ~ inv_glucose + (1|patient/p6e_seg_index), data = point_pb)
stats_2 <- extract_coef(model_2, "inv_glucose", one_sided = TRUE, expected_sign = 1)

# --- Model 3: LPR ~ 1/Glucose in Type Nb (POINT-LEVEL → Panel B) ---
model_3 <- lmer(lpr ~ inv_glucose + (1|patient/p6e_seg_index), data = point_nb)
stats_3 <- extract_coef(model_3, "inv_glucose", one_sided = TRUE, expected_sign = -1)

# Apply Benjamini-Hochberg FDR correction for Models 2 & 3 (family of 2 tests for hyperbolic hypothesis)
p_adj_23 <- p.adjust(c(stats_2$p_value, stats_3$p_value), method = "BH")
stats_2$p_value <- p_adj_23[1]
stats_3$p_value <- p_adj_23[2]

# --- Model 4: LPR ~ m (SEGMENT-LEVEL → Panel C) ---
model_4 <- lmer(median_lpr ~ m + (1|patient), data = segment_df, weights = n_points,
                control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
stats_4 <- extract_coef(model_4, "m", one_sided = TRUE, expected_sign = 1)

# --- Model 5: b ~ m (SEGMENT-LEVEL, no panel) ---
model_5 <- lmer(b ~ m + (1|patient), data = segment_df, weights = n_points,
                control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
stats_5 <- extract_coef(model_5, "m", one_sided = FALSE)

# --- Model 6: m ~ Glucose (SEGMENT-LEVEL → Panel D) ---
model_6 <- lmer(m ~ glucose + (1|patient), data = segment_df, weights = n_points,
                control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
stats_6 <- extract_coef(model_6, "glucose", one_sided = TRUE, expected_sign = -1)

# --- Model 7: LPR ~ Glucose (SEGMENT-LEVEL → Panel D) ---
model_7 <- lmer(median_lpr ~ glucose + (1|patient), data = segment_df, weights = n_points,
                control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
stats_7 <- extract_coef(model_7, "glucose", one_sided = FALSE)

p_text_1 <- format_p_text(stats_1$p_value)
p_text_2 <- format_p_text(stats_2$p_value)
p_text_3 <- format_p_text(stats_3$p_value)
p_text_4 <- format_p_text(stats_4$p_value)
p_text_5 <- format_p_text(stats_5$p_value)
p_text_6 <- format_p_text(stats_6$p_value)
p_text_7 <- format_p_text(stats_7$p_value)

# Save model summaries
sink(file.path(OUTPUT_DIR, "_1_model_summaries.txt"))
cat("=======================================================================\n")
cat("Model 1: Pyruvate (mM) ~ Glucose  [No panel]  (POINT-LEVEL)\n")
cat("Establishes pyruvate–glucose association; reported in text.\n")
cat(sprintf("Residual SD = %.5f mM (%.1f μM).\n", resid_sd_pyr, resid_sd_pyr * 1000))
cat(sprintf("N points = %d, N segments = %d, N patients = %d\n",
            n_points_total, n_segments, n_distinct(point_df$patient)))
cat("=======================================================================\n")
print(summary(model_1))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 2: LPR ~ 1/Glucose in Type Pb  [Panel A]  (POINT-LEVEL)\n")
cat("Within-segment: LPR–glucose association when b > 0.\n")
cat("Modelled as lpr ~ inv_glucose; plotted on glucose axis (hyperbolic shape).\n")
cat(sprintf("N points = %d, N segments = %d, N patients = %d\n",
            n_points_pb, n_seg_pb, n_pat_pb))
cat(sprintf("Benjamini-Hochberg FDR corrected p-value (n=2): %g\n", stats_2$p_value))
cat("=======================================================================\n")
print(summary(model_2))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 3: LPR ~ 1/Glucose in Type Nb  [Panel B]  (POINT-LEVEL)\n")
cat("Within-segment: LPR–glucose association when b < 0.\n")
cat("Modelled as lpr ~ inv_glucose; plotted on glucose axis (hyperbolic shape).\n")
cat(sprintf("N points = %d, N segments = %d, N patients = %d\n",
            n_points_nb, n_seg_nb, n_pat_nb))
cat(sprintf("Benjamini-Hochberg FDR corrected p-value (n=2): %g\n", stats_3$p_value))
cat("=======================================================================\n")
print(summary(model_3))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 4: Median LPR ~ m  [Panel C]  (SEGMENT-LEVEL)\n")
cat("How faithfully does LPR track the underlying gradient?\n")
cat("Theoretical β = 1 if LPR perfectly tracked m.\n")
cat(sprintf("N segments = %d, N patients = %d\n", n_segments, n_patients))
cat("=======================================================================\n")
print(summary(model_4))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 5: b ~ m  [No panel]  (SEGMENT-LEVEL)\n")
cat("Exploratory: quantifies the m–b coupling.\n")
cat(sprintf("N segments = %d, N patients = %d\n", n_segments, n_patients))
cat("=======================================================================\n")
print(summary(model_5))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 6: m ~ Glucose  [Panel D]  (SEGMENT-LEVEL)\n")
cat("Between-segment: is substrate availability associated with gradient?\n")
cat(sprintf("N segments = %d, N patients = %d\n", n_segments, n_patients))
cat("=======================================================================\n")
print(summary(model_6))
cat("\n\n")
cat("=======================================================================\n")
cat("Model 7: Median LPR ~ Glucose  [Panel D]  (SEGMENT-LEVEL)\n")
cat("Between-segment: does segment-level LPR track glucose?\n")
cat(sprintf("N segments = %d, N patients = %d\n", n_segments, n_patients))
cat("=======================================================================\n")
print(summary(model_7))
sink()

# ==============================================================================
# Visualization (4-Panel Figure)
# ==============================================================================

common_theme <- theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
    plot.tag.position = c(0.02, 0.99),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = margin(t = 3, r = 3, b = 0, l = 3, unit = "mm")
  )

# Viewports (point-level for A–B, segment-level for C–D)
glucose_range_pb <- quantile(point_pb$glucose, c(0.025, 0.975), na.rm = TRUE)
glucose_range_nb <- quantile(point_nb$glucose, c(0.025, 0.975), na.rm = TRUE)
glucose_range_seg <- quantile(segment_df$glucose, c(0.025, 0.975), na.rm = TRUE)
m_range <- quantile(segment_df$m, c(0.025, 0.975), na.rm = TRUE)

glucose_seq_pb <- seq(glucose_range_pb[1], glucose_range_pb[2], length.out = 100)
glucose_seq_nb <- seq(glucose_range_nb[1], glucose_range_nb[2], length.out = 100)
glucose_seq_seg <- seq(glucose_range_seg[1], glucose_range_seg[2], length.out = 100)
m_seq <- seq(m_range[1], m_range[2], length.out = 100)

# --- Panel A: LPR ~ 1/Glucose within Pb (Model 2, plotted on glucose axis) ---
pred_a <- data.frame(inv_glucose = safe_divide(1, glucose_seq_pb))
res_a <- get_lmer_pred_se(model_2, pred_a)
pred_a$glucose <- glucose_seq_pb
pred_a$y <- res_a$fit
pred_a$ymin <- pred_a$y - 1.96 * res_a$se
pred_a$ymax <- pred_a$y + 1.96 * res_a$se

ylim_a <- range(c(pred_a$ymin[1], pred_a$ymax[1],
                  pred_a$ymin[nrow(pred_a)], pred_a$ymax[nrow(pred_a)]))

p_a <- ggplot() +
  geom_ribbon(data = pred_a, aes(x = glucose, ymin = ymin, ymax = ymax),
              fill = "#44AA99", alpha = 0.2) +
  geom_line(data = pred_a, aes(x = glucose, y = y),
            color = "#44AA99", linewidth = 1) +
  coord_cartesian(xlim = glucose_range_pb, ylim = ylim_a) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(tag = "A", x = "Glucose (mM)", y = "LPR") +
  common_theme

# --- Panel B: LPR ~ 1/Glucose within Nb (Model 3, plotted on glucose axis) ---
pred_b <- data.frame(inv_glucose = safe_divide(1, glucose_seq_nb))
res_b <- get_lmer_pred_se(model_3, pred_b)
pred_b$glucose <- glucose_seq_nb
pred_b$y <- res_b$fit
pred_b$ymin <- pred_b$y - 1.96 * res_b$se
pred_b$ymax <- pred_b$y + 1.96 * res_b$se

ylim_b <- range(c(pred_b$ymin[1], pred_b$ymax[1],
                  pred_b$ymin[nrow(pred_b)], pred_b$ymax[nrow(pred_b)]))

p_b <- ggplot() +
  geom_ribbon(data = pred_b, aes(x = glucose, ymin = ymin, ymax = ymax),
              fill = "#DDCC77", alpha = 0.2) +
  geom_line(data = pred_b, aes(x = glucose, y = y),
            color = "#DDCC77", linewidth = 1) +
  coord_cartesian(xlim = glucose_range_nb, ylim = ylim_b) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(tag = "B", x = "Glucose (mM)", y = "LPR") +
  common_theme

# --- Panel C: LPR ~ m (Model 4, segment-level) ---
pred_c <- data.frame(m = m_seq)
res_c <- get_lmer_pred_se(model_4, pred_c)
pred_c$y <- res_c$fit
pred_c$ymin <- pred_c$y - 1.96 * res_c$se
pred_c$ymax <- pred_c$y + 1.96 * res_c$se

ref_c <- data.frame(m = m_range, lpr = m_range)

ylim_c <- range(c(pred_c$ymin[1], pred_c$ymax[1],
                  pred_c$ymin[nrow(pred_c)], pred_c$ymax[nrow(pred_c)],
                  m_range))

p_c <- ggplot() +
  geom_line(data = ref_c, aes(x = m, y = lpr),
            linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_ribbon(data = pred_c, aes(x = m, ymin = ymin, ymax = ymax),
              fill = "#CC3311", alpha = 0.2) +
  geom_line(data = pred_c, aes(x = m, y = y),
            color = "#CC3311", linewidth = 1) +
  coord_cartesian(xlim = m_range, ylim = ylim_c) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(tag = "C", x = expression(italic(m)), y = "LPR") +
  common_theme

# --- Panel D: m AND LPR ~ Glucose (Models 6 & 7, segment-level) ---
pred_d_m <- data.frame(glucose = glucose_seq_seg)
res_d_m <- get_lmer_pred_se(model_6, pred_d_m)
pred_d_m$y <- res_d_m$fit
pred_d_m$ymin <- pred_d_m$y - 1.96 * res_d_m$se
pred_d_m$ymax <- pred_d_m$y + 1.96 * res_d_m$se

pred_d_lpr <- data.frame(glucose = glucose_seq_seg)
res_d_lpr <- get_lmer_pred_se(model_7, pred_d_lpr)
pred_d_lpr$y <- res_d_lpr$fit
pred_d_lpr$ymin <- pred_d_lpr$y - 1.96 * res_d_lpr$se
pred_d_lpr$ymax <- pred_d_lpr$y + 1.96 * res_d_lpr$se

ylim_d <- range(c(pred_d_m$ymin[1], pred_d_m$ymax[1],
                  pred_d_m$ymin[nrow(pred_d_m)], pred_d_m$ymax[nrow(pred_d_m)],
                  pred_d_lpr$ymin[1], pred_d_lpr$ymax[1],
                  pred_d_lpr$ymin[nrow(pred_d_lpr)], pred_d_lpr$ymax[nrow(pred_d_lpr)]))

p_d <- ggplot() +
  geom_ribbon(data = pred_d_lpr, aes(x = glucose, ymin = ymin, ymax = ymax),
              fill = "#E69F00", alpha = 0.2) +
  geom_line(data = pred_d_lpr, aes(x = glucose, y = y),
            color = "#E69F00", linewidth = 1) +
  geom_ribbon(data = pred_d_m, aes(x = glucose, ymin = ymin, ymax = ymax),
              fill = "#0072B2", alpha = 0.2) +
  geom_line(data = pred_d_m, aes(x = glucose, y = y),
            color = "#0072B2", linewidth = 1) +
  coord_cartesian(xlim = glucose_range_seg, ylim = ylim_d) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(tag = "D", x = "Glucose (mM)", y = "Value") +
  common_theme

# Assemble 4-panel figure
fig <- p_a + p_b + p_c + p_d +
  plot_layout(nrow = 1)

save_plot_portable(
  filename = file.path(OUTPUT_DIR, "1_figure.png"),
  plot = fig,
  width = 185, height = 55, units = "mm", dpi = 300
)

# ==============================================================================
# Summary Statistics CSV
# ==============================================================================

summary_stats <- data.frame(
  Metric = c(
    "N_points", "N_segments", "N_patients",
    "N_points_Pb", "N_segments_Pb", "N_patients_Pb",
    "N_points_Nb", "N_segments_Nb", "N_patients_Nb",
    "Model_1_beta_glucose_mM_per_mM", "Model_1_beta_glucose_uM_per_mM",
    "Model_1_SE_uM_per_mM", "Model_1_p_value", "Model_1_resid_SD_mM",
    "Model_2_beta_inv_glucose_Pb", "Model_2_SE_Pb", "Model_2_p_value_bh_Pb",
    "Model_3_beta_inv_glucose_Nb", "Model_3_SE_Nb", "Model_3_p_value_bh_Nb",
    "Model_4_beta_m", "Model_4_SE", "Model_4_p_value",
    "Model_5_beta_m_bm", "Model_5_SE_bm", "Model_5_p_value_bm",
    "Model_6_beta_glucose", "Model_6_SE", "Model_6_p_value",
    "Model_7_beta_glucose_lpr", "Model_7_SE_lpr", "Model_7_p_value_lpr"
  ),
  Value = c(
    n_points_total, n_segments, n_patients,
    n_points_pb, n_seg_pb, n_pat_pb,
    n_points_nb, n_seg_nb, n_pat_nb,
    stats_1$beta, stats_1$beta * 1000,
    stats_1$se * 1000, stats_1$p_value, resid_sd_pyr,
    stats_2$beta, stats_2$se, stats_2$p_value,
    stats_3$beta, stats_3$se, stats_3$p_value,
    stats_4$beta, stats_4$se, stats_4$p_value,
    stats_5$beta, stats_5$se, stats_5$p_value,
    stats_6$beta, stats_6$se, stats_6$p_value,
    stats_7$beta, stats_7$se, stats_7$p_value
  )
)

write_csv(summary_stats, file.path(OUTPUT_DIR, "_1_stats.csv"))

# ==============================================================================
# Completion Log
# ==============================================================================
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
write_csv(
  data.frame(execution_time_seconds = execution_time),
  file.path(OUTPUT_DIR, "__execution_time.csv")
)

cat("\n\n8_substrate_lpr_fidelity.R complete.\n\n")

# ==============================================================================
# Figure Legend
# ==============================================================================
legends <- list(
  list(
    target = "Figure 1 (1_figure.png)",
    caption = sprintf(
      paste0(
        "Glucose associations with LPR and gradient (m). ",
        "(A) Within-segment LPR vs glucose for Type Pb segments (b > 0; %d segments, %d patients; ",
        "β_1/Glucose = %.2f mM, %s). ",
        "(B) Within-segment LPR vs glucose for Type Nb segments (b < 0; %d segments, %d patients; ",
        "β_1/Glucose = %.2f mM, %s). ",
        "Within-segment models use 1/Glucose as predictor to capture the hyperbolic LPR–glucose ",
        "relationship (LPR = m + b/P); predictions are plotted on the glucose x-axis. ",
        "P-values for within-segment models (A, B) are Benjamini-Hochberg FDR corrected (n=2). ",
        "(C) Segment-level median LPR vs gradient m (β = %.3f, %s); dashed black line indicates the ",
        "theoretical β = 1 relationship expected if LPR perfectly tracked the gradient. ",
        "(D) Segment-level gradient (m, blue; β = %.3f mM⁻¹, %s) and median LPR (orange; β = %.3f mM⁻¹, %s) ",
        "vs median glucose. ",
        "Lines = LME fixed effects; shaded = 95%% CI. ",
        "A–B: nested random intercepts (segment within patient). ",
        "C–D: patient random intercepts, weighted by segment size. ",
        "Sample: %d segments, %d patients."
      ),
      n_seg_pb, n_pat_pb, stats_2$beta, p_text_2,
      n_seg_nb, n_pat_nb, stats_3$beta, p_text_3,
      stats_4$beta, p_text_4,
      stats_6$beta, p_text_6,
      stats_7$beta, p_text_7,
      n_segments, n_patients
    ),
    abbreviations = "b, linear intercept; BH, Benjamini-Hochberg; CI, confidence interval; FDR, false discovery rate; LME, linear mixed-effects model; LPR, lactate/pyruvate ratio; m, linear gradient; Nb, negative intercept (b < 0); Pb, positive intercept (b > 0)"
  ),
  list(
    target = "Text-only results (no panel)",
    caption = sprintf(
      paste0(
        "Model 1 (pyruvate ~ glucose, point-level): β = %.2f μM/mM (%s); ",
        "residual SD = %.1f μM. %d datapoints, %d segments, %d patients. ",
        "Model 5 (b ~ m, segment-level): β = %.4f mM (%s). ",
        "%d segments, %d patients."
      ),
      stats_1$beta * 1000, p_text_1,
      resid_sd_pyr * 1000, n_points_total, n_segments, n_patients,
      stats_5$beta, p_text_5,
      n_segments, n_patients
    ),
    abbreviations = "b, linear intercept; LPR, lactate/pyruvate ratio; m, linear gradient; SD, standard deviation"
  )
)

writeLines(build_notes(legends, title = "8_substrate_lpr_fidelity"), file.path(OUTPUT_DIR, "0_notes.txt"))
