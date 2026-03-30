# ==============================================================================
# Script: 9_outcome.R
# Manuscript relevance: 3.6, Fig. 6
# ==============================================================================
# PURPOSE:
#   Test whether segment-level LP dynamics parameters differ systematically
#   between patients with different functional outcomes (GOSE). Examines
#   whether worse outcomes associate with higher LPR, steeper gradients (m),
#   or different intercepts (b).
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#
# OUTPUT:
#   - 1_output/9_outcome/__execution_time.csv: Runtime log
#   - 1_output/9_outcome/_1_emm_by_gose.csv: Estimated marginal means by GOSE (Analysis 1)
#   - 1_output/9_outcome/_1_model_summaries.txt: Full lmer summaries (Analysis 1)
#   - 1_output/9_outcome/_1_stats.csv: Model contrast statistics (Analysis 1)
#   - 1_output/9_outcome/0_notes.txt: Figure legend
#   - 1_output/9_outcome/2_figure.png: 3-panel EMM visualization (Analysis 2)
#
# DATA FILTERS:
#   - Retained segments only: p6e_m > 0 (positive gradient)
#   - GOSE: Valid values 1–8, collapsed to 4 levels:
#       Death (1), Severe (2–4), Moderate (5–6), Good (7–8)
#
# ANALYSES:
# ------------------------------------------------------------------------------
# (1) LME Models: Parameter ~ GOSE (_1_emm_by_gose.csv, _1_model_summaries.txt, _1_stats.csv)
#     Three LME models testing outcome associations:
#
#     (1a) LME: LPR ~ GOSE [VALIDATION]
#          - Test if segment LPR is higher in worse outcomes
#          Unit:       Segment-level (LPR = segment median of point-level values)
#          Model:      lpr ~ gose_factor + (1|patient)
#          Weighting:  n_points (observations per segment)
#          Hypothesis: Negative estimate = higher in worse outcomes (one-sided)
#          Test:       Linear contrast via emmeans (coefficients: -3, -1, +1, +3)
#          MCC:        None (independent hypothesis)
#
#     (1b) LME: m ~ GOSE [PRIMARY]
#          - Test if segment gradient is higher in worse outcomes
#          Unit:       Segment-level (m = single value per segment from linear fit)
#          Model:      m ~ gose_factor + (1|patient)
#          Weighting:  n_points
#          Hypothesis: Negative estimate = higher in worse outcomes (one-sided)
#          Test:       Linear contrast via emmeans
#          MCC:        None (primary hypothesis)
#
#     (1c) LME: b ~ GOSE [EXPLORATORY]
#          - Explore intercept–outcome relationship
#          Unit:       Segment-level (b = single value per segment from linear fit)
#          Model:      b ~ gose_factor + (1|patient)
#          Weighting:  n_points
#          Hypothesis: Two-sided (exploratory)
#          Test:       Linear contrast via emmeans
#          MCC:        None (exploratory)
#
# (2) Figure: EMMs by GOSE (2_figure.png)
#     - Visualize group means with 95% CIs
#     Panels:     (A) LPR, (B) m, (C) b
#     Elements:   Point estimates ± 95% CI by GOSE category
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 9_outcome.R\n")

options(dplyr.show_progress = FALSE, readr.show_col_types = FALSE, readr.show_progress = FALSE)
emm_options(lmer.df = "satterthwaite")  # Explicit: emmeans uses Satterthwaite (not KR)

# ==============================================================================
# Input / Output
# ==============================================================================
INPUT_FILE <- here::here("1_output", "2_linear_segmentation", "1_results.csv")

if (!file.exists(INPUT_FILE)) {
  cat("\n\nError: Results file not found. Run 2_linear_segmentation.py first.\n")
  stop()
}

OUTPUT_DIR <- here::here("1_output", "9_outcome")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

save_plot_portable <- function(filename, plot, ...) {
  if (isTRUE(capabilities("cairo"))) {
    ggsave(filename = filename, plot = plot, type = "cairo", ...)
  } else {
    warning(sprintf("Cairo backend unavailable; saving %s with the default device.", basename(filename)))
    ggsave(filename = filename, plot = plot, ...)
  }
}

# EMM plot by GOSE category (point + CI, trimmed to union of CIs with 5% padding)
plot_emm_by_gose <- function(emm_df, y_lab, colors) {
  # Calculate union of CIs for y-axis limits with 5% padding
  y_min <- min(emm_df$lower.CL)
  y_max <- max(emm_df$upper.CL)
  y_range <- y_max - y_min
  y_padding <- y_range * 0.05
  
  ggplot(emm_df, aes(x = gose_factor, y = emmean, color = gose_factor)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                  width = 0.15, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    coord_cartesian(ylim = c(y_min - y_padding, y_max + y_padding)) +
    labs(x = "GOSE Category", y = y_lab) +
    theme_minimal() +
    theme(
      text = element_text(size = 8, colour = "black"),
      axis.text = element_text(size = 6, colour = "black"),
      axis.title = element_text(size = 7, colour = "black"),
      plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
      plot.tag.position = c(0.02, 0.99),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 3, r = 3, b = 0, l = 3, unit = "mm")
    )
}

# ==============================================================================
# Data Loading and Preparation
# ==============================================================================
start_time <- Sys.time()

# Filter: Retain only segments with positive gradient (p6e_m > 0)
df_retained <- readr::read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  filter(is.finite(p6e_m), p6e_m > 0) %>%
  as_tibble()

if (nrow(df_retained) == 0) {
  cat("\n\nError: No retained segment data found.\n")
  stop()
}

# Aggregate to Segment Level and Prepare GOSE Factors
segment_stats_gose <- df_retained %>%
  group_by(patient, p6e_seg_index) %>%
  summarise(
    seg_median_LPR = median(lpr, na.rm = TRUE),
    seg_m = first(p6e_m),
    seg_b = first(p6e_b),
    gose = first(gose),
    n_points = n(),
    .groups = "drop"
  ) %>%
  arrange(patient, p6e_seg_index) %>%
  filter(!is.na(gose), gose >= 1, gose <= 8) %>%
  mutate(
    # Collapse 8-level GOSE into 4 clinical categories
    gose_factor = factor(
      case_when(
        gose == 1 ~ "Death",
        gose >= 2 & gose <= 4 ~ "Severe",
        gose >= 5 & gose <= 6 ~ "Moderate",
        gose >= 7 & gose <= 8 ~ "Good"
      ),
      levels = c("Death", "Severe", "Moderate", "Good")
    )
  )

# Sample sizes by GOSE category (patient level)
patient_stats <- segment_stats_gose %>%
  group_by(patient) %>%
  summarise(
    gose_category = first(gose_factor),
    .groups = "drop"
  )

n_patients <- nrow(patient_stats)
n_segments <- nrow(segment_stats_gose)
gose_counts <- patient_stats %>% count(gose_category) %>% pull(n, name = gose_category)
get_gose_count <- function(vec, key) {
  if (is.null(vec) || length(vec) == 0 || !(key %in% names(vec))) {
    return(0L)
  }
  val <- vec[[key]]
  if (is.na(val)) 0L else as.integer(val)
}
n_death <- get_gose_count(gose_counts, "Death")
n_severe <- get_gose_count(gose_counts, "Severe")
n_moderate <- get_gose_count(gose_counts, "Moderate")
n_good <- get_gose_count(gose_counts, "Good")

# Segment counts by GOSE category
seg_death <- sum(segment_stats_gose$gose_factor == "Death")
seg_severe <- sum(segment_stats_gose$gose_factor == "Severe")
seg_moderate <- sum(segment_stats_gose$gose_factor == "Moderate")
seg_good <- sum(segment_stats_gose$gose_factor == "Good")

# Colors for 4-level GOSE (gradient from worst to best)
colors_gose <- c(
  "Death" = "#D55E00",     # Red-orange
  "Severe" = "#E69F00",    # Orange
  "Moderate" = "#56B4E9",  # Light blue
  "Good" = "#009E73"       # Green
)

# ==============================================================================
# LME Models: Parameter ~ GOSE
# ==============================================================================

lme_m <- lmer(seg_m ~ gose_factor + (1|patient), data = segment_stats_gose,
              weights = n_points,
              control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
lme_lpr <- lmer(seg_median_LPR ~ gose_factor + (1|patient), data = segment_stats_gose,
                weights = n_points,
                control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
lme_b <- lmer(seg_b ~ gose_factor + (1|patient), data = segment_stats_gose,
              weights = n_points,
              control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))

# Extract EMMs (actual group means) from each model
emm_m_by_gose <- as.data.frame(emmeans(lme_m, ~ gose_factor))
emm_lpr_by_gose <- as.data.frame(emmeans(lme_lpr, ~ gose_factor))
emm_b_by_gose <- as.data.frame(emmeans(lme_b, ~ gose_factor))

# Linear contrast: tests if means show a linear trend from Death → Good
# Coefficients [-3, -1, +1, +3] are the standard linear contrast for 4 equally-spaced levels
# They sum to 0 and are equally spaced, which tests for a linear (monotonic) trend
# Negative estimate = means DECREASE from Death → Good = parameter HIGHER in worse outcomes
linear_contrast <- list("linear_trend" = c(-3, -1, 1, 3))

# Extract linear trend p-value from emmeans contrast
extract_linear_contrast <- function(model, one_sided_negative = TRUE) {
  emm <- emmeans(model, ~ gose_factor)
  contr <- contrast(emm, linear_contrast)
  contr_summary <- summary(contr)
  
  estimate <- contr_summary$estimate
  se <- contr_summary$SE
  t_val <- contr_summary$t.ratio
  p_two <- contr_summary$p.value
  
  if (one_sided_negative) {
    # One-sided: H1 is estimate < 0 (worse outcome → higher parameter)
    # With our contrast coding, negative estimate means Death > Good
    p_val <- ifelse(t_val < 0, p_two / 2, 1 - p_two / 2)
  } else {
    p_val <- p_two
  }
  
  list(estimate = estimate, se = se, t_value = t_val, p_value = p_val)
}

trend_stats_m <- extract_linear_contrast(lme_m, one_sided_negative = TRUE)
trend_stats_lpr <- extract_linear_contrast(lme_lpr, one_sided_negative = TRUE)
trend_stats_b <- extract_linear_contrast(lme_b, one_sided_negative = FALSE)  # Exploratory

# Save outcome results (linear contrast test)
outcome_results <- tibble(
  Variable = c("LPR", "m", "b"),
  Role = c("Validation", "Primary", "Exploratory"),
  Linear_contrast_estimate = c(trend_stats_lpr$estimate, trend_stats_m$estimate, trend_stats_b$estimate),
  SE = c(trend_stats_lpr$se, trend_stats_m$se, trend_stats_b$se),
  t_value = c(trend_stats_lpr$t_value, trend_stats_m$t_value, trend_stats_b$t_value),
  p_value = c(trend_stats_lpr$p_value, trend_stats_m$p_value, trend_stats_b$p_value),
  p_type = c("one-sided", "one-sided", "two-sided"),
  interpretation = c(
    "Negative estimate: higher LPR in worse outcome groups",
    "Negative estimate: higher m in worse outcome groups", 
    "Two-sided test for b"
  )
)
write_csv(outcome_results, file.path(OUTPUT_DIR, "_1_stats.csv"))

# Save EMMs with confidence intervals
emm_results <- bind_rows(
  emm_lpr_by_gose %>% mutate(variable = "LPR"),
  emm_m_by_gose %>% mutate(variable = "m"),
  emm_b_by_gose %>% mutate(variable = "b")
) %>%
  select(variable, gose_factor, emmean, SE, lower.CL, upper.CL)
write_csv(emm_results, file.path(OUTPUT_DIR, "_1_emm_by_gose.csv"))

# ==============================================================================
# Figure: EMMs by GOSE Category (Analysis 2)
# ==============================================================================

# Panel A: EMM of LPR by GOSE category
p_lpr_gose <- plot_emm_by_gose(emm_lpr_by_gose, "EMM LPR", colors_gose)

# Panel B: EMM of m by GOSE category
p_m_gose <- plot_emm_by_gose(emm_m_by_gose, 
                              expression(paste("EMM Gradient (", italic("m"), ")")),
                              colors_gose)

# Panel C: EMM of b by GOSE category
p_b_gose <- plot_emm_by_gose(emm_b_by_gose, 
                              expression(paste("EMM Intercept (", italic("b"), ", mM)")),
                              colors_gose)

# Assemble multi-panel figure
fig_outcome <- (p_lpr_gose | p_m_gose | p_b_gose) + 
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
    plot.tag.position = c(0.02, 0.99),
    legend.position = "none"
  )

save_plot_portable(
  filename = file.path(OUTPUT_DIR, "2_figure.png"),
  plot = fig_outcome,
  width = 185, height = 55, units = "mm", dpi = 300
)


# ==============================================================================
# Save Full Model Summaries
# ==============================================================================
sink(file.path(OUTPUT_DIR, "_1_model_summaries.txt"))
cat("=== LME Models: Parameter ~ GOSE (Outcome Association Analysis) ===\n\n")
cat("STATISTICAL APPROACH:\n")
cat("Segment-level LME with patient random intercepts — a standard approach for\n")
cat("hierarchical data where multiple segments are nested within patients.\n\n")
cat("Model: parameter ~ gose_factor + (1|patient), weighted by n_points\n")
cat("  - gose_factor is UNORDERED → each group gets its own mean (actual EMMs)\n")
cat("  - Patient random intercepts account for within-patient correlation\n")
cat("  - Linear contrast via emmeans tests monotonic trend across outcome groups\n")
cat("  - Negative estimate = higher parameter in worse outcome groups\n")
cat("  - LPR/m: one-sided tests; b: two-sided test (exploratory)\n\n")
cat("INTERPRETATION:\n")
cat("This tests ASSOCIATION: do segment-level parameters differ systematically\n")
cat("between patients with different functional outcomes?\n\n")
cat(sprintf("Sample sizes: %d patients (Death %d; Severe %d; Moderate %d; Good %d)\n",
            n_patients, n_death, n_severe, n_moderate, n_good))
cat(sprintf("               %d segments (Death %d; Severe %d; Moderate %d; Good %d)\n\n",
            n_segments, seg_death, seg_severe, seg_moderate, seg_good))
cat("--- m ~ gose_factor [Primary] ---\n")
print(summary(lme_m))
cat("\n--- Linear contrast for m ---\n")
print(contrast(emmeans(lme_m, ~ gose_factor), linear_contrast))
cat("\n--- LPR ~ gose_factor [Validation] ---\n")
print(summary(lme_lpr))
cat("\n--- Linear contrast for LPR ---\n")
print(contrast(emmeans(lme_lpr, ~ gose_factor), linear_contrast))
cat("\n--- b ~ gose_factor [Exploratory] ---\n")
print(summary(lme_b))
cat("\n--- Linear contrast for b ---\n")
print(contrast(emmeans(lme_b, ~ gose_factor), linear_contrast))
cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("=== EMMs by GOSE Category (actual group means from same models) ===\n\n")
cat("--- m EMMs by GOSE ---\n")
print(emm_m_by_gose)
cat("\n--- LPR EMMs by GOSE ---\n")
print(emm_lpr_by_gose)
cat("\n--- b EMMs by GOSE ---\n")
print(emm_b_by_gose)
sink()

# ==============================================================================
# Completion Log
# ==============================================================================
end_time <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

log_df <- data.frame(execution_time_seconds = execution_time_seconds)
write_csv(log_df, file.path(OUTPUT_DIR, "__execution_time.csv"))

cat("\n\n9_outcome.R complete.\n\n")

# ==============================================================================
# Build figure legend (only output for notes)
# ==============================================================================

# Format p-values for legend
format_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
p_text_lpr <- format_p(trend_stats_lpr$p_value)
p_text_m <- format_p(trend_stats_m$p_value)
p_text_b <- format_p(trend_stats_b$p_value)

legends <- list(
  list(
    target = "Figure 1 (2_figure.png)",
    caption = sprintf(
      "Estimated marginal means of segment-level LP parameters by functional outcome. (A) LPR, (B) gradient (m), (C) intercept (b) across GOSE categories: Death (n=%d patients, %d segments), Severe (n=%d, %d), Moderate (n=%d, %d), Good (n=%d, %d). Points = EMMs from mixed-effects models; error bars = 95%% CI. Linear contrasts tested monotonic trends (coefficients: -3, -1, +1, +3): LPR β = %.2f (%s), m β = %.2f (%s), b β = %.2f (%s).",
      n_death, seg_death, n_severe, seg_severe, n_moderate, seg_moderate, n_good, seg_good,
      trend_stats_lpr$estimate, p_text_lpr, trend_stats_m$estimate, p_text_m, trend_stats_b$estimate, p_text_b
    ),
    abbreviations = "CI, confidence interval; EMM, estimated marginal mean; GOSE, Glasgow Outcome Scale-Extended; LPR, lactate/pyruvate ratio"
  )
)

writeLines(build_notes(legends, title = "9_outcome"), file.path(OUTPUT_DIR, "0_notes.txt"))
