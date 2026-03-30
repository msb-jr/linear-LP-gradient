# ==============================================================================
# Script: 7_never_always_LPR25.R
# Manuscript relevance: 3.4, Table S3
# ==============================================================================
# PURPOSE:
#   Test whether persistently elevated LPR patients exhibit distinct linear LP
#   dynamics (steeper gradients, higher intercepts, more Type Pb segments) via
#   comparison of two patient cohorts: 
#   - "Always LPR > 25": Patients whose minimum LPR exceeds 25 for all datapoints
#   - "Never LPR > 25": Patients whose maximum LPR never exceeds 25
#   Differences are expected due to LPR = m + b/P (from L = mP + b). 
#   Aims to elucidate the polarization of patients' LPR about the threshold of 25,
#   as demonstrated in our group's previously published work:
#   doi: 10.1177/0271678X211042112 ;  doi: 10.1371/journal.pone.0331310
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#
# OUTPUT:
#   - 1_output/7_never_always_LPR25/__execution_time.csv: Runtime log
#   - 1_output/7_never_always_LPR25/0_notes.txt: Table and figure legends
#   - 1_output/7_never_always_LPR25/1_stats.csv: Cohort comparison table (Analysis 1)
#   - 1_output/7_never_always_LPR25/2_figure.png: Cohort comparison figure (Analysis 2)
#
# DATA FILTERS:
#   - Cohort assignment uses ALL datapoints before discarding segments with m ≤ 0
#   - Statistical analysis: Retained segments only (p6e_m > 0), except for (1a)
#
# ANALYSES:
# ------------------------------------------------------------------------------
# (1) Linear Parameter Comparison (1_stats.csv)
#
#     (1a) Pre-filter Percent Type ZNm (m ≤ 0)
#          - Report per-patient % of datapoints in Type Zm (m = 0) or Nm (m < 0)
#            segments before filtering out Zm and Nm (ZNm) segments
#          - Similar % supports comparability of Always vs Never LPR > 25 cohorts on 
#            the basis of retained segments
#          Approach:   Patient % → Cohort median (IQR)
#
#     (1b) Linear Gradient (m) Comparison
#          - Compare segment m between Always/Never cohorts
#          Approach:   Segment → Patient (weighted mean, w = n_points) → Cohort median (IQR)
#          Hypothesis: Always > Never (one-sided)
#          Test:       Wilcoxon rank-sum (Mann-Whitney U)
#          MCC:        Benjamini-Hochberg across 3 tests (m, b, % Type Pb)
#
#     (1c) Linear Intercept (b) Comparison
#          - Analogous to 1b above
#
#     (1d) Percent Type Pb (b > 0) Comparison
#          - Compare per-patient % of datapoints in Type Pb (b > 0) segments
#          Approach:   Patient % → Cohort median (IQR)
#          Hypothesis: Always > Never (one-sided)
#          Test:       Wilcoxon rank-sum
#          MCC:        BH across 3 tests (m, b, % Type Pb)
#
# (2) Cohort Visualization (2_figure.png)
#     - Visualize hyperbolic LPR~pyruvate relationship (from L = mP +b) by cohort
#     Panels:     (A) "Always LPR > 25", (B) "Never LPR > 25"
#     Elements:   - Main curve: Cohort median of patient-level weighted means
#                 - Spaghetti line opacity ∝ segment weight
#                 - Colors: Pb (teal), Nb (gold), Zb (gray)
#                 - References: y = 0 (red dotted), y = m_median (brown dashed);
#                   may not be visible depending on viewport
#     Viewport:   95% quantile cropping (2.5th–97.5th percentile) per cohort
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(colorspace)
  library(patchwork)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 7_never_always_LPR25.R\n")

# ==============================================================================
# Input / Output
# ==============================================================================
start_time <- Sys.time()

input_file <- here::here("1_output", "2_linear_segmentation", "1_results.csv")

if (!file.exists(input_file)) {
  cat("\n\nError: Results file not found. Run 2_linear_segmentation.py first.\n")
  stop()
}

output_dir <- here::here("1_output", "7_never_always_LPR25")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

safe_divide <- function(numerator, denominator, tol = 1e-12, fill = NA_real_) {
  result <- numerator / denominator
  invalid <- !is.finite(result) | !is.finite(denominator) | abs(denominator) <= tol
  result[invalid] <- fill
  result
}

# ==============================================================================
# Data Loading and Pre-processing
# ==============================================================================
df_raw <- read.csv(input_file)

# Identify patient cohorts by summarizing ALL LPR values per patient
# Rationale: Cohort definitions capture the observed LPR range
patient_lpr_summary <- df_raw %>%
  group_by(patient) %>%
  summarise(
    max_lpr = max(lpr, na.rm = TRUE),
    min_lpr = min(lpr, na.rm = TRUE),
    .groups = 'drop'
  )

# Filter to get the patient IDs for each cohort (based on all datapoints)
always_high_patients <- patient_lpr_summary %>%
  filter(min_lpr > 25) %>%
  pull(patient)

never_high_patients <- patient_lpr_summary %>%
  filter(max_lpr <= 25) %>%
  pull(patient)

if (length(always_high_patients) == 0 || length(never_high_patients) == 0) {
  cat("\n\nError: One or both cohorts are empty after applying the LPR > 25 thresholds.\n")
  stop()
}

# Prepare base dataframe for downstream analyses using RETAINED segments only
# Retained segments = segments with p6e gradient m > 0 (strictly positive).
# Segments with m ≤ 0 are discarded as they lack the expected LDH-driven LP relationship.
df_prepared <- df_raw %>%
  dplyr::filter(
    is.finite(p6e_m),
    p6e_m > 0  # Retained segments only
  ) %>%
  mutate(
    patient = factor(patient),
    segment = factor(paste(patient, p6e_seg_index, sep = "_")),
    # Stratify by intercept:
    # - neg: b < 0 (Type Nb)
    # - pos: b > 0 (Type Pb)
    # - zero: b = 0 (Type Zb)
    grp = factor(case_when(
      p6e_b < 0 ~ "neg",
      p6e_b > 0 ~ "pos",
      TRUE ~ "zero"
    ))
  )

# ==============================================================================
# Aggregation and Hyperbolic Visualization
# ==============================================================================

# Prepare data for each cohort (segments nested within patients)
df_always_high <- df_prepared %>% filter(patient %in% always_high_patients)
df_never_high  <- df_prepared %>% filter(patient %in% never_high_patients)

# Also prepare pre-filter dataframes for ZNm calculation
df_raw_always <- df_raw %>% filter(patient %in% always_high_patients)
df_raw_never <- df_raw %>% filter(patient %in% never_high_patients)

# --- % Type Pb: cohort median (IQR) of patient-level values ---
# Calculate per-patient proportion of datapoints in Type Pb segments (b > 0).
# Summarize as cohort median (IQR); displayed as percent.
patient_props_always <- df_always_high %>%
  group_by(patient) %>%
  summarise(prop_type_pb = mean(p6e_b > 0, na.rm = TRUE), .groups = 'drop')

patient_props_never <- df_never_high %>%
  group_by(patient) %>%
  summarise(prop_type_pb = mean(p6e_b > 0, na.rm = TRUE), .groups = 'drop')

# Cohort median (IQR) of patient-level proportions (converted to percent at display)
summary_prop_always <- tibble(
  med = median(patient_props_always$prop_type_pb, na.rm = TRUE),
  q1 = quantile(patient_props_always$prop_type_pb, 0.25, na.rm = TRUE),
  q3 = quantile(patient_props_always$prop_type_pb, 0.75, na.rm = TRUE)
)

summary_prop_never <- tibble(
  med = median(patient_props_never$prop_type_pb, na.rm = TRUE),
  q1 = quantile(patient_props_never$prop_type_pb, 0.25, na.rm = TRUE),
  q3 = quantile(patient_props_never$prop_type_pb, 0.75, na.rm = TRUE)
)

# Wilcoxon rank-sum test comparing patient-level values between cohorts
# Hypothesis: Always > Never (one-sided) for m, b, and % Type Pb
get_wilcox_p <- function(x, y, alternative = "two.sided") {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative = alternative, exact = FALSE)$p.value)
}
p_prop_wilcox <- get_wilcox_p(patient_props_always$prop_type_pb, patient_props_never$prop_type_pb, alternative = "greater")

# Per-patient percent of datapoints in Type ZNm segments (m ≤ 0), pre-filter
# Similar rates support comparability: filtering affects both cohorts equally
patient_pct_znm_always <- df_raw_always %>%
  group_by(patient) %>%
  summarise(pct_znm = mean(p6e_m <= 0, na.rm = TRUE) * 100, .groups = 'drop')

patient_pct_znm_never <- df_raw_never %>%
  group_by(patient) %>%
  summarise(pct_znm = mean(p6e_m <= 0, na.rm = TRUE) * 100, .groups = 'drop')

summary_znm_always <- tibble(
  med = median(patient_pct_znm_always$pct_znm, na.rm = TRUE),
  q1 = quantile(patient_pct_znm_always$pct_znm, 0.25, na.rm = TRUE),
  q3 = quantile(patient_pct_znm_always$pct_znm, 0.75, na.rm = TRUE)
)

summary_znm_never <- tibble(
  med = median(patient_pct_znm_never$pct_znm, na.rm = TRUE),
  q1 = quantile(patient_pct_znm_never$pct_znm, 0.25, na.rm = TRUE),
  q3 = quantile(patient_pct_znm_never$pct_znm, 0.75, na.rm = TRUE)
)

# Helper: compute patient-level weighted averages (weighted by segment size)
# Rationale: Coefficients should reflect the population of datapoints in each cohort.
# Segments with more datapoints contribute more to each patient's typical coefficient,
# as they better correspond to that patient's observed LPR.
compute_patient_level_weighted <- function(df) {
  # First get unique m,b per segment and count datapoints
  segment_stats <- df %>%
    group_by(patient, segment) %>%
    summarise(
      n_points = n(),
      seg_m = first(p6e_m),  # Same for all points in segment
      seg_b = first(p6e_b),
      .groups = 'drop'
    )
  
  # Then compute weighted means by patient
  segment_stats %>%
    group_by(patient) %>%
    summarise(
      wmean_b = safe_divide(sum(n_points * seg_b), sum(n_points)),
      wmean_m = safe_divide(sum(n_points * seg_m), sum(n_points)),
      .groups = 'drop'
    )
}

# Helper: summarize across patients (cohort median/IQR of patient-level weighted means)
summarize_patient_wmeans <- function(patient_wmeans_df) {
  tibble(
    q1_b = quantile(patient_wmeans_df$wmean_b, 0.25, na.rm = TRUE),
    cohort_median_b = median(patient_wmeans_df$wmean_b, na.rm = TRUE),
    q3_b = quantile(patient_wmeans_df$wmean_b, 0.75, na.rm = TRUE),
    q1_m = quantile(patient_wmeans_df$wmean_m, 0.25, na.rm = TRUE),
    cohort_median_m = median(patient_wmeans_df$wmean_m, na.rm = TRUE),
    q3_m = quantile(patient_wmeans_df$wmean_m, 0.75, na.rm = TRUE)
  )
}

patient_wmeans_always <- compute_patient_level_weighted(df_always_high)
patient_wmeans_never  <- compute_patient_level_weighted(df_never_high)

summary_always <- summarize_patient_wmeans(patient_wmeans_always)
summary_never  <- summarize_patient_wmeans(patient_wmeans_never)

# --- Create Hyperbolic Plots (cohort median of patient wmeans fit + spaghetti) ---

# Nonparametric comparisons (Wilcoxon rank-sum) of patient-level weighted means between cohorts
# Hypothesis: "Always LPR > 25" patients have HIGHER m and b than "Never LPR > 25" patients
pb_wilcox_raw <- get_wilcox_p(patient_wmeans_always$wmean_b, patient_wmeans_never$wmean_b, alternative = "greater")
p_m_wilcox_raw <- get_wilcox_p(patient_wmeans_always$wmean_m, patient_wmeans_never$wmean_m, alternative = "greater")

# Apply Benjamini-Hochberg (FDR) correction for multiple testing (3 comparisons: m, b, and % Type Pb)
p_values_raw <- c(p_m_wilcox_raw, pb_wilcox_raw, p_prop_wilcox)
p_values_adjusted <- p.adjust(p_values_raw, method = "BH")
p_m_wilcox <- p_values_adjusted[1]
pb_wilcox <- p_values_adjusted[2]
p_prop_wilcox_adj <- p_values_adjusted[3]

generate_hyperbolic_plot_data <- function(df, patient_wmeans_df) {
  # Cohort median of patient-level weighted means
  m_med <- median(patient_wmeans_df$wmean_m, na.rm = TRUE)
  b_med <- median(patient_wmeans_df$wmean_b, na.rm = TRUE)
  
  # Segment-level weights for alpha scaling (share of patient datapoints)
  segment_weights <- df %>%
    group_by(patient, segment) %>%
    summarise(n_points = n(), .groups = "drop") %>%
    group_by(patient) %>%
    mutate(segment_alpha = safe_divide(n_points, sum(n_points))) %>%
    ungroup()
  
  # Population grid for fit (in mM for calculations)
  grid <- data.frame(
    pyruvate = seq(min(df$pyruvate, na.rm = TRUE), max(df$pyruvate, na.rm = TRUE), length.out = 200)
  )
  grid$fit_lpr <- m_med + safe_divide(b_med, grid$pyruvate)
  # Convert pyruvate to μM for display (1 mM = 1000 μM)
  grid$pyruvate_um <- grid$pyruvate * 1000
  
  # Spaghetti data: segment-specific ranges only (in mM for calculations)
  spaghetti_df <- do.call(rbind, lapply(unique(df$segment), function(seg_id) {
    seg_data <- df[df$segment == seg_id, , drop = FALSE]
    seg_grid <- data.frame(pyruvate = seq(min(seg_data$pyruvate, na.rm = TRUE), max(seg_data$pyruvate, na.rm = TRUE), length.out = 100))
    seg_grid$segment <- seg_id
    # Use each segment's p6e coefficients for hyperbolic prediction
    seg_coeffs <- seg_data %>% distinct(segment, .keep_all = TRUE)
    seg_grid$y2_pred <- seg_coeffs$p6e_m + safe_divide(seg_coeffs$p6e_b, seg_grid$pyruvate)
    # Convert pyruvate to μM for display (1 mM = 1000 μM)
    seg_grid$pyruvate_um <- seg_grid$pyruvate * 1000
    seg_grid
  }))
  spaghetti_df <- spaghetti_df %>%
    dplyr::left_join(segment_weights %>% select(segment, segment_alpha), by = "segment") %>%
    dplyr::left_join(dplyr::distinct(df, segment, grp), by = "segment")
  
  list(grid = grid, spaghetti = spaghetti_df, m_med = m_med)
}

plot_data_always <- generate_hyperbolic_plot_data(df_always_high, patient_wmeans_always)
plot_data_never  <- generate_hyperbolic_plot_data(df_never_high, patient_wmeans_never)

# --- Compute 95% viewport limits (middle 95% of data) for each cohort ---
# Viewports use 2.5th–97.5th percentile cropping to focus on typical data.

# Always cohort 95% viewport
viewport_always <- list(
  xlim = quantile(df_always_high$pyruvate, c(0.025, 0.975), na.rm = TRUE) * 1000,  # μM
  ylim = quantile(df_always_high$lpr, c(0.025, 0.975), na.rm = TRUE)
)

# Never cohort 95% viewport
viewport_never <- list(
  xlim = quantile(df_never_high$pyruvate, c(0.025, 0.975), na.rm = TRUE) * 1000,  # μM
  ylim = quantile(df_never_high$lpr, c(0.025, 0.975), na.rm = TRUE)
)

brown_color <- "#996633"

p_hyper_always <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.6) +
  geom_hline(yintercept = plot_data_always$m_med, linetype = "dashed", color = darken(brown_color, 0.2), linewidth = 0.6, alpha = 0.8) +
  geom_line(
    data = plot_data_always$spaghetti,
    aes(x = pyruvate_um, y = y2_pred, group = segment, color = grp, alpha = segment_alpha),
    linewidth = 0.3
  ) +
  geom_line(data = plot_data_always$grid, aes(x = pyruvate_um, y = fit_lpr), colour = darken(brown_color, 0.2), linewidth = 0.8) +
  coord_cartesian(xlim = viewport_always$xlim, ylim = viewport_always$ylim) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  labs(tag = "A", title = NULL, x = "Pyruvate (μM)", y = NULL) +
  scale_color_manual(values = c(pos = "#44AA99", neg = "#DDCC77", zero = "#888888")) +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
    plot.tag.position = c(0.02, 0.99),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 3, r = 3, b = 0, l = 3, unit = "mm")
  )

p_hyper_never <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.6) +
  geom_hline(yintercept = plot_data_never$m_med, linetype = "dashed", color = darken(brown_color, 0.2), linewidth = 0.6, alpha = 0.8) +
  geom_line(
    data = plot_data_never$spaghetti,
    aes(x = pyruvate_um, y = y2_pred, group = segment, color = grp, alpha = segment_alpha),
    linewidth = 0.3
  ) +
  geom_line(data = plot_data_never$grid, aes(x = pyruvate_um, y = fit_lpr), colour = darken(brown_color, 0.2), linewidth = 0.8) +
  coord_cartesian(xlim = viewport_never$xlim, ylim = viewport_never$ylim) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
  labs(tag = "B", title = NULL, x = "Pyruvate (μM)", y = NULL) +
  scale_color_manual(values = c(pos = "#44AA99", neg = "#DDCC77", zero = "#888888")) +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
    plot.tag.position = c(0.02, 0.99),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 3, r = 3, b = 0, l = 3, unit = "mm")
  )

combined_hyper_fig <- (p_hyper_always | p_hyper_never) + 
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(t = 3, r = 2, b = 0, l = 2, unit = "mm"))

# Add shared y-axis title for LPR
combined_hyper_fig <- wrap_elements(combined_hyper_fig) +
  labs(tag = "LPR") +
  theme(
    plot.tag = element_text(size = 7, colour = "black", angle = 90),
    plot.tag.position = "left"
  )

figure_file <- file.path(output_dir, "2_figure.png")
if (isTRUE(capabilities("cairo"))) {
  ggsave(
    filename = figure_file,
    plot = combined_hyper_fig,
    width = 90, height = 55, units = "mm", dpi = 300, type = "cairo"
  )
} else {
  warning("Cairo backend unavailable; saving 2_figure.png with the default device.")
  ggsave(
    filename = figure_file,
    plot = combined_hyper_fig,
    width = 90, height = 55, units = "mm", dpi = 300
  )
}

# ==============================================================================
# Create Summary Table
# ==============================================================================

# Construct cohort columns
always_col_label <- "Always LPR > 25"
never_col_label <- "Never LPR > 25"

# Strings for linear coefficients: cohort median (Q1, Q3) of patient-level weighted means
always_b_str <- sprintf("%.2f (%.2f, %.2f)", summary_always$cohort_median_b, summary_always$q1_b, summary_always$q3_b)
always_m_str <- sprintf("%.2f (%.2f, %.2f)", summary_always$cohort_median_m, summary_always$q1_m, summary_always$q3_m)
never_b_str  <- sprintf("%.2f (%.2f, %.2f)", summary_never$cohort_median_b, summary_never$q1_b, summary_never$q3_b)
never_m_str  <- sprintf("%.2f (%.2f, %.2f)", summary_never$cohort_median_m, summary_never$q1_m, summary_never$q3_m)

# Strings for percent of Type Pb datapoints per patient
always_prop_str <- sprintf("%.1f%% (%.1f%%, %.1f%%)", summary_prop_always$med * 100, summary_prop_always$q1 * 100, summary_prop_always$q3 * 100)
never_prop_str <- sprintf("%.1f%% (%.1f%%, %.1f%%)", summary_prop_never$med * 100, summary_prop_never$q1 * 100, summary_prop_never$q3 * 100)

# Strings for pre-filter percent of Type ZNm datapoints per patient
always_znm_str <- sprintf("%.1f%% (%.1f%%, %.1f%%)", summary_znm_always$med, summary_znm_always$q1, summary_znm_always$q3)
never_znm_str <- sprintf("%.1f%% (%.1f%%, %.1f%%)", summary_znm_never$med, summary_znm_never$q1, summary_znm_never$q3)

# Calculate effect sizes (difference in cohort medians: Always - Never)
delta_m <- summary_always$cohort_median_m - summary_never$cohort_median_m
delta_b <- summary_always$cohort_median_b - summary_never$cohort_median_b
delta_prop <- (summary_prop_always$med - summary_prop_never$med) * 100  # Convert to percentage points
delta_znm <- summary_znm_always$med - summary_znm_never$med  # Already in percentage points

# Create the summary tibble
summary_data <- tibble(
  parameter = c(
    "Pre-filter % Type ZNm [c]",
    "Linear Gradient (m) [d]", "Linear Intercept (b, mM) [d]",
    "% Type Pb"
  )
)
summary_data[[always_col_label]] <- c(
  always_znm_str,
  always_m_str,
  always_b_str,
  always_prop_str
)
summary_data[[never_col_label]] <- c(
  never_znm_str,
  never_m_str,
  never_b_str,
  never_prop_str
)
summary_data[["Δ [a]"]] <- c(
  sprintf("%.1f%%", delta_znm),
  sprintf("%.2f", delta_m),
  sprintf("%.2f", delta_b),
  sprintf("%.1f%%", delta_prop)
)
summary_data[["P-value [b]"]] <- c(
  "—",
  ifelse(p_m_wilcox < 0.001, "< 0.001", sprintf("%.3f", p_m_wilcox)),
  ifelse(pb_wilcox < 0.001, "< 0.001", sprintf("%.3f", pb_wilcox)),
  ifelse(p_prop_wilcox_adj < 0.001, "< 0.001", sprintf("%.3f", p_prop_wilcox_adj))
)

# Save raw summary data to CSV
write.csv(summary_data, file.path(output_dir, "1_stats.csv"), row.names = FALSE)

# ==============================================================================
# Demographic Statistics
# ==============================================================================
compute_demo_stats <- function(df_cohort) {
  n_patients <- dplyr::n_distinct(df_cohort$patient)
  n_segments <- dplyr::n_distinct(df_cohort$segment)
  
  spp <- df_cohort %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(n_segments = dplyr::n_distinct(segment), .groups = "drop")
  
  if (nrow(spp) > 0) {
    med_spp <- stats::median(spp$n_segments, na.rm = TRUE)
    q1_spp <- stats::quantile(spp$n_segments, 0.25, na.rm = TRUE)
    q3_spp <- stats::quantile(spp$n_segments, 0.75, na.rm = TRUE)
  } else {
    med_spp <- NA_real_
    q1_spp <- NA_real_
    q3_spp <- NA_real_
  }
  
  list(
    n_patients = n_patients,
    n_segments = n_segments,
    med_spp = med_spp,
    q1_spp = q1_spp,
    q3_spp = q3_spp
  )
}

# Compute stats for each cohort
demo_always <- compute_demo_stats(df_always_high)
demo_never <- compute_demo_stats(df_never_high)

# Pre-filter cohort sizes (for legend)
n_patients_always_all <- dplyr::n_distinct(df_raw_always$patient)
n_segments_always_all <- dplyr::n_distinct(paste(df_raw_always$patient, df_raw_always$p6e_seg_index))
n_patients_never_all <- dplyr::n_distinct(df_raw_never$patient)
n_segments_never_all <- dplyr::n_distinct(paste(df_raw_never$patient, df_raw_never$p6e_seg_index))

# ==============================================================================
# Log Execution Time
# ==============================================================================
end_time <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
time_df <- data.frame(
  execution_time_seconds = execution_time_seconds
)
write.csv(
  time_df,
  file.path(output_dir, "__execution_time.csv"),
  row.names = FALSE
)

# Build simple legends for table and figure
legends <- list(
  list(
    target = "Table 1 (1_stats.csv)",
    caption = sprintf(
      "Comparison of linear LP parameters between 'Always LPR > 25' and 'Never LPR > 25' cohorts. Values are cohort median (IQR) of patient-level values. Cohort sizes: Always pre-filter %d patients, %d segments; post-filter %d patients, %d segments. Never pre-filter %d patients, %d segments; post-filter %d patients, %d segments.",
      n_patients_always_all, n_segments_always_all, demo_always$n_patients, demo_always$n_segments,
      n_patients_never_all, n_segments_never_all, demo_never$n_patients, demo_never$n_segments
    ),
    footnotes = c(
      "[a] Δ = Always − Never.",
      "[b] P-values from one-sided Wilcoxon tests (H1: Always > Never) with BH correction.",
      "[c] Pre-filter only: per-patient percent of datapoints with m ≤ 0 (zero or negative gradient) before filtering. All other rows use filtered data (m > 0).",
      "[d] Patient-level coefficients weight segments by datapoint counts."
    ),
    abbreviations = "BH, Benjamini-Hochberg; IQR, interquartile range; LPR, lactate/pyruvate ratio; ZNm, zero or negative gradient (m ≤ 0); Pb, positive intercept (b > 0)"
  ),
  list(
    target = "Figure 1 (2_figure.png)",
    caption = sprintf(
      "Linear LP-estimated hyperbolic LPR-pyruvate relationships in cohorts defined by observed LPR. (A) 'Always LPR > 25' cohort (n = %d): patients whose minimum LPR exceeded 25 for all samples. (B) 'Never LPR > 25' cohort (n = %d): patients whose maximum LPR never exceeded 25. Spaghetti lines represent individual segments within patients (teal = Type Pb, gold = Type Nb, gray = Type Zb); opacity reflects segment weight. Brown curves represent the cohort median trend; brown dashed lines = cohort median gradient (m). The red dotted line marks the boundary of nonphysiological extrapolation (negative LPR).",
      demo_always$n_patients, demo_never$n_patients
    ),
    abbreviations = "LPR, lactate/pyruvate ratio; Nb, negative intercept; Pb, positive intercept; Zb, zero intercept"
  )
)

writeLines(build_notes(legends, title = "7_never_always_LPR25"), file.path(output_dir, "0_notes.txt"))
cat("\n\n7_never_always_LPR25.R complete.\n\n")
