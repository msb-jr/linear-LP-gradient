# ==============================================================================
# Script: 6_hyperbolic_LPR_error.R
# Manuscript relevance: 3.3.ii, Fig. 4, Table S2, Table S5
# ==============================================================================
# PURPOSE:
#   Quantify and visualize how the observed lactate/pyruvate ratio (LPR)
#   deviates from the linear gradient 'm' (hyperbolic asymptote) at different
#   pyruvate levels. Demonstrates implications of LPR = m + b/P (from L = mP + b),
#   in which the LPR deviates from m with changes in pyruvate when b ≠ 0.
#   Most substantial deviations at low pyruvate.
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#
# OUTPUT:
#   - 1_output/6_hyperbolic_LPR_error/__execution_time.csv: Runtime log
#   - 1_output/6_hyperbolic_LPR_error/0_notes.txt: Table and (general) figure legends
#   - 1_output/6_hyperbolic_LPR_error/1_stats.csv: Cohort medians/IQRs table (Analysis 1)
#   - 1_output/6_hyperbolic_LPR_error/2_exemplar/: Exemplar patient figs/legends (Analysis 2)
#   - 1_output/6_hyperbolic_LPR_error/3_all/: All patients (Analysis 3)
#
# DATA FILTERS (OVERVIEW):
#   - Retained segments only: p6e_m > 0 (positive gradient)
#   - Linear intercept stratification: Pb (b > 0), Nb (b < 0), Zb (b = 0)
#
# ANALYSES:
# ------------------------------------------------------------------------------
# (1) LPR Error Statistics (1_stats.csv)
#
#     (1a) Segment-Level LPR Error Metrics
#          For each segment compute:
#             - p_min, p_max: Pyruvate extremes
#             - rel_error_Pmin: 100 × (LPR_obs - m) / m at min pyruvate
#             - rel_error_Pmax: 100 × (LPR_obs - m) / m at max pyruvate
#             - LPR_range: max(lpr) - min(lpr)
#
#     (1b) Patient-Level Metrics (within each segment cohort)
#          Cohorts: unstratified, Type Pb only, Type Nb only.
#          For each cohort, group by patient and compute:
#            - Weighted mean (weights = segment n_points) of rel_error_Pmin/Pmax
#            - Patient medians of segment-level LPR_range, p_min, p_max
#
#     (1c) Cohort-Level Summaries (across patients in each cohort)
#            - Median (IQR) of the patient-level metrics from (1b)
#
# (2) Exemplar Patient Plots (2_exemplar/)
#     - Generate LPR error visualization for patients exemplifying theory
#     Selection:  ALL criteria must be met:
#                 (a) All segments have m > 0 (no discarded)
#                 (b) No Type Zb segments (only clear Pb/Nb)
#                 (c) Has both Type Pb AND Type Nb segments
#                 (d) All LPR values on 'correct' side of m per segment
#                 (e) |error| at Pmin ≥ 75th percentile (large error at low P)
#                 (f) |error| at Pmax ≤ 25th percentile (small error at high P)
#
# (3) All-Patient Plots (3_all/)
#     - Generate LPR error visualization for every patient
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(parallel)
  library(patchwork)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 6_hyperbolic_LPR_error.R\n")

options(dplyr.show_progress = FALSE)

# ==============================================================================
# Input / Output
# ==============================================================================
input_file <- here::here("1_output", "2_linear_segmentation", "1_results.csv")

if (!file.exists(input_file)) {
  cat("\n\nError: Results file not found. Run 2_linear_segmentation.py first.\n")
  stop()
}

output_dir <- here::here("1_output", "6_hyperbolic_LPR_error")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Helper Functions
# ==============================================================================
median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(list(q1 = NA_real_, median = NA_real_, q3 = NA_real_))
  }
  list(
    q1     = stats::quantile(x, 0.25, na.rm = TRUE, type = 7),
    median = stats::median(x, na.rm = TRUE),
    q3     = stats::quantile(x, 0.75, na.rm = TRUE, type = 7)
  )
}

format_interval <- function(med, q1, q3, digits = 1, suffix = "") {
  if (!all(is.finite(c(med, q1, q3)))) {
    return("N/A")
  }
  fmt <- paste0("%.", digits, "f")
  med_str <- sprintf(fmt, med)
  q1_str <- sprintf(fmt, q1)
  q3_str <- sprintf(fmt, q3)
  paste0(med_str, suffix, " (", q1_str, suffix, ", ", q3_str, suffix, ")")
}

safe_divide <- function(numerator, denominator, tol = 1e-12, fill = NA_real_) {
  result <- numerator / denominator
  invalid <- !is.finite(result) | !is.finite(denominator) | abs(denominator) <= tol
  result[invalid] <- fill
  result
}

# ==============================================================================
# Timing Start
# ==============================================================================
start_time <- Sys.time()

# ==============================================================================
# Data Loading and Pre-processing
# ==============================================================================
df_raw <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)

# Filter for retained segments (positive gradient) and valid data
df_prepared <- df_raw %>%
  dplyr::filter(
    is.finite(pyruvate),
    pyruvate > 0,
    is.finite(p6e_m),
    is.finite(p6e_b),
    p6e_m > 0  # Retained segments only
  ) %>%
  dplyr::mutate(
    patient = factor(patient),
    segment = factor(paste(patient, p6e_seg_index, sep = "_")),
    # Cohort type based on sign of p6e intercept
    seg_type = dplyr::case_when(
      p6e_b > 0 ~ "Pb",
      p6e_b < 0 ~ "Nb",
      TRUE      ~ "Zb"
    )
  )

if (nrow(df_prepared) == 0) {
  cat("\n\nError: No rows remain after filtering for retained segments (p6e_m > 0) and valid pyruvate/coefficients.\n")
  stop()
}

# For each segment, calculate LPR error metrics:
#   - Relative LPR error at Pmin and Pmax: 100 * (LPR_obs - m) / m
#     Positive values = LPR overestimates m; Negative = LPR underestimates m
#   - LPR range within segment
segment_metrics <- df_prepared %>%
  dplyr::group_by(patient, segment, seg_type) %>%
  dplyr::summarise(
    p6e_m = dplyr::first(p6e_m),
    p6e_b = dplyr::first(p6e_b),
    p6e_r = dplyr::first(p6e_r),
    n_points = n(),
    
    # Segment pyruvate extremes
    p_min = if (any(is.finite(pyruvate))) min(pyruvate, na.rm = TRUE) else NA_real_,
    p_max = if (any(is.finite(pyruvate))) max(pyruvate, na.rm = TRUE) else NA_real_,
    
    # Relative LPR Error at Pmin: 100 * (LPR_obs - m) / m
    rel_error_pmin = {
      m0 <- dplyr::first(p6e_m)
      
      if (is.finite(m0) && m0 != 0 && any(is.finite(pyruvate))) {
        idx_min <- which.min(pyruvate)
        lpr_val <- lpr[idx_min]
        
        if (is.finite(lpr_val)) safe_divide(100 * (lpr_val - m0), m0) else NA_real_
      } else {
        NA_real_
      }
    },
    
    # Relative LPR Error at Pmax: 100 * (LPR_obs - m) / m
    rel_error_pmax = {
      m0 <- dplyr::first(p6e_m)
      
      if (is.finite(m0) && m0 != 0 && any(is.finite(pyruvate))) {
        idx_max <- which.max(pyruvate)
        lpr_val <- lpr[idx_max]
        
        if (is.finite(lpr_val)) safe_divide(100 * (lpr_val - m0), m0) else NA_real_
      } else {
        NA_real_
      }
    },
    
    # LPR range within segment
    lpr_range = {
      if (any(is.finite(lpr))) {
        max(lpr, na.rm = TRUE) - min(lpr, na.rm = TRUE)
      } else {
        NA_real_
      }
    },
    
    .groups = "drop"
  )

# ==============================================================================
# Patient-Level Aggregation
# ==============================================================================
compute_patient_level <- function(seg_df) {
  if (nrow(seg_df) == 0) {
    return(dplyr::tibble(
      patient = character(0),
      wmean_rel_error_pmin = numeric(0),
      wmean_rel_error_pmax = numeric(0),
      median_lpr_range = numeric(0),
      median_p_min = numeric(0),
      median_p_max = numeric(0)
    ))
  }

  # Weighted mean for error metrics (parameter-derived); median for LPR range and pyruvate extremes (raw data)
  seg_df %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(
      wmean_rel_error_pmin = stats::weighted.mean(rel_error_pmin, w = n_points, na.rm = TRUE),
      wmean_rel_error_pmax = stats::weighted.mean(rel_error_pmax, w = n_points, na.rm = TRUE),
      median_lpr_range = stats::median(lpr_range, na.rm = TRUE),
      median_p_min = stats::median(p_min, na.rm = TRUE),
      median_p_max = stats::median(p_max, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(wmean_rel_error_pmin) | is.finite(wmean_rel_error_pmax)
    )
}

patient_all <- compute_patient_level(segment_metrics)
patient_Pb <- compute_patient_level(dplyr::filter(segment_metrics, seg_type == "Pb"))
patient_Nb <- compute_patient_level(dplyr::filter(segment_metrics, seg_type == "Nb"))

# ==============================================================================
# Cohort-Level Aggregation (median/IQR of patient-level values)
# ==============================================================================
summarise_cohort <- function(patient_df, group_code) {
  if (nrow(patient_df) == 0) {
    return(dplyr::tibble(
      group = group_code,
      q1_rel_err_pmin = NA_real_, median_rel_err_pmin = NA_real_, q3_rel_err_pmin = NA_real_,
      q1_rel_err_pmax = NA_real_, median_rel_err_pmax = NA_real_, q3_rel_err_pmax = NA_real_,
      q1_lpr_range = NA_real_, median_lpr_range = NA_real_, q3_lpr_range = NA_real_,
      q1_p_min = NA_real_, median_p_min = NA_real_, q3_p_min = NA_real_,
      q1_p_max = NA_real_, median_p_max = NA_real_, q3_p_max = NA_real_
    ))
  }

  s_pmin <- median_iqr(patient_df$wmean_rel_error_pmin)
  s_pmax <- median_iqr(patient_df$wmean_rel_error_pmax)
  s_range <- median_iqr(patient_df$median_lpr_range)
  s_p_min <- median_iqr(patient_df$median_p_min)
  s_p_max <- median_iqr(patient_df$median_p_max)

  dplyr::tibble(
    group = group_code,
    q1_rel_err_pmin = s_pmin$q1, median_rel_err_pmin = s_pmin$median, q3_rel_err_pmin = s_pmin$q3,
    q1_rel_err_pmax = s_pmax$q1, median_rel_err_pmax = s_pmax$median, q3_rel_err_pmax = s_pmax$q3,
    q1_lpr_range = s_range$q1, median_lpr_range = s_range$median, q3_lpr_range = s_range$q3,
    q1_p_min = s_p_min$q1, median_p_min = s_p_min$median, q3_p_min = s_p_min$q3,
    q1_p_max = s_p_max$q1, median_p_max = s_p_max$median, q3_p_max = s_p_max$q3
  )
}

cohort_all <- summarise_cohort(patient_all, "Unstratified")
cohort_Pb <- summarise_cohort(patient_Pb, "Pb")
cohort_Nb <- summarise_cohort(patient_Nb, "Nb")

all_cohort_stats <- dplyr::bind_rows(cohort_all, cohort_Pb, cohort_Nb) %>%
  dplyr::mutate(level = "Patient")

# ==============================================================================
# Prepare Summary Table (CSV + Notes TXT)
# ==============================================================================
summary_data <- all_cohort_stats %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    group_label = dplyr::case_when(
      group == "Unstratified" ~ "Unstratified",
      group == "Pb"          ~ "Type Pb (b > 0)",
      group == "Nb"          ~ "Type Nb (b < 0)",
      TRUE                    ~ group
    ),
    `Pmin (µM)` = format_interval(1000 * median_p_min, 1000 * q1_p_min, 1000 * q3_p_min, digits = 0, suffix = ""),
    `Pmax (µM)` = format_interval(1000 * median_p_max, 1000 * q1_p_max, 1000 * q3_p_max, digits = 0, suffix = ""),
    `Rel. LPR Error at Pmin (%)` = format_interval(median_rel_err_pmin, q1_rel_err_pmin, q3_rel_err_pmin, digits = 1, suffix = "%"),
    `Rel. LPR Error at Pmax (%)` = format_interval(median_rel_err_pmax, q1_rel_err_pmax, q3_rel_err_pmax, digits = 1, suffix = "%"),
    `LPR Range` = format_interval(median_lpr_range, q1_lpr_range, q3_lpr_range, digits = 1, suffix = "")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    dplyr::case_when(group == "Unstratified" ~ 1L, group == "Pb" ~ 2L, group == "Nb" ~ 3L, TRUE ~ 4L)
  )

summary_display <- summary_data %>%
  dplyr::select(
    `Stratification Group` = group_label,
    `Pmin (µM)`,
    `Rel. LPR Error at Pmin (%)`,
    `Pmax (µM)`,
    `Rel. LPR Error at Pmax (%)`,
    `LPR Range`
  )

summary_csv_path <- file.path(output_dir, "1_stats.csv")
write.csv(summary_display, summary_csv_path, row.names = FALSE)

# ==============================================================================
# Demographic Footnotes
# ==============================================================================
segments_all <- df_raw %>%
  dplyr::filter(is.finite(p6e_m)) %>%
  dplyr::mutate(
    seg_type = dplyr::case_when(
      p6e_b > 0 ~ "Pb",
      p6e_b < 0 ~ "Nb",
      TRUE      ~ "Zb"
    ),
    segment = paste(patient, p6e_seg_index, sep = "_")
  ) %>%
  dplyr::distinct(segment, .keep_all = TRUE)

segments_retained <- df_prepared %>%
  dplyr::distinct(segment, .keep_all = TRUE)

compute_demo_stats <- function(group_code) {
  if (group_code == "Unstratified") {
    seg_all_g <- segments_all
    seg_retained_g <- segments_retained
  } else if (group_code == "Pb") {
    seg_all_g <- dplyr::filter(segments_all, seg_type == "Pb")
    seg_retained_g <- dplyr::filter(segments_retained, seg_type == "Pb")
  } else if (group_code == "Nb") {
    seg_all_g <- dplyr::filter(segments_all, seg_type == "Nb")
    seg_retained_g <- dplyr::filter(segments_retained, seg_type == "Nb")
  } else {
    seg_all_g <- segments_all[0, ]
    seg_retained_g <- segments_retained[0, ]
  }

  n_patients_all  <- dplyr::n_distinct(seg_all_g$patient)
  n_segments_all  <- nrow(seg_all_g)
  n_patients_retained <- dplyr::n_distinct(seg_retained_g$patient)
  n_segments_retained <- nrow(seg_retained_g)

  pct_patients <- safe_divide(100 * n_patients_retained, n_patients_all)
  pct_segments <- safe_divide(100 * n_segments_retained, n_segments_all)

  spp_retained <- seg_retained_g %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(n_segments = dplyr::n(), .groups = "drop")

  if (nrow(spp_retained) > 0) {
    med_spp <- stats::median(spp_retained$n_segments, na.rm = TRUE)
    q1_spp  <- stats::quantile(spp_retained$n_segments, 0.25, na.rm = TRUE)
    q3_spp  <- stats::quantile(spp_retained$n_segments, 0.75, na.rm = TRUE)
  } else {
    med_spp <- NA_real_
    q1_spp  <- NA_real_
    q3_spp  <- NA_real_
  }

  list(
    n_patients_all  = n_patients_all,
    n_patients_retained = n_patients_retained,
    pct_patients    = pct_patients,
    n_segments_all  = n_segments_all,
    n_segments_retained = n_segments_retained,
    pct_segments    = pct_segments,
    med_spp         = med_spp,
    q1_spp          = q1_spp,
    q3_spp          = q3_spp
  )
}

demo_unstrat <- compute_demo_stats("Unstratified")
demo_Pb <- compute_demo_stats("Pb")
demo_Nb <- compute_demo_stats("Nb")

# ==============================================================================
# Exemplar Patient Plots
# ==============================================================================
# Selection criteria:
#   1. All segments have m > 0 (no discarded segments)
#   2. No Type Zb segments (b = 0) - only clear Pb/Nb examples
#   3. Has both Type Pb (b > 0) AND Type Nb (b < 0) segments
#   4. In each segment: ALL LPR values on 'correct' side of m
#   5. In each segment: |LPR error| at Pmin in top quartile (≥75th percentile)
#   6. In each segment: |LPR error| at Pmax in bottom quartile (≤25th percentile)
#
# Segments with m ≤ 0 are excluded from exemplar and plots_all outputs.

# Plot aesthetics
FIGURE_WIDTH_MM <- 185.0
FIGURE_HEIGHT_MM <- 55.0
MAX_SEG_COLS <- 3

save_plot_portable <- function(filename, plot, ...) {
  if (isTRUE(capabilities("cairo"))) {
    ggsave(filename = filename, plot = plot, type = "cairo", ...)
  } else {
    warning(sprintf("Cairo backend unavailable; saving %s with the default device.", basename(filename)))
    ggsave(filename = filename, plot = plot, ...)
  }
}

common_theme <- theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

# Colors
COLOR_PB <- "#44AA99"   # Teal for Type Pb (b > 0)
COLOR_NB <- "#DDCC77"   # Gold for Type Nb (b < 0)
COLOR_ZB <- "#888888"   # Gray for Type Zb (b = 0)
COLOR_PB_FILL <- "#44AA9933"  # Teal with alpha for ribbon
COLOR_NB_FILL <- "#DDCC7733"  # Gold with alpha for ribbon
COLOR_ZB_FILL <- "#88888833"  # Gray with alpha for ribbon

# Create output directory for exemplar plots
exemplar_dir <- file.path(output_dir, "2_exemplar")
if (!dir.exists(exemplar_dir)) {
  dir.create(exemplar_dir, recursive = TRUE)
}

# Compute segment-level validation metrics for exemplar selection
compute_segment_validation <- function(seg_df) {
  m <- seg_df$p6e_m[1]
  b <- seg_df$p6e_b[1]
  r_val <- seg_df$p6e_r[1]
  
  if (!is.finite(m) || !is.finite(b) || m <= 0) {
    return(list(valid = FALSE, reason = "m <= 0 or non-finite coefficients"))
  }
  
  seg_type <- if (b > 0) "Pb" else if (b < 0) "Nb" else "Zb"
  
  lpr_vals <- seg_df$lpr[is.finite(seg_df$lpr) & is.finite(seg_df$pyruvate) & seg_df$pyruvate > 0]
  pyr_vals <- seg_df$pyruvate[is.finite(seg_df$lpr) & is.finite(seg_df$pyruvate) & seg_df$pyruvate > 0]
  
  if (length(lpr_vals) < 2) {
    return(list(valid = FALSE, reason = "fewer than 2 valid LPR points"))
  }
  
  # Check ALL LPR values on correct side of m
  if (seg_type == "Pb") {
    all_correct_side <- all(lpr_vals > m)
  } else if (seg_type == "Nb") {
    all_correct_side <- all(lpr_vals < m)
  } else {
    all_correct_side <- FALSE
  }
  
  # Compute all |errors|
  all_errors <- safe_divide(100 * (lpr_vals - m), m)
  all_abs_errors <- abs(all_errors)
  
  idx_pmin <- which.min(pyr_vals)
  idx_pmax <- which.max(pyr_vals)
  
  abs_error_pmin <- all_abs_errors[idx_pmin]
  abs_error_pmax <- all_abs_errors[idx_pmax]
  
  # Percentile-based criterion
  q75 <- quantile(all_abs_errors, 0.75)
  q25 <- quantile(all_abs_errors, 0.25)
  
  pmin_has_max_error <- abs_error_pmin >= q75
  pmax_has_min_error <- abs_error_pmax <= q25
  
  list(
    valid = TRUE,
    m = m,
    b = b,
    r = r_val,
    seg_type = seg_type,
    lpr_min = min(lpr_vals),
    lpr_max = max(lpr_vals),
    p_min = min(pyr_vals),
    p_max = max(pyr_vals),
    error_pmin = all_errors[idx_pmin],
    error_pmax = all_errors[idx_pmax],
    all_correct_side = all_correct_side,
    pmin_has_max_error = pmin_has_max_error,
    pmax_has_min_error = pmax_has_min_error,
    n_points = length(lpr_vals)
  )
}

# Validate patient for exemplar selection
validate_exemplar_patient <- function(patient_id, df_raw) {
  df_patient <- df_raw %>%
    filter(patient == patient_id, is.finite(p6e_m), is.finite(p6e_b))
  
  if (nrow(df_patient) == 0) {
    return(list(valid = FALSE, reason = "No data"))
  }
  
  segments <- unique(df_patient$p6e_seg_index)
  if (length(segments) < 2) {
    return(list(valid = FALSE, reason = "Fewer than 2 segments"))
  }
  
  seg_metrics <- list()
  has_pb <- FALSE
  has_nb <- FALSE
  
  for (seg_idx in segments) {
    seg_df <- df_patient %>% filter(p6e_seg_index == seg_idx)
    metrics <- compute_segment_validation(seg_df)
    
    # Criterion 1: m > 0 and sufficient data
    if (!metrics$valid) {
      return(list(valid = FALSE, reason = sprintf("Segment %d: %s", seg_idx, metrics$reason)))
    }
    
    seg_type <- metrics$seg_type
    
    # Criterion 2: No Type Zb segments (only clear Pb/Nb examples)
    if (seg_type == "Zb") {
      return(list(valid = FALSE, reason = sprintf("Segment %d is Type Zb (b = 0)", seg_idx)))
    }
    
    if (seg_type == "Pb") has_pb <- TRUE
    if (seg_type == "Nb") has_nb <- TRUE
    
    # Criterion 4: All LPR on 'correct' side
    if (!metrics$all_correct_side) {
      return(list(valid = FALSE, reason = sprintf("Segment %d: not all LPR on correct side", seg_idx)))
    }
    
    # Criterion 5: Pmin in top quartile
    if (!metrics$pmin_has_max_error) {
      return(list(valid = FALSE, reason = sprintf("Segment %d: Pmin not in top quartile", seg_idx)))
    }
    
    # Criterion 6: Pmax in bottom quartile
    if (!metrics$pmax_has_min_error) {
      return(list(valid = FALSE, reason = sprintf("Segment %d: Pmax not in bottom quartile", seg_idx)))
    }
    
    metrics$seg_idx <- seg_idx
    seg_metrics[[length(seg_metrics) + 1]] <- metrics
  }
  
  # Criterion 3: Must have both types
  if (!has_pb) {
    return(list(valid = FALSE, reason = "No Type Pb segments"))
  }
  if (!has_nb) {
    return(list(valid = FALSE, reason = "No Type Nb segments"))
  }
  
  list(valid = TRUE, reason = "Valid", seg_metrics = seg_metrics, df_patient = df_patient)
}

# ==============================================================================
# Shared Plot Building Functions
# ==============================================================================

#' Build LPR error plot from prepared data
#' @param plot_df Data frame with columns: time_since_injury, lpr, p6e_seg_index, seg_type, m
#' @return ggplot object
build_lpr_error_plot <- function(plot_df) {
  # Build m-line segments for geom_segment (with seg_type for coloring)
  m_lines <- plot_df %>%
    group_by(p6e_seg_index) %>%
    summarise(
      t_min = min(time_since_injury, na.rm = TRUE),
      t_max = max(time_since_injury, na.rm = TRUE),
      m = first(m),
      seg_type = first(seg_type),
      .groups = "drop"
    ) %>%
    filter(is.finite(m))
  
  # Build polygon data for shading between LPR line and m line
  # For each segment: trace LPR points left-to-right, then m line right-to-left
  polygon_df <- plot_df %>%
    filter(!is.na(seg_type), is.finite(m), is.finite(lpr)) %>%
    arrange(p6e_seg_index, time_since_injury) %>%
    group_by(p6e_seg_index) %>%
    group_modify(~ {
      seg_data <- .x
      m_val <- seg_data$m[1]
      seg_type_val <- seg_data$seg_type[1]
      # LPR path (left to right)
      lpr_path <- data.frame(
        x = seg_data$time_since_injury,
        y = seg_data$lpr
      )
      # m line path (right to left, to close the polygon)
      m_path <- data.frame(
        x = rev(seg_data$time_since_injury),
        y = rep(m_val, nrow(seg_data))
      )
      # Combine into closed polygon
      bind_rows(lpr_path, m_path) %>%
        mutate(seg_type = seg_type_val, m = m_val)
    }) %>%
    ungroup()
  
  ggplot(plot_df, aes(x = time_since_injury)) +
    # Shaded polygon between LPR line and m line (error visualization)
    geom_polygon(
      data = filter(polygon_df, seg_type == "Pb"),
      aes(x = x, y = y, group = p6e_seg_index),
      fill = COLOR_PB_FILL, color = NA
    ) +
    geom_polygon(
      data = filter(polygon_df, seg_type == "Nb"),
      aes(x = x, y = y, group = p6e_seg_index),
      fill = COLOR_NB_FILL, color = NA
    ) +
    geom_polygon(
      data = filter(polygon_df, seg_type == "Zb"),
      aes(x = x, y = y, group = p6e_seg_index),
      fill = COLOR_ZB_FILL, color = NA
    ) +
    # Colored dashed m lines (matching segment type)
    geom_segment(
      data = filter(m_lines, seg_type == "Pb"),
      aes(x = t_min, xend = t_max, y = m, yend = m),
      color = COLOR_PB, linetype = "dashed", linewidth = 0.7
    ) +
    geom_segment(
      data = filter(m_lines, seg_type == "Nb"),
      aes(x = t_min, xend = t_max, y = m, yend = m),
      color = COLOR_NB, linetype = "dashed", linewidth = 0.7
    ) +
    geom_segment(
      data = filter(m_lines, seg_type == "Zb"),
      aes(x = t_min, xend = t_max, y = m, yend = m),
      color = COLOR_ZB, linetype = "dashed", linewidth = 0.7
    ) +
    # LPR lines connecting points within each segment
    geom_line(
      aes(y = lpr, color = seg_type, group = p6e_seg_index),
      linewidth = 0.5, alpha = 0.8
    ) +
    # LPR points with shapes: circles for Pb, triangles for Nb, squares for Zb
    geom_point(
      aes(y = lpr, color = seg_type, shape = seg_type),
      size = 1.5, alpha = 0.9
    ) +
    scale_color_manual(
      values = c("Pb" = COLOR_PB, "Nb" = COLOR_NB, "Zb" = COLOR_ZB),
      na.value = "gray50"
    ) +
    scale_shape_manual(
      values = c("Pb" = 16, "Nb" = 17, "Zb" = 15)  # 16 = circle, 17 = triangle, 15 = square
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(x = "Time Since Injury (hours)", y = "LPR") +
    common_theme
}

#' Build a single lactate vs pyruvate scatter panel for one segment
#' @param seg_data Data frame for one segment (must include pyruvate, lactate, p6e_m, p6e_b)
#' @param seg_type_val Character: "Pb", "Nb", or "Zb"
#' @param m_val Numeric gradient
#' @param b_val Numeric intercept
#' @return ggplot object or NULL
build_segment_panel <- function(seg_data, seg_type_val, m_val, b_val) {
  color <- if (seg_type_val == "Pb") COLOR_PB else if (seg_type_val == "Nb") COLOR_NB else COLOR_ZB
  shape_val <- if (seg_type_val == "Pb") 16 else if (seg_type_val == "Nb") 17 else 15

  seg_data <- seg_data %>%
    filter(is.finite(pyruvate), is.finite(lactate), pyruvate > 0)

  if (nrow(seg_data) == 0) return(NULL)

  x_range <- range(seg_data$pyruvate, na.rm = TRUE)
  line_df <- data.frame(
    pyruvate_um = x_range * 1000,
    lactate = m_val * x_range + b_val
  )

  ggplot(seg_data, aes(x = pyruvate * 1000, y = lactate)) +
    geom_point(color = color, shape = shape_val, size = 1.5, alpha = 0.8) +
    geom_line(data = line_df, aes(x = pyruvate_um, y = lactate),
              color = color, linewidth = 0.7) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(x = "Pyruvate (\u00b5M)", y = "Lactate (mM)") +
    common_theme
}

#' Build a combined lactate vs pyruvate panel with all retained segments overlaid
#' @param df_patient Patient data (must include p6e_seg_index, p6e_m, p6e_b, pyruvate, lactate)
#' @return ggplot object or NULL
build_combined_segment_panel <- function(df_patient) {
  seg_data_all <- df_patient %>%
    filter(
      is.finite(p6e_m), p6e_m > 0, is.finite(p6e_b),
      is.finite(pyruvate), is.finite(lactate), pyruvate > 0
    ) %>%
    mutate(
      seg_type = case_when(
        p6e_b > 0 ~ "Pb",
        p6e_b < 0 ~ "Nb",
        TRUE ~ "Zb"
      )
    )

  if (nrow(seg_data_all) == 0) return(NULL)

  line_data <- seg_data_all %>%
    group_by(p6e_seg_index, seg_type) %>%
    summarise(
      m = first(p6e_m),
      b = first(p6e_b),
      x_min = min(pyruvate),
      x_max = max(pyruvate),
      .groups = "drop"
    ) %>%
    mutate(
      y_min = m * x_min + b,
      y_max = m * x_max + b
    )

  ggplot(seg_data_all, aes(x = pyruvate * 1000, y = lactate, color = seg_type, shape = seg_type)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_segment(
      data = line_data,
      aes(x = x_min * 1000, xend = x_max * 1000, y = y_min, yend = y_max, color = seg_type),
      linewidth = 0.7, inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c("Pb" = COLOR_PB, "Nb" = COLOR_NB, "Zb" = COLOR_ZB),
      na.value = "gray50"
    ) +
    scale_shape_manual(values = c("Pb" = 16, "Nb" = 17, "Zb" = 15)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(x = "Pyruvate (\u00b5M)", y = "Lactate (mM)") +
    common_theme
}

#' Build per-segment lactate vs pyruvate panels for a patient
#' @param df_patient Patient data (must include p6e_seg_index, p6e_m, p6e_b, pyruvate, lactate)
#' @return List of ggplot objects (one per retained segment)
build_patient_segment_panels <- function(df_patient) {
  seg_indices <- df_patient %>%
    filter(is.finite(p6e_m), p6e_m > 0, is.finite(p6e_b)) %>%
    distinct(p6e_seg_index) %>%
    arrange(p6e_seg_index) %>%
    pull(p6e_seg_index)

  panels <- list()
  for (idx in seg_indices) {
    seg_data <- df_patient %>% filter(p6e_seg_index == idx)
    m_val <- seg_data$p6e_m[1]
    b_val <- seg_data$p6e_b[1]
    seg_type_val <- if (b_val > 0) "Pb" else if (b_val < 0) "Nb" else "Zb"
    p <- build_segment_panel(seg_data, seg_type_val, m_val, b_val)
    if (!is.null(p)) panels[[length(panels) + 1]] <- p
  }
  panels
}

#' Combine per-segment panels (top rows) with LPR vs time panel (bottom)
#' Panels are labeled A, B, C, ... and segment panels stretch to fill each row.
#' @param seg_panels List of ggplot objects from build_patient_segment_panels
#' @param lpr_panel ggplot from build_lpr_error_plot
#' @return List with $plot (patchwork or ggplot) and $fig_height (mm)
assemble_figure <- function(seg_panels, lpr_panel) {
  n_segs <- length(seg_panels)
  if (n_segs == 0) {
    return(list(plot = lpr_panel, fig_height = FIGURE_HEIGHT_MM))
  }

  n_seg_rows <- ceiling(n_segs / MAX_SEG_COLS)

  # Balanced row sizes: e.g. 4 panels → 2+2 not 3+1
  base_per_row <- n_segs %/% n_seg_rows
  extra <- n_segs %% n_seg_rows
  row_sizes <- rep(base_per_row, n_seg_rows)
  if (extra > 0) row_sizes[seq_len(extra)] <- row_sizes[seq_len(extra)] + 1

  # Build each segment row with | (panels fill row width equally)
  seg_rows <- list()
  panel_idx <- 1
  for (r in seq_len(n_seg_rows)) {
    n_in_row <- row_sizes[r]
    seg_rows[[r]] <- Reduce(`|`, seg_panels[panel_idx:(panel_idx + n_in_row - 1)])
    panel_idx <- panel_idx + n_in_row
  }

  # Stack segment rows + LPR panel with /
  combined <- Reduce(`/`, c(seg_rows, list(lpr_panel)))

  combined <- combined +
    patchwork::plot_layout(heights = rep(1, n_seg_rows + 1)) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(size = 9, face = "bold", colour = "black",
                              hjust = 1, vjust = -1),
      plot.tag.position = c(0.02, 0.99),
      plot.margin = margin(t = 3, r = 3, b = 2, l = 3, unit = "mm")
    )

  fig_height <- (n_seg_rows + 1) * FIGURE_HEIGHT_MM

  list(plot = combined, fig_height = fig_height)
}

# Create plot for one exemplar patient (uses shared build_lpr_error_plot)
create_exemplar_plot <- function(patient_id, df_patient, seg_metrics) {
  seg_info <- data.frame(
    seg_idx = sapply(seg_metrics, function(x) x$seg_idx),
    seg_type = sapply(seg_metrics, function(x) x$seg_type),
    m = sapply(seg_metrics, function(x) x$m),
    r = sapply(seg_metrics, function(x) x$r)
  )

  plot_df <- df_patient %>%
    filter(is.finite(lpr), is.finite(time_since_injury)) %>%
    left_join(seg_info, by = c("p6e_seg_index" = "seg_idx")) %>%
    arrange(time_since_injury)

  seg_panels <- build_patient_segment_panels(df_patient)
  combined_panel <- build_combined_segment_panel(df_patient)
  if (!is.null(combined_panel)) seg_panels <- c(list(combined_panel), seg_panels)
  lpr_panel <- build_lpr_error_plot(plot_df)
  assemble_figure(seg_panels, lpr_panel)
}

# Write metadata file for exemplar patient
write_exemplar_metadata <- function(patient_id, df_patient, seg_metrics, txt_path) {
  n_points <- nrow(df_patient %>% filter(is.finite(lpr)))
  t_range <- range(df_patient$time_since_injury, na.rm = TRUE)
  n_pb <- sum(sapply(seg_metrics, function(x) x$seg_type == "Pb"))
  n_nb <- sum(sapply(seg_metrics, function(x) x$seg_type == "Nb"))
  n_zb <- sum(sapply(seg_metrics, function(x) x$seg_type == "Zb"))
  
  lines <- c(
    sprintf("Patient: %s", patient_id),
    paste(rep("=", 60), collapse = ""),
    "",
    "Hyperbolic LPR Error Exemplar",
    paste(rep("-", 40), collapse = ""),
    "",
    sprintf("Total points: %d", n_points),
    sprintf("Time span: %.1f–%.1f hours", t_range[1], t_range[2]),
    sprintf("Segments: %d total (%d Type Pb, %d Type Nb, %d Type Zb)", length(seg_metrics), n_pb, n_nb, n_zb),
    "",
    "Figure elements:",
    "  - LPR: Teal circles = Type Pb (b > 0), Gold triangles = Type Nb (b < 0), Gray squares = Type Zb (b = 0; none shown in exemplars)",
    "  - Solid lines connect LPR points within each segment",
    "  - Dashed lines: Segment gradient (m) values, colored by segment type",
    "  - Shaded regions: LPR error (deviation from m)",
    "",
    "Segment Details",
    paste(rep("-", 40), collapse = ""),
    ""
  )
  
  for (metrics in seg_metrics) {
    seg_idx <- metrics$seg_idx
    seg_type <- metrics$seg_type
    m <- metrics$m
    b <- metrics$b
    r <- metrics$r
    p_min_um <- metrics$p_min * 1000
    p_max_um <- metrics$p_max * 1000
    lpr_min <- metrics$lpr_min
    lpr_max <- metrics$lpr_max
    err_pmin <- metrics$error_pmin
    err_pmax <- metrics$error_pmax
    n_pts <- metrics$n_points
    
    lines <- c(lines,
      sprintf("Segment %d (Type %s):", seg_idx, seg_type),
      sprintf("  n = %d points", n_pts),
      sprintf("  m = %.4f, b = %.4f mM", m, b),
      if (is.finite(r)) sprintf("  r = %.3f", r) else NULL,
      sprintf("  Pyruvate range: %.0f–%.0f µM", p_min_um, p_max_um),
      sprintf("  LPR range: %.2f–%.2f (m = %.2f)", lpr_min, lpr_max, m),
      if (seg_type == "Pb") {
        sprintf("    → Criterion: ALL LPR > m? min LPR (%.2f) > m (%.2f)? %s", lpr_min, m, lpr_min > m)
      } else {
        sprintf("    → Criterion: ALL LPR < m? max LPR (%.2f) < m (%.2f)? %s", lpr_max, m, lpr_max < m)
      },
      sprintf("  LPR error at Pmin: %+.1f%%", err_pmin),
      sprintf("  LPR error at Pmax: %+.1f%%", err_pmax),
      ""
    )
  }
  
  lines <- c(lines,
    paste(rep("=", 60), collapse = ""),
    "Selection criteria met:",
    "  ✓ All segments have m > 0",
    "  ✓ No Type Zb segments (b = 0)",
    "  ✓ Has both Type Pb and Type Nb segments",
    "  ✓ ALL LPR values on correct side of m in each segment",
    "  ✓ |LPR error| at Pmin in top quartile in each segment",
    "  ✓ |LPR error| at Pmax in bottom quartile in each segment"
  )
  
  writeLines(lines, txt_path)
}

# Find and generate exemplar plots
all_patients <- unique(df_raw$patient)
qualifying_patients <- c()

for (pid in all_patients) {
  result <- validate_exemplar_patient(pid, df_raw)
  if (result$valid) {
    qualifying_patients <- c(qualifying_patients, pid)
    
    # Create plot
    fig <- create_exemplar_plot(pid, result$df_patient, result$seg_metrics)
    png_path <- file.path(exemplar_dir, sprintf("%s.png", pid))
    save_plot_portable(
      png_path, fig$plot,
      width = FIGURE_WIDTH_MM, height = fig$fig_height, units = "mm", dpi = 300
    )
    
    # Write metadata
    txt_path <- file.path(exemplar_dir, sprintf("%s_notes.txt", pid))
    write_exemplar_metadata(pid, result$df_patient, result$seg_metrics, txt_path)
  }
}

n_qualifying <- length(qualifying_patients)

# Write summary file
summary_lines <- c(
  "Hyperbolic LPR Error Exemplar Plots",
  paste(rep("=", 60), collapse = ""),
  "",
  sprintf("Total patients screened: %d", length(all_patients)),
  sprintf("Qualifying patients: %d", n_qualifying),
  "",
  "Selection criteria:",
  "  1. All segments have m > 0 (no discarded segments)",
  "  2. No Type Zb segments (b = 0) - only clear Pb/Nb examples",
  "  3. Has both Type Pb (b > 0) AND Type Nb (b < 0) segments",
  "  4. In each segment: ALL LPR values on correct side of m",
  "     - Type Pb (b > 0): ALL LPR > m",
  "     - Type Nb (b < 0): ALL LPR < m",
  "  5. In each segment: |LPR error| at Pmin in top quartile (>=75th pctl)",
  "  6. In each segment: |LPR error| at Pmax in bottom quartile (<=25th pctl)",
  "",
  "Note: Segments with m <= 0 are excluded from ALL plots (consistent with",
  "      'retained segments' definition used throughout the pipeline).",
  "",
  "Qualifying patient IDs:",
  paste0("  ", qualifying_patients)
)
writeLines(summary_lines, file.path(exemplar_dir, "_summary.txt"))

# ==============================================================================
# All Patient Plots (plots_all)
# ==============================================================================
# Generate plots for ALL patients with at least one retained segment (m > 0)

plots_all_dir <- file.path(output_dir, "3_all")
if (!dir.exists(plots_all_dir)) {
  dir.create(plots_all_dir, recursive = TRUE)
}

# Function to create plot for any patient (uses shared build_lpr_error_plot)
create_patient_plot <- function(patient_id, df_raw) {
  df_patient <- df_raw %>%
    filter(patient == patient_id, is.finite(p6e_m), is.finite(p6e_b), p6e_m > 0)
  
  if (nrow(df_patient) == 0) {
    return(NULL)
  }
  
  # Build segment info
  seg_info <- df_patient %>%
    group_by(p6e_seg_index) %>%
    summarise(
      m = first(p6e_m),
      b = first(p6e_b),
      r = first(p6e_r),
      seg_type = case_when(
        first(p6e_b) > 0 ~ "Pb",
        first(p6e_b) < 0 ~ "Nb",
        TRUE ~ "Zb"
      ),
      .groups = "drop"
    ) %>%
    rename(seg_idx = p6e_seg_index)
  
  # Prepare plot data
  plot_df <- df_patient %>%
    filter(is.finite(lpr), is.finite(time_since_injury)) %>%
    left_join(seg_info, by = c("p6e_seg_index" = "seg_idx")) %>%
    arrange(time_since_injury)
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  seg_panels <- build_patient_segment_panels(df_patient)
  combined_panel <- build_combined_segment_panel(df_patient)
  if (!is.null(combined_panel)) seg_panels <- c(list(combined_panel), seg_panels)
  lpr_panel <- build_lpr_error_plot(plot_df)
  fig <- assemble_figure(seg_panels, lpr_panel)

  list(plot = fig$plot, fig_height = fig$fig_height, seg_info = seg_info, df_patient = df_patient)
}

# Write metadata file for any patient (plots_all)
write_patient_metadata <- function(patient_id, df_patient, seg_info, txt_path) {
  n_points <- nrow(df_patient %>% filter(is.finite(lpr)))
  t_range <- range(df_patient$time_since_injury, na.rm = TRUE)
  n_pb <- sum(seg_info$seg_type == "Pb")
  n_nb <- sum(seg_info$seg_type == "Nb")
  n_zb <- sum(seg_info$seg_type == "Zb")
  n_segs <- nrow(seg_info)
  
  lines <- c(
    sprintf("Patient: %s", patient_id),
    paste(rep("=", 60), collapse = ""),
    "",
    "Hyperbolic LPR Error Plot",
    paste(rep("-", 40), collapse = ""),
    "",
    sprintf("Total points: %d", n_points),
    sprintf("Time span: %.1f–%.1f hours", t_range[1], t_range[2]),
    sprintf("Segments: %d total (%d Type Pb, %d Type Nb, %d Type Zb)", n_segs, n_pb, n_nb, n_zb),
    "",
    "Figure elements:",
    "  - LPR: Teal circles = Type Pb (b > 0), Gold triangles = Type Nb (b < 0), Gray squares = Type Zb (b = 0)",
    "  - Solid lines connect LPR points within each segment",
    "  - Dashed lines: Segment gradient (m) values, colored by segment type",
    "  - Shaded regions: LPR error (deviation from m)",
    "",
    "Segment Details",
    paste(rep("-", 40), collapse = ""),
    ""
  )
  
  for (i in seq_len(nrow(seg_info))) {
    seg_idx <- seg_info$seg_idx[i]
    seg_type <- seg_info$seg_type[i]
    m <- seg_info$m[i]
    b <- seg_info$b[i]
    r <- seg_info$r[i]
    
    seg_data <- df_patient %>% 
      filter(p6e_seg_index == seg_idx, is.finite(lpr), is.finite(pyruvate), pyruvate > 0)
    
    n_pts <- nrow(seg_data)
    if (n_pts > 0) {
      p_min_um <- min(seg_data$pyruvate) * 1000
      p_max_um <- max(seg_data$pyruvate) * 1000
      lpr_min <- min(seg_data$lpr)
      lpr_max <- max(seg_data$lpr)
      
      # Compute errors at Pmin and Pmax
      idx_pmin <- which.min(seg_data$pyruvate)
      idx_pmax <- which.max(seg_data$pyruvate)
      err_pmin <- safe_divide(100 * (seg_data$lpr[idx_pmin] - m), m)
      err_pmax <- safe_divide(100 * (seg_data$lpr[idx_pmax] - m), m)
      
      lines <- c(lines,
        sprintf("Segment %d (Type %s):", seg_idx, seg_type),
        sprintf("  n = %d points", n_pts),
        sprintf("  m = %.4f, b = %.4f mM", m, b),
        if (is.finite(r)) sprintf("  r = %.3f", r) else NULL,
        sprintf("  Pyruvate range: %.0f–%.0f µM", p_min_um, p_max_um),
        sprintf("  LPR range: %.2f–%.2f (m = %.2f)", lpr_min, lpr_max, m),
        sprintf("  LPR error at Pmin: %+.1f%%", err_pmin),
        sprintf("  LPR error at Pmax: %+.1f%%", err_pmax),
        ""
      )
    }
  }
  
  writeLines(lines, txt_path)
}

# Get all patients with retained segments
patients_with_segments <- df_raw %>%
  filter(is.finite(p6e_m), p6e_m > 0) %>%
  distinct(patient) %>%
  pull(patient)

# Create worker factory for parallel processing (explicit variable passing via closure)
# This avoids relying on global variable access in forked processes
make_patient_plot_worker <- function(df_raw, output_dir, fig_width) {
  force(df_raw)
  force(output_dir)
  force(fig_width)

  function(pid) {
    result <- create_patient_plot(pid, df_raw)

    if (!is.null(result)) {
      png_path <- file.path(output_dir, sprintf("%s.png", pid))
      save_plot_portable(
        png_path, result$plot,
        width = fig_width, height = result$fig_height, units = "mm", dpi = 300
      )

      txt_path <- file.path(output_dir, sprintf("%s_notes.txt", pid))
      write_patient_metadata(pid, result$df_patient, result$seg_info, txt_path)

      return(TRUE)
    }
    return(FALSE)
  }
}

# Create worker with explicit dependencies
process_patient_plot <- make_patient_plot_worker(
  df_raw = df_raw,
  output_dir = plots_all_dir,
  fig_width = FIGURE_WIDTH_MM
)

# Parallel processing
# mclapply uses forking (Unix/macOS only); forced to 1 core on Windows
n_cores <- if (.Platform$OS.type == "windows") 1L else max(1L, detectCores() - 1L)

results <- mclapply(patients_with_segments, process_patient_plot, mc.cores = n_cores)

# ==============================================================================
# Execution Time Log
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

cat("\n\n6_hyperbolic_LPR_error.R complete.\n\n")

# Build simple legends for table and figures
legends <- list(
  list(
    target = "Table 1 (1_stats.csv)",
    caption = sprintf(
      paste(
        "LPR error statistics quantifying deviation from the linear gradient (m) across pyruvate concentrations.",
        "Values are cohort medians (IQR) of patient-level metrics.",
        "Patient metrics: segment-weighted means (weights = n_points) of relative LPR error at segment Pmin and Pmax; patient medians of segment-level Pmin, Pmax, and LPR range.",
        "Relative LPR error = 100 × (LPR - m) / m at segment Pmin and Pmax.",
        "Cohort sizes: Unstratified = %d patients, %d segments; Type Pb (b > 0) = %d patients, %d segments; Type Nb (b < 0) = %d patients, %d segments.",
        sep = " "
      ),
      demo_unstrat$n_patients_retained, demo_unstrat$n_segments_retained,
      demo_Pb$n_patients_retained, demo_Pb$n_segments_retained,
      demo_Nb$n_patients_retained, demo_Nb$n_segments_retained
    ),
    footnotes = character(),
    abbreviations = "IQR, interquartile range; LPR, lactate/pyruvate ratio; Pmax, maximum pyruvate; Pmin, minimum pyruvate; Pb, positive intercept; Nb, negative intercept"
  ),
  list(
    target = "Figure (2_exemplar/*.png)",
    caption = sprintf(
      paste(
        "Exemplar patient LPR error plots.",
        "Top row(s): per-segment lactate vs pyruvate scatter with OLS regression lines, colored teal (Type Pb, b > 0) or gold (Type Nb, b < 0).",
        "Bottom panel: LPR versus time since injury (hours); dashed lines = gradient asymptote (m); shaded regions = deviation of LPR from m.",
        "Teal circles = Type Pb; gold triangles = Type Nb; Zb (b = 0) segments absent by exemplar criteria.",
        "%d patients met all exemplar criteria.",
        sep = " "
      ),
      n_qualifying
    ),
    abbreviations = "LPR, lactate/pyruvate ratio; Nb, negative intercept; OLS, ordinary least squares; Pb, positive intercept; Zb, zero intercept"
  ),
  list(
    target = "Figure (3_all/*.png)",
    caption = paste(
      "Individual patient LPR error plots for all patients with at least one retained segment.",
      "Top row(s): per-segment lactate vs pyruvate scatter with OLS regression lines, colored by segment type.",
      "Bottom panel: LPR versus time since injury (hours); same visual encoding as exemplar plots.",
      "Gray = Type Zb segments (b = 0) when present."
    ),
    abbreviations = "LPR, lactate/pyruvate ratio; Nb, negative intercept; OLS, ordinary least squares; Pb, positive intercept; Zb, zero intercept"
  )
)

writeLines(build_notes(legends, title = "6_hyperbolic_LPR_error"), file.path(output_dir, "0_notes.txt"))
