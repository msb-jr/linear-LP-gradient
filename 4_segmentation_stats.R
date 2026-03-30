# ==============================================================================
# Script: 4_segmentation_stats.R
# Manuscript relevance: 3.1, Table S1, Table S4
# ==============================================================================
# PURPOSE:
#   Summarize the input to / output of the linear segmentation pipeline (p6e)
#   and the yield after applying the positive-gradient (m > 0) retention rule.
#
# INPUT:
#   - 1_output/1_pre-processing/1_data.csv: Pre-processed data (p6e in)
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data (p6e out)
#
# OUTPUT:
#   - 1_output/4_segmentation_stats/__execution_time.csv: Runtime log
#   - 1_output/4_segmentation_stats/0_notes.txt: Table footnotes
#   - 1_output/4_segmentation_stats/1_stats.csv: Summary table
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 4_segmentation_stats.R\n")

# ==============================================================================
# Input / Output
# ==============================================================================
start_time <- Sys.time()

output_dir <- here::here("1_output", "4_segmentation_stats")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results_input_file <- here::here("1_output", "2_linear_segmentation", "1_results.csv")
source_input_file <- here::here("1_output", "1_pre-processing", "1_data.csv")

if (!file.exists(results_input_file)) {
  cat("\n\nError: Results file not found. Run 2_linear_segmentation.py first.\n")
  stop()
}

# ==============================================================================
# Load and Pre-process Data
# ==============================================================================

# Load result data
point_data <- read.csv(results_input_file)

# Load source data
if (!file.exists(source_input_file)) {
  cat("\n\nError: Source file not found:", source_input_file, "\n")
  stop()
}
source_point_data <- read.csv(source_input_file)

# --- Source Data Stats ---
datapoints_per_patient <- source_point_data %>%
  group_by(patient) %>%
  summarise(
    n_datapoints = n(),
    data_span_days = (max(time_since_injury, na.rm = TRUE) - min(time_since_injury, na.rm = TRUE)) / 24,
    .groups = 'drop'
  )

n_total_patients <- nrow(datapoints_per_patient)
n_total_datapoints <- sum(datapoints_per_patient$n_datapoints)

# ==============================================================================
# Helper Functions
# ==============================================================================

format_median_iqr_simple <- function(df, var_name, decimals = 0) {
  vals <- df[[var_name]]
  if (all(is.na(vals))) return(list(Median=NA, IQR=NA, Range=NA))
  
  q1 <- quantile(vals, 0.25, na.rm = TRUE)
  med <- median(vals, na.rm = TRUE)
  q3 <- quantile(vals, 0.75, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  fmt <- paste0("%.", decimals, "f")
  list(
    Median = sprintf(fmt, med),
    IQR = paste0(sprintf(fmt, q1), "–", sprintf(fmt, q3)),
    Range = paste0(sprintf(fmt, min_val), "–", sprintf(fmt, max_val))
  )
}

format_cohort_median_of_patient_medians <- function(segment_df, var_name, decimals = 0) {
  patient_meds <- segment_df %>%
    group_by(patient) %>%
    summarise(patient_median = median(.data[[var_name]], na.rm = TRUE), .groups = 'drop')
  
  if (nrow(patient_meds) == 0) return(list(Median=NA, IQR=NA, Range=NA))
  
  vals <- patient_meds$patient_median
  q1 <- quantile(vals, 0.25, na.rm = TRUE)
  med <- median(vals, na.rm = TRUE)
  q3 <- quantile(vals, 0.75, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  fmt <- paste0("%.", decimals, "f")
  list(
    Median = sprintf(fmt, med),
    IQR = paste0(sprintf(fmt, q1), "–", sprintf(fmt, q3)),
    Range = paste0(sprintf(fmt, min_val), "–", sprintf(fmt, max_val))
  )
}

safe_divide <- function(numerator, denominator, tol = 1e-12, fill = NA_real_) {
  result <- numerator / denominator
  invalid <- !is.finite(result) | !is.finite(denominator) | abs(denominator) <= tol
  result[invalid] <- fill
  result
}


# ==============================================================================
# Generate Segmentation Statistics
# ==============================================================================

# --- Segment-Level Stats ---
# Group by patient and p6e_seg_index
segment_data <- point_data %>%
  group_by(patient, p6e_seg_index) %>%
  summarise(
    m = first(p6e_m),
    n_datapoints = n(),
    data_span = max(time_since_injury, na.rm = TRUE) - min(time_since_injury, na.rm = TRUE),
    .groups = 'drop'
  )

n_total_segments <- nrow(segment_data)

# Retained Segments (Positive Gradient)
pos_gradient_segments <- segment_data %>% filter(m > 0)
n_retained_segments <- nrow(pos_gradient_segments)
pct_retained_segments <- safe_divide(100 * n_retained_segments, n_total_segments)

# Datapoints in retained segments
retained_datapoints <- point_data %>%
  semi_join(pos_gradient_segments, by = c("patient", "p6e_seg_index"))
n_retained_datapoints <- nrow(retained_datapoints)
pct_retained_datapoints <- safe_divide(100 * n_retained_datapoints, n_total_datapoints)

# Retained Patients
n_patients_retained <- n_distinct(pos_gradient_segments$patient)
pct_patients_retained <- safe_divide(100 * n_patients_retained, n_total_patients)

# Segments per retained patient
retained_segments_per_patient <- pos_gradient_segments %>%
  group_by(patient) %>%
  summarise(n_retained_segments = n(), .groups = 'drop')

# --- Build Table ---

formatted_stats <- data.frame(
  Section = character(), Metric = character(), Value = character(),
  Median = character(), IQR = character(), Range = character(),
  stringsAsFactors = FALSE
)

# Input Data Section
dpp <- format_median_iqr_simple(datapoints_per_patient, "n_datapoints", 0)
dspan <- format_median_iqr_simple(datapoints_per_patient, "data_span_days", 1)

formatted_stats <- rbind(formatted_stats, data.frame(
  Section = "Input Data",
  Metric = c("Datapoints", "Patients", "Datapoints per patient", "Data span per patient (days)"),
  Value = c(as.character(n_total_datapoints), as.character(n_total_patients), NA, NA),
  Median = c(NA, NA, dpp$Median, dspan$Median),
  IQR = c(NA, NA, dpp$IQR, dspan$IQR),
  Range = c(NA, NA, dpp$Range, dspan$Range),
  stringsAsFactors = FALSE
))

# Results Section
stats_dps <- format_cohort_median_of_patient_medians(pos_gradient_segments, "n_datapoints", 0)
stats_span <- format_cohort_median_of_patient_medians(pos_gradient_segments, "data_span", 1)
stats_spp <- format_median_iqr_simple(retained_segments_per_patient, "n_retained_segments", 0)

formatted_stats <- rbind(formatted_stats, data.frame(
  Section = "Segmentation Results",
  Metric = c(
    "Segments [a]",
    "Retained segments [b]",
    "Datapoints across retained segments [c]",
    "Datapoints per retained segment [d]",
    "Data span per retained segment (hours) [d]",
    "Retained patients [e]",
    "Retained segments per retained patient"
  ),
  Value = c(
    as.character(n_total_segments),
    paste0(n_retained_segments, " (", round(pct_retained_segments, 1), "%)"),
    paste0(n_retained_datapoints, " (", round(pct_retained_datapoints, 1), "%)"),
    NA,
    NA,
    paste0(n_patients_retained, " (", round(pct_patients_retained, 1), "%)"),
    NA
  ),
  Median = c(NA, NA, NA, stats_dps$Median, stats_span$Median, NA, stats_spp$Median),
  IQR = c(NA, NA, NA, stats_dps$IQR, stats_span$IQR, NA, stats_spp$IQR),
  Range = c(NA, NA, NA, stats_dps$Range, stats_span$Range, NA, stats_spp$Range),
  stringsAsFactors = FALSE
))

# Clean up output table for CSV
csv_out <- formatted_stats %>%
  rename(section=Section, metric=Metric, `count (percent)`=Value, median=Median, iqr=IQR, range=Range) %>%
  select(section, metric, `count (percent)`, median, iqr, range)

# Save CSV
csv_path <- file.path(output_dir, "1_stats.csv")
write.csv(csv_out, csv_path, row.names = FALSE)

cat("\n\n4_segmentation_stats.R complete.\n\n")


# ==============================================================================
# Execution Time Log
# ==============================================================================
end_time <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
time_df <- data.frame(execution_time_seconds = execution_time_seconds)
write.csv(time_df, file.path(output_dir, "__execution_time.csv"), row.names = FALSE)

# Build table legend with footnotes
legends <- list(
  list(
    target = "Table 1 (1_stats.csv)",
    footnotes = c(
      "[a] Segments from PELT + Elbow workflow (Script 2).",
      "[b] Retention requires positive gradient (p6e_m > 0).",
      "[c] Datapoints in segments with m > 0; percentage relative to total input datapoints.",
      "[d] Medians: cohort median (IQR) of patient-level medians.",
      "[e] Retained patients have at least one retained segment."
    ),
    abbreviations = "IQR, interquartile range; PELT, Pruned Exact Linear Time"
  )
)

note_path <- file.path(output_dir, "0_notes.txt")
writeLines(build_notes(legends, title = "4_segmentation_stats"), note_path)
