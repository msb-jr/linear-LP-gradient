# ==============================================================================
# Script: 5_linear_LP_models.R
# Manuscript relevance: 3.2, 3.3.i, Fig. 3
# ==============================================================================
# PURPOSE:
#   Analyze lactate–pyruvate segmentation output to understand how segment-level
#   linear fits differ from patient-level behavior, test the significance of the
#   linear intercept term (b), visualize retained segments' linear LP coefficients
#   along with their hyperbolic LPR–P implications, and generate descriptive
#   statistics of the coefficients.
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#
# OUTPUT:
#   - 1_output/5_linear_LP_models/__execution_time.csv: Runtime log
#   - 1_output/5_linear_LP_models/_1_fisher_stats.csv: Fisher-z LME results (Analysis 1)
#   - 1_output/5_linear_LP_models/_2_aic_stats.csv: AIC comparison results (Analysis 2)
#   - 1_output/5_linear_LP_models/0_notes.txt: Figure and table legends
#   - 1_output/5_linear_LP_models/3_figure.png: 4×2 panel visualization (Analysis 3)
#   - 1_output/5_linear_LP_models/4_coefficients.csv: Coefficient summary table (Analysis 4)
#   - 1_output/5_linear_LP_models/5_param_variation.csv: Parameter variation table (Analysis 5)
#
# DATA FILTERS (OVERVIEW):
#   - Analysis 1: All segments
#   - Analyses 2–4: Retained segments only (p6e_m > 0)
#   - Linear intercept stratification:
#       Type Pb: b > 0 (positive intercept)
#       Type Nb: b < 0 (negative intercept)
#       Type Zb: b = 0 (included in unstratified only)
#
# ANALYSES:
# ------------------------------------------------------------------------------
# (1) Patient vs Segment Fisher-z LME (1_fisher_stats.csv)
#     - Test whether segment-level correlations (p6e_r) exceed
#       those at the patient-level (patient_r)
#     Approach:   Weighted LME on Fisher's z-transformed Pearson correlations
#     Model:      z ~ type_factor + (1|patient), type_factor ∈ {Patient, Segment}
#     Weighting:  max(n_points - 3, 1)      # max(_,1) is defensive
#     Hypothesis: β₁ > 0 (Segment mean > Patient mean; one-sided)
#     Test:       Satterthwaite t-test
#     MCC:        None (independent hypothesis)
#
# (2) Origin-Constrained vs Unconstrained AIC LME (2_aic_stats.csv)
#     - Test whether allowing linear intercept (b) improves segment fits
#     Approach:   Fit OLS (through origin vs unconstrained), compute ΔAIC
#     Model:      delta_aic ~ 1 + (1|patient)
#     Weighting:  None (ΔAIC normalizes for segment size)
#     Hypothesis: ΔAIC > 0 (unconstrained preferred; one-sided)
#     Test:       Satterthwaite t-test
#     MCC:        None (distinct from Analysis 1)
#
# (3) Multi-Panel Visualization (3_figure.png)
#     - Visualize retained segment behavior (L vs P → LPR vs P) by intercept sign
#     Panels:     A–B: All retained data (unstratified)
#                 C–D: Type Pb only
#                 E–F: Type Nb only
#                 G–H: Pb ∩ Nb x-range, y from model predictions
#     Elements:   - Trend lines: cohort median of patient-level segment-weighted means
#                 - Spaghetti line opacity ∝ segment weight
#                 - References: y = 0 (red dotted), y = m (dashed);
#                   may not be visible depending on viewport
#     Viewport:   95% quantile cropping (2.5th–97.5th percentile)
#
# (4) Coefficient Summary Table (4_coefficients.csv)
#     - Report cohort-level summaries of m, b, R² via aggregation:
#        (4a) Segment-Level: m, b, R² per retained segment
#        (4b) Patient-Level: segment-weighted means (weights = n_points)
#        (4c) Cohort-Level: median (IQR) of patient-level metrics;
#             Cohorts: unstratified, Type Pb only, Type Nb only
#
# (5) Within-Patient Parameter Variation (5_param_variation.csv)
#     (5a) Patient-Level: for patients with ≥2 retained segments,
#          within-patient variation = max(param) − min(param) for m and b
#     (5b) Cohort-Level: median (IQR) of patient-level variation
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(colorspace)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 5_linear_LP_models.R\n")

options(dplyr.show_progress = FALSE)
emm_options(lmer.df = "satterthwaite")

# ==============================================================================
# Input / Output
# ==============================================================================
start_time <- Sys.time()

input_file <- here::here("1_output", "2_linear_segmentation", "1_results.csv")

if (!file.exists(input_file)) {
    cat("\n\nError: Results file not found at:", input_file, "\nRun 2_linear_segmentation.py first.\n")
    stop()
}

output_dir <- here::here("1_output", "5_linear_LP_models")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Fisher's z transformation: r -> z
fisher_z <- function(r) {
    # Clamp to avoid infinite values at r = ±1
    r_clamped <- pmax(pmin(r, 0.9999), -0.9999)
    0.5 * log((1 + r_clamped) / (1 - r_clamped))
}

#' Inverse Fisher's z transformation: z -> r
fisher_z_inv <- function(z) {
    (exp(2 * z) - 1) / (exp(2 * z) + 1)
}

#' Safely divide while guarding against zero or non-finite denominators
safe_divide <- function(numerator, denominator, tol = 1e-12) {
    result <- numerator / denominator
    invalid <- !is.finite(result) | !is.finite(denominator) | abs(denominator) <= tol
    result[invalid] <- NA_real_
    result
}

# ==============================================================================
# Data Loading
# ==============================================================================
df_raw <- read.csv(input_file, stringsAsFactors = FALSE)

# Create segment identifier
df_raw <- df_raw %>%
    mutate(
        segment = paste(patient, p6e_seg_index, sep = "_"),
        patient = factor(patient)
    )

# ==============================================================================
# Analysis 1: Fisher's z Comparison (Patient vs Segment Fits) - All Segments
# ==============================================================================

# Get one row per segment (for segment-level stats)
segments_all <- df_raw %>%
    group_by(segment) %>%
    summarise(
        patient = first(patient),
        p6e_r = first(p6e_r),
        patient_r = first(patient_r),
        n_points_segment = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(p6e_r) & !is.na(patient_r))

# Get patient-level stats (one row per patient)
patient_stats <- df_raw %>%
    group_by(patient) %>%
    summarise(
        patient_r = first(patient_r),
        n_points_patient = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(patient_r))

# Build long-format data for LME
# Each patient contributes:
#   - 1 row for patient-level fit (Type = 0)
#   - K rows for segment-level fits (Type = 1)

# Patient-level rows
patient_rows <- patient_stats %>%
    mutate(
        type = 0,  # 0 = patient-level
        r = patient_r,
        n_points = n_points_patient,
        z = fisher_z(r),
        weight = pmax(n_points - 3, 1)  # Fisher-z variance weight (n_points - 3), floored at 1
    ) %>%
    select(patient, type, r, z, n_points, weight)

# Segment-level rows
segment_rows <- segments_all %>%
    mutate(
        type = 1,  # 1 = segment-level
        r = p6e_r,
        n_points = n_points_segment,
        z = fisher_z(r),
        weight = pmax(n_points - 3, 1)  # Same Fisher-z weight for segment-level rows
    ) %>%
    select(patient, type, r, z, n_points, weight)

# Combine
lme_data <- bind_rows(patient_rows, segment_rows) %>%
    mutate(
        patient = factor(patient),
        type_factor = factor(type, levels = c(0, 1), labels = c("Patient", "Segment"))
    )

# Fit weighted LME: z ~ type_factor + (1 | patient)
# Type = "Patient" is reference; Type = "Segment" is the effect of interest
# Singular fit (ICC ≈ 0) is expected from variance at segment-level
# (within-patient) rather than patient-level, indicating linear LP dynamics
# are governed by local physiology.
fisher_lme <- tryCatch({
    lmer(z ~ type_factor + (1 | patient), data = lme_data, weights = weight, REML = FALSE,
         control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
}, error = function(e) {
    cat("\n\nError: LME fitting failed:", e$message, "\n")
    NULL
})

fisher_lme_results <- NULL
if (!is.null(fisher_lme)) {
    # Extract fixed effects
    fe <- fixef(fisher_lme)
    fe_summary <- summary(fisher_lme)$coefficients
    
    # Confidence intervals for fixed effects
    fe_ci <- tryCatch(confint(fisher_lme, parm = "beta_", method = "Wald"), 
                      error = function(e) NULL)
    
    # The factor coefficient is named "type_factorSegment" (Segment vs Patient reference)
    coef_name <- "type_factorSegment"
    
    # Calculate one-sided p-value for segment > patient effect (Satterthwaite t-test)
    t_val <- fe_summary[coef_name, "t value"]
    p_two_sided <- fe_summary[coef_name, "Pr(>|t|)"]
    p_value_one_sided <- ifelse(t_val > 0, p_two_sided / 2, 1 - p_two_sided / 2)
    
    # Estimated Marginal Means
    emm <- emmeans(fisher_lme, ~ type_factor)
    emm_df <- as.data.frame(emm)
    
    # Back-transform EMMs from z to r
    emm_df$r_estimate <- fisher_z_inv(emm_df$emmean)
    emm_df$r_lower <- fisher_z_inv(emm_df$lower.CL)
    emm_df$r_upper <- fisher_z_inv(emm_df$upper.CL)
    
    # Extract ICC from random effects
    vc <- as.data.frame(VarCorr(fisher_lme))
    var_patient <- vc$vcov[vc$grp == "patient"]
    var_residual <- vc$vcov[vc$grp == "Residual"]
    icc <- safe_divide(var_patient, var_patient + var_residual)
    
    fisher_lme_results <- list(
        b0 = fe["(Intercept)"],
        b1 = fe[coef_name],
        b1_se = fe_summary[coef_name, "Std. Error"],
        b1_ci_lower = if(!is.null(fe_ci)) fe_ci[coef_name, 1] else NA,
        b1_ci_upper = if(!is.null(fe_ci)) fe_ci[coef_name, 2] else NA,
        p_value_one_sided = p_value_one_sided,
        emm_patient_z = emm_df$emmean[emm_df$type_factor == "Patient"],
        emm_segment_z = emm_df$emmean[emm_df$type_factor == "Segment"],
        emm_patient_r = emm_df$r_estimate[emm_df$type_factor == "Patient"],
        emm_patient_r_ci = c(emm_df$r_lower[emm_df$type_factor == "Patient"],
                            emm_df$r_upper[emm_df$type_factor == "Patient"]),
        emm_segment_r = emm_df$r_estimate[emm_df$type_factor == "Segment"],
        emm_segment_r_ci = c(emm_df$r_lower[emm_df$type_factor == "Segment"],
                            emm_df$r_upper[emm_df$type_factor == "Segment"]),
        icc = icc,
        n_patients = n_distinct(lme_data$patient),
        n_segments = sum(lme_data$type == 1)
    )
}

# Save Fisher's z comparison results
if (!is.null(fisher_lme_results)) {
    fisher_summary <- data.frame(
        metric = c(
            "B1 (Segment effect on Fisher's z)",
            "B1 Standard Error",
            "B1 95% CI Lower",
            "B1 95% CI Upper",
            "p-value (one-sided, Segment > Patient)",
            "EMM Patient-level r",
            "EMM Patient-level r 95% CI Lower",
            "EMM Patient-level r 95% CI Upper",
            "EMM Segment-level r",
            "EMM Segment-level r 95% CI Lower",
            "EMM Segment-level r 95% CI Upper",
            "ICC (Patient)",
            "N Patients",
            "N Segments"
        ),
        value = c(
            fisher_lme_results$b1,
            fisher_lme_results$b1_se,
            fisher_lme_results$b1_ci_lower,
            fisher_lme_results$b1_ci_upper,
            fisher_lme_results$p_value_one_sided,
            fisher_lme_results$emm_patient_r,
            fisher_lme_results$emm_patient_r_ci[1],
            fisher_lme_results$emm_patient_r_ci[2],
            fisher_lme_results$emm_segment_r,
            fisher_lme_results$emm_segment_r_ci[1],
            fisher_lme_results$emm_segment_r_ci[2],
            fisher_lme_results$icc,
            fisher_lme_results$n_patients,
            fisher_lme_results$n_segments
        ),
        stringsAsFactors = FALSE
    )
    write.csv(fisher_summary, file.path(output_dir, "_1_fisher_stats.csv"), row.names = FALSE)
}

# ==============================================================================
# Analysis 2: AIC Comparison (Constrained vs Unconstrained) - Retained Segments
# ==============================================================================

# Retained segments: those with positive gradient (p6e_m > 0)
df_retained <- df_raw %>%
    filter(p6e_m > 0)  # All weighting/aggregation in Analyses 2–4 operates only on this subset

n_retained_segments <- n_distinct(df_retained$segment)
n_retained_patients <- n_distinct(df_retained$patient)

# For each retained segment, fit both models and compute ∆AIC
aic_results <- df_retained %>%
    group_by(patient, segment) %>%
    group_modify(~{
        dat <- .x
        n_pts <- nrow(dat)
        
        # Need at least 3 points for meaningful comparison
        if (n_pts < 3) {
            return(tibble(
                n_points = n_pts,
                aic_constrained = NA_real_,
                aic_unconstrained = NA_real_,
                delta_aic = NA_real_,
                m_constrained = NA_real_,
                m_unconstrained = NA_real_,
                b_unconstrained = NA_real_
            ))
        }
        
        # Model A: Constrained (lactate = m * pyruvate, no intercept)
        fit_constrained <- tryCatch(
            lm(lactate ~ pyruvate - 1, data = dat),
            error = function(e) NULL
        )
        
        # Model B: Unconstrained (lactate = m * pyruvate + b)
        fit_unconstrained <- tryCatch(
            lm(lactate ~ pyruvate, data = dat),
            error = function(e) NULL
        )
        
        if (!is.null(fit_constrained) && !is.null(fit_unconstrained)) {
            aic_c <- AIC(fit_constrained)
            aic_u <- AIC(fit_unconstrained)
            
            tibble(
                n_points = n_pts,
                aic_constrained = aic_c,
                aic_unconstrained = aic_u,
                delta_aic = aic_c - aic_u,  # Positive = unconstrained is better
                m_constrained = coef(fit_constrained)[1],
                m_unconstrained = coef(fit_unconstrained)["pyruvate"],
                b_unconstrained = coef(fit_unconstrained)["(Intercept)"]
            )
        } else {
            tibble(
                n_points = n_pts,
                aic_constrained = NA_real_,
                aic_unconstrained = NA_real_,
                delta_aic = NA_real_,
                m_constrained = NA_real_,
                m_unconstrained = NA_real_,
                b_unconstrained = NA_real_
            )
        }
    }) %>%
    ungroup() %>%
    filter(!is.na(delta_aic))

# Fit LME: ∆AIC ~ 1 + (1 | patient)
# Singular fit (ICC ≈ 0) is expected because of ∆AIC variation at segment-level
# All retained segments contribute one row (no additional weights beyond that)
aic_lme <- tryCatch({
    lmer(delta_aic ~ 1 + (1 | patient), data = aic_results, REML = FALSE,
         control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
}, error = function(e) {
    cat("\n\nError: AIC LME fitting failed:", e$message, "\n")
    NULL
})

aic_lme_results <- NULL
if (!is.null(aic_lme)) {
    # Extract fixed effect (intercept = typical ∆AIC)
    fe <- fixef(aic_lme)
    fe_summary <- summary(aic_lme)$coefficients
    
    # Confidence interval
    fe_ci <- tryCatch(confint(aic_lme, parm = "beta_", method = "Wald"), 
                      error = function(e) NULL)
    
    # One-sided p-value for intercept (testing if ∆AIC > 0) via Satterthwaite t-test
    t_val <- fe_summary["(Intercept)", "t value"]
    p_two_sided <- fe_summary["(Intercept)", "Pr(>|t|)"]
    p_value_one_sided <- ifelse(t_val > 0, p_two_sided / 2, 1 - p_two_sided / 2)
    
    # Percentage of segments with ∆AIC > 2
    pct_aic_gt2 <- 100 * mean(aic_results$delta_aic > 2, na.rm = TRUE)
    
    # Extract ICC from random effects
    vc_aic <- as.data.frame(VarCorr(aic_lme))
    var_patient_aic <- vc_aic$vcov[vc_aic$grp == "patient"]
    var_residual_aic <- vc_aic$vcov[vc_aic$grp == "Residual"]
    icc_aic <- safe_divide(var_patient_aic, var_patient_aic + var_residual_aic)
    
    aic_lme_results <- list(
        fixed_effect = fe["(Intercept)"],
        fixed_effect_se = fe_summary["(Intercept)", "Std. Error"],
        ci_lower = if(!is.null(fe_ci)) fe_ci["(Intercept)", 1] else NA,
        ci_upper = if(!is.null(fe_ci)) fe_ci["(Intercept)", 2] else NA,
        p_value_one_sided = p_value_one_sided,
        pct_delta_aic_gt2 = pct_aic_gt2,
        n_segments = nrow(aic_results),
        n_patients = n_distinct(aic_results$patient),
        median_delta_aic = median(aic_results$delta_aic),
        iqr_delta_aic = IQR(aic_results$delta_aic),
        icc = icc_aic
    )
}

# Save AIC comparison results
if (!is.null(aic_lme_results)) {
    aic_summary <- data.frame(
        metric = c(
            "Fixed Effect (Typical ∆AIC)",
            "Fixed Effect SE",
            "95% CI Lower",
            "95% CI Upper",
            "p-value (one-sided, ∆AIC > 0)",
            "% Segments with ∆AIC > 2",
            "Median ∆AIC",
            "IQR ∆AIC",
            "ICC (Patient)",
            "N Segments",
            "N Patients"
        ),
        value = c(
            aic_lme_results$fixed_effect,
            aic_lme_results$fixed_effect_se,
            aic_lme_results$ci_lower,
            aic_lme_results$ci_upper,
            aic_lme_results$p_value_one_sided,
            aic_lme_results$pct_delta_aic_gt2,
            aic_lme_results$median_delta_aic,
            aic_lme_results$iqr_delta_aic,
            aic_lme_results$icc,
            aic_lme_results$n_segments,
            aic_lme_results$n_patients
        ),
        stringsAsFactors = FALSE
    )
    write.csv(aic_summary, file.path(output_dir, "_2_aic_stats.csv"), row.names = FALSE)
}

# ==============================================================================
# Analysis 3: Visualization (Retained Segments Only) - 4x2 Multi-Panel Figure
# ==============================================================================

# Prepare retained data for visualization with intercept-based stratification
df_vis <- df_retained %>%
    mutate(
        segment = factor(segment),
        # Stratify by intercept sign:
        # - neg: b < 0 (Type Nb)
        # - pos: b > 0 (Type Pb)
        # - zero: b = 0 (Type Zb) - excluded from stratified groups
        grp = case_when(
            p6e_b < 0 ~ "neg",
            p6e_b > 0 ~ "pos",
            TRUE ~ "zero"
        ),
        p6e_r2 = p6e_r^2,
        lpr = safe_divide(lactate, pyruvate)  # Compute LPR for quantile limits
    )

# --- Compute 95% viewport limits (middle 95% of data) for panels A–F ---
# Unstratified limits (for A, B)
unstrat_xlim <- quantile(df_vis$pyruvate, c(0.025, 0.975), na.rm = TRUE) * 1000  # μM
unstrat_y1_lim <- quantile(df_vis$lactate, c(0.025, 0.975), na.rm = TRUE)
unstrat_y2_lim <- quantile(df_vis$lpr, c(0.025, 0.975), na.rm = TRUE)

# Stratified limits (for C, D, E, F)
pos_data <- df_vis %>% filter(grp == "pos")
neg_data <- df_vis %>% filter(grp == "neg")

# Segment and patient counts per stratum (for Table 1 caption)
n_segments_pos <- n_distinct(pos_data$segment)
n_patients_pos <- n_distinct(pos_data$patient)
n_segments_neg <- n_distinct(neg_data$segment)
n_patients_neg <- n_distinct(neg_data$patient)

viewport_limits <- list()
if (nrow(pos_data) > 0) {
    viewport_limits[["pos"]] <- list(
        xlim = quantile(pos_data$pyruvate, c(0.025, 0.975), na.rm = TRUE) * 1000,
        y1_lim = quantile(pos_data$lactate, c(0.025, 0.975), na.rm = TRUE),
        y2_lim = quantile(pos_data$lpr, c(0.025, 0.975), na.rm = TRUE)
    )
}
if (nrow(neg_data) > 0) {
    viewport_limits[["neg"]] <- list(
        xlim = quantile(neg_data$pyruvate, c(0.025, 0.975), na.rm = TRUE) * 1000,
        y1_lim = quantile(neg_data$lactate, c(0.025, 0.975), na.rm = TRUE),
        y2_lim = quantile(neg_data$lpr, c(0.025, 0.975), na.rm = TRUE)
    )
}

# Intersection of group-specific 95% x-limits (for panels G–H)
# Returns the x-range where BOTH cohorts have credible data coverage
intersect_group_xlimits <- function(groups) {
    limits <- lapply(groups, function(g) {
        if (!is.null(viewport_limits[[g]])) viewport_limits[[g]][["xlim"]] else NULL
    })
    limits <- Filter(Negate(is.null), limits)
    if (length(limits) == 0) return(NULL)
    # Intersection: max of lower bounds, min of upper bounds
    lower <- max(sapply(limits, `[`, 1), na.rm = TRUE)
    upper <- min(sapply(limits, `[`, 2), na.rm = TRUE)
    if (lower >= upper) return(NULL)  # No valid intersection
    c(lower, upper)
}

# Segment-level datapoint counts and patient-level totals (used for weights + alpha)
segment_sizes <- df_retained %>%
    count(patient, segment, name = "segment_n_points")

patient_sizes <- segment_sizes %>%
    group_by(patient) %>%
    summarise(patient_n_points = sum(segment_n_points), .groups = "drop")

segment_weights <- segment_sizes %>%
    left_join(patient_sizes, by = "patient") %>%
    mutate(segment_alpha = safe_divide(segment_n_points, patient_n_points))

segment_counts_simple <- segment_sizes %>%
    select(segment, n_points = segment_n_points)

# Get segment-level coefficients (one row per segment)
segment_coeffs <- df_vis %>%
    distinct(segment, .keep_all = TRUE) %>%
    select(patient, segment, grp, p6e_m, p6e_b, p6e_r, p6e_r2) %>%
    left_join(segment_weights %>% select(segment, segment_alpha), by = "segment")

# --- Aggregation: Patient-level weighted means, then cohort-level medians ---
# Patient-level weighted means for each group (weights = segment n_points)
patient_weighted_stratified <- segment_coeffs %>%
    left_join(segment_counts_simple, by = "segment") %>%
    group_by(grp, patient) %>%
    summarise(
        wmean_m = weighted.mean(p6e_m, w = n_points, na.rm = TRUE),  # Heavier weight for larger segments
        wmean_b = weighted.mean(p6e_b, w = n_points, na.rm = TRUE),
        wmean_r2 = weighted.mean(p6e_r2, w = n_points, na.rm = TRUE),
        .groups = "drop"
    )

# Cohort-level medians for each group (unweighted; each patient counts once)
# Quartiles computed before medians to keep quantile summaries independent
segment_summary_stats <- patient_weighted_stratified %>%
    group_by(grp) %>%
    summarise(
        q1_m = quantile(wmean_m, 0.25, na.rm = TRUE),
        q3_m = quantile(wmean_m, 0.75, na.rm = TRUE),
        q1_b = quantile(wmean_b, 0.25, na.rm = TRUE),
        q3_b = quantile(wmean_b, 0.75, na.rm = TRUE),
        q1_r2 = quantile(wmean_r2, 0.25, na.rm = TRUE),
        q3_r2 = quantile(wmean_r2, 0.75, na.rm = TRUE),
        median_m = median(wmean_m, na.rm = TRUE),
        median_b = median(wmean_b, na.rm = TRUE),
        median_r2 = median(wmean_r2, na.rm = TRUE),
        .groups = "drop"
    )

# --- Patient-level weighted means (unstratified) ---
patient_weighted_unstrat <- segment_coeffs %>%
    left_join(segment_counts_simple, by = "segment") %>%
    group_by(patient) %>%
    summarise(
        wmean_m = weighted.mean(p6e_m, w = n_points, na.rm = TRUE),  # Within-patient weights = segment size
        wmean_b = weighted.mean(p6e_b, w = n_points, na.rm = TRUE),
        wmean_r2 = weighted.mean(p6e_r2, w = n_points, na.rm = TRUE),
        .groups = "drop"
    )

# Quartiles computed before medians to keep quantile summaries independent
cohort_medians <- patient_weighted_unstrat %>%
    summarise(  # No additional weights here; treat patients equally for cohort summary
        q1_m = quantile(wmean_m, 0.25, na.rm = TRUE),
        q3_m = quantile(wmean_m, 0.75, na.rm = TRUE),
        q1_b = quantile(wmean_b, 0.25, na.rm = TRUE),
        q3_b = quantile(wmean_b, 0.75, na.rm = TRUE),
        q1_r2 = quantile(wmean_r2, 0.25, na.rm = TRUE),
        q3_r2 = quantile(wmean_r2, 0.75, na.rm = TRUE),
        m = median(wmean_m, na.rm = TRUE),
        b = median(wmean_b, na.rm = TRUE),
        r2 = median(wmean_r2, na.rm = TRUE)
    )

# --- Create Stratified Plots (Rows 2 & 3: Type Pb and Type Nb) ---
make_placeholder_plot <- function(message) {
    ggplot() +
        annotate("text", x = 0, y = 0, label = message, hjust = 0.5, vjust = 0.5, size = 3) +
        xlim(-1, 1) +
        ylim(-1, 1) +
        theme_void() +
        theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
}

plots_y1 <- list()
plots_y2 <- list()

for (g in c("pos", "neg")) {
    g_stats <- filter(segment_summary_stats, grp == g)
    cohort_df <- df_vis[df_vis$grp == g, ]
    
    if (nrow(cohort_df) == 0 || nrow(g_stats) == 0) next
    
    # Grid for trend line
    grid <- data.frame(
        x_val = seq(min(cohort_df$pyruvate, na.rm = TRUE),
                    max(cohort_df$pyruvate, na.rm = TRUE), length = 200)
    )
    grid$x_val_um <- grid$x_val * 1000
    grid$fit_y1 <- g_stats$median_b + g_stats$median_m * grid$x_val
    grid$fit_y2 <- g_stats$median_m + safe_divide(g_stats$median_b, grid$x_val)
    
    # Spaghetti lines for this group
    seg_ids <- unique(as.character(cohort_df$segment))
    spaghetti_df <- do.call(rbind, lapply(seg_ids, function(seg_id) {
        seg_data <- df_vis %>% filter(segment == seg_id)
        seg_coef <- segment_coeffs %>% filter(segment == seg_id)
        if (nrow(seg_data) < 2 || nrow(seg_coef) == 0) return(NULL)
        
        seg_grid <- data.frame(
            x_val = seq(min(seg_data$pyruvate, na.rm = TRUE),
                        max(seg_data$pyruvate, na.rm = TRUE), length = 100)
        )
        seg_grid$x_val_um <- seg_grid$x_val * 1000
        seg_grid$segment <- seg_id
        seg_grid$y1_pred <- seg_coef$p6e_b + seg_coef$p6e_m * seg_grid$x_val
        seg_grid$y2_pred <- seg_coef$p6e_m + safe_divide(seg_coef$p6e_b, seg_grid$x_val)
        seg_grid$segment_alpha <- seg_coef$segment_alpha
        seg_grid
    }))
    
    line_color <- if (g == "pos") darken("#44AA99", 0.2) else darken("#DDCC77", 0.2)
    
    # Get 95% viewport limits for this group
    g_viewport <- viewport_limits[[g]]
    
    p_y1 <- ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
        geom_line(data = spaghetti_df, aes(x = x_val_um, y = y1_pred, group = segment, alpha = segment_alpha),
                  colour = "grey75", linewidth = 0.2) +
        geom_line(data = grid, aes(x = x_val_um, y = fit_y1), colour = line_color, linewidth = 1.4) +
        labs(x = NULL, y = "Lactate (mM)") +
        coord_cartesian(xlim = g_viewport$xlim, ylim = g_viewport$y1_lim) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        theme_minimal() +
        theme(
            text = element_text(size = 8, colour = "black"),
            axis.text = element_text(size = 6, colour = "black"),
            axis.title = element_text(size = 7, colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        ) +
        scale_alpha_identity()
    plots_y1[[g]] <- p_y1
    
    p_y2 <- ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
        geom_line(data = spaghetti_df, aes(x = x_val_um, y = y2_pred, group = segment, alpha = segment_alpha),
                  colour = "grey75", linewidth = 0.2) +
        geom_line(data = grid, aes(x = x_val_um, y = fit_y2), colour = line_color, linewidth = 1.4) +
        geom_hline(yintercept = g_stats$median_m, linetype = "dashed", 
                   color = line_color, linewidth = 0.6, alpha = 0.8) +
        labs(x = NULL, y = "LPR") +
        coord_cartesian(xlim = g_viewport$xlim, ylim = g_viewport$y2_lim) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        theme_minimal() +
        theme(
            text = element_text(size = 8, colour = "black"),
            axis.text = element_text(size = 6, colour = "black"),
            axis.title = element_text(size = 7, colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        ) +
        scale_alpha_identity()
    plots_y2[[g]] <- p_y2
}

# Ensure stratified plots exist even if no segments are available
if (is.null(plots_y1[["pos"]]) || is.null(plots_y2[["pos"]])) {
    placeholder <- make_placeholder_plot("No Type Pb segments available")
    if (is.null(plots_y1[["pos"]])) {
        plots_y1[["pos"]] <- placeholder
    }
    if (is.null(plots_y2[["pos"]])) {
        plots_y2[["pos"]] <- placeholder
    }
}

if (is.null(plots_y1[["neg"]]) || is.null(plots_y2[["neg"]])) {
    placeholder <- make_placeholder_plot("No Type Nb segments available")
    if (is.null(plots_y1[["neg"]])) {
        plots_y1[["neg"]] <- placeholder
    }
    if (is.null(plots_y2[["neg"]])) {
        plots_y2[["neg"]] <- placeholder
    }
}

# --- Create Unstratified Plots (Row 1: A, B) ---
unstrat_grid <- data.frame(
    x_val = seq(min(df_vis$pyruvate, na.rm = TRUE),
                max(df_vis$pyruvate, na.rm = TRUE), length = 200)
)
unstrat_grid$x_val_um <- unstrat_grid$x_val * 1000
unstrat_grid$fit_y1 <- cohort_medians$b + cohort_medians$m * unstrat_grid$x_val
unstrat_grid$fit_y2 <- cohort_medians$m + safe_divide(cohort_medians$b, unstrat_grid$x_val)

# Spaghetti for unstratified (colored by group)
unstrat_spaghetti <- do.call(rbind, lapply(unique(as.character(segment_coeffs$segment)), function(seg_id) {
    seg_data <- df_vis %>% filter(segment == seg_id)
    seg_coef <- segment_coeffs %>% filter(segment == seg_id)
    if (nrow(seg_data) < 2 || nrow(seg_coef) == 0) return(NULL)
    
    seg_grid <- data.frame(
        x_val = seq(min(seg_data$pyruvate, na.rm = TRUE),
                    max(seg_data$pyruvate, na.rm = TRUE), length = 100)
    )
    seg_grid$x_val_um <- seg_grid$x_val * 1000
    seg_grid$segment <- seg_id
    seg_grid$grp <- seg_coef$grp
    seg_grid$y1_pred <- seg_coef$p6e_b + seg_coef$p6e_m * seg_grid$x_val
    seg_grid$y2_pred <- seg_coef$p6e_m + safe_divide(seg_coef$p6e_b, seg_grid$x_val)
    seg_grid$segment_alpha <- seg_coef$segment_alpha
    seg_grid
}))

# Guard: if unstrat_spaghetti is NULL or empty, use placeholders
if (is.null(unstrat_spaghetti) || nrow(unstrat_spaghetti) == 0) {
    p_y1_unstrat <- make_placeholder_plot("No retained segments available")
    p_y2_unstrat <- make_placeholder_plot("No retained segments available")
} else {
    p_y1_unstrat <- ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
        geom_line(data = filter(unstrat_spaghetti, grp == "zero"), 
                  aes(x = x_val_um, y = y1_pred, group = segment, alpha = segment_alpha),
                  colour = "#888888", linewidth = 0.2) +
        geom_line(data = filter(unstrat_spaghetti, grp == "pos"), 
                  aes(x = x_val_um, y = y1_pred, group = segment, alpha = segment_alpha),
                  colour = "#44AA99", linewidth = 0.2) +
        geom_line(data = filter(unstrat_spaghetti, grp == "neg"), 
                  aes(x = x_val_um, y = y1_pred, group = segment, alpha = segment_alpha),
                  colour = "#DDCC77", linewidth = 0.2) +
        geom_line(data = unstrat_grid, aes(x = x_val_um, y = fit_y1), colour = darken("#996633", 0.2), linewidth = 1.4) +
        labs(x = NULL, y = "Lactate (mM)") +
        coord_cartesian(xlim = unstrat_xlim, ylim = unstrat_y1_lim) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        ) +
        scale_alpha_identity()

    p_y2_unstrat <- ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
        geom_line(data = filter(unstrat_spaghetti, grp == "zero"), 
                  aes(x = x_val_um, y = y2_pred, group = segment, alpha = segment_alpha),
                  colour = "#888888", linewidth = 0.2) +
        geom_line(data = filter(unstrat_spaghetti, grp == "pos"), 
                  aes(x = x_val_um, y = y2_pred, group = segment, alpha = segment_alpha),
                  colour = "#44AA99", linewidth = 0.2) +
        geom_line(data = filter(unstrat_spaghetti, grp == "neg"), 
                  aes(x = x_val_um, y = y2_pred, group = segment, alpha = segment_alpha),
                  colour = "#DDCC77", linewidth = 0.2) +
        geom_line(data = unstrat_grid, aes(x = x_val_um, y = fit_y2), colour = darken("#996633", 0.2), linewidth = 1.4) +
        geom_hline(yintercept = cohort_medians$m, linetype = "dashed", color = darken("#996633", 0.2), linewidth = 0.6, alpha = 0.8) +
        labs(x = NULL, y = "LPR") +
        coord_cartesian(xlim = unstrat_xlim, ylim = unstrat_y2_lim) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        ) +
        scale_alpha_identity()
}

# --- Create Combined Plots (Row 4: G, H) ---
# Compute SHARED x-axis as intersection of group-specific 95% limits
# Intersection = x-range where BOTH cohorts have credible data coverage
combined_xlim <- intersect_group_xlimits(c("pos", "neg"))

# If no valid intersection exists, use placeholders
if (is.null(combined_xlim)) {
    p_y1_combined <- make_placeholder_plot("No valid x-axis intersection for combined LP panel")
    p_y2_combined <- make_placeholder_plot("No valid x-axis intersection for combined LPR panel")
} else {
    # y-limits: derived from model predictions at x-endpoints
    # Convert x from μM to mM for model evaluation
    x_endpoints_mM <- combined_xlim / 1000
    y1_endpoints <- c()
    y2_endpoints <- c()

    for (g in c("pos", "neg")) {
        g_stats <- filter(segment_summary_stats, grp == g)
        if (nrow(g_stats) == 0) next
        
        # Panel G: Lactate = m * pyruvate + b
        y1_at_xmin <- g_stats$median_b + g_stats$median_m * x_endpoints_mM[1]
        y1_at_xmax <- g_stats$median_b + g_stats$median_m * x_endpoints_mM[2]
        y1_endpoints <- c(y1_endpoints, y1_at_xmin, y1_at_xmax)
        
        # Panel H: LPR = m + b / pyruvate
        y2_at_xmin <- g_stats$median_m + safe_divide(g_stats$median_b, x_endpoints_mM[1])
        y2_at_xmax <- g_stats$median_m + safe_divide(g_stats$median_b, x_endpoints_mM[2])
        y2_endpoints <- c(y2_endpoints, y2_at_xmin, y2_at_xmax)
    }

    # Set y-limits from model predictions (default ggplot2 expansion will add padding)
    if (length(y1_endpoints) > 0) {
        combined_y1_lim <- range(y1_endpoints, na.rm = TRUE)
    } else {
        combined_y1_lim <- c(0, 8)
    }

    # Include asymptote lines in Panel H y-limits
    stats_pos_plot <- segment_summary_stats %>% filter(grp == "pos")
    stats_neg_plot <- segment_summary_stats %>% filter(grp == "neg")
    m_pos <- if (nrow(stats_pos_plot) > 0) stats_pos_plot$median_m[1] else NA_real_
    m_neg <- if (nrow(stats_neg_plot) > 0) stats_neg_plot$median_m[1] else NA_real_

    if (length(y2_endpoints) > 0) {
        y2_all <- c(y2_endpoints, m_pos, m_neg)  # Include asymptote lines
        combined_y2_lim <- range(y2_all, na.rm = TRUE)
    } else {
        combined_y2_lim <- c(0, 50)
    }

    # Generate grids over INTERSECTION range (not each cohort's full domain)
    # Y-limits derived from curve endpoints; ggplot2 default expansion adds padding
    combined_grids <- lapply(c("pos", "neg"), function(g) {
        g_stats <- filter(segment_summary_stats, grp == g)
        if (nrow(g_stats) == 0) return(NULL)
        
        # Grid spans the intersection x-range (in mM)
        grid <- data.frame(
            x_val = seq(x_endpoints_mM[1], x_endpoints_mM[2], length = 200)
        )
        grid$x_val_um <- grid$x_val * 1000
        grid$fit_y1 <- g_stats$median_b + g_stats$median_m * grid$x_val
        grid$fit_y2 <- g_stats$median_m + safe_divide(g_stats$median_b, grid$x_val)
        grid$group <- ifelse(g == "pos", "Type 'Pb' Segments", "Type 'Nb' Segments")
        grid
    })
    combined_grid_list <- Filter(Negate(is.null), combined_grids)
    combined_grid <- if (length(combined_grid_list) > 0) bind_rows(combined_grid_list) else NULL

    if (!is.null(combined_grid) && nrow(combined_grid) > 0) {
        p_y1_combined <- ggplot(combined_grid, aes(x = x_val_um, group = group, colour = group)) +
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
            geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
            geom_line(aes(y = fit_y1), linewidth = 1.1) +
            scale_color_manual(values = c("Type 'Pb' Segments" = darken("#44AA99", 0.2), 
                                          "Type 'Nb' Segments" = darken("#DDCC77", 0.2))) +
            labs(x = "Pyruvate (μM)", y = "Lactate (mM)") +
            coord_cartesian(xlim = combined_xlim, ylim = combined_y1_lim) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
            theme_minimal() +
            theme(
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)
            )
        
        p_y2_combined <- ggplot(combined_grid, aes(x = x_val_um, group = group, colour = group)) +
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "red", alpha = 0.1) +
            geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
            geom_line(aes(y = fit_y2), linewidth = 1.1) +
            # Asymptote lines clipped to x-range
            annotate("segment", x = combined_xlim[1], xend = combined_xlim[2], y = m_pos, yend = m_pos,
                     linetype = "dashed", color = darken("#44AA99", 0.2), linewidth = 0.6, alpha = 0.8) +
            annotate("segment", x = combined_xlim[1], xend = combined_xlim[2], y = m_neg, yend = m_neg,
                     linetype = "dashed", color = darken("#DDCC77", 0.2), linewidth = 0.6, alpha = 0.8) +
            scale_color_manual(values = c("Type 'Pb' Segments" = darken("#44AA99", 0.2), 
                                          "Type 'Nb' Segments" = darken("#DDCC77", 0.2))) +
            labs(x = "Pyruvate (μM)", y = "LPR") +
            coord_cartesian(xlim = combined_xlim, ylim = combined_y2_lim) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
            theme_minimal() +
            theme(
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)
            )
    } else {
        p_y1_combined <- make_placeholder_plot("Insufficient data for combined LP panel")
        p_y2_combined <- make_placeholder_plot("Insufficient data for combined LPR panel")
    }
}

# --- Add Panel Labels and Assemble 4x2 Figure ---
p_y1_unstrat <- p_y1_unstrat + labs(tag = "A")
p_y2_unstrat <- p_y2_unstrat + labs(tag = "B")
plots_y1[["pos"]] <- plots_y1[["pos"]] + labs(tag = "C")
plots_y2[["pos"]] <- plots_y2[["pos"]] + labs(tag = "D")
plots_y1[["neg"]] <- plots_y1[["neg"]] + labs(tag = "E")
plots_y2[["neg"]] <- plots_y2[["neg"]] + labs(tag = "F")
p_y1_combined <- p_y1_combined + labs(tag = "G")
p_y2_combined <- p_y2_combined + labs(tag = "H")

figure <- (p_y1_unstrat | p_y2_unstrat) /
          (plots_y1[["pos"]] | plots_y2[["pos"]]) /
          (plots_y1[["neg"]] | plots_y2[["neg"]]) /
          (p_y1_combined | p_y2_combined) &
    theme(
        text = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 7, colour = "black"),
        plot.tag = element_text(size = 9, face = "bold", colour = "black", hjust = 1, vjust = -1),
        plot.tag.position = c(0.02, 0.99),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 5, b = 2, l = 5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    )

figure_path <- file.path(output_dir, "3_figure.png")
if (isTRUE(capabilities("cairo"))) {
    ggsave(
        filename = figure_path,
        plot = figure, width = 185, height = 193, units = "mm", dpi = 300, type = "cairo"
    )
} else {
    warning("Cairo backend unavailable; saving 3_figure.png with the default device.")
    ggsave(
        filename = figure_path,
        plot = figure, width = 185, height = 193, units = "mm", dpi = 300
    )
}

# ==============================================================================
# Analysis 4: Summary Table (Retained Segments Only)
# ==============================================================================

# Helper to format median (IQR)
format_med_iqr <- function(med, q1, q3) {
    sprintf("%.2f (%.2f – %.2f)", med, q1, q3)
}

# Build summary table (unstratified + stratified rows)
summary_rows <- list()

# Unstratified
summary_rows[["Unstratified"]] <- data.frame(
    group = "Unstratified",
    m = format_med_iqr(cohort_medians$m, cohort_medians$q1_m, cohort_medians$q3_m),
    b = format_med_iqr(cohort_medians$b, cohort_medians$q1_b, cohort_medians$q3_b),
    r2 = format_med_iqr(cohort_medians$r2, cohort_medians$q1_r2, cohort_medians$q3_r2),
    stringsAsFactors = FALSE
)

# Type Pb (positive intercept)
stats_pos <- segment_summary_stats %>% filter(grp == "pos")
if (nrow(stats_pos) > 0) {
    summary_rows[["Type Pb"]] <- data.frame(
        group = "Type Pb (b > 0)",
        m = format_med_iqr(stats_pos$median_m, stats_pos$q1_m, stats_pos$q3_m),
        b = format_med_iqr(stats_pos$median_b, stats_pos$q1_b, stats_pos$q3_b),
        r2 = format_med_iqr(stats_pos$median_r2, stats_pos$q1_r2, stats_pos$q3_r2),
        stringsAsFactors = FALSE
    )
}

# Type Nb (negative intercept)
stats_neg <- segment_summary_stats %>% filter(grp == "neg")
if (nrow(stats_neg) > 0) {
    summary_rows[["Type Nb"]] <- data.frame(
        group = "Type Nb (b < 0)",
        m = format_med_iqr(stats_neg$median_m, stats_neg$q1_m, stats_neg$q3_m),
        b = format_med_iqr(stats_neg$median_b, stats_neg$q1_b, stats_neg$q3_b),
        r2 = format_med_iqr(stats_neg$median_r2, stats_neg$q1_r2, stats_neg$q3_r2),
        stringsAsFactors = FALSE
    )
}

summary_table <- bind_rows(summary_rows)
names(summary_table) <- c("Group", "Gradient (m)", "Intercept (b, mM)", "R²")

write.csv(summary_table, file.path(output_dir, "4_coefficients.csv"), row.names = FALSE)

# ==============================================================================
# Analysis 5: Within-Patient Parameter Variation
# ==============================================================================
# For patients with multiple segments, calculate how much m and b vary within
# each patient via max(parameter) - min(parameter)

# Get segment-level coefficients with patient info
segment_coeffs_for_variation <- df_retained %>%
    group_by(patient, segment) %>%
    summarise(
        seg_m = first(p6e_m),
        seg_b = first(p6e_b),
        n_points = n(),
        .groups = "drop"
    )

# Compute within-patient variation (for patients with ≥2 segments)
patient_param_variation <- segment_coeffs_for_variation %>%
    group_by(patient) %>%
    filter(n() > 1) %>%
    summarise(
        variation_m = max(seg_m, na.rm = TRUE) - min(seg_m, na.rm = TRUE),
        variation_b = max(seg_b, na.rm = TRUE) - min(seg_b, na.rm = TRUE),
        n_segments = n(),
        .groups = "drop"
    )

n_patients_multi_seg <- nrow(patient_param_variation)

# Count patients with exactly 1 retained segment (variation not estimable)
n_patients_single_seg <- n_retained_patients - n_patients_multi_seg

# Summarize within-patient variation (median and IQR across patients)
# Store raw values for notes
param_variation_m <- list(
    median = median(patient_param_variation$variation_m, na.rm = TRUE),
    q1 = quantile(patient_param_variation$variation_m, 0.25, na.rm = TRUE),
    q3 = quantile(patient_param_variation$variation_m, 0.75, na.rm = TRUE)
)
param_variation_b <- list(
    median = median(patient_param_variation$variation_b, na.rm = TRUE),
    q1 = quantile(patient_param_variation$variation_b, 0.25, na.rm = TRUE),
    q3 = quantile(patient_param_variation$variation_b, 0.75, na.rm = TRUE)
)

# Format for CSV output
param_variation_summary <- data.frame(
    variable = c("m", "b (mM)"),
    `within-patient variation [a]` = c(
        sprintf("%.2f (%.2f–%.2f)", param_variation_m$median, param_variation_m$q1, param_variation_m$q3),
        sprintf("%.2f (%.2f–%.2f)", param_variation_b$median, param_variation_b$q1, param_variation_b$q3)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

write.csv(param_variation_summary, file.path(output_dir, "5_param_variation.csv"), row.names = FALSE)

# ==============================================================================
# Execution Time
# ==============================================================================
end_time <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
time_df <- data.frame(execution_time_seconds = execution_time_seconds)
write.csv(time_df, file.path(output_dir, "__execution_time.csv"), row.names = FALSE)

# Build legends for figure + coefficient table + param variation table
legends <- list(
  list(
    target = "Figure 1 (3_figure.png)",
    caption = sprintf(
      "Left panels: Lactate vs Pyruvate; Right panels: LPR vs Pyruvate. Dashed horizontal lines indicate gradient asymptotes (m). Red shading denotes non-physiological negative values (may not be seen in plots depending on viewport). (A-B) Unstratified: all %d retained segments from %d patients showing spaghetti plots of segment-level fits colored by intercept type (teal = Pb, gold = Nb, gray = Zb) with cohort median trend (brown). (C-D) Type Pb segments (b > 0): %d segments, %d patients; positive intercepts produce decreasing LPR hyperbolae. (E-F) Type Nb segments (b < 0): %d segments, %d patients; negative intercepts produce increasing LPR hyperbolae. (G-H) Combined overlay comparing Type Pb and Type Nb cohort medians.",
      n_retained_segments, n_retained_patients,
      n_segments_pos, n_patients_pos,
      n_segments_neg, n_patients_neg
    ),
    abbreviations = "LPR, lactate/pyruvate ratio; Nb, negative intercept; Pb, positive intercept; Zb, zero intercept"
  ),
  list(
    target = "Table 1 (4_coefficients.csv)",
    caption = sprintf(
      "Linear model coefficients for retained segments (m > 0). Values are cohort median (IQR) of patient-level segment-weighted means. Unstratified: %d segments, %d patients; Type Pb: %d segments, %d patients; Type Nb: %d segments, %d patients.",
      n_retained_segments, n_retained_patients,
      n_segments_pos, n_patients_pos,
      n_segments_neg, n_patients_neg
    ),
    abbreviations = "b, lactate–pyruvate linear intercept; IQR, interquartile range; m, lactate–pyruvate linear gradient; Nb, negative intercept; Pb, positive intercept; R², within-segment fit"
  ),
  list(
    target = "Table 2 (5_param_variation.csv)",
    caption = sprintf(
      "Within-patient variation of linear model parameters. Values are cohort median (IQR) of patient-level variation. Based on %d patients with ≥2 retained segments (%d patients with only 1 retained segment excluded).",
      n_patients_multi_seg, n_patients_single_seg
    ),
    footnotes = c(
      "[a] Variation = max − min across segments for each patient."
    ),
    abbreviations = "b, lactate–pyruvate linear intercept; IQR, interquartile range; m, lactate–pyruvate linear gradient"
  )
)

writeLines(build_notes(legends, title = "5_linear_LP_models"), file.path(output_dir, "0_notes.txt"))

cat("\n\n5_linear_LP_models.R complete.\n\n")
