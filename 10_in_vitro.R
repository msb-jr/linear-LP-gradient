# ==============================================================================
# Script: 10_in_vitro.R
# Manuscript relevance: 3.8, Fig. 7
# ==============================================================================
# PURPOSE:
#   Validate the linear lactate–pyruvate framework in vitro by testing whether
#   electron transport chain (ETC) inhibition via rotenone increases the
#   lactate–pyruvate gradient (m), providing mechanistic evidence linking the gradient
#   to mitochondrial dysfunction.
#
# INPUT:
#   - 0_input/in_vitro.csv: Pooled timecourse data (control vs rotenone)
#     Columns: type (wt|rot), time (hours), lactate (mM), pyruvate (mM)
#
# OUTPUT:
#   - 1_output/10_in_vitro/__execution_time.csv: Runtime log
#   - 1_output/10_in_vitro/_1_model_summaries.txt: OLS + ANCOVA summaries
#   - 1_output/10_in_vitro/_1_stats.csv: Regression and ANCOVA statistics
#   - 1_output/10_in_vitro/0_notes.txt: Figure legend
#   - 1_output/10_in_vitro/1_figure.png: Scatterplot with regression fits
#
# ANALYSES:
# ------------------------------------------------------------------------------
# (1) Per-group OLS regressions (_1_stats.csv, _1_model_summaries.txt)
#     - Separate OLS fits: lactate ~ pyruvate for control and rotenone
#     - Reports: gradient, intercept, R², Pearson r, p-value for each group
#     Test:       Standard OLS F-test
#
# (2) ANCOVA (_1_stats.csv, _1_model_summaries.txt)
#     - Tests whether rotenone shifts the gradient and/or intercept
#     Model:      lactate ~ pyruvate * type
#     - pyruvate:type interaction tests gradient difference
#     - type main effect tests intercept difference
#     Test:       Standard F-test
#
# (3) Figure (1_figure.png)
#     - Single-panel scatterplot: lactate vs pyruvate (µM on x-axis)
#     - Points and regression lines coloured by group
#     - 95% CI ribbons around regression fits
#     Dimensions: 90 mm width (single column)
#
# MULTIPLE COMPARISONS CORRECTION:
#   Not applied. The ANCOVA interaction term (gradient shift) tests the primary, 
#   a priori physiological hypothesis. The main effect (intercept shift) is a 
#   distinct, exploratory test of a different parameter. Per-group OLS fits are 
#   descriptive validations. All tests address independent components of the model.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(here::here("_shared", "notes.R"))

cat("\n\nStarting 10_in_vitro.R\n")

start_time <- Sys.time()

save_plot_portable <- function(filename, plot, ...) {
  if (isTRUE(capabilities("cairo"))) {
    ggsave(filename = filename, plot = plot, type = "cairo", ...)
  } else {
    warning(sprintf("Cairo backend unavailable; saving %s with the default device.", basename(filename)))
    ggsave(filename = filename, plot = plot, ...)
  }
}

format_p_text <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

# ==============================================================================
# Input / Output
# ==============================================================================
INPUT_FILE <- here::here("0_input", "in_vitro.csv")

if (!file.exists(INPUT_FILE)) {
  cat("\n\nError: In vitro data file not found at 0_input/in_vitro.csv\n")
  stop()
}

OUTPUT_DIR <- here::here("1_output", "10_in_vitro")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Data Loading
# ==============================================================================

df <- read_csv(INPUT_FILE, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    type = factor(type, levels = c("wt", "rot"), labels = c("Control", "Rotenone")),
    pyruvate_um = pyruvate * 1000
  )

df_ctrl <- df %>% filter(type == "Control")
df_rot <- df %>% filter(type == "Rotenone")

n_ctrl <- nrow(df_ctrl)
n_rot <- nrow(df_rot)

# ==============================================================================
# Per-Group OLS Regressions
# ==============================================================================

ols_ctrl <- lm(lactate ~ pyruvate, data = df_ctrl)
ols_rot <- lm(lactate ~ pyruvate, data = df_rot)

extract_ols <- function(model, data) {
  s <- summary(model)
  cor_test <- cor.test(data$pyruvate, data$lactate)
  list(
    gradient = coef(model)[["pyruvate"]],
    intercept = coef(model)[["(Intercept)"]],
    r = cor_test$estimate,
    r2 = s$r.squared,
    p_value = cor_test$p.value
  )
}

ols_ctrl_stats <- extract_ols(ols_ctrl, df_ctrl)
ols_rot_stats <- extract_ols(ols_rot, df_rot)

# ==============================================================================
# ANCOVA: lactate ~ pyruvate * type
# ==============================================================================

ancova_model <- lm(lactate ~ pyruvate * type, data = df)
ancova_summary <- summary(ancova_model)
ancova_anova <- anova(ancova_model)

ancova_interaction_p <- ancova_anova["pyruvate:type", "Pr(>F)"]
ancova_type_p <- ancova_anova["type", "Pr(>F)"]
ancova_dof <- ancova_model$df.residual

p_text_interaction <- format_p_text(ancova_interaction_p)
p_text_type <- format_p_text(ancova_type_p)

# ==============================================================================
# Save Model Summaries
# ==============================================================================

sink(file.path(OUTPUT_DIR, "_1_model_summaries.txt"))
cat("=======================================================================\n")
cat("Per-Group OLS: Control (lactate ~ pyruvate)\n")
cat(sprintf("N = %d\n", n_ctrl))
cat("=======================================================================\n")
print(summary(ols_ctrl))
cat("\n\n")
cat("=======================================================================\n")
cat("Per-Group OLS: Rotenone (lactate ~ pyruvate)\n")
cat(sprintf("N = %d\n", n_rot))
cat("=======================================================================\n")
print(summary(ols_rot))
cat("\n\n")
cat("=======================================================================\n")
cat("ANCOVA: lactate ~ pyruvate * type\n")
cat("Interaction (pyruvate:type) tests gradient difference.\n")
cat("Main effect (type) tests intercept difference.\n")
cat(sprintf("N = %d (Control %d, Rotenone %d), residual df = %d\n",
            n_ctrl + n_rot, n_ctrl, n_rot, ancova_dof))
cat("=======================================================================\n")
print(ancova_summary)
cat("\n\nANOVA table:\n")
print(ancova_anova)
sink()

# ==============================================================================
# Visualization
# ==============================================================================

common_theme <- theme_minimal() +
  theme(
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "mm")
  )

group_colors <- c("Control" = "#0072B2", "Rotenone" = "#CC3311")

fig <- ggplot(df, aes(x = pyruvate_um, y = lactate, colour = type, fill = type)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(x = expression(paste("Pyruvate (", mu, "M)")),
       y = "Lactate (mM)") +
  common_theme

save_plot_portable(
  filename = file.path(OUTPUT_DIR, "1_figure.png"),
  plot = fig,
  width = 90, height = 70, units = "mm", dpi = 300
)

# ==============================================================================
# Summary Statistics CSV
# ==============================================================================

summary_stats <- data.frame(
  Metric = c(
    "N_control", "N_rotenone",
    "Control_gradient", "Control_intercept_mM", "Control_r", "Control_R2", "Control_p",
    "Rotenone_gradient", "Rotenone_intercept_mM", "Rotenone_r", "Rotenone_R2", "Rotenone_p",
    "ANCOVA_gradient_shift_p", "ANCOVA_intercept_shift_p", "ANCOVA_residual_df"
  ),
  Value = c(
    n_ctrl, n_rot,
    ols_ctrl_stats$gradient, ols_ctrl_stats$intercept, ols_ctrl_stats$r,
    ols_ctrl_stats$r2, ols_ctrl_stats$p_value,
    ols_rot_stats$gradient, ols_rot_stats$intercept, ols_rot_stats$r,
    ols_rot_stats$r2, ols_rot_stats$p_value,
    ancova_interaction_p, ancova_type_p, ancova_dof
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

cat("\n\n10_in_vitro.R complete.\n\n")

# ==============================================================================
# Figure Legend
# ==============================================================================
legends <- list(
  list(
    target = "Figure 1 (1_figure.png)",
    caption = sprintf(
      paste0(
        "In vitro lactate–pyruvate relationship under ETC inhibition. ",
        "Scatterplot of lactate vs pyruvate for control (blue; n = %d; gradient = %.1f, intercept = %.3f mM, R² = %.2f) ",
        "and rotenone-treated (red; n = %d; gradient = %.1f, intercept = %.3f mM, R² = %.2f) conditions. ",
        "ANCOVA: gradient shift %s; intercept shift %s (residual df = %d). ",
        "Data pooled from multiple timecourse experiments with sampling between 1 and 24 hours. ",
        "Lines = OLS regression fits; shaded = 95%% CI."
      ),
      n_ctrl, ols_ctrl_stats$gradient, ols_ctrl_stats$intercept, ols_ctrl_stats$r2,
      n_rot, ols_rot_stats$gradient, ols_rot_stats$intercept, ols_rot_stats$r2,
      p_text_interaction, p_text_type, ancova_dof
    ),
    abbreviations = "CI, confidence interval; ETC, electron transport chain; OLS, ordinary least squares"
  )
)

writeLines(build_notes(legends, title = "10_in_vitro"), file.path(OUTPUT_DIR, "0_notes.txt"))
