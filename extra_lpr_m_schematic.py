# ==============================================================================
# Script: extra_lpr_m_schematic.py
# Manuscript relevance: Fig. 1
# ==============================================================================
# PURPOSE:
#   Generate figures illustrating how within-patient variation
#   in the linear gradient (m) manifests in lactate-pyruvate scatter plots
#   and in the LPR time series, using type-specific (Nb / Pb) parameters.
#   These figures will be inserted into the Fig. 1 CMD schematic.
#
# INPUT:
#   - None (generates synthetic data from representative parameters)
#
# OUTPUT:
#   - 1_output/extra_lpr_m_schematic/0_notes.txt: Figure descriptions
#   - 1_output/extra_lpr_m_schematic/1_figure.png: LPR vs time (both segments)
#   - 1_output/extra_lpr_m_schematic/2_figure.png: Gold triangles (Type Nb)
#   - 1_output/extra_lpr_m_schematic/3_figure.png: Teal circles (Type Pb)
#
# PARAMETERS (derived from Cambridge cohort):
#   Segment 1 — Type Nb (b < 0), gold triangles:
#     m = 29.04, b = -0.55 mM, R² = 0.85, Pyr 75–170 μM
#   Segment 2 — Type Pb (b > 0), teal circles:
#     m = 17.97, b = 0.87 mM, R² = 0.63, Pyr 61–171 μM
#   Sources: 5_linear_LP_models/4_coefficients.csv
#            6_hyperbolic_LPR_error/1_stats.csv
#
#   Shared:
#   - n = 30 per segment (median datapoints per retained segment)
#     Source: 4_segmentation_stats/1_stats.csv
#   - Time span: 41.4 h per segment (median), hours 2–85 total
#     Sources: 4_segmentation_stats/1_stats.csv
#              1_pre-processing.py: FILTER_EARLY
#   - Minimum time spacing: 1 hour between consecutive points
# ==============================================================================

import os
os.environ['MPLBACKEND'] = 'Agg'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from _shared.notes import build_notes

print("\n\nStarting extra_lpr_m_schematic.py")

# ==============================================================================
# Input / Output
# ==============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "1_output" / "extra_lpr_m_schematic"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ==============================================================================
# Configuration / Toggles
# ==============================================================================
PLOT_DPI = 300
FONTSIZE_PANEL_LABEL = 9
FONTSIZE_LABELS = 7
FONTSIZE_TICKS = 6
BASE_MARKER_SIZE = 10
LINE_WIDTH = 1.5

# Colors and markers (consistent with 6_hyperbolic_LPR_error.R)
COLOR_PB = '#44AA99'   # Teal for Type Pb (b > 0)
COLOR_NB = '#DDCC77'   # Gold for Type Nb (b < 0)
MARKER_PB = 'o'        # Circle for Type Pb
MARKER_NB = '^'        # Triangle for Type Nb


def darken_color(color_hex, factor=0.7):
    """Darkens a hex color by reducing brightness."""
    rgb = mcolors.hex2color(color_hex)
    hsv = mcolors.rgb_to_hsv(np.array(rgb))
    hsv[2] = max(0, hsv[2] * factor)
    return mcolors.hsv_to_rgb(hsv)


COLOR_PB_DARK = darken_color(COLOR_PB)
COLOR_NB_DARK = darken_color(COLOR_NB)

# --- General parameters ---
SEED = 42
N_POINTS = 30
TIME_START_H = 2
TIME_END_H = 85
TIME_MID_H = (TIME_START_H + TIME_END_H) / 2
MIN_TIME_GAP_H = 1.0

# --- Type Nb (b < 0) — Segment 1, gold ---
# Source: 5_linear_LP_models/4_coefficients.csv (median row "Type Nb")
M_NB = 29.04
B_NB = -0.55         # median
R2_NB = 0.85
# Source: 6_hyperbolic_LPR_error/1_stats.csv (Type Nb median Pmin/Pmax)
PYR_MIN_NB_UM = 75
PYR_MAX_NB_UM = 170

# --- Type Pb (b > 0) — Segment 2, teal ---
# Source: 5_linear_LP_models/4_coefficients.csv (median row "Type Pb")
M_PB = 17.97
B_PB = 0.87          # median
R2_PB = 0.63
# Source: 6_hyperbolic_LPR_error/1_stats.csv (Type Pb median Pmin/Pmax)
PYR_MIN_PB_UM = 61
PYR_MAX_PB_UM = 171

# Segment prevalence (from 5_linear_LP_models/0_notes.txt)
N_SEG_PB = 967
N_SEG_NB = 372
N_SEG_TOTAL = 1339
PCT_PB = N_SEG_PB / N_SEG_TOTAL * 100   # 72.2%
PCT_NB = N_SEG_NB / N_SEG_TOTAL * 100   # 27.8%


# ==============================================================================
# Helper Functions
# ==============================================================================
def apply_axis_style(ax):
    """Apply consistent styling to axes."""
    ax.tick_params(axis='both', which='major', labelsize=FONTSIZE_TICKS, colors='black')
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)


def generate_spaced_times(n, t_start, t_end, min_gap, rng):
    """Generate n sorted time points in [t_start, t_end] with ≥ min_gap spacing."""
    slack = (t_end - t_start) - (n - 1) * min_gap
    raw = np.sort(rng.uniform(0, 1, n))
    raw = (raw - raw[0]) / (raw[-1] - raw[0]) * slack
    times = t_start + raw + np.arange(n) * min_gap
    return times


def generate_segment_data(m, n, pyr_min_um, pyr_max_um, b, r2_target, rng):
    """Generate synthetic L vs P data with exact R² for OLS regression.

    Noise is constructed orthogonal to the regressor so that
    the OLS fit recovers the true gradient/intercept and the
    coefficient of determination equals r2_target exactly.
    """
    pyr_mm = rng.uniform(pyr_min_um / 1000, pyr_max_um / 1000, n)

    y_pred = m * pyr_mm + b

    noise = rng.standard_normal(n)
    noise = noise - noise.mean()
    pyr_c = pyr_mm - pyr_mm.mean()
    noise = noise - (np.dot(noise, pyr_c) / np.dot(pyr_c, pyr_c)) * pyr_c

    ss_reg = np.sum((y_pred - y_pred.mean()) ** 2)
    ss_noise = np.sum(noise ** 2)
    scale = np.sqrt(ss_reg * (1 - r2_target) / (r2_target * ss_noise)) if ss_noise > 0 else 0

    lactate = y_pred + scale * noise
    pyr_um = pyr_mm * 1000
    lpr = lactate / pyr_mm

    return pyr_um, pyr_mm, lactate, lpr


# ==============================================================================
# Data Generation
# ==============================================================================
rng = np.random.default_rng(SEED)

# Segment 1: Type Nb (gold)
pyr_um_1, pyr_mm_1, lac_1, lpr_1 = generate_segment_data(
    M_NB, N_POINTS, PYR_MIN_NB_UM, PYR_MAX_NB_UM, B_NB, R2_NB, rng)

# Segment 2: Type Pb (teal)
pyr_um_2, pyr_mm_2, lac_2, lpr_2 = generate_segment_data(
    M_PB, N_POINTS, PYR_MIN_PB_UM, PYR_MAX_PB_UM, B_PB, R2_PB, rng)

# Time points with ≥ 1 h spacing
time_1 = generate_spaced_times(N_POINTS, TIME_START_H, TIME_MID_H, MIN_TIME_GAP_H, rng)
time_2 = generate_spaced_times(N_POINTS, TIME_MID_H, TIME_END_H, MIN_TIME_GAP_H, rng)

# Shared axis limits for scatter panels
all_pyr_um = np.concatenate([pyr_um_1, pyr_um_2])
all_lac = np.concatenate([lac_1, lac_2])
pyr_pad = (all_pyr_um.max() - all_pyr_um.min()) * 0.08
lac_pad = (all_lac.max() - all_lac.min()) * 0.08
xlim_scatter = (all_pyr_um.min() - pyr_pad, all_pyr_um.max() + pyr_pad)
ylim_scatter = (all_lac.min() - lac_pad, all_lac.max() + lac_pad)

# ==============================================================================
# Figure 1: LPR vs Time
# ==============================================================================
fig1, ax1 = plt.subplots(figsize=(3.0, 1.5))

ax1.scatter(time_1, lpr_1, color=COLOR_NB, marker=MARKER_NB,
            s=BASE_MARKER_SIZE, alpha=0.7, zorder=2)
ax1.scatter(time_2, lpr_2, color=COLOR_PB, marker=MARKER_PB,
            s=BASE_MARKER_SIZE, alpha=0.7, zorder=2)

ax1.set_xlabel("Time (hours)", fontsize=FONTSIZE_LABELS)
ax1.set_ylabel("LPR", fontsize=FONTSIZE_LABELS)
ax1.grid(True, linestyle='-', alpha=0.7)
apply_axis_style(ax1)

fig1.tight_layout(pad=0.3)
fig1.savefig(OUTPUT_DIR / "1_figure.png", dpi=PLOT_DPI)
plt.close(fig1)

# ==============================================================================
# Figure 2: Type Nb — Gold triangles (m = 29.04)
# ==============================================================================
fig2, ax2 = plt.subplots(figsize=(1.5, 1.5))

ax2.scatter(pyr_um_1, lac_1, color=COLOR_NB, marker=MARKER_NB,
            s=BASE_MARKER_SIZE, alpha=0.7, zorder=2)

x_line = np.array([pyr_mm_1.min(), pyr_mm_1.max()])
y_line = M_NB * x_line + B_NB
ax2.plot(x_line * 1000, y_line, color=COLOR_NB_DARK, linewidth=LINE_WIDTH, zorder=3)

ax2.set_xlabel("Pyruvate (\u03bcM)", fontsize=FONTSIZE_LABELS)
ax2.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
ax2.set_xlim(xlim_scatter)
ax2.set_ylim(ylim_scatter)
ax2.grid(True, linestyle='-', alpha=0.7)
apply_axis_style(ax2)

fig2.tight_layout(pad=0.3)
fig2.savefig(OUTPUT_DIR / "2_figure.png", dpi=PLOT_DPI)
plt.close(fig2)

# ==============================================================================
# Figure 3: Type Pb — Teal circles (m = 17.97)
# ==============================================================================
fig3, ax3 = plt.subplots(figsize=(1.5, 1.5))

ax3.scatter(pyr_um_2, lac_2, color=COLOR_PB, marker=MARKER_PB,
            s=BASE_MARKER_SIZE, alpha=0.7, zorder=2)

x_line = np.array([pyr_mm_2.min(), pyr_mm_2.max()])
y_line = M_PB * x_line + B_PB
ax3.plot(x_line * 1000, y_line, color=COLOR_PB_DARK, linewidth=LINE_WIDTH, zorder=3)

ax3.set_xlabel("Pyruvate (\u03bcM)", fontsize=FONTSIZE_LABELS)
ax3.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
ax3.set_xlim(xlim_scatter)
ax3.set_ylim(ylim_scatter)
ax3.grid(True, linestyle='-', alpha=0.7)
apply_axis_style(ax3)

fig3.tight_layout(pad=0.3)
fig3.savefig(OUTPUT_DIR / "3_figure.png", dpi=PLOT_DPI)
plt.close(fig3)

# ==============================================================================
# Notes
# ==============================================================================
legends = [
    {
        "target": "Figure 1 (1_figure.png)",
        "caption": (
            f"LPR vs time for the two synthetic segments shown in Figures 2\u20133. "
            f"Gold triangles (Type Nb): m = {M_NB}, LPR range "
            f"{lpr_1.min():.1f}\u2013{lpr_1.max():.1f}. "
            f"Teal circles (Type Pb): m = {M_PB}, LPR range "
            f"{lpr_2.min():.1f}\u2013{lpr_2.max():.1f}. "
            f"Time span: {TIME_START_H}\u2013{TIME_END_H} hours, representing two "
            f"contiguous segments of approximately "
            f"{(TIME_END_H - TIME_START_H) / 2:.1f} hours each "
            f"(minimum {MIN_TIME_GAP_H:.0f} h spacing between consecutive points)."
        ),
        "abbreviations": "LPR, lactate/pyruvate ratio; m, gradient; Nb, negative intercept; Pb, positive intercept",
    },
    {
        "target": "Figure 2 (2_figure.png)",
        "caption": (
            f"Lactate vs pyruvate scatter plot for a synthetic Type Nb segment "
            f"(gold triangles). m = {M_NB}, b = {B_NB} mM, "
            f"R\u00b2 = {R2_NB}, n = {N_POINTS}. "
            f"Pyruvate range: {pyr_um_1.min():.0f}\u2013{pyr_um_1.max():.0f} \u03bcM. "
            f"Type Nb segments comprise {PCT_NB:.1f}% of retained segments "
            f"({N_SEG_NB}/{N_SEG_TOTAL})."
        ),
        "abbreviations": (
            "L, lactate (mM); LPR, lactate/pyruvate ratio; m, gradient; "
            "mM, millimolar; Nb, negative intercept; P, pyruvate; "
            "Pb, positive intercept; \u03bcM, micromolar"
        ),
    },
    {
        "target": "Figure 3 (3_figure.png)",
        "caption": (
            f"Lactate vs pyruvate scatter plot for a synthetic Type Pb segment "
            f"(teal circles). m = {M_PB}, b = {B_PB} mM, "
            f"R\u00b2 = {R2_PB}, n = {N_POINTS}. "
            f"Pyruvate range: {pyr_um_2.min():.0f}\u2013{pyr_um_2.max():.0f} \u03bcM. "
            f"Type Pb segments comprise {PCT_PB:.1f}% of retained segments "
            f"({N_SEG_PB}/{N_SEG_TOTAL})."
        ),
        "abbreviations": (
            "L, lactate (mM); LPR, lactate/pyruvate ratio; m, gradient; "
            "mM, millimolar; Nb, negative intercept; P, pyruvate; "
            "Pb, positive intercept; \u03bcM, micromolar"
        ),
    },
    {
        "target": "Parameter sources",
        "caption": (
            f"All parameters are representative of the Cambridge TBI cohort. "
            f"Type Nb: median m = {M_NB}, median b = {B_NB} mM, "
            f"median R\u00b2 = {R2_NB}, pyruvate sampled from "
            f"{PYR_MIN_NB_UM}\u2013{PYR_MAX_NB_UM} \u03bcM "
            f"(observed: {pyr_um_1.min():.0f}\u2013{pyr_um_1.max():.0f} \u03bcM). "
            f"Type Pb: median m = {M_PB}, median b = {B_PB} mM, "
            f"median R\u00b2 = {R2_PB}, pyruvate sampled from "
            f"{PYR_MIN_PB_UM}\u2013{PYR_MAX_PB_UM} \u03bcM "
            f"(observed: {pyr_um_2.min():.0f}\u2013{pyr_um_2.max():.0f} \u03bcM). "
            f"(from 5_linear_LP_models.R \u2192 4_coefficients.csv and "
            f"6_hyperbolic_LPR_error.R \u2192 1_stats.csv). "
            f"Segment prevalence: Type Pb {PCT_PB:.1f}% ({N_SEG_PB}/{N_SEG_TOTAL}), "
            f"Type Nb {PCT_NB:.1f}% ({N_SEG_NB}/{N_SEG_TOTAL}) "
            f"(from 5_linear_LP_models.R \u2192 0_notes.txt). "
            f"Median datapoints per retained segment = {N_POINTS} "
            f"(from 4_segmentation_stats.R \u2192 1_stats.csv). "
            f"Time span {TIME_END_H - TIME_START_H:.0f} h for two segments based "
            f"on median 41.4 h per retained segment "
            f"(from 4_segmentation_stats.R \u2192 1_stats.csv). "
            f"Random seed = {SEED}."
        ),
    },
]

notes_path = OUTPUT_DIR / "0_notes.txt"
notes_path.write_text(
    build_notes(legends, title="extra_lpr_m_schematic"), encoding="utf-8"
)

print("\nextra_lpr_m_schematic.py complete.\n")
