# ==============================================================================
# Script: extra_linear_to_hyperbolic.py
# Manuscript relevance: 2.3.i, Fig. S1
# ==============================================================================
# PURPOSE:
#   Generate a conceptual figure illustrating the mathematical transformation
#   from linear (L = mP + b) to hyperbolic (LPR = m + b/P) representation.
#   Demonstrates that LPR ≠ m whenever b ≠ 0, with error magnitude inversely
#   proportional to pyruvate concentration.
#
# INPUT:
#   - None (generates synthetic data from mathematical relationships)
#
# OUTPUT:
#   - 1_output/extra_linear_to_hyperbolic/0_notes.txt: Figure description
#   - 1_output/extra_linear_to_hyperbolic/1_figure.png: 3-row conceptual figure
#
# MATHEMATICAL FRAMEWORK:
#   Linear model:     L = m·P + b
#   Hyperbolic model: LPR = L/P = m + b/P
#   Key insight:      LPR approaches m as P → ∞; deviates as P → 0
#
# FIGURE PANELS:
# ------------------------------------------------------------------------------
# Row 1: Effect of intercept (b) at fixed gradient (m = 25)
#   (A) Linear: L vs P for b ∈ {-1, -0.25, 0, +0.25, +1}
#   (B) Hyperbolic: LPR vs P showing the effect of the sign and magnitude of b
#
# Row 2: Effect of gradient (m) at fixed positive intercept (b = +1)
#   (C) Linear: L vs P for m ∈ {15, 25, 35}
#   (D) Hyperbolic: LPR vs P showing deviation above m at low P
#
# Row 3: Effect of gradient (m) at fixed negative intercept (b = -1)
#   (E) Linear: L vs P for m ∈ {15, 25, 35}
#   (F) Hyperbolic: LPR vs P showing deviation below m at low P
#
# VISUALIZATION PARAMETERS:
#   - Colors: Teal shades (b > 0), Gold shades (b < 0), Black (b = 0)
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from _shared.notes import build_notes


# ==============================================================================
# Variables / Toggles
# ==============================================================================
# --- Figure Size ---
WIDTH_MM = 185  # two columns
PANEL_HEIGHT_MM = 55  # Panel height in mm
N_ROWS = 3  # Number of rows in this figure
MM_TO_INCH = 1 / 25.4
FIG_WIDTH_IN = WIDTH_MM * MM_TO_INCH
FIG_HEIGHT_IN = (PANEL_HEIGHT_MM * N_ROWS) * MM_TO_INCH  # Panel height per row
# --- Font Sizes ---
FONT_PANEL_LABEL = 9
FONT_AXIS_LABEL = 7
FONT_LEGEND = 6
FONT_TICKS = 6

# --- Line Width ---
LINE_WIDTH = 1.5  # Adjusted for smaller figure size (in points)

# Pyruvate range for the plots (in mM for calculations, converted to μM for display)
PYRUVATE_RANGE_MM = np.linspace(0.001, 0.6, 500)
PYRUVATE_RANGE_UM = PYRUVATE_RANGE_MM * 1000  # Convert to μM for display

# Parameters for the plots
M_CONSTANT = 25
B_CONSTANT_POS = 1
B_CONSTANT_NEG = -1

M_VALUES_TO_TEST = [35, 25, 15] # Reversed for legend order
B_VALUES_TO_TEST = [1, 0.25, 0, -0.25, -1]

# --- Color palette definitions ---
# Colors for positive, zero/neutral, and negative intercepts
COLOR_B_POS = "#2D6B5D"  # Dark Teal (from gradient 35)
COLOR_B_POS_LIGHT = "#3A9684" # Medium Teal (from gradient 25)
COLOR_B_ZERO = "black"   # Black
COLOR_B_NEG_LIGHT = "#BDA955" # Medium Gold (from gradient 25)
COLOR_B_NEG = "#A89752"  # Dark Gold (from gradient 35)

# Thematic shades for varying gradient 'm' - highest m is darkest
TEAL_SHADES = ["#2D6B5D", "#3A9684", "#6DB8B4"] # Dark, Medium, Light
GOLD_SHADES = ["#A89752", "#BDA955", "#D0BD58"] # Dark, Medium, Light

# Reference line colors
NEG_REF_LINE_COLOR = 'red'
NEG_REF_AREA_COLOR = 'red' 

print("\n\nStarting extra_linear_to_hyperbolic.py")

# ==============================================================================
# Output
# ==============================================================================
OUTPUT_FILENAME = "1_figure.png"

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "1_output" / "extra_linear_to_hyperbolic"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
output_path = OUTPUT_DIR / OUTPUT_FILENAME

# ==============================================================================
# Definitions / Functions
# ==============================================================================
def calculate_lactate(pyruvate, m, b):
    """Calculates lactate from the linear model."""
    return m * pyruvate + b

def calculate_lpr(pyruvate, m, b):
    """Calculates LPR from the hyperbolic model."""
    return m + b / pyruvate

def apply_common_aesthetics(ax):
    """Applies annotations for axis baselines."""
    ax.axhline(0, color=NEG_REF_LINE_COLOR, linestyle=':', linewidth=1, alpha=0.7)
    ax.axhspan(ax.get_ylim()[0], 0, color=NEG_REF_AREA_COLOR, alpha=0.05, zorder=0)
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.8)

# ==============================================================================
# Processing / Visualization
# ==============================================================================

plt.style.use('seaborn-v0_8-whitegrid')

# Create a 3x2 subplot figure
fig, axes = plt.subplots(3, 2, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

# --- Panel A & B: Effect of varying Intercept (b) ---
ax_A = axes[0, 0]
ax_B = axes[0, 1]
color_map_b = {1: COLOR_B_POS, 0.25: COLOR_B_POS_LIGHT, 0: COLOR_B_ZERO, -0.25: COLOR_B_NEG_LIGHT, -1: COLOR_B_NEG}

for b_val in B_VALUES_TO_TEST:
    lactate = calculate_lactate(PYRUVATE_RANGE_MM, M_CONSTANT, b_val)
    lpr = calculate_lpr(PYRUVATE_RANGE_MM, M_CONSTANT, b_val)
    
    # Set zorder to draw b=0 line behind others
    z_order = 1 if b_val == 0 else 2
    
    ax_A.plot(PYRUVATE_RANGE_UM, lactate, label=f'b = {b_val}', color=color_map_b[b_val], linewidth=LINE_WIDTH, zorder=z_order)
    ax_B.plot(PYRUVATE_RANGE_UM, lpr, label=f'b = {b_val}', color=color_map_b[b_val], linewidth=LINE_WIDTH, zorder=z_order)

ax_A.text(-0.1, 1.15, 'A', transform=ax_A.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_A.set_xlabel('')
ax_A.set_ylabel('Lactate (mM)', fontsize=FONT_AXIS_LABEL)
legend_A = ax_A.legend(title=f"m = {M_CONSTANT}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black')
legend_A.get_frame().set_linewidth(0.6)
ax_A.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_A.set_xlim(0, 150)  # Display range of 0–150 μM (0.15 mM)
ax_A.set_ylim(-0.2, 4)
apply_common_aesthetics(ax_A)


ax_B.text(-0.1, 1.15, 'B', transform=ax_B.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_B.set_xlabel('')
ax_B.set_ylabel('LPR', fontsize=FONT_AXIS_LABEL)
legend_B = ax_B.legend(title=f"m = {M_CONSTANT}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black', loc='upper right')
legend_B.get_frame().set_linewidth(0.6)
ax_B.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_B.set_ylim(-10, 60)
apply_common_aesthetics(ax_B)


# --- Panel C & D: Effect of varying Gradient (m) with b = +1 (Teal Shades) ---
ax_C = axes[1, 0]
ax_D = axes[1, 1]
color_map_m_pos = dict(zip(M_VALUES_TO_TEST, TEAL_SHADES))

for m_val in M_VALUES_TO_TEST:
    lactate = calculate_lactate(PYRUVATE_RANGE_MM, m_val, B_CONSTANT_POS)
    lpr = calculate_lpr(PYRUVATE_RANGE_MM, m_val, B_CONSTANT_POS)
    
    ax_C.plot(PYRUVATE_RANGE_UM, lactate, label=f'm = {m_val}', color=color_map_m_pos[m_val], linewidth=LINE_WIDTH)
    ax_D.plot(PYRUVATE_RANGE_UM, lpr, label=f'm = {m_val}', color=color_map_m_pos[m_val], linewidth=LINE_WIDTH)

ax_C.text(-0.1, 1.15, 'C', transform=ax_C.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_C.set_xlabel('')
ax_C.set_ylabel('Lactate (mM)', fontsize=FONT_AXIS_LABEL)
legend_C = ax_C.legend(title=f"b = {B_CONSTANT_POS}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black')
legend_C.get_frame().set_linewidth(0.6)
ax_C.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_C.set_ylim(0, 23)
apply_common_aesthetics(ax_C)


ax_D.text(-0.1, 1.15, 'D', transform=ax_D.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_D.set_xlabel('')
ax_D.set_ylabel('LPR', fontsize=FONT_AXIS_LABEL)
legend_D = ax_D.legend(title=f"b = {B_CONSTANT_POS}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black')
legend_D.get_frame().set_linewidth(0.6)
ax_D.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_D.set_ylim(10, 60)
# Horizontal dashed lines at each m value (asymptotes) in matching colors
for m_val in M_VALUES_TO_TEST:
    ax_D.axhline(m_val, color=color_map_m_pos[m_val], linestyle='--', linewidth=0.6, alpha=0.8, zorder=1)
apply_common_aesthetics(ax_D)


# --- Panel E & F: Effect of varying Gradient (m) with b = -1 (Gold Shades) ---
ax_E = axes[2, 0]
ax_F = axes[2, 1]
color_map_m_neg = dict(zip(M_VALUES_TO_TEST, GOLD_SHADES))

for m_val in M_VALUES_TO_TEST:
    lactate = calculate_lactate(PYRUVATE_RANGE_MM, m_val, B_CONSTANT_NEG)
    lpr = calculate_lpr(PYRUVATE_RANGE_MM, m_val, B_CONSTANT_NEG)
    
    ax_E.plot(PYRUVATE_RANGE_UM, lactate, label=f'm = {m_val}', color=color_map_m_neg[m_val], linewidth=LINE_WIDTH)
    ax_F.plot(PYRUVATE_RANGE_UM, lpr, label=f'm = {m_val}', color=color_map_m_neg[m_val], linewidth=LINE_WIDTH)

ax_E.text(-0.1, 1.15, 'E', transform=ax_E.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_E.set_xlabel('Pyruvate (μM)', fontsize=FONT_AXIS_LABEL)
ax_E.set_ylabel('Lactate (mM)', fontsize=FONT_AXIS_LABEL)
legend_E = ax_E.legend(title=f"b = {B_CONSTANT_NEG}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black')
legend_E.get_frame().set_linewidth(0.6)
ax_E.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_E.set_ylim(-2, 21)
apply_common_aesthetics(ax_E)


ax_F.text(-0.1, 1.15, 'F', transform=ax_F.transAxes, size=FONT_PANEL_LABEL, weight='bold', va='top', ha='left')
ax_F.set_xlabel('Pyruvate (μM)', fontsize=FONT_AXIS_LABEL)
ax_F.set_ylabel('LPR', fontsize=FONT_AXIS_LABEL)
legend_F = ax_F.legend(title=f"b = {B_CONSTANT_NEG}", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND, frameon=True, framealpha=1.0, facecolor='white', edgecolor='black')
legend_F.get_frame().set_linewidth(0.6)
ax_F.tick_params(axis='both', which='major', labelsize=FONT_TICKS)
ax_F.set_ylim(-10, 40)
# Horizontal dashed lines at each m value (asymptotes) in matching colors
for m_val in M_VALUES_TO_TEST:
    ax_F.axhline(m_val, color=color_map_m_neg[m_val], linestyle='--', linewidth=0.6, alpha=0.8, zorder=1)
apply_common_aesthetics(ax_F)

# Adjust layout and save the figure
plt.tight_layout(pad=1.0, h_pad=1.0, w_pad=2.0)
plt.savefig(output_path, dpi=300)
plt.close(fig)

# ==============================================================================
# Notes Generation
# ==============================================================================
# Build figure legend (only output for notes)
legends = [
    {
        "target": "Figure 1 (1_figure.png)",
        "caption": (
            f"Linear to hyperbolic transformation. Left: L = mP + b. Right: LPR = m + b/P. "
            f"(A-B) Varying intercept (b = {', '.join(str(b) for b in B_VALUES_TO_TEST)} mM) at m = {M_CONSTANT}: "
            f"positive b → decreasing hyperbola (teal), negative b → increasing (gold), b = 0 → constant LPR = m (black). "
            f"(C-D) Varying gradient (m = {', '.join(str(m) for m in M_VALUES_TO_TEST)}) at b = {B_CONSTANT_POS} mM. "
            f"(E-F) Varying gradient at b = {B_CONSTANT_NEG} mM. "
            f"Dashed lines = asymptotes at m. Red shading = non-physiological negative values."
        ),
        "abbreviations": "L, lactate; LPR, lactate/pyruvate ratio; mM, millimolar; P, pyruvate; μM, micromolar",
    },
]

notes_path = OUTPUT_DIR / "0_notes.txt"
notes_path.write_text(build_notes(legends, title="extra_linear_to_hyperbolic"), encoding="utf-8")

print("\nextra_linear_to_hyperbolic.py complete.\n")

