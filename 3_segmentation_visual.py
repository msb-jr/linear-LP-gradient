# ==============================================================================
# Script: 3_segmentation_visual.py
# Manuscript relevance: Fig. 2, Fig. S2, Fig. S3
# ==============================================================================
# PURPOSE:
#   Generate per-patient figures visualizing the PELT segmentation results.
#
# INPUT:
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data
#   - 0_input/distinct_patients.csv (optional): Patients of particular interest [a]
#
# OUTPUT:
#   - 1_output/3_segmentation_visual/__execution_time.csv: Runtime log
#   - 1_output/3_segmentation_visual/0_notes.txt: Multi-patient figure legend
#   - 1_output/3_segmentation_visual/1_figure.png: Multi-patient figure [b]
#   - 1_output/3_segmentation_visual/plots_all/: Per-patient figures + legends
#   - 1_output/3_segmentation_visual/plots_discarded/: Patients with m≤0 segments
#   - 1_output/3_segmentation_visual/plots_distinct/: Patients from optional CSV
#   - 1_output/3_segmentation_visual/plots_no_change/: Patients with 0 breakpoints
#
# [a] Example distinct_patients.csv: patient,multi
#                                    a7b9s94,1
#                                    ak94j20,0
#                                    cs9d8fj,1
#
# [b] multi-patient figure setup in `Configuration / Toggles` (below)
# ==============================================================================

import sys
import os
import time
import pandas as pd
import numpy as np

# Set matplotlib backend before importing pyplot to prevent GUI issues
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
import shutil
import multiprocessing as mp

# ==============================================================================
# Input / Output
# ==============================================================================
from pathlib import Path

from _shared.notes import build_notes

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_FILE = SCRIPT_DIR / "1_output" / "2_linear_segmentation" / "1_results.csv"
OUTPUT_DIR = SCRIPT_DIR / "1_output" / "3_segmentation_visual"
OUTPUT_ALL_PATIENTS = OUTPUT_DIR / "plots_all"
OUTPUT_NO_BREAKPOINTS = OUTPUT_DIR / "plots_no_change"
OUTPUT_HAS_DISCARDED = OUTPUT_DIR / "plots_discarded"
OUTPUT_DISTINCT = OUTPUT_DIR / "plots_distinct"

# Optional distinct patients file
DISTINCT_PATIENTS_FILE = SCRIPT_DIR / "0_input" / "distinct_patients.csv"
# [a] See example structure in header

# ==============================================================================
# Configuration / Toggles
# ==============================================================================
PLOT_DPI = 300
FONTSIZE_PANEL_LABEL = 9
FONTSIZE_LABELS = 7
FONTSIZE_TICKS = 6
BASE_MARKER_SIZE = 10
AXIS_PADDING_FACTOR = 0.05

# Layout
MM_TO_INCH = 1.0 / 25.4
DOUBLE_COLUMN_WIDTH_MM = 185.0
DOUBLE_COLUMN_WIDTH_INCH = DOUBLE_COLUMN_WIDTH_MM * MM_TO_INCH
PANEL_HEIGHT_MM = 55.0
PANEL_HEIGHT = PANEL_HEIGHT_MM * MM_TO_INCH  # Convert to inches
MAX_COLS_DISCARDED = 5  # Max columns when has discarded segments (A, B, C, D, E)
MAX_COLS_RETAINED = 3   # Max columns when all retained (A, B, C)

# Colors
OKABE_ITO_COLORS = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                    '#0072B2', '#D55E00', '#CC79A7', '#000000']
CMAP_TIME = 'cividis'
DISCARDED_MARKER = 'x'
NEUTRAL_GREY = 'dimgrey'

# Parallel processing (-1 = all cores; Windows forced to 1 — fork not supported)
N_JOBS = 1 if sys.platform == 'win32' else -1

# ------------------------------------------------------------------------------
# [b] MULTI-PATIENT FIGURE CONFIGURATION:
# ------------------------------------------------------------------------------
MAKE_MULTI_PATIENT_FIGURE = True
# Silently skipped if distinct_patients.csv missing or misconfigured
# Set False to skip panel figure entirely

PANEL_LAYOUT = (3, 4)
# (rows, cols) grid
# Do not need to fill all panels. If too many patients selected, overflow ignored.

PANEL_COLUMN = 'multi'
# Inclusion column name
# Name of a column in 0_input/distinct_patients.csv. (example in header)
# Patients with non-zero values in this column are included.
# ------------------------------------------------------------------------------

# ==============================================================================
# Helper Functions
# ==============================================================================
def darken_color(color_hex, factor=0.7):
    """Darkens a hex color by reducing brightness."""
    if isinstance(color_hex, (list, tuple)):
        rgb = color_hex
    else:
        rgb = mcolors.hex2color(color_hex)
    hsv = mcolors.rgb_to_hsv(rgb)
    hsv[2] = max(0, hsv[2] * factor)
    return mcolors.hsv_to_rgb(hsv)

OKABE_ITO_COLORS_DARKER = [darken_color(c, factor=0.7) for c in OKABE_ITO_COLORS]


def shorten_method_name(method_str):
    """Shortens Elbow method description for display."""
    if not isinstance(method_str, str):
        return str(method_str)
    if "Largest Normalized Cost Drop Fallback" in method_str:
        return "LgNormDrop FB"
    elif "Minimum Cost Fallback" in method_str:
        return "MinCost FB"
    elif "Few Points" in method_str:
        return "Few Pts"
    elif "kneed" in method_str:
        return "Kneed"
    return method_str


def get_axis_limits(pyruvate_series, lactate_series, padding_factor):
    """Returns axis limits for pyruvate (μM) and lactate (mM)."""
    if pyruvate_series.empty or lactate_series.empty:
        return (0, 1000), (0, 1)

    pyr = pd.to_numeric(pyruvate_series, errors='coerce').dropna()
    lac = pd.to_numeric(lactate_series, errors='coerce').dropna()

    if pyr.empty or lac.empty:
        return (0, 1000), (0, 1)

    # Convert pyruvate from mM to μM
    pyr_min_um = pyr.min() * 1000
    pyr_max_um = pyr.max() * 1000
    lac_min, lac_max = lac.min(), lac.max()

    if pd.isna(pyr_min_um) or pd.isna(pyr_max_um) or pd.isna(lac_min) or pd.isna(lac_max):
        return (0, 1000), (0, 1)

    pyr_range = pyr_max_um - pyr_min_um if pyr_max_um > pyr_min_um else 200
    lac_range = lac_max - lac_min if lac_max > lac_min else 0.2

    xlim = (pyr_min_um - pyr_range * padding_factor, pyr_max_um + pyr_range * padding_factor)
    ylim = (lac_min - lac_range * padding_factor, lac_max + lac_range * padding_factor)
    return xlim, ylim


def apply_axis_style(ax):
    """Apply consistent styling to axes."""
    ax.tick_params(axis='both', which='major', labelsize=FONTSIZE_TICKS, colors='black')
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)


def get_panel_label(index):
    """Convert index (0, 1, 2, ...) to panel label (A, B, C, ...)."""
    return chr(ord('A') + index)


# ==============================================================================
# Metadata Text File Generation
# ==============================================================================
def write_metadata_to_file(metadata, txt_path, patient_id):
    """Write metadata dictionary to a text file describing the figure panels."""
    try:
        with open(txt_path, 'w') as f:
            f.write(f"Patient: {patient_id}\n")
            f.write("=" * 50 + "\n\n")

            # Panel descriptions
            f.write("Panel descriptions\n")
            f.write("------------------\n\n")

            # Panel A: Raw data by time
            raw_n = metadata.get('raw_n', 'N/A')
            t_min = metadata.get('raw_time_min', None)
            t_max = metadata.get('raw_time_max', None)
            if isinstance(t_min, (int, float)) and isinstance(t_max, (int, float)):
                span_text = f"from day {t_min:.1f} (blue) to day {t_max:.1f} (yellow)"
            else:
                span_text = "time span not available"
            f.write(f"A) Raw lactate vs pyruvate across time, colored by days since injury using the 'cividis' colormap.\n")
            f.write(f"   n = {raw_n}; time span {span_text}.\n\n")

            # Panel B: PELT penalty optimization
            opt_bkps = metadata.get('optimal_n_bkps', 'N/A')
            opt_cost = metadata.get('optimal_cost', 'N/A')
            elbow_method = metadata.get('elbow_method', 'N/A')
            f.write("B) PELT penalty optimization curve (segmentation cost vs. number of breakpoints).\n")
            f.write(f"   Optimal number of breakpoints: {opt_bkps}; optimal cost: {opt_cost}; Elbow method: {elbow_method}.\n\n")

            # Panel C: All segments with OLS lines
            num_segments = metadata.get('num_segments', 0)
            has_discarded = metadata.get('has_discarded', False)
            f.write("C) All PELT segments with OLS regression lines.\n")
            f.write(f"   Number of segments: {num_segments}.\n")
            if has_discarded:
                f.write("   Solid lines indicate retained segments (m > 0); dashed lines indicate discarded segments (m ≤ 0).\n\n")
            else:
                f.write("   All segments shown with solid lines (all retained, m > 0).\n\n")

            # Panels D & E (only if has_discarded)
            if has_discarded:
                discarded_n = metadata.get('discarded_n', 0)
                discarded_pct = metadata.get('discarded_pct', 0)
                retained_n = metadata.get('retained_n', 0)
                retained_pct = metadata.get('retained_pct', 0)
                
                f.write("D) Discarded points from segments with non-positive gradient (m ≤ 0).\n")
                f.write(f"   Total discarded: {discarded_n} ({discarded_pct:.1f}%).\n")
                f.write("   Retained points shown faded in background; discarded points shown as X markers colored by segment.\n\n")

                f.write("E) Retained segments only (m > 0) with OLS regression lines.\n")
                f.write(f"   Total retained: {retained_n} ({retained_pct:.1f}%).\n\n")

            # Segment details
            segments_info = metadata.get('segments', [])
            if segments_info:
                f.write("Segment-specific details\n")
                f.write("------------------------\n\n")

                for seg_info in segments_info:
                    seg_num = seg_info.get('seg_num', '?')
                    panel_label = seg_info.get('panel_label', '?')
                    n_points = seg_info.get('n_points', 'N/A')
                    pct_points = seg_info.get('pct_points', 0)
                    t_min_seg = seg_info.get('time_min', None)
                    t_max_seg = seg_info.get('time_max', None)
                    m = seg_info.get('m', None)
                    b = seg_info.get('b', None)
                    r = seg_info.get('r', None)
                    is_discarded = seg_info.get('is_discarded', False)
                    status = "DISCARDED" if is_discarded else "RETAINED"

                    f.write(f"Segment {seg_num} (panel {panel_label}) - {status}:\n")
                    
                    if isinstance(t_min_seg, (int, float)) and isinstance(t_max_seg, (int, float)):
                        f.write(f"  Time span: days {t_min_seg:.1f}–{t_max_seg:.1f}.\n")
                    
                    f.write(f"  Points: n = {n_points} ({pct_points:.1f}% of total).\n")
                    
                    if isinstance(m, (int, float)) and isinstance(b, (int, float)):
                        f.write(f"  OLS fit (L = m * P + b), with lactate (L) and pyruvate (P) in mM:\n")
                        f.write(f"    L = {m:.4f} * P + {b:.4f}\n")
                        if isinstance(r, (int, float)):
                            f.write(f"    Pearson r = {r:.4f}\n")
                    
                    f.write("\n")

    except Exception as e:
        print(f"\nError: Writing metadata for patient {patient_id}: {e}")


# ==============================================================================
# Plot Functions
# ==============================================================================
def plot_raw_by_time(ax, df, xlim, ylim, label="A"):
    """Panel A: Raw data colored by time."""
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    t_days = pd.to_numeric(df['time_since_injury'], errors='coerce') / 24.0
    pyruvate_um = df['pyruvate'] * 1000
    ax.scatter(pyruvate_um, df['lactate'], c=t_days, cmap=CMAP_TIME,
               s=BASE_MARKER_SIZE, alpha=0.7)

    ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def plot_elbow_curve(ax, df_history, elbow_method, label="B"):
    """Panel B: PELT penalty optimization (Elbow) curve."""
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df_history.empty:
        ax.text(0.5, 0.5, "No optimization data", ha="center", va="center", transform=ax.transAxes)
        return

    x = df_history['n_bkps']
    y = df_history['min_cost']
    optimal = df_history['optimal']

    ax.plot(x, y, marker='o', linestyle='-', color='darkcyan', alpha=0.7, zorder=1)

    # Highlight optimal point
    knee_data = df_history[optimal == 1]
    if not knee_data.empty:
        knee_x = knee_data['n_bkps'].iloc[0]
        knee_y = knee_data['min_cost'].iloc[0]
        # Get axis limits for bounded lines
        x_min, _ = ax.get_xlim()
        y_min, _ = ax.get_ylim()
        # Vertical line from bottom up to point (not above)
        ax.plot([knee_x, knee_x], [y_min, knee_y], linestyle='--', color='r', alpha=0.7)
        # Horizontal line from left to point (not beyond)
        ax.plot([x_min, knee_x], [knee_y, knee_y], linestyle='--', color='r', alpha=0.7)
        ax.plot(knee_x, knee_y, 'ro', markersize=10, markeredgecolor='black', markeredgewidth=1.5)

    ax.set_xlabel("Number of Breakpoints", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Segmentation Cost", fontsize=FONTSIZE_LABELS)

    x_unique = df_history['n_bkps'].unique()
    if len(x_unique) <= 15:
        ax.set_xticks(sorted(x_unique))
        ax.tick_params(axis='x', labelrotation=45)
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=10))

    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def plot_all_segments_with_lines(ax, df, segment_col, xlim, ylim, label="C"):
    """
    Panel C: All segments with OLS lines.
    Dashed lines for segments with m ≤ 0.
    """
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    segments = sorted(df[segment_col].unique())

    for idx, seg_num in enumerate(segments):
        seg_data = df[df[segment_col] == seg_num]
        if seg_data.empty:
            continue

        color_idx = idx % len(OKABE_ITO_COLORS)
        color = OKABE_ITO_COLORS[color_idx]
        dark_color = OKABE_ITO_COLORS_DARKER[color_idx]

        # Plot points
        pyruvate_um = seg_data['pyruvate'] * 1000
        ax.scatter(pyruvate_um, seg_data['lactate'], color=color,
                   s=BASE_MARKER_SIZE, alpha=0.7)

        # Plot OLS line
        m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
        b = seg_data['p6e_b'].iloc[0] if 'p6e_b' in seg_data.columns else np.nan

        if pd.notna(m) and pd.notna(b) and len(seg_data) >= 2:
            x_vals = seg_data['pyruvate'].values
            x_range = np.array([x_vals.min(), x_vals.max()])
            y_pred = m * x_range + b
            x_range_um = x_range * 1000
            linestyle = '--' if m <= 0 else '-'
            ax.plot(x_range_um, y_pred, color=dark_color, linestyle=linestyle, linewidth=1.5)

    ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def plot_discarded_points(ax, df, segment_col, xlim, ylim, label="D"):
    """
    Panel D: Discarded points (from segments with m ≤ 0).
    Shows retained points faded in background, discarded points with X markers colored by segment.
    """
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    segments = sorted(df[segment_col].unique())

    # Identify retained vs discarded segments
    retained_mask = pd.Series(False, index=df.index)
    discarded_mask = pd.Series(False, index=df.index)

    for seg_num in segments:
        seg_data = df[df[segment_col] == seg_num]
        m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
        if pd.notna(m) and m > 0:
            retained_mask |= (df[segment_col] == seg_num)
        else:
            discarded_mask |= (df[segment_col] == seg_num)

    # Plot retained points faded in background
    retained_df = df[retained_mask]
    if not retained_df.empty:
        pyruvate_um = retained_df['pyruvate'] * 1000
        ax.scatter(pyruvate_um, retained_df['lactate'], color=NEUTRAL_GREY,
                   marker='o', s=BASE_MARKER_SIZE, alpha=0.3, zorder=0)

    # Plot discarded points with X markers, colored by segment
    discarded_df = df[discarded_mask]
    if not discarded_df.empty:
        for idx, seg_num in enumerate(segments):
            seg_data = discarded_df[discarded_df[segment_col] == seg_num]
            if seg_data.empty:
                continue

            color_idx = idx % len(OKABE_ITO_COLORS)
            color = OKABE_ITO_COLORS[color_idx]

            pyruvate_um = seg_data['pyruvate'] * 1000
            ax.scatter(pyruvate_um, seg_data['lactate'], color=color,
                       marker=DISCARDED_MARKER, s=BASE_MARKER_SIZE * 1.5,
                       alpha=0.8, linewidths=1, zorder=1)
    else:
        ax.text(0.5, 0.5, "No discarded points", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def plot_retained_segments(ax, df, segment_col, xlim, ylim, label="E"):
    """
    Panel E: Only retained segments (m > 0) with OLS lines.
    """
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    segments = sorted(df[segment_col].unique())
    has_retained = False

    for idx, seg_num in enumerate(segments):
        seg_data = df[df[segment_col] == seg_num]
        if seg_data.empty:
            continue

        m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
        b = seg_data['p6e_b'].iloc[0] if 'p6e_b' in seg_data.columns else np.nan

        # Skip segments with non-positive gradient
        if pd.isna(m) or m <= 0:
            continue

        has_retained = True
        color_idx = idx % len(OKABE_ITO_COLORS)
        color = OKABE_ITO_COLORS[color_idx]
        dark_color = OKABE_ITO_COLORS_DARKER[color_idx]

        # Plot points
        pyruvate_um = seg_data['pyruvate'] * 1000
        ax.scatter(pyruvate_um, seg_data['lactate'], color=color,
                   s=BASE_MARKER_SIZE, alpha=0.7)

        # Plot OLS line
        if pd.notna(b) and len(seg_data) >= 2:
            x_vals = seg_data['pyruvate'].values
            x_range = np.array([x_vals.min(), x_vals.max()])
            y_pred = m * x_range + b
            x_range_um = x_range * 1000
            ax.plot(x_range_um, y_pred, color=dark_color, linestyle='-', linewidth=1.5)

    if not has_retained:
        ax.text(0.5, 0.5, "No retained segments", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def plot_single_segment(ax, df, seg_num, seg_index, xlim, ylim, label):
    """
    Plot for an individual segment panel: shows segment data with OLS line.
    Discarded segments (m ≤ 0) shown with X markers and dashed line.
    """
    ax.text(-0.1, 1.15, label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    seg_data = df[df['p6e_seg_index'] == seg_num]
    if seg_data.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
    b = seg_data['p6e_b'].iloc[0] if 'p6e_b' in seg_data.columns else np.nan
    is_discarded = pd.isna(m) or m <= 0

    color_idx = seg_index % len(OKABE_ITO_COLORS)
    color = OKABE_ITO_COLORS[color_idx]
    dark_color = OKABE_ITO_COLORS_DARKER[color_idx]

    # Plot points
    pyruvate_um = seg_data['pyruvate'] * 1000
    marker = DISCARDED_MARKER if is_discarded else 'o'
    ax.scatter(pyruvate_um, seg_data['lactate'], color=color,
               marker=marker, s=BASE_MARKER_SIZE, alpha=0.7)

    # Plot OLS line
    if pd.notna(m) and pd.notna(b) and len(seg_data) >= 2:
        x_vals = seg_data['pyruvate'].values
        x_range = np.array([x_vals.min(), x_vals.max()])
        y_pred = m * x_range + b
        x_range_um = x_range * 1000
        linestyle = '--' if is_discarded else '-'
        ax.plot(x_range_um, y_pred, color=dark_color, linestyle=linestyle, linewidth=1.5)

    ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


# ==============================================================================
# Panel Plot Function for Multi-patient Figure
# ==============================================================================
def plot_patient_final_for_panel(ax, df_patient, panel_label, show_xlabel=True, show_ylabel=True):
    """
    Create the plot for one patient's panel in a multi-patient figure.
    Discarded segments shown as gray X markers with dashed lines.
    """
    # Add panel label
    ax.text(-0.1, 1.15, panel_label, transform=ax.transAxes,
            fontsize=FONTSIZE_PANEL_LABEL, fontweight='bold', va='top', ha='left')

    if df_patient.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    # Get axis limits
    xlim, ylim = get_axis_limits(df_patient['pyruvate'], df_patient['lactate'], AXIS_PADDING_FACTOR)

    # Plot each segment
    segments = sorted(df_patient['p6e_seg_index'].unique())
    for seg_idx, seg_num in enumerate(segments):
        seg_data = df_patient[df_patient['p6e_seg_index'] == seg_num]
        m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
        b = seg_data['p6e_b'].iloc[0] if 'p6e_b' in seg_data.columns else np.nan
        is_discarded = pd.isna(m) or m <= 0
        
        color_idx = seg_idx % len(OKABE_ITO_COLORS)
        color = OKABE_ITO_COLORS[color_idx]
        dark_color = OKABE_ITO_COLORS_DARKER[color_idx]
        
        pyruvate_um = seg_data['pyruvate'] * 1000
        
        if is_discarded:
            # Discarded: gray X markers, dashed line
            ax.scatter(pyruvate_um, seg_data['lactate'], color=NEUTRAL_GREY,
                       marker=DISCARDED_MARKER, s=BASE_MARKER_SIZE * 1.2,
                       alpha=0.6, linewidths=0.8, zorder=1)
            if pd.notna(m) and pd.notna(b) and len(seg_data) >= 2:
                x_vals = seg_data['pyruvate'].values
                x_range = np.array([x_vals.min(), x_vals.max()])
                y_pred = m * x_range + b
                x_range_um = x_range * 1000
                ax.plot(x_range_um, y_pred, color=NEUTRAL_GREY, linestyle='--', linewidth=1.2, alpha=0.6)
        else:
            # Retained: colored circles, solid line
            ax.scatter(pyruvate_um, seg_data['lactate'], color=color,
                       s=BASE_MARKER_SIZE, alpha=0.7, zorder=2)
            if pd.notna(m) and pd.notna(b) and len(seg_data) >= 2:
                x_vals = seg_data['pyruvate'].values
                x_range = np.array([x_vals.min(), x_vals.max()])
                y_pred = m * x_range + b
                x_range_um = x_range * 1000
                ax.plot(x_range_um, y_pred, color=dark_color, linestyle='-', linewidth=1.5, zorder=3)

    # Axis labels (conditional)
    if show_xlabel:
        ax.set_xlabel("Pyruvate (μM)", fontsize=FONTSIZE_LABELS)
    if show_ylabel:
        ax.set_ylabel("Lactate (mM)", fontsize=FONTSIZE_LABELS)
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(axis='both', which='major', labelsize=FONTSIZE_TICKS)
    ax.grid(True, linestyle='-', alpha=0.7)
    apply_axis_style(ax)


def create_segmentation_panel(patient_ids, layout, df_results, output_path):
    """
    Creates a multi-panel figure showing the final segmentation for multiple patients.
    Each patient occupies one panel (labeled A, B, C, ...).
    Sorted by number of segments (ascending).
    Returns the number of patients plotted.
    """
    rows, cols = layout
    
    # Calculate segment counts for sorting
    segment_counts = {}
    for pid in patient_ids:
        df_patient = df_results[df_results['patient'] == pid]
        if not df_patient.empty:
            segment_counts[pid] = df_patient['p6e_seg_index'].nunique()
        else:
            segment_counts[pid] = 0
    
    # Sort by number of segments (ascending)
    sorted_patient_ids = sorted(patient_ids, key=lambda pid: segment_counts.get(pid, 0))
    
    # Figure setup
    fig_width = DOUBLE_COLUMN_WIDTH_INCH
    fig_height = PANEL_HEIGHT * rows  # Standard panel height per row
    
    fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), 
                             constrained_layout=False, squeeze=False)
    axes_flat = axes.flatten()
    
    # Track number of patients plotted for legend
    n_patients_plotted = 0
    
    for i, patient_id in enumerate(sorted_patient_ids):
        if i >= len(axes_flat):
            break
        
        ax = axes_flat[i]
        panel_label = get_panel_label(i)
        
        # Determine axis label visibility
        is_leftmost = (i % cols == 0)
        is_bottom = (i + cols >= len(sorted_patient_ids)) or (i // cols == rows - 1)
        
        df_patient = df_results[df_results['patient'] == patient_id].copy()
        
        if df_patient.empty:
            ax.text(0.5, 0.5, f"Patient {patient_id}\nNo data", ha="center", va="center", 
                    transform=ax.transAxes, fontsize=FONTSIZE_TICKS)
            ax.axis('off')
            continue
        
        # Plot panel
        plot_patient_final_for_panel(
            ax, df_patient, panel_label,
            show_xlabel=is_bottom, show_ylabel=is_leftmost
        )
        n_patients_plotted += 1
    
    # Hide unused subplots
    for j in range(len(sorted_patient_ids), len(axes_flat)):
        axes_flat[j].axis('off')
    
    # Adjust layout
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.subplots_adjust(hspace=0.35, wspace=0.25)
    
    # Save figure
    try:
        fig.savefig(output_path, dpi=PLOT_DPI)
    except Exception as e:
        print(f"\nError: Saving panel: {e}")
    finally:
        plt.close(fig)
    
    return n_patients_plotted


def build_panel_figure_legend(n_patients: int) -> list:
    """Build legend list for the main panel figure."""
    return [
        {
            "target": "Figure 1 (1_figure.png)",
            "caption": (
                f"Lactate-pyruvate segmentation results for {n_patients} representative patients. "
                f"Each panel shows one patient's segment-level scatter plots with OLS regression lines. "
                f"Colored circles and solid lines indicate retained segments (m > 0); "
                f"gray X markers and dashed lines indicate discarded segments (m ≤ 0). "
                f"Segment colors follow the Okabe-Ito palette; panels are ordered by number of segments."
            ),
            "abbreviations": "OLS, ordinary least squares; m, gradient; b, intercept",
        },
    ]


# ==============================================================================
# Main Visualization Function
# ==============================================================================
def create_patient_visualization(patient_id, df_patient, df_history, output_path):
    """
    Creates a multi-panel figure for a single patient.
    Row 1: Summary panels (A, B, C, [D, E if has discarded])
    Row 2+: Individual segment panels (F, G, H, ...)
    Also writes a metadata text file with the same name as the PNG.
    """
    if df_patient.empty:
        return False

    # Initialize metadata dictionary
    metadata = {}

    # Get axis limits
    xlim, ylim = get_axis_limits(df_patient['pyruvate'], df_patient['lactate'], AXIS_PADDING_FACTOR)

    # Raw data stats
    total_points = len(df_patient)
    metadata['raw_n'] = total_points
    t_hours = pd.to_numeric(df_patient['time_since_injury'], errors='coerce').dropna()
    if not t_hours.empty:
        t_days = t_hours / 24.0
        metadata['raw_time_min'] = float(t_days.min())
        metadata['raw_time_max'] = float(t_days.max())

    # Get segments and check for discarded
    segments = sorted(df_patient['p6e_seg_index'].unique())
    num_segments = len(segments)
    metadata['num_segments'] = num_segments

    has_discarded = False
    discarded_points = 0
    retained_points = 0
    segments_info = []

    for seg_idx, seg_num in enumerate(segments):
        seg_data = df_patient[df_patient['p6e_seg_index'] == seg_num]
        m = seg_data['p6e_m'].iloc[0] if 'p6e_m' in seg_data.columns else np.nan
        b = seg_data['p6e_b'].iloc[0] if 'p6e_b' in seg_data.columns else np.nan
        r = seg_data['p6e_r'].iloc[0] if 'p6e_r' in seg_data.columns else np.nan
        
        is_seg_discarded = pd.isna(m) or m <= 0
        if is_seg_discarded:
            has_discarded = True
            discarded_points += len(seg_data)
        else:
            retained_points += len(seg_data)

        # Time span for segment
        t_seg_hours = pd.to_numeric(seg_data['time_since_injury'], errors='coerce').dropna()
        t_min_seg = float(t_seg_hours.min() / 24.0) if not t_seg_hours.empty else None
        t_max_seg = float(t_seg_hours.max() / 24.0) if not t_seg_hours.empty else None

        segments_info.append({
            'seg_num': int(seg_num),
            'seg_idx': seg_idx,
            'n_points': len(seg_data),
            'pct_points': (len(seg_data) / total_points * 100) if total_points > 0 else 0,
            'time_min': t_min_seg,
            'time_max': t_max_seg,
            'm': float(m) if pd.notna(m) else None,
            'b': float(b) if pd.notna(b) else None,
            'r': float(r) if pd.notna(r) else None,
            'is_discarded': is_seg_discarded,
            'panel_label': None  # Will be set later
        })

    metadata['has_discarded'] = has_discarded
    metadata['discarded_n'] = discarded_points
    metadata['discarded_pct'] = (discarded_points / total_points * 100) if total_points > 0 else 0
    metadata['retained_n'] = retained_points
    metadata['retained_pct'] = (retained_points / total_points * 100) if total_points > 0 else 0

    # Get elbow method and optimal breakpoints
    elbow_method = "Unknown"
    if 'p6e_elbow_method' in df_patient.columns:
        method_vals = df_patient['p6e_elbow_method'].dropna()
        if not method_vals.empty:
            elbow_method = shorten_method_name(method_vals.iloc[0])
    metadata['elbow_method'] = elbow_method

    # Get optimal breakpoints from history
    if not df_history.empty:
        optimal_row = df_history[df_history['optimal'] == 1]
        if not optimal_row.empty:
            metadata['optimal_n_bkps'] = int(optimal_row['n_bkps'].iloc[0])
            metadata['optimal_cost'] = f"{optimal_row['min_cost'].iloc[0]:.2f}"
        else:
            metadata['optimal_n_bkps'] = 'N/A'
            metadata['optimal_cost'] = 'N/A'
    else:
        metadata['optimal_n_bkps'] = 'N/A'
        metadata['optimal_cost'] = 'N/A'

    # Determine layout based on whether patient has discarded segments
    # - With discarded: 5 cols (A, B, C, D, E in row 1, segments wrap at 5)
    # - All retained: 3 cols (A, B, C in row 1, segments wrap at 3)
    max_cols = MAX_COLS_DISCARDED if has_discarded else MAX_COLS_RETAINED

    # Row 2+: Individual segments (up to max_cols per row)
    segment_rows = (num_segments + max_cols - 1) // max_cols
    total_rows = 1 + segment_rows

    # Figure dimensions
    fig_width = DOUBLE_COLUMN_WIDTH_INCH
    fig_height = PANEL_HEIGHT * total_rows

    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.04, hspace=0.08, wspace=0.1)
    gs = GridSpec(total_rows, max_cols, figure=fig)

    # ========== Row 1: Summary panels ==========
    panel_idx = 0

    # Panel A: Raw data by time
    ax_a = fig.add_subplot(gs[0, 0])
    plot_raw_by_time(ax_a, df_patient, xlim, ylim, label=get_panel_label(panel_idx))
    panel_idx += 1

    # Panel B: Elbow curve
    ax_b = fig.add_subplot(gs[0, 1])
    plot_elbow_curve(ax_b, df_history, elbow_method, label=get_panel_label(panel_idx))
    panel_idx += 1

    # Panel C: All segments with OLS lines
    ax_c = fig.add_subplot(gs[0, 2])
    plot_all_segments_with_lines(ax_c, df_patient, 'p6e_seg_index', xlim, ylim,
                                 label=get_panel_label(panel_idx))
    panel_idx += 1

    if has_discarded:
        # Panel D: Discarded points
        ax_d = fig.add_subplot(gs[0, 3])
        plot_discarded_points(ax_d, df_patient, 'p6e_seg_index', xlim, ylim,
                              label=get_panel_label(panel_idx))
        panel_idx += 1

        # Panel E: Retained segments only
        ax_e = fig.add_subplot(gs[0, 4])
        plot_retained_segments(ax_e, df_patient, 'p6e_seg_index', xlim, ylim,
                               label=get_panel_label(panel_idx))
        panel_idx += 1

    # ========== Row 2+: Individual segment panels ==========
    for seg_idx, seg_num in enumerate(segments):
        row = 1 + seg_idx // max_cols
        col = seg_idx % max_cols

        ax_seg = fig.add_subplot(gs[row, col])
        current_label = get_panel_label(panel_idx)
        plot_single_segment(ax_seg, df_patient, seg_num, seg_idx, xlim, ylim,
                            label=current_label)
        
        # Update segment info with panel label
        segments_info[seg_idx]['panel_label'] = current_label
        panel_idx += 1

    # Turn off unused cells in segment rows
    last_seg_row = 1 + (num_segments - 1) // max_cols
    last_seg_col = (num_segments - 1) % max_cols
    for row in range(1, total_rows):
        for col in range(max_cols):
            if row > last_seg_row or (row == last_seg_row and col > last_seg_col):
                ax_empty = fig.add_subplot(gs[row, col])
                ax_empty.axis('off')

    # Store segments info in metadata
    metadata['segments'] = segments_info

    txt_output_dir = output_path.parent

    try:
        fig.savefig(output_path, dpi=PLOT_DPI)
        # Write metadata text file using patient-specific naming
        txt_path = txt_output_dir / f"{patient_id}_notes.txt"
        write_metadata_to_file(metadata, txt_path, patient_id)
        
    except Exception as e:
        print(f"\nError: Saving plot for patient {patient_id}: {e}")
        return False
    finally:
        plt.close(fig)

    return has_discarded


def visualization_worker(args):
    """Worker function for parallel processing."""
    patient_id, df_patient, df_history, output_path = args
    return patient_id, create_patient_visualization(patient_id, df_patient, df_history, output_path)


def generate_visualization_tasks(patient_ids, df_results, df_history_all, output_base_dir):
    """Generator that yields task arguments one by one to save memory."""
    for pid in patient_ids:
        df_patient = df_results[df_results['patient'] == pid].copy()
        df_history = df_history_all[df_history_all['patient'] == pid].copy()
        output_path = output_base_dir / f"{pid}.png"
        yield (pid, df_patient, df_history, output_path)


# ==============================================================================
# Main
# ==============================================================================
def main():
    print("\n\nStarting 3_segmentation_visual.py")
    
    # Initialize panel metadata (will be populated if panel figure is generated)
    panel_metadata = None

    # Load results
    if not RESULTS_FILE.exists():
        print(f"\nError: Results file not found at {RESULTS_FILE}")
        sys.exit(1)

    try:
        df_results = pd.read_csv(RESULTS_FILE)
        df_results['patient'] = df_results['patient'].astype(str)
    except Exception as e:
        print(f"\nError: Loading results: {e}")
        sys.exit(1)

    # Create output directories
    OUTPUT_ALL_PATIENTS.mkdir(parents=True, exist_ok=True)
    OUTPUT_NO_BREAKPOINTS.mkdir(parents=True, exist_ok=True)
    OUTPUT_HAS_DISCARDED.mkdir(parents=True, exist_ok=True)
    OUTPUT_DISTINCT.mkdir(parents=True, exist_ok=True)

    # Load distinct patient IDs if file exists
    distinct_patient_ids = set()
    if DISTINCT_PATIENTS_FILE.exists():
        try:
            df_distinct = pd.read_csv(DISTINCT_PATIENTS_FILE)
            if 'patient' in df_distinct.columns:
                distinct_patient_ids = set(df_distinct['patient'].astype(str).unique())
        except Exception:
            pass

    # Build elbow history from p6e_cost_history
    df_history_all = build_elbow_history(df_results)

    # Generate panel figure if enabled
    if MAKE_MULTI_PATIENT_FIGURE and DISTINCT_PATIENTS_FILE.exists():
        try:
            df_panel = pd.read_csv(DISTINCT_PATIENTS_FILE)
            if 'patient' in df_panel.columns and PANEL_COLUMN in df_panel.columns:
                df_panel[PANEL_COLUMN] = pd.to_numeric(df_panel[PANEL_COLUMN], errors='coerce').fillna(0)
                panel_patient_ids = df_panel[df_panel[PANEL_COLUMN] != 0]['patient'].astype(str).tolist()
                if panel_patient_ids:
                    max_patients = PANEL_LAYOUT[0] * PANEL_LAYOUT[1]
                    if len(panel_patient_ids) > max_patients:
                        panel_patient_ids = panel_patient_ids[:max_patients]
                    figure_output_path = OUTPUT_DIR / "1_figure.png"
                    panel_metadata = create_segmentation_panel(panel_patient_ids, PANEL_LAYOUT, df_results, figure_output_path)
        except Exception:
            panel_metadata = None

    # Prepare tasks generator
    patient_ids = sorted(df_results['patient'].unique())

    # Run visualization in parallel using imap for memory efficiency
    results = []
    num_processes = N_JOBS if N_JOBS != -1 else None
    
    task_gen = generate_visualization_tasks(patient_ids, df_results, df_history_all, OUTPUT_ALL_PATIENTS)

    with mp.Pool(processes=num_processes) as pool:
        # imap consumes the generator lazily, keeping memory usage low
        results = list(pool.imap(visualization_worker, task_gen))

    # Categorize outputs
    for pid, has_discarded in results:
        source_path = OUTPUT_ALL_PATIENTS / f"{pid}.png"
        source_txt_path = OUTPUT_ALL_PATIENTS / f"{pid}_notes.txt"
        if not source_path.exists():
            continue

        df_patient = df_results[df_results['patient'] == pid]

        # Check number of breakpoints
        total_segs = pd.to_numeric(df_patient['p6e_total_segs'], errors='coerce').dropna()
        num_bkps = int(total_segs.iloc[0]) - 1 if not total_segs.empty else -1

        try:
            if num_bkps == 0:
                shutil.copy2(source_path, OUTPUT_NO_BREAKPOINTS / f"{pid}.png")
                if source_txt_path.exists():
                    shutil.copy2(source_txt_path, OUTPUT_NO_BREAKPOINTS / f"{pid}_notes.txt")
            if has_discarded:
                shutil.copy2(source_path, OUTPUT_HAS_DISCARDED / f"{pid}.png")
                if source_txt_path.exists():
                    shutil.copy2(source_txt_path, OUTPUT_HAS_DISCARDED / f"{pid}_notes.txt")
            
            # Copy distinct patients
            if pid in distinct_patient_ids:
                shutil.copy2(source_path, OUTPUT_DISTINCT / f"{pid}.png")
                if source_txt_path.exists():
                    shutil.copy2(source_txt_path, OUTPUT_DISTINCT / f"{pid}_notes.txt")
        except Exception as e:
            print(f"\nError: Copying plot for patient {pid}: {e}")

    print("\n3_segmentation_visual.py complete.\n")
    
    return panel_metadata


def build_elbow_history(df_results):
    """Parse p6e_cost_history to build elbow optimization history."""
    rows = []
    for pid, grp in df_results.groupby('patient'):
        # Get selected number of breakpoints
        selected_n_bkps = None
        if 'p6e_total_segs' in grp.columns:
            total_segs = pd.to_numeric(grp['p6e_total_segs'], errors='coerce').dropna()
            if not total_segs.empty:
                selected_n_bkps = int(total_segs.iloc[0]) - 1

        # Parse history string
        if 'p6e_cost_history' not in grp.columns:
            continue
        hist_series = grp['p6e_cost_history'].dropna()
        if hist_series.empty:
            continue

        hist_str = str(hist_series.iloc[0]).strip()
        if not hist_str:
            continue

        # Remove header if present
        if hist_str.startswith('axes='):
            semicolon_idx = hist_str.find(';')
            if semicolon_idx != -1:
                hist_str = hist_str[semicolon_idx + 1:]

        # Parse tuples
        for tok in [t.strip() for t in hist_str.split(';') if t.strip()]:
            if not (tok.startswith('(') and tok.endswith(')')):
                continue
            inner = tok[1:-1]
            parts = inner.split(',')
            if len(parts) >= 2:
                try:
                    k = int(float(parts[0]))
                    cost = float(parts[1])
                    rows.append({
                        'patient': pid,
                        'n_bkps': k,
                        'min_cost': cost,
                        'optimal': 1 if (selected_n_bkps is not None and k == selected_n_bkps) else 0
                    })
                except Exception:
                    pass

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['patient', 'n_bkps', 'min_cost', 'optimal'])


if __name__ == "__main__":
    # Multiprocessing start method (analogous to R's mclapply behavior):
    # - macOS/Linux: use 'fork' (copy-on-write, efficient)
    # - Windows: defaults to 'spawn' (no fork support; slower but functional)
    if sys.platform != 'win32':
        mp.set_start_method('fork', force=True)

    start_time = time.time()
    panel_metadata = main()
    execution_time = time.time() - start_time
    
    # Save execution time
    pd.DataFrame({'execution_time_seconds': [execution_time]}).to_csv(
        OUTPUT_DIR / "__execution_time.csv", index=False
    )
    
    # Write figure legend notes
    if panel_metadata is not None:
        legends = build_panel_figure_legend(panel_metadata)
        notes_path = OUTPUT_DIR / "0_notes.txt"
        notes_path.write_text(build_notes(legends, title="3_segmentation_visual"), encoding="utf-8")
