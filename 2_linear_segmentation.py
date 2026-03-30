# ==============================================================================
# Script: 2_linear_segmentation.py
# Manuscript relevance: 2.3.ii
# ==============================================================================
# PURPOSE:
#   Detect changepoints in lactate–pyruvate linear regression using PELT (Pruned
#   Exact Linear Time), optimize segment count via the Elbow method, and fit
#   OLS regression to both patient-level and segment-level data to extract
#   linear model parameters (gradient m, intercept b, Pearson r).
#
# INPUT:
#   - 1_output/1_pre-processing/1_data.csv: Pre-processed data
#
# OUTPUT:
#   - 1_output/2_linear_segmentation/__execution_time.csv: Runtime log
#   - 1_output/2_linear_segmentation/__warning.csv: Created if datapoints dropped
#   - 1_output/2_linear_segmentation/1_results.csv: Segment-annotated data with
#       patient_* (patient-level) and p6e_* (segment-level) columns.
#
# PIPELINE:
# ------------------------------------------------------------------------------
# (1) Patient-Level OLS (patient_* columns)
#     - Fit a single linear model (L ~ P) across each patient's datapoints
#     Approach:   OLS regression per patient: Lactate = m·Pyruvate + b
#
# (2) PELT Changepoint Detection
#     - Identify changepoints in the L–P linear regression over time
#     Approach:   PELT algorithm with CostLinear cost function, performing a
#                 recursive search over a range of penalties for unique solutions
#     Constraint: Minimum segment length = 6 datapoints
#
# (3) Elbow Optimization
#     - Select optimal number of breakpoints (k)
#     Approach:   Find elbow point on a cost-complexity curve using Kneedle
#                 algorithm with fallbacks.
#
# (4) Segment-Level OLS (p6e_* columns)
#     - Fit linear models within each Elbow-optimized segment
#     Approach:   OLS regression per segment: Lactate = m·Pyruvate + b
# ==============================================================================

import os
import sys
import time
import warnings
import multiprocessing
from typing import Dict, List, Tuple

# Set thread limits before importing numerical libraries
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import ruptures as rpt
import kneed
from sklearn.linear_model import LinearRegression
from scipy import stats


# ==============================================================================
# Input / Output
# ==============================================================================
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR / "1_output" / "1_pre-processing" / "1_data.csv"
OUTPUT_DIR = SCRIPT_DIR / "1_output" / "2_linear_segmentation"
OUTPUT_FILE = OUTPUT_DIR / "1_results.csv"
TIME_FILE = OUTPUT_DIR / "__execution_time.csv"
WARNING_FILE = OUTPUT_DIR / "__warning.csv"

# ==============================================================================
# Configuration / Toggles
# ==============================================================================
MIN_SEG_LEN = 6  # Minimum points per segment
N_JOBS_PARALLEL = 1 if sys.platform == 'win32' else multiprocessing.cpu_count()

# ==============================================================================
# Utilities
# ==============================================================================
def _safe_numeric(series: pd.Series) -> pd.Series:
    """Convert series to numeric, coercing errors to NaN."""
    return pd.to_numeric(series, errors="coerce")

def compute_ols_stats(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Fit OLS (y ~ x) and return (m, b, r, sse).
    
    Returns:
        m: gradient
        b: intercept
        r: Pearson correlation coefficient
        sse: sum of squared errors (for PELT cost function)
    """
    if len(x) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan, np.nan
    
    # Remove non-finite values
    mask = np.isfinite(x) & np.isfinite(y)
    x_clean = x[mask]
    y_clean = y[mask]
    
    if len(x_clean) < 2:
        return np.nan, np.nan, np.nan, np.nan
    
    try:
        # Fit OLS
        X = x_clean.reshape(-1, 1)
        model = LinearRegression().fit(X, y_clean)
        m = float(model.coef_[0])
        b = float(model.intercept_)
        
        # Calculate Pearson r
        r, _ = stats.pearsonr(x_clean, y_clean)
        r = float(r)
        
        # Calculate SSE (for PELT cost function)
        y_pred = model.predict(X)
        sse = float(np.sum((y_clean - y_pred) ** 2))
        
        return m, b, r, sse
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

# ==============================================================================
# PELT Helpers
# ==============================================================================
def compute_segment_sse(segment_pyruvate: np.ndarray, segment_lactate: np.ndarray) -> float:
    """Compute SSE for a segment (used by PELT cost function)."""
    _, _, _, sse = compute_ols_stats(segment_pyruvate, segment_lactate)
    return sse if np.isfinite(sse) else np.inf

def _normalize_breakpoints(bkps_raw, n_points: int) -> List[int]:
    """Ensure valid, sorted breakpoints ending at n_points."""
    if not isinstance(bkps_raw, list):
        return [n_points]
    valid = {b for b in bkps_raw if b <= n_points}
    valid.add(n_points)
    return sorted(valid)

def _discover_solutions_with_pelt(signal: np.ndarray,
                                  min_seg_len: int,
                                  pen_upper: float) -> Dict[Tuple[int, ...], float]:
    """
    Run PELT over penalty range [0, pen_upper] using recursive search.
    Returns map: breakpoints_tuple -> generating_penalty.

    Requires pen_upper equal to the unsegmented SSE so zero-breakpoint solutions remain accessible.
    """
    algo = rpt.Pelt(model="linear", min_size=min_seg_len).fit(signal)

    pen_low = 0.0
    pen_high = pen_upper

    discovered_bkps_map: Dict[Tuple[int, ...], float] = {}
    queue: List[Tuple[float, float, Tuple[int, ...], Tuple[int, ...]]] = []

    # Endpoints
    try:
        bkps_low_raw = algo.predict(pen=pen_low)
        bkps_low_tuple = tuple(_normalize_breakpoints(bkps_low_raw, signal.shape[0]))
        discovered_bkps_map[bkps_low_tuple] = pen_low
    except Exception:
        bkps_low_tuple = tuple()

    try:
        bkps_high_raw = algo.predict(pen=pen_high)
        bkps_high_tuple = tuple(_normalize_breakpoints(bkps_high_raw, signal.shape[0]))
        discovered_bkps_map[bkps_high_tuple] = pen_high
    except Exception:
        bkps_high_tuple = tuple()

    if bkps_low_tuple and bkps_high_tuple and bkps_low_tuple != bkps_high_tuple:
        queue.append((pen_low, pen_high, bkps_low_tuple, bkps_high_tuple))

    # Recursive bisection
    while queue:
        p1, p2, bkps1, bkps2 = queue.pop(0)
        if abs(p2 - p1) < 1e-9:
            continue
        p_mid = p1 + (p2 - p1) / 2.0
        if not (p1 < p_mid < p2):
            continue
        try:
            bkps_mid_raw = algo.predict(pen=p_mid)
            bkps_mid_tuple = tuple(_normalize_breakpoints(bkps_mid_raw, signal.shape[0]))
            if bkps_mid_tuple not in discovered_bkps_map:
                discovered_bkps_map[bkps_mid_tuple] = p_mid
        except Exception:
            eps = 1e-9
            left = max(p1, p_mid - eps)
            right = min(p2, p_mid + eps)
            if left < p_mid:
                queue.append((p1, left, bkps1, bkps1))
            if right > p_mid:
                queue.append((right, p2, bkps2, bkps2))
            continue

        if bkps1 != bkps_mid_tuple:
            queue.append((p1, p_mid, bkps1, bkps_mid_tuple))
        if bkps_mid_tuple != bkps2:
            queue.append((p_mid, p2, bkps_mid_tuple, bkps2))

    return discovered_bkps_map

def _process_single_patient_pelt(patient_df: pd.DataFrame,
                                 patient_id: str,
                                 min_seg_len: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run PELT for one patient.
    Returns (stats_df, unfiltered_df).
    """
    patient_df = patient_df.sort_values(by="time_since_injury").copy()
    lact = patient_df["lactate"].values
    pyru = patient_df["pyruvate"].values
    n_points = len(lact)

    if n_points < max(2, min_seg_len):
        # Handle patients without enough datapoints for segmentation
        sse = compute_segment_sse(pyru, lact)
        stats_df = pd.DataFrame([{
            "patient": patient_id,
            "num_breakpoints": 0,
            "cost_of_segmentation": sse if np.isfinite(sse) else np.nan,
            "generating_penalty": 0.0,
            "pen_upper": sse if np.isfinite(sse) and sse > 0 else 1.0
        }])
        tmp = patient_df.copy()
        tmp["solution_id"] = 0
        tmp["segment_number"] = 1
        return stats_df, tmp

    # PELT search
    # CostLinear signal layout: [dependent_var, regressor_1, ..., regressor_p]
    # Lactate = m * pyruvate + b → [lactate, pyruvate, ones] where ones encodes the intercept
    signal = np.column_stack((lact, pyru, np.ones(n_points)))
    
    # Use global SSE so penalty search covers zero-breakpoint solutions
    # Penalties >= global_SSE prevent additional changepoints
    global_sse = compute_segment_sse(pyru, lact)
    pen_upper = global_sse if np.isfinite(global_sse) and global_sse > 0 else 1.0
    
    discovered = _discover_solutions_with_pelt(signal, min_seg_len=min_seg_len, pen_upper=pen_upper)

    solutions_rows = []
    unfiltered_solution_rows = []

    solution_id_counter = 0
    for bkps_tuple, gen_pen in sorted(discovered.items(), key=lambda kv: (len(kv[0]) - 1, kv[1], kv[0])):
        if not bkps_tuple or bkps_tuple[-1] != n_points:
            continue
        boundaries = [0] + list(bkps_tuple)
        total_cost = 0.0
        seg_details = []
        valid = True
        for i in range(len(boundaries) - 1):
            start = boundaries[i]
            end = boundaries[i + 1]
            if start >= end:
                valid = False
                break
            sse = compute_segment_sse(pyru[start:end], lact[start:end])
            if not np.isfinite(sse):
                valid = False
                break
            total_cost += sse
            seg_details.append((start, end))

        current_solution_id = solution_id_counter
        solutions_rows.append({
            "patient": patient_id,
            "solution_id": current_solution_id,
            "num_breakpoints": len(bkps_tuple) - 1,
            "cost_of_segmentation": total_cost if valid else np.nan,
            "generating_penalty": gen_pen,
            "pen_upper": pen_upper
        })

        if valid and seg_details:
            tmp_rows = []
            for seg_idx, (start, end) in enumerate(seg_details, start=1):
                seg_df = patient_df.iloc[start:end].copy()
                seg_df["solution_id"] = current_solution_id
                seg_df["segment_number"] = seg_idx
                tmp_rows.append(seg_df)
            if tmp_rows:
                unfiltered_solution_rows.append(pd.concat(tmp_rows, ignore_index=False))
        solution_id_counter += 1

    stats_df = pd.DataFrame(solutions_rows)
    unfiltered_df = pd.concat(unfiltered_solution_rows, ignore_index=True) if unfiltered_solution_rows else pd.DataFrame()
    return stats_df, unfiltered_df

def _pelt_worker(args_tuple):
    """Worker for multiprocessing."""
    pid, pdf, min_seg_len = args_tuple
    return _process_single_patient_pelt(pdf, pid, min_seg_len)

def run_pelt_processing(input_df: pd.DataFrame, min_seg_len: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run PELT on all patients in parallel."""
    processed = input_df.copy()
    processed.columns = [c.lower().strip() for c in processed.columns]

    required_cols = ["patient", "lactate", "pyruvate", "time_since_injury"]
    if not all(col in processed.columns for col in required_cols):
        missing = [c for c in required_cols if c not in processed.columns]
        raise ValueError(f"Input must contain: {', '.join(required_cols)}. Missing: {', '.join(missing)}")

    processed["lactate"] = _safe_numeric(processed["lactate"])
    processed["pyruvate"] = _safe_numeric(processed["pyruvate"])
    processed["time_since_injury"] = _safe_numeric(processed["time_since_injury"])
    processed = processed.dropna(subset=["patient", "time_since_injury", "lactate", "pyruvate"])
    processed["patient"] = processed["patient"].astype(str)

    patient_ids = processed["patient"].unique()
    tasks = [(pid, processed[processed["patient"] == pid].copy(), min_seg_len) for pid in patient_ids]

    pelt_stats_all = []
    pelt_unfiltered_all = []

    if tasks:
        with multiprocessing.Pool(processes=min(N_JOBS_PARALLEL, len(tasks))) as pool:
            results = pool.map(_pelt_worker, tasks)
        for (_, _, _), (stats_df, unfiltered_df) in zip(tasks, results):
            if stats_df is not None and not stats_df.empty:
                pelt_stats_all.append(stats_df)
            if unfiltered_df is not None and not unfiltered_df.empty:
                pelt_unfiltered_all.append(unfiltered_df)

    pelt_stats_df = pd.concat(pelt_stats_all, ignore_index=True) if pelt_stats_all else pd.DataFrame()
    pelt_unfiltered_df = pd.concat(pelt_unfiltered_all, ignore_index=True) if pelt_unfiltered_all else pd.DataFrame()
    return pelt_stats_df, pelt_unfiltered_df

# ==============================================================================
# Elbow Method
# ==============================================================================
def _elbow_for_patient(x: np.ndarray, y: np.ndarray) -> Tuple[int, float, str]:
    """Find elbow using kneed with fallbacks."""
    n = len(x)
    if n == 0:
        return (None, None, "None")
    if n == 1:
        return (int(x[0]), float(y[0]), "Few Points (1)")
    if n == 2:
        idx = 0 if y[0] <= y[1] else 1
        return (int(x[idx]), float(y[idx]), "Few Points (2)")

    # Normalize y to [0, 1]
    y_min, y_max = np.min(y), np.max(y)
    y_norm = (y - y_min) / (y_max - y_min) if (y_max - y_min) > 0 else np.zeros_like(y)

    try:
        kl = kneed.KneeLocator(x, y_norm, S=1.0, curve="convex", direction="decreasing", interp_method="interp1d", online=False)
        if kl.knee is not None:
            elbow_val = int(kl.knee)
            elbow_idx = np.where(x == elbow_val)[0][0]
            return (elbow_val, float(y[elbow_idx]), "kneed (S=1.0, normalized)")
    except Exception:
        pass

    # Fallback 1: Largest normalized drop
    delta_k = x[1:] - x[:-1]
    raw_drops = y[:-1] - y[1:]
    valid = np.where(delta_k > 0)[0]
    if len(valid) > 0:
        normalized_drops = raw_drops[valid] / delta_k[valid]
        total_range = (np.max(y) - np.min(y)) if len(y) > 0 else 0.0
        if len(normalized_drops) > 0 and np.any(normalized_drops > 1e-9):
            best_drop_idx = valid[np.argmax(normalized_drops)]
            significant = (total_range > 1e-9 and raw_drops[best_drop_idx] > (0.01 * total_range)) or \
                          (total_range <= 1e-9 and raw_drops[best_drop_idx] > 1e-9)
            if significant and (best_drop_idx + 1) < len(x):
                return (int(x[best_drop_idx + 1]), float(y[best_drop_idx + 1]), "Largest Normalized Cost Drop Fallback")

    # Fallback 2: Minimum cost
    if len(y) > 0:
        min_cost_val = np.min(y)
        idxs = np.where(np.isclose(y, min_cost_val))[0]
        final_idx = idxs[np.argmin(x[idxs])]
        return (int(x[final_idx]), float(y[final_idx]), "Minimum Cost Fallback")

    return (None, None, "None")

def run_elbow_processing(pelt_stats_df: pd.DataFrame,
                         pelt_unfiltered_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Select optimal segmentation per patient using Elbow method."""
    if pelt_stats_df.empty or pelt_unfiltered_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    curves = []
    min_bkps_per_patient = {}  # Track minimum breakpoints discovered per patient
    for pid, grp in pelt_stats_df.groupby("patient"):
        grp_sorted = grp.sort_values(by=["num_breakpoints", "cost_of_segmentation", "generating_penalty"])
        best_per_k = grp_sorted.groupby("num_breakpoints", as_index=False).first()
        curves.append((pid, best_per_k))
        # Track minimum breakpoints discovered on cost curve
        min_bkps_per_patient[pid] = int(grp["num_breakpoints"].min())

    elbow_rows = []
    history_rows = []
    wanted_pairs = []
    for pid, best_df in curves:
        x = best_df["num_breakpoints"].values
        y = best_df["cost_of_segmentation"].values.astype(float)
        p = best_df["generating_penalty"].values
        elbow_x, elbow_y, method = _elbow_for_patient(x, y)

        if elbow_x is None:
            elbow_x = 0
            if (best_df["num_breakpoints"] == 0).any():
                elbow_y = float(best_df.loc[best_df["num_breakpoints"] == 0, "cost_of_segmentation"].iloc[0])
            else:
                elbow_y = float(np.min(y)) if len(y) > 0 else np.nan
            method = "None -> default K=0"

        try:
            chosen_row = best_df.loc[best_df["num_breakpoints"] == elbow_x].iloc[0]
            chosen_solution_id = int(chosen_row["solution_id"])
            chosen_penalty = float(chosen_row["generating_penalty"])
        except Exception:
            chosen_solution_id = None
            chosen_penalty = np.nan

        elbow_rows.append({
            "patient": pid,
            "elbow_num_breakpoints": int(elbow_x),
            "elbow_cost": float(elbow_y) if elbow_y is not None else np.nan,
            "elbow_method": method,
            "elbow_solution_id": chosen_solution_id,
            "elbow_penalty": chosen_penalty,
            "min_bkps_discovered": min_bkps_per_patient.get(pid, np.nan),
        })

        for i in range(len(x)):
            history_rows.append({
                "patient": pid,
                "n_bkps": int(x[i]),
                "min_cost": float(y[i]) if pd.notna(y[i]) else np.nan,
                "min_penalty": float(p[i]) if i < len(p) and pd.notna(p[i]) else np.nan,
                "optimal": 1 if int(x[i]) == int(elbow_x) else 0
            })

        if chosen_solution_id is not None:
            wanted_pairs.append((pid, chosen_solution_id))

    elbow_summary_df = pd.DataFrame(elbow_rows)
    elbow_history_df = pd.DataFrame(history_rows)

    if not wanted_pairs:
        return pd.DataFrame(), elbow_summary_df, elbow_history_df

    key_df = pd.DataFrame(wanted_pairs, columns=["patient", "solution_id"]).drop_duplicates()
    elbow_optimized_df = pelt_unfiltered_df.merge(key_df, on=["patient", "solution_id"], how="inner")
    return elbow_optimized_df, elbow_summary_df, elbow_history_df

# ==============================================================================
# OLS Fitting
# ==============================================================================
def process_segment_ols(segment_df: pd.DataFrame) -> pd.DataFrame:
    """Fit OLS to a single segment."""
    warnings.filterwarnings("ignore")
    x = segment_df["pyruvate"].values
    y = segment_df["lactate"].values
    
    m, b, r, _ = compute_ols_stats(x, y)
    
    df_out = segment_df.copy()
    df_out["p6e_m"] = m
    df_out["p6e_b"] = b
    df_out["p6e_r"] = r
    return df_out

def _ols_worker(task):
    """Worker for parallel OLS processing."""
    _name, group = task
    return process_segment_ols(group)

def run_ols_processing(elbow_optimized_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate segment-level OLS fits in parallel."""
    if elbow_optimized_df.empty:
        return pd.DataFrame()
    
    segment_groups = elbow_optimized_df.groupby(["patient", "segment_number"])
    tasks = [(name, group.copy()) for name, group in segment_groups]

    if not tasks:
        return pd.DataFrame()

    with multiprocessing.Pool(processes=min(N_JOBS_PARALLEL, len(tasks))) as pool:
        results = pool.map(_ols_worker, tasks)

    if not results:
        return pd.DataFrame()
    return pd.concat(results, ignore_index=True)

def compute_patient_level_fits(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Compute global patient-level OLS fits."""
    results = []
    for pid, group in raw_df.groupby("patient"):
        x = group["pyruvate"].values
        y = group["lactate"].values
        
        m, b, r, _ = compute_ols_stats(x, y)
        
        results.append({
            "patient": pid,
            "patient_m": m,
            "patient_b": b,
            "patient_r": r
        })
    return pd.DataFrame(results)

# ==============================================================================
# Main
# ==============================================================================
def main():
    print("\n\nStarting 2_linear_segmentation.py")
    main_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    if not INPUT_FILE.exists():
        print(f"\nError: Input file not found at {INPUT_FILE}")
        sys.exit(1)

    # Load
    raw = pd.read_csv(INPUT_FILE)
    raw.columns = [c.lower().strip() for c in raw.columns]

    required = ["patient", "lactate", "pyruvate", "time_since_injury"]
    if not all(c in raw.columns for c in required):
        missing = [c for c in required if c not in raw.columns]
        print(f"\nError: Missing columns: {', '.join(missing)}")
        sys.exit(1)

    raw["time_since_injury"] = _safe_numeric(raw["time_since_injury"])
    raw["lactate"] = _safe_numeric(raw["lactate"])
    raw["pyruvate"] = _safe_numeric(raw["pyruvate"])
    raw = raw.dropna(subset=["patient", "lactate", "pyruvate", "time_since_injury"])
    raw["patient"] = raw["patient"].astype(str)

    pelt_stats_df, pelt_unfiltered_df = run_pelt_processing(raw, MIN_SEG_LEN)

    elbow_optimized_df, elbow_summary_df, elbow_history_df = run_elbow_processing(pelt_stats_df, pelt_unfiltered_df)
    
    if elbow_optimized_df.empty:
        print("\nNo data after elbow selection; exiting.")
        execution_time = time.time() - main_start
        pd.DataFrame({"execution_time_seconds": [execution_time]}).to_csv(TIME_FILE, index=False)
        return

    final_df = run_ols_processing(elbow_optimized_df)
    
    patient_stats = compute_patient_level_fits(raw)
    final_df = final_df.merge(patient_stats, on="patient", how="left")
    
    # Merge elbow metadata
    if not elbow_summary_df.empty:
        final_df = final_df.merge(
            elbow_summary_df[["patient", "elbow_num_breakpoints", "elbow_cost", "elbow_method", "elbow_penalty", "min_bkps_discovered"]],
            on="patient",
            how="left"
        )
        final_df["p6e_total_segs"] = pd.to_numeric(final_df["elbow_num_breakpoints"], errors="coerce") + 1
        final_df.drop(columns=["elbow_num_breakpoints"], inplace=True)
    
    # Rename columns
    rename_map = {
        "solution_id": "p6e_cost_index",
        "elbow_penalty": "p6e_penalty",
        "elbow_cost": "p6e_cost",
        "segment_number": "p6e_seg_index",
        "elbow_method": "p6e_elbow_method",
        "min_bkps_discovered": "p6e_min_bkps_discovered",
    }
    for old_c, new_c in rename_map.items():
        if old_c in final_df.columns:
            final_df.rename(columns={old_c: new_c}, inplace=True)

    # Build elbow history string
    elbow_hist_map = {}
    if not elbow_history_df.empty:
        try:
            for pid, grp in elbow_history_df.groupby("patient"):
                triples = grp.sort_values("n_bkps")[["n_bkps", "min_cost", "min_penalty"]].values
                header = "axes=(n_bkps,cost,penalty);"
                series = ";".join([f"({int(nb)},{float(mc)},{float(mp)})" for nb, mc, mp in triples])
                elbow_hist_map[pid] = header + series
        except Exception:
            pass
    final_df["p6e_cost_history"] = final_df["patient"].map(elbow_hist_map) if elbow_hist_map else np.nan
    
    # Column order
    ordered_cols_spec = [
        "patient",
        "time_since_injury",
        "lpr",
        "lactate",
        "pyruvate",
        "glucose",
        "gose",
        "p6e_seg_index",
        "p6e_total_segs",
        "p6e_min_bkps_discovered",
        "p6e_cost_index",
        "p6e_cost",
        "p6e_penalty",
        "p6e_cost_history",
        "p6e_elbow_method",
        "p6e_m",
        "p6e_b",
        "p6e_r",
        "patient_m",
        "patient_b",
        "patient_r",
    ]
    
    valid_cols = [c for c in ordered_cols_spec if c in final_df.columns]
    final_df = final_df[valid_cols]

    # Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    execution_time = time.time() - main_start
    pd.DataFrame({"execution_time_seconds": [execution_time]}).to_csv(TIME_FILE, index=False)

    # Check for datapoint mismatches
    input_counts = raw.groupby("patient").size().rename("input_datapoints")
    output_counts = final_df.groupby("patient").size().rename("output_datapoints")
    merged = input_counts.to_frame().join(output_counts, how="left").fillna(0).astype(int).reset_index()
    mismatches = merged[merged["input_datapoints"] != merged["output_datapoints"]]
    if not mismatches.empty:
        mismatches.to_csv(WARNING_FILE, index=False)

    print("\n2_linear_segmentation.py complete.\n")

if __name__ == "__main__":
    main()
