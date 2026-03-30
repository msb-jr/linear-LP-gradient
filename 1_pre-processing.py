# ==============================================================================
# Script: 1_pre-processing.py
# Manuscript relevance: 2.2, SM2.ii
# ==============================================================================
# PURPOSE:
#   Convert the original data into a cleaned dataset for all downstream analyses.
#
# INPUT:
#   - 0_input/original_data.csv: See README.md for expected columns and units.
#                                Note: INPUT PYRUVATE NEEDS TO BE IN UNITS OF μM.
#
# OUTPUT:
#   - 1_output/1_pre-processing/__execution_time.csv: Runtime log
#   - 1_output/1_pre-processing/1_data.csv: Pre-processed data
#
# PREPROCESSING STEPS:
# ------------------------------------------------------------------------------
# (1) Column Selection (FILTER_COLUMNS = True)
#     - Only keep columns required by downstream scripts
#     Columns:    patient, time_since_injury, lactate, pyruvate, glucose, gose
#
# (2) Early Sample Removal (FILTER_EARLY = True)
#     - Remove unreliable baseline readings from each patient's data
#     Approach:   Drop samples where time_since_injury ≤ (patient's minimum time + 2 hours)
#
# (3) Cerebral Microdialysis Concentration Filtering (FILTER_ISCUS_RANGE = True)
#     - Null out implausible values (outside ranges from Iscus Flex 2017 manual)
#     Ranges:     lactate 0.1–12 mM, pyruvate 10–1500 μM, glucose 0.1–25 mM
#
# (4) Missing Data Removal (FILTER_NULLS = True)
#     - Ensure complete cases for key variables
#     Variables:   patient, time_since_injury, lactate, pyruvate
#
# (5) Low-Sample Patient Exclusion (FILTER_LOW_SAMPLE_PATIENTS = True)
#     - Ensure sufficient data for meaningful segmentation by Script 2
#     Minimum:   24 samples per patient
#
# (6) Pyruvate Unit Conversion & LPR Computation (No toggle)
#     - pyruvate from μM to mM; lpr = lactate/pyruvate
# ==============================================================================

import time
import pandas as pd
from pathlib import Path


# ==============================================================================
# Input / Output
# ==============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR / "0_input" / "original_data.csv"
OUTPUT_DIR = SCRIPT_DIR / "1_output" / "1_pre-processing"
OUTPUT_FILE = OUTPUT_DIR / "1_data.csv"

# ==============================================================================
# Configuration / Toggles
# ==============================================================================
# Pre-processing toggles (DEFAULT = True)
FILTER_COLUMNS = True              # Keep only essential columns
FILTER_EARLY = True                # Remove early samples (first 2 hours)
FILTER_ISCUS_RANGE = True          # Set out-of-range values to NA (Iscus Flex 2017 manual)
FILTER_NULLS = True                # Remove rows with nulls key columns
FILTER_LOW_SAMPLE_PATIENTS = True  # Remove patients with < 24 samples


# Detection ranges from Iscus Flex 2017 manual
# (values outside these ranges are set to NA)
LACTATE_RANGE = (0.1, 12)      # mM
PYRUVATE_RANGE = (10, 1500)    # μM (before downstream conversion to mM)
GLUCOSE_RANGE = (0.1, 25)      # mM


def main() -> None:
    """
    Pre-process the original dataset and write a cleaned CSV for downstream analyses.
    """
    print("\n\nStarting 1_pre-processing.py")
    start_time = time.time()

    # ==========================================================================
    # Pre-processing
    # ==========================================================================
    data = pd.read_csv(INPUT_FILE)

    # Keep only relevant columns
    if FILTER_COLUMNS:
        columns_to_keep = [
            'patient',
            'time_since_injury',
            'lactate',
            'pyruvate',
            'glucose',
            'gose',
        ]
        data = data[columns_to_keep]

    # Filter out early samples for each patient
    if FILTER_EARLY:
        # Get the minimum time_since_injury for each patient
        min_time_per_patient = data.groupby('patient')['time_since_injury'].transform('min')
        # Keep only samples that are more than 2 hours after first sample
        data = data[data['time_since_injury'] > (min_time_per_patient + 2)]

    # Apply detection range filters (set out-of-range values to NA, does not remove rows)
    # Based on Iscus Flex 2017 manual specifications
    if FILTER_ISCUS_RANGE:
        # Lactate: 0.1 - 12 mM
        data.loc[
            (data['lactate'] < LACTATE_RANGE[0]) | (data['lactate'] > LACTATE_RANGE[1]),
            'lactate'
        ] = None
        # Pyruvate: 10 - 1500 μM (still in μM at this point)
        data.loc[
            (data['pyruvate'] < PYRUVATE_RANGE[0]) | (data['pyruvate'] > PYRUVATE_RANGE[1]),
            'pyruvate'
        ] = None
        # Glucose: 0.1 - 25 mM
        data.loc[
            (data['glucose'] < GLUCOSE_RANGE[0]) | (data['glucose'] > GLUCOSE_RANGE[1]),
            'glucose'
        ] = None

    # Keep only rows with non-null values for key columns
    if FILTER_NULLS:
        columns_to_check = ['patient', 'time_since_injury', 'lactate', 'pyruvate']
        data = data.dropna(subset=columns_to_check)

    # Filter out patients with fewer than 24 samples
    if FILTER_LOW_SAMPLE_PATIENTS:
        sample_counts = data.groupby('patient').size()
        patients_to_keep = sample_counts[sample_counts >= 24].index
        data = data[data['patient'].isin(patients_to_keep)]

    # Convert pyruvate from µM to mM (ASSUMES INPUT DATA HAS PYRUVATE IN µM)
    data['pyruvate'] = data['pyruvate'] * 0.001

    # Calculate LPR from lactate and pyruvate (both now in mM)
    # Insert column after time_since_injury
    safe_pyruvate = data['pyruvate'].replace(0, pd.NA)
    lpr_values = data['lactate'].divide(safe_pyruvate)
    if 'lpr' in data.columns:
        data = data.drop(columns=['lpr'])
    time_index = data.columns.get_loc('time_since_injury')
    data.insert(time_index + 1, 'lpr', lpr_values)

    # Save the pre-processed data
    # Ensure the output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUTPUT_FILE, index=False)

    # Save execution time
    execution_time = time.time() - start_time
    execution_time_df = pd.DataFrame({
        'execution_time_seconds': [execution_time]
    })
    execution_time_file = OUTPUT_DIR / '__execution_time.csv'
    execution_time_df.to_csv(execution_time_file, index=False)

    print("\n1_pre-processing.py complete.\n")


if __name__ == "__main__":
    main()