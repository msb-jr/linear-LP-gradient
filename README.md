# Linear LP Gradient Analysis

## Citation

If you use this code, please cite:

> Baker MS, Heihre JM, Kulinski T, et al. The microdialysis-derived lactate–pyruvate gradient indicates anaerobic activity after brain injury. *Manuscript under review*.

---

## Requirements

- **Python 3.9+** ([python.org](https://www.python.org/downloads/))
- **R 4.0+** ([r-project.org](https://cran.r-project.org/))

Verify installation: `python3 --version` and `R --version`

---

## Download

Navigate to where you want the project, then:

```bash
git clone https://github.com/msb-jr/linear-LP-gradient.git
cd linear-LP-gradient
```

Or click **Code → Download ZIP** on GitHub and extract.

---

## Setup

### Python Environment

```bash
python3 -m venv venv              # Create (first time only)
source venv/bin/activate          # Activate (Mac/Linux)
# venv\Scripts\activate.bat       # Activate (Windows)
pip install -r requirements.txt   # Install packages
```

### R Environment

```bash
Rscript -e "renv::restore()"
```

---

## Pipeline

### Input Data

Place your data file at:

```
0_input/original_data.csv
```

**Expected content:** CSV with the following columns<sup>**a**</sup>:

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `patient` | string or numeric | — | Unique patient identifier |
| `time_since_injury` | numeric | hours | Time elapsed since injury |
| `lactate` | numeric | mM | Cerebral microdialysis lactate |
| `pyruvate` | numeric | **μM** | Cerebral microdialysis pyruvate |
| `glucose`<sup>**b**</sup> | numeric | mM | Cerebral microdialysis glucose |
| `gose`<sup>**c**</sup> | integer (1-8) | — | Glasgow Outcome Scale-Extended |

**Note:** `lpr` (lactate/pyruvate ratio) is computed automatically by Script 1 from `lactate` and `pyruvate`. 
<br><sup>**a**</sup> Columns can contain missing values, they will be handled automatically as needed. 
<br><sup>**b**</sup> Used in Script 8
<br><sup>**c**</sup> Used in Script 9

Additionally, Script 10 requires a separate in vitro dataset at `0_input/in_vitro.csv` with columns: `type` (wt|rot), `time` (hours), `lactate` (mM), `pyruvate` (mM).

---

### Scripts

Run scripts 1 → 2 in order, then can run 3+. See each script header for full details.

| Script | Description | Manuscript relevance |
|--------|-------------|----------------------|
| `1_pre-processing.py` | Data pre-processing | 2.2, SM2.ii |
| `2_linear_segmentation.py` | PELT linear segmentation + annotation | 2.3.ii |
| `3_segmentation_visual.py` | Per-patient segmentation plots | Fig. 2, Fig. S2, Fig. S3 |
| `4_segmentation_stats.R` | Segmentation yield summary | 3.1, Table S1, Table S4 |
| `5_linear_LP_models.R` | Linear model summaries | 3.2, 3.3.i, Fig. 3 |
| `6_hyperbolic_LPR_error.R` | LPR error analysis | 3.3.ii, Fig. 4, Table S2, Table S5 |
| `7_never_always_LPR25.R` | LPR cohort comparisons | 3.4, Table S3 |
| `8_substrate_lpr_fidelity.R` | Substrate delivery and LPR fidelity | 3.5, Fig. 5 |
| `9_outcome.R` | Outcome associations (GOSE) | 3.6, Fig. 6 |
| `10_in_vitro.R` | In vitro ETC inhibition validation | 3.8, Fig. 7 |
| `extra_linear_to_hyperbolic.py` | Conceptual linear→hyperbolic relationship | 2.3.i, Fig. S1 |
| `extra_lpr_m_schematic.py` | CMD Schematic with LPR timecourse & linear regressions | Fig. 1 |

### Run

```bash
# Individual
python 1_pre-processing.py
python 2_linear_segmentation.py
# ... etc

# All at once
python 1_pre-processing.py && \
python 2_linear_segmentation.py && \
python 3_segmentation_visual.py && \
Rscript 4_segmentation_stats.R && \
Rscript 5_linear_LP_models.R && \
Rscript 6_hyperbolic_LPR_error.R && \
Rscript 7_never_always_LPR25.R && \
Rscript 8_substrate_lpr_fidelity.R && \
Rscript 9_outcome.R && \
Rscript 10_in_vitro.R && \
python extra_linear_to_hyperbolic.py && \
python extra_lpr_m_schematic.py
```

Outputs are saved to `1_output/<script_name>/`.

---

## Cleaning Up

**Remove project environments:**
- Python: delete the `venv/` folder
- R: delete the `renv/library/` folder