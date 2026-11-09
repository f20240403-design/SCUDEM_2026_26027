# Overnight RMSE study: final values

Run on the night of 4–5 November 2026: 3.6 hours on 10 MATLAB workers, from 00:17 to 03:52. Every value here is final and supersedes the first pass in [`../rmse_first_pass/`](../rmse_first_pass/README.md). The report section "Final validation with per-child records" presents these results.

These outputs also feed the rest of the project:
- the report's growth tables and arm table;
- the closed-form figures in [`../figures/`](../figures/README.md), which read `results/full/stage0.mat`;
- `../scudem_problem_a_variables.xlsx`.

## Folder layout

```
final_study/
├── data/      downloaded raw datasets, one folder each (read as data only; nothing executed)
│   ├── ghidini/   casefile_clinical_notes.xlsx   (Zenodo 12805856)
│   ├── snyder/    individuals.csv                     (NIST AnthroKids, Snyder 1977 individual records)
│   ├── nhanes/    bmx_1999.xpt, demo_1999.xpt, p_bmx_2017_2020.xpt, p_demo_2017_2020.xpt   (CDC NHANES)
│   └── imprint/   imprint_subject01.zip                       (Zenodo 23023603, subject 1)
├── scripts/   data preparation (Python) and the study (MATLAB)
└── results/
    ├── episodes_with_children.csv, snyder_individuals_tidy.csv, nhanes_children_tidy.csv   (prepared inputs)
    ├── full/    the 4-hour run: every final number, table, figure and checkpoint
    ├── smoke/   2-minute end-to-end test with tiny sample sizes (not for reporting)
    ├── bench/   timing benchmark of stage S0 (not for reporting)
    └── *_console.txt   raw MATLAB console output of the runs
```

## Re-running

Run from the `final_study/` folder. The data-preparation steps take seconds:

```
python -I scripts/build_episodes.py "data/ghidini/casefile_clinical_notes.xlsx" results/episodes_with_children.csv
```

```
python -I scripts/tidy_growth.py data/snyder/individuals.csv data/nhanes results
```

The study also reads `../ghidini_demographics.csv` from the project folder, for the ages at amputation used to fill in the missing ones. Then run it. `smoke` takes about 2 minutes and `full` about 3.6 hours on 10 workers. The study needs MATLAB R2024b with the Parallel Computing Toolbox:

```
matlab -batch "cd scripts; overnight_rmse('full')"
```

Every stage is checkpointed to `results/<mode>/stageN.mat`. A re-run loads the finished stages and only redoes what is missing, so redoing the outputs alone takes about 1 minute. To redo a stage from scratch, delete its `stageN.mat` and every later one.

| Script | Role |
|---|---|
| `scripts/build_episodes.py` | Raw casefile workbook → socket episodes with child IDs, dates, years since first visit, age at amputation |
| `scripts/tidy_growth.py` | Snyder and NHANES individual records → tidy tables in cm and kg |
| `scripts/inspect_xlsx.py` | Helper: print the sheets and first rows of a workbook |
| `scripts/overnight_rmse.m` | The study: stages S0–S8, checkpoints, tables, figures, `summary.txt` |
| `scripts/ov_fit_growth.m` | Logistic and two-phase growth fits (weighted least squares on yearly means, K capped at 1.10 × max mean) |
| `scripts/ov_curve.m` | Growth-curve lookup with a fast numerical inverse |
| `scripts/ov_simulate.m` | Monte Carlo time to intervention and cause for every episode (child-level spread, competing criteria) |
| `scripts/ov_obs.m` | Observed Kaplan–Meier and Aalen–Johansen curves |
| `scripts/ov_metrics.m` | RMSE, R², cause curves, C-index and IPCW integrated Brier score |

## 1. Datasets from the literature survey (`literature_datasets_survey.docx`)

| # | Dataset | Status | Used? | Why |
|---|---|---|---|---|
| 1 | Ghidini casefiles, Zenodo 12805856 | Open; raw workbook downloaded | **Yes** (core) | 814 dated visits with child IDs → per-episode fitting ages and child-level resampling |
| 2 | Ghidini growth cohort, Zenodo 14187212 | Restricted (access on request) | No | Would give ρ and residual girth; needs Iseul's access request |
| 3 | OpenLimbTT (GitHub) | Open | No | Adult residual-limb shape model; no growth or timing information |
| 4 | IMPRINT, Zenodo 23023603 | Open; subject 1 downloaded (28 MB) | Cross-check | Soft-tissue and liner moduli (0.25 and 0.6 MPa) bound the liner stiffness k_s; lifting stays at most about 1 mm |
| 5 | Residual-limb scans, Dryad | Open (sign-in for bulk download) | No | Adult limb geometry; no constant of the model depends on it |
| 6 | Smartphone photogrammetry (Sensors 2026) | Supplementary tables | No | Measurement method, not model input |
| 7 | Kinematic + EMG dataset (Sci Data 2025) | Open | No | Muscle activity; the model is quasi-static and uses masses only |
| 8 | MeganePro MDS1 | Open | No | Grasp EMG; no constant of the model |
| 9 | Ninapro DB1–DB3 | Open (academic) | No | Grasp EMG; no constant of the model |
| 10 | WAY-EEG-GAL | Open | Already used | Object masses 0.165 / 0.33 / 0.66 kg (m_o) |
| 11 | Otmani 2023 code and data (GitLab) | Sign-in required | No | Its finding (adult tables scaled by length) is already used for ℓ_p |
| 12 | NHANES | Open; 1999–2000 and 2017–2020 body measures downloaded | **Yes** (variant) | Modern calf and thigh circumference and upper-leg length (ages 8+): independent growth cross-check |
| — | Snyder 1977 individual records (NIST AnthroKids) | Open; downloaded | **Yes** (core) | 3,895 children: all growth curves refitted from individuals, including thigh and whole-leg length, with child-level bootstrap |

Each download is in its own folder under `data/`. The files are only read as data; nothing downloaded was executed.

## 2. Data preparation

- `scripts/build_episodes.py` rebuilds the socket episodes from the raw workbook, using the same rules as the data pack. It reproduces the original 273 episodes and 155 events exactly, now with child IDs, dates and years since each child's first visit. Output: `results/episodes_with_children.csv`.
- `scripts/tidy_growth.py` converts the individual records to tidy tables:
  - **Snyder:** mm → cm, 0.1 kg → kg, age from thousandths of a year, zeros dropped. The age-10 means match the published table.
  - **NHANES:** both survey cycles, ages 2–19.
  - **Outputs:** `results/snyder_individuals_tidy.csv`, `results/nhanes_children_tidy.csv`.

## 3. The study (`scripts/overnight_rmse.m`)

**Model.** The competing-risks leg model has five constants:
- Γ, the tightness tolerance;
- δ₀/ρ, the clearance;
- ΔL_max, the length tolerance;
- σ, the child-level spread;
- the delay from amputation to first visit.

Each socket's fitting age is its child's age at amputation plus the delay plus the years since the child's first visit. Sex mixes are level-specific. The fit target is the stacked RMSE of the Kaplan–Meier survival curve and the tight and short cumulative-incidence curves. Monte Carlo uses 400 draws per episode during fitting and 4,000 for the final evaluation.

| Stage | What | Size |
|---|---|---|
| S0 | Growth fits from individual records, logistic and two-phase; child-level growth bootstrap | 14 curves × 2,000 resamples |
| S1 | Global parameter sweep | 967,680 parameter sets (63 min) |
| S2 | Multi-start Nelder–Mead until RMSE_new/RMSE_old ≥ 0.9999 | Converged in 2 rounds |
| S3 | Child-level bootstrap: refit on resampled children; noise floor; held-out metrics | 5,250 refits (80 min) |
| S4 | Leave-one-child-out cross-validation | 32 folds × 3 starts |
| S5 | Model variants, each multi-start with bootstrap | 10 variants × 160–400 refits |
| S6 | Growth-curve uncertainty propagated into the constants | 1,500 refits |
| S7 | Arm prediction: median child with CI, predictive spread, Sobol indices | 3,000 draws; 200,000 samples per age; 131,072 Sobol base samples |
| S8 | Final evaluation (4,000 draws per episode), tables, figures | — |

## 4. Final results

**Constants**: estimate, then the 95% interval over 5,250 refits on resampled children.

| Constant | Estimate | 95% CI | Note |
|---|---|---|---|
| Γ (tightness tolerance) | 0.255 | 0.147–0.260 | Only bounded below: tightness almost never binds |
| δ₀/ρ (leg clearance, mean child) | 1.60 cm | 1.16–3.04 | 1.31 cm for the median child |
| ΔL_max (leg length tolerance, mean child) | 2.26 cm | 1.50–4.51 | 1.85 cm for the median child, inside the clinical 1–2 cm |
| σ (child-to-child spread) | 0.63 | 0.43–0.96 | |
| delay (amputation to first visit) | ≈ 0 yr | 0–0.0002 | Children are fitted soon after amputation |

**Performance**

| Set | RMSE survival | R² | RMSE tight / short | C-index | IBS (null) |
|---|---|---|---|---|---|
| TT calibration | 0.0150 | 0.997 | 0.014 / 0.008 | 0.58 | 0.187 (0.177) |
| TT leave-one-child-out | 0.0152 | 0.997 | 0.014 / 0.009 | 0.57 | 0.190 (0.177) |
| KD+TF held out | 0.0228 | 0.992 | 0.143 / 0.131 | 0.64 | 0.160 (0.166) |
| Child-level noise floor (median) | 0.033 | — | 0.030 / 0.024 | — | — |

**What the results mean**
- **Population level:** excellent. The fit is below the noise floor, and leaving each child out gives the same RMSE, so it is not overfitted.
- **Individual children:** prediction is weak. C-index is about 0.58, and the Brier score is no better than giving every socket the population curve.
- **Leg mechanism:** first problems are 65% distal clearance, 34% length and 1% tightness. The casefiles' "tight" problems are therefore mostly the limb reaching the end of the socket.
- **Variants:**
  - Per-episode ages are what give the model its ranking ability: C-index 0.58–0.64, against 0.51–0.52 with the stand-in ages.
  - The NHANES growth curves, a per-socket spread and an unbounded Γ all change nothing material.
  - Thigh and whole-leg curves for the knee and transfemoral levels predict them far too early (RMSE 0.17).
- **Arm prediction:**
  - The median child fitted at 4–10 years needs intervention after 9.8–12.3 months (boys) or 10.7–15.8 months (girls). At age 10 the 95% CI is 7.6–15.8 months (boys) and 9.7–20.7 (girls).
  - **Sobol indices** (share of the variance in T\*): the arm ΔL_max explains about 75%, δ₀ about 34% and ρ about 11%. Γ, the masses, k_s and growth uncertainty each explain less than 2%.

## 5. Files

- `results/full/summary.txt`: every number in one place
- `results/full/final_constants.csv`, `final_metrics.csv`, `variants.csv`, `noise_floor_cluster_bootstrap.csv`, `growth_fits_individuals.csv`, `arm_prediction.csv`, `arm_sobol.csv`, `refinement_history.csv`
- `results/full/fig_*.png`: the final fit, bootstrap distributions, variants, arm prediction, Sobol indices, convergence and leave-one-child-out figures
- `results/full/stage0.mat` … `stage7.mat`, `final.mat`: checkpoints. Re-running `overnight_rmse('full')` (see Re-running) loads them and only redoes the outputs, in about 1 minute. `stage0.mat` holds the growth fits (`P`), their 2,000 bootstrap draws (`PB`) and the NHANES fits (`PN`).
- `results/full/progress.log`: timestamped log of the run
- `results/smoke/`: the 2-minute end-to-end test run with tiny sample sizes; not for reporting
- `results/full/fig_final_fit.png`, `fig_bootstrap_constants.png`, `fig_arm_prediction.png`, `fig_arm_sobol.png`, `fig_variants.png` and `fig_loco.png` are copied to `../figures/` as `final_fit.png`, `final_bootstrap_constants.png`, `final_arm_prediction.png`, `final_arm_sobol.png`, `final_variants.png` and `final_leave_one_child_out.png` (the report includes the first four)
