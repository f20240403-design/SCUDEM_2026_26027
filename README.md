# A growth–statics model for paediatric prosthetic socket replacement

SCUDEM XI 2026, Problem A. When does a child's below-elbow (transradial) prosthetic socket need adjusting or replacing as the residual limb grows?

The model has four layers:
- a logistic growth ODE for limb length and circumference;
- a socket-fit map (tightness and distal clearance);
- planar lift statics, whose curl load pistons the socket;
- a closed-form three-criterion threshold that gives the time to intervention T\*.

Paediatric casefiles exist only for the leg, so the model is **validated on the leg and used to predict for the arm**.

**Report:** [`scudem_problem_a_report.tex`](scudem_problem_a_report.tex); compiled version: [`scudem_problem_a_report.pdf`](scudem_problem_a_report.pdf).

## Headline results

- **Leg validation** (Cambodian casefiles, 273 socket episodes, 59 children):
  - transtibial calibration: RMSE 0.0150 (R² 0.997);
  - leaving each child out in turn: 0.0152;
  - held-out knee and transfemoral levels: 0.0228 (R² 0.992);
  - all are below the data's own noise floor of 0.033.
- **Final leg constants** (95% CIs from 5,250 refits on resampled children):

| Constant | Estimate | 95% CI |
|---|---|---|
| Clearance δ₀/ρ | 1.60 cm | 1.16–3.04 |
| Length tolerance ΔL_max | 2.26 cm | 1.50–4.51 |
| Child-to-child spread σ | 0.63 | 0.43–0.96 |
| Tightness tolerance Γ | 0.255 | only bounded below: ≳ 0.15 |

- **Arm prediction:** the median child fitted at ages 4–10 needs intervention after 9.8–12.3 months (boys) or 10.7–15.8 months (girls). Distal clearance is the criterion that fires first.
- **Arm uncertainty:** the arm length tolerance and the built-in clearance explain most of it (Sobol total indices about 0.75 and 0.34).

## Repository layout

| Path | Contents |
|---|---|
| `scudem_problem_a_report.tex` / `.pdf` | Final report (pdfLaTeX; compiled on Overleaf) |
| `scudem_problem_a_overleaf.zip` | The tex plus its 11 figures, ready to upload to Overleaf (rebuilt by `build_overleaf_zip.py`) |
| `scudem_problem_a_variables.xlsx` | Every variable and constant with its value, type and source |
| [`final_study/`](final_study/README.md) | **Final values.** Per-child RMSE study: raw data, preparation scripts, MATLAB study (3.6 h on 10 workers), all results |
| [`rmse_first_pass/`](rmse_first_pass/README.md) | First-pass calibration with stand-in ages (superseded, kept for the report's first-pass section) |
| [`figures/`](figures/README.md) | Every figure in the report, the closed-form figure scripts and the heat maps |
| `fit_growth.py` | First-pass logistic fits to the Snyder age-group table |
| `ghidini_socket_episodes.csv`, `ghidini_demographics.csv`, `snyder1977_growth_by_age_sex.csv` | Data pack used by the first pass (described in `data_pack_readme.pdf`) |
| `work_log.md` | Chronological log of the work and decisions |
| `literature_datasets_survey.docx` | Survey of related papers and datasets |
| `scudem_problem_a_working_derivation.pdf` | Earlier 8-page working derivation |
| `ghidini_2025_prosthetist_perspectives.pdf` | Ghidini et al. 2025, the main clinical reference |
| `rename_map.csv` | Old → new file names from the tidy-up for this repository |

## Reproducing

Requirements: MATLAB R2024b (the Parallel Computing Toolbox for the final study), and Python 3 with `pandas` and `openpyxl`.

```
python build_overleaf_zip.py
```

```
matlab -batch "cd figures; make_figures; arm_sensitivity; make_heatmaps"
```

- **The final study** is in [`final_study/README.md`](final_study/README.md): data preparation, then `overnight_rmse('full')`, which resumes from checkpoints.
- **The first pass** is in [`rmse_first_pass/README.md`](rmse_first_pass/README.md).

## Data sources

| Data | Source | Used for |
|---|---|---|
| Paediatric lower-limb casefiles | Ghidini & Edgar, Zenodo [10.5281/zenodo.12805856](https://doi.org/10.5281/zenodo.12805856) | Leg validation (episodes, events, ages) |
| Children's anthropometry | Snyder et al. 1977, individual records via NIST AnthroKids | All growth curves |
| NHANES 1999–2000, 2017–2020 body measures | CDC/NCHS | Growth cross-check |
| IMPRINT finite-element dataset | Zenodo [10.5281/zenodo.23023603](https://doi.org/10.5281/zenodo.23023603) | Liner stiffness cross-check |

Third-party datasets and the included paper remain under their original licences and terms of use; see each source for details.
