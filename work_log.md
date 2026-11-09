# SCUDEM XI 2026, Problem A: data and analysis work log

Work done on Wednesday, 4 November 2026. Covers the data pulled, how each table was built, the growth-curve fits, and what is still open.

## Summary

| Step | Status | Output file |
|---|---|---|
| Pull the Ghidini casefile data (validation) | Done | `ghidini_socket_episodes.csv`, `ghidini_demographics.csv` |
| Pull the Snyder 1977 child anthropometry (growth) | Done | `snyder1977_growth_by_age_sex.csv` |
| Fit the logistic growth ODE, with RMSE | Done | `growth_fits_logistic.csv`, `fit_growth.py` |
| Ghidini growth cohort (Zenodo 14187212) | Locked; needs an access request from Iseul | — |
| Otmani segment data (GitLab) | Locked behind a sign-in; not needed | — |
| NHANES | Not pulled; Snyder covers everything NHANES lacks | — |
| Calibration and validation against the casefiles | Not started; can start now | — |
| Mechanical parameters from literature | With the mechanical teammate | — |

---

## 1. Validation data: Ghidini casefiles

**Source:** Ghidini & Edgar, *Data related to 62 Clinical Casefiles of Children with Lower Limb Loss/Absence from Cambodia*, Zenodo, doi:10.5281/zenodo.12805856 (open). Raw file: `casefile_clinical_notes.xlsx`, with four sheets: Demographics, TT (transtibial), KD (knee disarticulation) and TF (transfemoral).

**What the raw data holds**
- About 820 dated visits across 62 children.
- Each visit records the reason for attending, the clinical evaluation (category and subcategory) and the intervention.
- The demographics give limb-loss level, cause and age at amputation. There are **no dates of birth and no fitting ages**.

**How the episodes were built** (`ghidini_socket_episodes.csv`)
1. An episode starts at every visit where a new prosthesis was provided.
2. **Event** = the first later visit coded "Growth" with "socket tight" and/or "short".
   - "Socket tight" corresponds to the model's tightness criterion.
   - "Short" corresponds to the model's length-discrepancy criterion.
3. **Censored** = the episode ended at the next new prosthesis for a non-growth reason, or at the child's last visit.
4. Zero-length episodes (two prostheses on the same day) were dropped.

**Columns:** `level`, `days`, `months`, `event` (1 = growth problem, 0 = censored), `cause` (tight, short, tight+short), `action` (adjustment or new_prosthesis).

**Results**

| | All levels | Transtibial only |
|---|---|---|
| Episodes | 273 | 188 |
| Growth events | 155 | 102 |
| Kaplan–Meier median time to first growth problem | **11.8 months** | 11.8 months |
| Share still free of growth problems at 12 months | 48% | 48% |

- Event causes: tight 73, short 51, both 31.
- Fixed by: new prosthesis 98, adjustment 57.

**Cross-check against the paper** ([PMC12572713](https://pmc.ncbi.nlm.nih.gov/articles/PMC12572713)): the paper reports adjustments every 8.3 ± 6.6 months and replacements every 16.4 ± 9.9 months. Our episodes measure something different (time from a new prosthesis to its first growth problem), so the numbers aren't expected to match.

**Limitations**
- No fitting ages, so predictions can only be compared as a distribution, not child by child.
- Distal clearance isn't recorded as its own problem.
- This is leg data. Use it for the leg validation run only.

---

## 2. Growth data: Snyder 1977

**Source:** Snyder et al. (1977), *Anthropometry of Infants, Children, and Youths to Age 18 for Product Safety Design*, UMTRI. Individual-level data from NIST AnthroKids: https://math.nist.gov/~SRessler/anthrokids/orig77/individuals.csv (3,900 children, 123 columns).

**Processing**
- Converted raw values from mm to cm, and from 0.1 kg to kg.
- Sex code 1 = male. Confirmed from stature at ages 16–18: code 1 averages 175.9 cm, code 2 averages 162.8 cm.
- Zeros (not measured) were excluded.
- Grouped by sex and whole year of age (2–18), giving n, mean and SD per group.

**Measurements in `snyder1977_growth_by_age_sex.csv`**

| `measure` | Snyder variable | Model use |
|---|---|---|
| `forearm_length_radiale_stylion` | RADIALE-STYLION LENGTH | Intact forearm ℓ_fa(a); residual length L = ρ·ℓ_fa |
| `elbow_hand_length` | ELBOW-HAND LENGTH | ℓ_h at fitting (replaces a literature value) |
| `forearm_circumference` | FOREARM CIRCUMFERENCE | C(a), arm |
| `upper_arm_length_shoulder_elbow` | SHOULDER-ELBOW LENGTH | ℓ_ua(a), shoulder extension |
| `body_mass` | WEIGHT | M_b(a), strength proxy |
| `shank_length_tibiale_minus_sphyrion` | TIBIALE HEIGHT − SPHYRION HEIGHT | Shank length, leg validation |
| `calf_circumference` | CALF CIRCUMFERENCE | Circumference, leg validation |

**Notes**
- Radius and shank lengths were measured on a subset of about 1,240 children, so those groups are small (n ≈ 7–53). Other measurements have n ≈ 30–157 per group.
- Age-18 groups are small.
- These are intact limbs of US children from the 1970s. Residual limbs need the ratio ρ, and the Cambodian population differs. State both as assumptions.

---

## 3. Growth ODE fits (no literature values needed)

**Model:** the closed-form solution of the logistic ODE dL/da = rL(1 − L/K),

$$L(a) = \frac{K}{1 + e^{-r(a - a_{50})}}$$

where a₅₀ is the age at half the adult size.

**Method** (`fit_growth.py`)
- Nonlinear least squares on the age-sex means.
- Each mean is weighted by its standard error, SD/√n.
- Adult size K is bounded between the largest observed mean and 10% above it.
- Weighted RMSE = √(Σ wᵢ residualᵢ²), with wᵢ = nᵢ / Σn.

**Why K is bounded:** the data stops at age 18, before some measurements level off. Unbounded, the male forearm length fit gave K ≈ 33 cm, which is implausibly large, and the male circumference fits gave impossible values. **State this cap as an assumption in the write-up.**

**Results** (`growth_fits_logistic.csv`)

| Measure | Sex | K (cm) | r (1/yr) | a₅₀ (yr) | Weighted RMSE (cm) | RMSE % of K | Largest residual (cm) | K at cap? |
|---|---|---|---|---|---|---|---|---|
| Forearm length | M | 29.70 | 0.156 ± 0.007 | 4.19 | 0.39 | 1.3 | 0.81 | yes |
| Forearm length | F | 25.85 | 0.166 ± 0.010 | 2.67 | 0.42 | 1.6 | 1.31 | yes |
| Forearm circumference | M | 31.02 | 0.108 ± 0.005 | 3.56 | 0.88 | 2.8 | 2.36 | yes |
| Forearm circumference | F | 25.52 | 0.125 ± 0.007 | −0.01 | 0.37 | 1.4 | 0.65 | yes |
| Shank length | M | 43.18 | 0.181 ± 0.007 | 4.74 | 0.31 | 0.7 | 1.27 | no |
| Shank length | F | 38.28 | 0.215 ± 0.009 | 3.38 | 0.46 | 1.2 | 1.36 | no |
| Calf circumference | M | 41.14 | 0.126 ± 0.005 | 3.63 | 0.78 | 1.9 | 1.76 | yes |
| Calf circumference | F | 37.29 | 0.144 ± 0.006 | 2.15 | 0.58 | 1.6 | 1.08 | yes |

**Findings**
1. The logistic ODE fits within **0.3–0.9 cm** (at most 2.8% of adult size) for every curve.
2. The worst fit is **male forearm circumference**, which misses by up to 2.4 cm around puberty. This is the growth spurt that the model's Remark 1 already flags, so present it as evidence for that remark.
3. Six of the eight adult sizes sit at the cap. The rate r and the age a₅₀ partly trade off against K, so treat each fitted curve as a whole rather than reading meaning into the individual parameters.

---

## 4. What's next

### Calibration and validation against the casefiles (can start now)

This compares the predicted distribution of T\* with the casefile Kaplan–Meier curve.

1. **Calibrate** ε_max, δ₀ and ΔL_max with the leg growth fits, minimising the MSE between the predicted and observed survival curves on a monthly grid.
2. **Hold out data:** calibrate on the transtibial episodes (188), and report the RMSE on the transfemoral and knee-disarticulation episodes (85).
3. **Independent check:** if the calibrated values fall within the mechanical teammate's literature ranges, the model reproduces real clinical data with physically plausible tolerances.
4. **Fitting age a₀:** use the age-at-amputation distribution from the demographics (mean 10.1 years, SD 4.5) as a stand-in, and say so.
5. **Pistoning:** set u_lift = 0 for the first pass, then add it once k_s and m_p arrive.

### Waiting on others

| Item | Who | Notes |
|---|---|---|
| ε_max, δ₀, ΔL_max (arm and leg), s_C, k_s, m_p, with citations | Mechanical teammate | ℓ_h and c_h are no longer needed; they come from Snyder |
| Access to the Ghidini growth cohort (Zenodo 14187212) | Iseul | Sign in to Zenodo and click "Request access". The paper (doi:10.1097/PXR.0000000000000483) may give the ratio ρ |

---

## Files in this pack

| File | Contents |
|---|---|
| `README.md` | Data pack guide |
| `ghidini_socket_episodes.csv` | 273 socket episodes for Kaplan–Meier validation |
| `ghidini_demographics.csv` | 62 casefiles: level, cause, age at amputation |
| `snyder1977_growth_by_age_sex.csv` | Growth table by age and sex, 7 measurements |
| `growth_fits_logistic.csv` | Fitted logistic parameters and RMSE |
| `fit_growth.py` | Script that reproduces the growth fits |
| `work_log.md` | This file |

## Sources
- [Ghidini casefile dataset, Zenodo 12805856](https://zenodo.org/records/12805856)
- [Ghidini growth cohort, Zenodo 14187212 (restricted)](https://zenodo.org/records/14187212)
- [Ghidini et al. 2025, Advances in Rehabilitation Science and Practice (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12572713)
- [NIST AnthroKids, Snyder 1977 data](https://math.nist.gov/~SRessler/anthrokids/)
- [Otmani et al. 2023, International Biomechanics](https://www.tandfonline.com/doi/full/10.1080/23335432.2023.2268686)

---

## 5. Overnight final study (4–5 November 2026)

- **Raw data downloaded.** The raw Ghidini casefiles (Zenodo 12805856), the Snyder individual records (NIST AnthroKids), NHANES 1999–2000 and 2017–2020 body measures, and one IMPRINT subject (Zenodo 23023603). All 12 datasets in `literature_datasets_survey.docx` were reviewed; `final_study/README.md` says which were used and why.
- **Episodes rebuilt with child IDs.** Same rules as before: 273 episodes and 155 events, reproduced exactly. Per-episode fitting ages are now available for 196 episodes.
- **Overnight study.** `final_study/scripts/overnight_rmse.m` ran for 3.6 h on 10 workers:
  - a 967,680-set sweep;
  - 5,250 refits on resampled children;
  - leave-one-child-out cross-validation;
  - 10 model variants;
  - 1,500 refits on resampled growth curves;
  - the arm prediction with Sobol indices.
- **Final leg constants** (95% CI):
  - Γ = 0.255 (0.147–0.260), only bounded below;
  - δ₀/ρ = 1.60 cm (1.16–3.04);
  - ΔL_max = 2.26 cm (1.50–4.51), median child 1.85 cm;
  - σ = 0.63 (0.43–0.96);
  - delay ≈ 0.
- **Final performance:**
  - transtibial RMSE 0.0150 (R² 0.997);
  - leave-one-child-out 0.0152;
  - held-out knee and transfemoral 0.0228 (R² 0.992);
  - noise floor 0.033;
  - C-index 0.58 / 0.64.
- **Final arm prediction:** the median child fitted at 4–10 years needs intervention after 9.8–12.3 months (boys) or 10.7–15.8 (girls).
- **Updated:** `scudem_problem_a_report.tex` (new Section "Final validation with per-child records"), `scudem_problem_a_variables.xlsx`, `scudem_problem_a_overleaf.zip`.

### 5.1 Recheck after the overnight study (5 November 2026)
- The report and workbook now use the final growth fits from the Snyder individual records, not the first-pass age-group fits (differences at most 0.3 cm in K and 0.08 yr in a50).
- `figures/make_figures.m` and `figures/arm_sensitivity.m` now use those fits and the final Γ = 0.255. The headline, tornado, spurt and growth figures were redrawn.
- First-pass leftovers in the tex were corrected: tightness cut-off ages 16.2 / 11.0 yr, implied ε_max ≥ 0.08–0.10, leg-informed scenario ρ ≈ 0.94, spurt-effect numbers. The Overleaf zip was rebuilt.
- The workbook's Growth fits sheet has the final logistic values, child-bootstrap SDs, n, and a new two-phase table.
