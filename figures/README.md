# Figures for the report

This folder holds every figure that `../scudem_problem_a_report.tex` includes, plus the two MATLAB scripts that draw the closed-form (deterministic) figures. The tex reads its images from here (`\graphicspath{{figures/}}`), and `../scudem_problem_a_overleaf.zip` packs this folder together with the tex.

**Status (5 November 2026):** all figures use the final values. Both scripts use the final growth fits from the Snyder individual records and the final leg calibration (Γ = 0.255; band 0.147–0.260). The source of both is `../final_study/results/full/`.

## Run

```
matlab -batch "make_figures"
```

```
matlab -batch "arm_sensitivity"
```

Run them from this folder. Each takes under a minute. The console output is saved in `make_figures_log.txt` and `arm_sensitivity_log.txt`, regenerated on 5 November 2026. The numbers quoted in the report text come from these logs.

## Figures in the report

| File | Report figure | Made by | What it shows |
|---|---|---|---|
| `growth_fits.png` | `fig:fits` | `make_figures.m` | Snyder age-group means ± SE, with the logistic curves fitted to the individual records (forearm length and circumference; shank length and calf circumference) |
| `arm_headline.png` | `fig:headline` | `make_figures.m` | Arm T\* against fitting age for the mean-tolerance child. Central line (clearance fires first at every age) and corner band with every input at its low or high end |
| `arm_tornado.png` | `fig:tornado` | `arm_sensitivity.m` | One-at-a-time sensitivity of T\* at a₀ = 10 (boys 19.3, girls 25.1 months) |
| `arm_growth_spurt.png` | `fig:spurt` | `arm_sensitivity.m` | Single logistic against two-phase (growth spurt) curves |
| `first_pass_rho_sweep.png` | `fig:rho` | `make_figures.m` (redrawn from `../rmse_first_pass/results/missing_rho_sweep.csv`) | First pass: best R² as ρ varies, for four clearances |
| `first_pass_km_vs_model.png` | `fig:km` | `../rmse_first_pass/run_rmse_calibration.m` (copy of `../rmse_first_pass/results/fig1_km_vs_model.png`) | First pass: Kaplan–Meier curves against the timing model, TT calibration and KD+TF held out |
| `first_pass_competing_risks.png` | `fig:cr` | `../rmse_first_pass/address_drawbacks.m` (copy of `results/fig6_competing_risks.png`) | First pass: cumulative incidence of tight and short problems |
| `final_fit.png` | `fig:finalfit` | `../final_study/scripts/overnight_rmse.m` (copy of `fig_final_fit.png`) | Final fit: survival and cause curves with child-bootstrap bands |
| `final_bootstrap_constants.png` | `fig:bootconst` | overnight (copy) | Final constants over 5,250 child-level and 1,500 growth-curve bootstrap refits |
| `final_arm_prediction.png` | `fig:armfinal` | overnight (copy) | Final arm prediction: median child with 95% CI, 10–90% band over children, mean-tolerance line |
| `final_arm_sobol.png` | `fig:sobol` | overnight (copy) | Total Sobol indices of arm T\* at fitting ages 6, 10 and 14 |

The figures prefixed `final_` are copies of `../final_study/results/full/fig_*.png` (`fig_loco.png` became `final_leave_one_child_out.png`). If the final study is re-run, copy them again. `final_variants.png` and `final_leave_one_child_out.png` are kept here for slides but are not included in the report.

## Scripts

**`make_figures.m`**
- Growth-fit figure and closed-form arm prediction (`arm_headline.png`), with all edge cases: a target above the adult size never fires, and δ₀ ≤ u_lift fires at fitting.
- Prints the numbers the report quotes:
  - u_lift range and the k_s needed for 0.5 cm of pistoning;
  - tightness cut-off ages (16.2 / 11.0 yr at Γ = 0.255);
  - the mean-tolerance T\* table (matches the "Mean tol." column of the report's arm table);
  - the effect of the lift at a₀ = 10 (20.2 → 19.3 → 18.7 months);
  - the non-monotone length branch below a₅₀;
  - socket moments, axial load and shoulder index by age (lift statics and shoulder extension).
- Redraws the first-pass ρ sweep.

**`arm_sensitivity.m`** (reads `../final_study/results/full/stage0.mat` for the growth fits and their bootstrap SDs)
- **S1:** tornado at a₀ = 10.
- **S2:** child-to-child spread at the central values (σ = 0.63, an independent draw per criterion). Superseded in the report by the overnight arm table.
- **S3:** growth-spurt comparison using the final two-phase fits.
- **S4:** leg-informed scenario, with ρ = 1.5 / 1.60 ≈ 0.94. Clearance still fires first, and T\* rises by 1–2 months.
- Saves `arm_sensitivity.mat`.

## Heat maps (`make_heatmaps.m`, not yet in the report)

```
matlab -batch "make_heatmaps"
```

Takes about 1 minute. Reads `../final_study/results/full/stage0.mat`, `stage1.mat` and `stage3.mat`. The log is in `make_heatmaps_log.txt`. Colour follows magnitude in a single hue: orange for misfit (RMSE), blue for months and shares.

| File | What it shows | How to read it |
|---|---|---|
| `heat_leg_fit_landscape.png` | Leg calibration: stacked RMSE over δ₀/ρ × ΔL_max, Γ × δ₀/ρ and σ × δ₀/ρ, from the 967,680-set sweep. Each cell is the best fit over the constants not shown. | Lighter is a better fit. Dot: final estimate. Solid line: RMSE 0.015. Dashed line: child-level noise floor (0.032). δ₀/ρ, ΔL_max and σ form a compact best region. For Γ, everything from about 0.03 up to the 0.26 bound lies inside the noise floor, which is why Γ is only bounded below. |
| `heat_arm_age_clearance.png` | Arm T\* (mean-tolerance child) over fitting age × built-in clearance δ₀, boys and girls | White contours in months. Above δ₀ ≈ 2.06 cm the length tolerance (2 cm) fires first, so more clearance stops helping. Black bar: central design for ages 4–10. |
| `heat_arm_design.png` | Arm socket design at a₀ = 10: δ₀ × ΔL_max (ρ = 1) and δ₀ × ρ (ΔL_max = 2 cm) | Dashed: clearance and length tie (Proposition `prop:race`). T\* rises with δ₀ only until it reaches that line, then only ΔL_max matters. Dotted: the leg-implied ratio δ₀/ρ = 1.60 cm passes next to the central design. |
| `heat_arm_risk.png` | Share of children who have needed intervention by month t, by fitting age. Central design; children (σ), calibrated constants (5,250 bootstrap fits) and growth curves all vary. | Read across at 1 year: boys fitted at 4 / 10 / 16: 65% / 48% / 17%; girls: 58% / 32% / 7%. White contours at 10, 25, 50, 75 and 90%. |

## Central arm inputs used by both scripts

| Input | Central | Band (low / high end) |
|---|---|---|
| Γ (tightness tolerance) | 0.255 | 0.147 / 0.260 |
| δ₀ (distal clearance) | 1.5 cm | 1 / 3 cm |
| ρ (residual / intact length) | 1 | 1 / 0.5 |
| ΔL_max (arm length tolerance) | 2 cm | 1 / 3 cm |
| m_p, m_o | 0.24, 0.33 kg | heaviest 0.24 + 0.66; lightest 0.13 + 0.165 |
| k_s (liner stiffness) | 75 N/cm | tornado: 25–150 |
| u_lift (curl pistoning) | 0.065 cm | 0.033–0.102 cm |

The arm values of δ₀, ρ and ΔL_max are unmeasured. They are set by assumption and carried as ranges. The Sobol analysis shows they explain almost all of the uncertainty in the arm prediction.
