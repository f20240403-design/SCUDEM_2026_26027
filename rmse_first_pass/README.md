# RMSE against the Ghidini casefiles: first pass (MATLAB)

> **Status: first pass, superseded.** This folder is the first leg validation, run on 4 November 2026 with MATLAB R2024b. It used the data-pack episode table, with ages at amputation standing in for fitting ages. The **final** constants and metrics come from the overnight per-child study in [`../final_study/`](../final_study/README.md), which re-processed the raw casefiles with child identifiers and real per-episode fitting ages. Where the numbers differ, the overnight values are the ones in the report.
>
> The report keeps this pass as its "Validation against clinical data" section, and still uses three of its figures: `fig1_km_vs_model.png`, `fig6_competing_risks.png` and the ρ sweep. Their copies in `../figures/` are named `first_pass_km_vs_model.png`, `first_pass_competing_risks.png` and `first_pass_rho_sweep.png`.

| | First pass (this folder) | Final (`../final_study/`) |
|---|---|---|
| Fitting ages | Ages at amputation (n = 45), jittered | Per episode: age at amputation + delay + years since first visit (196 of 273 episodes) |
| Constants | Γ = 0.121, λ = 1.95 cm, σ = 0.77 | Γ = 0.255, δ₀/ρ = 1.60 cm, ΔL_max = 2.26 cm, σ = 0.63, delay ≈ 0, all with 95% CIs |
| TT RMSE / R² | 0.0165 / 0.9963 | 0.0150 / 0.997 (leaving each child out: 0.0152) |
| KD+TF held out | 0.0293 / 0.9870 | 0.0228 / 0.992 |
| Noise floor | 0.0324 (episode bootstrap) | 0.033 (child-level bootstrap) |

## Scripts and run order

Run each from this folder:

```
matlab -batch "run_rmse_calibration"
```

```
matlab -batch "sweep_missing_constants"
```

```
matlab -batch "address_drawbacks"
```

| Script | Run time | Log | Purpose |
|---|---|---|---|
| `run_rmse_calibration.m` | ~3 min | `run_log.txt` | Main calibration of (Γ, λ, σ) on TT, held-out test on KD+TF, noise floor, fitting-age sensitivity |
| `sweep_missing_constants.m` | ~5 min | `sweep_log.txt` | Sweep of ρ × δ₀ × leg ΔL_max; 7-constant flexible model pushed toward R² = 0.9999 |
| `address_drawbacks.m` | ~10 min | `drawbacks_log.txt` | D1 competing risks, D2 statistical benchmarks, D3 two-phase growth fits |

Helper functions:

| File | Purpose |
|---|---|
| `model_tstar.m` | T\* per simulated socket for the 3-constant timing model (Γ, λ, σ) |
| `model_tstar_cause.m` | Competing-risks version (Γ, δ₀/ρ, ΔL_max, σ) that also returns the criterion that fired |
| `model_tstar_flex.m` | 7-constant version (two spreads, a clinic-visit delay, a fitting-age shift) |
| `km_curve.m` | Kaplan–Meier estimator on a monthly grid |
| `aj_cif.m` | Aalen–Johansen cumulative incidence of "tight" and "short" problems |
| `literature_values.m` | First-pass leg growth fits (Snyder age-group table, `../fit_growth.py`) and literature ranges, with citations |

Inputs: `../ghidini_socket_episodes.csv`, `../ghidini_demographics.csv`, `../snyder1977_growth_by_age_sex.csv`.

## Method

1. **Observed curve.** The Kaplan–Meier estimate of the months from a new prosthesis to its first growth problem. It is evaluated monthly up to the last month with at least 10 sockets at risk: 25 months for TT and 18 for KD+TF.
2. **Model curve.** For 20,000 simulated sockets:
   - draw a fitting age (age at amputation ± 0.5 yr) and a sex (72.6% boys);
   - compute T\* from the closed-form threshold ages, using the Snyder shank-length and calf-circumference fits, with u_lift = 0;
   - the model curve is the share of sockets with T\* > t.
3. **Unknown constants.** Literature ranges come first. The rest are found by a grid sweep (17,160 sets), then repeated Nelder–Mead refinement until RMSE_new / RMSE_old ≥ 0.9999.
4. **Calibrate on TT (188 episodes); test on KD+TF (85 episodes)** with the TT constants unchanged.
5. **Noise floor.** The RMSE between bootstrap resamples of the TT curve and the curve itself: the best agreement any model can show with 188 episodes.

## Results

### Calibration (`run_rmse_calibration.m`; `results/rmse_summary.csv`, `results/calibrated_constants.csv`)

| Constant | Calibrated | Literature range |
|---|---|---|
| Γ = (1+s_C)(1+ε_max) − 1 | 0.121 | 0.005–0.26 |
| λ = min(δ₀/ρ, ΔL_max) | 1.95 cm | 0.5–2.0 cm |
| σ (child-to-child spread, log-SD) | 0.77 | none |

| Comparison | RMSE | R² |
|---|---|---|
| TT, calibration | 0.0165 | 0.9963 |
| KD+TF, held out | 0.0293 | 0.9870 |
| TT with σ = 0 | 0.0926 | 0.8839 |
| KD+TF refitted (reference only) | 0.0188 | 0.9947 |
| Noise floor, TT bootstrap (median) | 0.0324 | 0.9858 |

- The refinement converged at iteration 3 (`convergence_history.csv`, `fig3_convergence.png`).
- Only Γ is identifiable, not s_C and ε_max separately. With u_lift = 0, clearance and length are crossings of the same shank curve, so only λ matters for timing.
- Shifting the stand-in ages by 0–3 years moves Γ within 0.07–0.18 and λ within 1.1–1.95 cm at almost the same RMSE (`sensitivity_fitting_age.csv`, `fig2_rmse_surface.png`). The constants are weakly identified individually.

### Missing constants and the R² = 0.9999 target (`sweep_missing_constants.m`)

- **ρ is not identifiable.** The best R² stays between 0.9962 and 0.9964 for every ρ (`missing_rho_sweep.csv`, `fig4_rho_sweep.png`).
- **Only the ratio δ₀/ρ is constrained.** It needed to be at least about 1.7 cm with stand-in ages; the final value is 1.60 cm.
- **ℓ_p and the arm ΔL_max** do not enter the leg model, so no leg fit can determine them.
- **Freeing every constant** (7 constants; `flex_convergence.csv`, `flex_vs_base.csv`, `fig5_flex_vs_base.png`) raises R² only from 0.9963 to 0.9967. It does so with implausible values: Γ at the literature ceiling, and fitting ages moved 4.5 years later.
- **R² = 0.9999 would need RMSE ≤ 0.0027,** about a twelfth of the noise floor, so it is not attainable.

### Drawbacks (`address_drawbacks.m`)

**D1 Competing risks** (`competing_risks_constants.csv`, `fig6_competing_risks.png`). The model is fitted to the Aalen–Johansen incidence of "tight" (tightness or clearance) and "short" (length) problems, with tight+short events split 50/50.

| | Γ | δ₀/ρ | ΔL_max | σ |
|---|---|---|---|---|
| Fitted value | 0.124 | 4.6 cm | 4.6 cm | 0.98 |

- **TT fit:** survival RMSE 0.020, tight 0.016, short 0.011.
- **Attribution of first problems:** tightness 35%, clearance 32%, length 33%.
- **Constants as means:** the median-child lengths are 2.9 cm. That implausible value came from the stand-in ages; the final study gives 1.85 cm.
- **Coding of "both" events:** assigning them wholly to one side keeps Γ within 0.10–0.14.

**D2 Benchmarks** (`benchmarks.csv`):

| Model | Constants | RMSE TT | RMSE KD+TF |
|---|---|---|---|
| Exponential | 1 | 0.060 | 0.079 |
| Weibull | 2 | 0.035 | 0.048 |
| Growth model (timing fit) | 3 | 0.016 | 0.029 |

The growth model stays inside the data's 95% band in 25 of 26 months.

**D3 Two-phase growth** (`bilogistic_fits.csv`, `bilogistic_params.mat`): fitted on the age-group table. Superseded by the individual-record fits in `../final_study/results/full/stage0.mat`, which are in the report's two-phase table.

## Results files

| File | Contents |
|---|---|
| `rmse_summary.csv`, `calibrated_constants.csv` | Main RMSE table; constants with literature check |
| `convergence_history.csv`, `sweep_grid.csv` | Refinement history; RMSE of every grid point |
| `sensitivity_fitting_age.csv` | Stand-in age shifts |
| `survival_curves.csv` | Kaplan–Meier and model curves, numbers at risk |
| `missing_rho_sweep.csv`, `flex_convergence.csv`, `flex_vs_base.csv` | Missing-constant sweep and flexible model |
| `competing_risks_constants.csv`, `benchmarks.csv` | D1 and D2 |
| `bilogistic_fits.csv`, `bilogistic_params.mat` | D3 (first pass) |
| `fig1_km_vs_model.png` … `fig6_competing_risks.png` | KM vs model, RMSE surface, convergence, ρ sweep, flexible vs base, competing risks |
| `calibration_workspace.mat` | MATLAB workspace |

## Literature sources for the ranges

- Ghidini et al. 2025, *Adv. Rehabil. Sci. Pract.* 14, doi:10.1177/27536351251384354 (sock layers, distal space, 1 cm ankle block)
- Ghidini et al. 2025 poster, MSF Science Portal (volume up to 60%, shortening up to 2 cm)
- Volume changes of about 1% are detectable: PMC7368584
- A 5-ply sock adds about 4% volume: PMC12614321
- A 5% volume gain is clinically unacceptable: US Patent 11,395,748
- Pistoning of about 4 mm per 30 N (k_s): Gholizadeh et al. 2011, doi:10.1177/0309364611423130
- Child terminal device of 130 g: doi:10.3389/frobt.2022.1000159
