"""Fit logistic growth curves L(a) = K / (1 + exp(-r (a - a50))) to Snyder 1977 age-sex means.
Weighted by n (sigma = SD/sqrt(n)); adult size K bounded to [max mean, 1.10 x max mean].
Input: snyder1977_growth_by_age_sex.csv   Output: growth_fits_logistic.csv"""
import csv, numpy as np
from scipy.optimize import curve_fit
rows = list(csv.DictReader(open('snyder1977_growth_by_age_sex.csv')))
def logi(a, K, r, a50): return K / (1 + np.exp(-r * (a - a50)))
out = []
for meas in ['forearm_length_radiale_stylion', 'forearm_circumference',
             'shank_length_tibiale_minus_sphyrion', 'calf_circumference']:
    for sex in ['male', 'female']:
        d = [r for r in rows if r['measure'] == meas and r['sex'] == sex]
        a = np.array([float(r['age_mid_yr']) for r in d]); y = np.array([float(r['mean']) for r in d])
        n = np.array([int(r['n']) for r in d]); sd = np.array([float(r['sd']) for r in d])
        ymax = y.max()
        p, cov = curve_fit(logi, a, y, p0=[ymax * 1.02, 0.3, 5], sigma=sd / np.sqrt(n), absolute_sigma=True,
                           bounds=([ymax, 0.01, -20], [1.10 * ymax, 2, 20]), maxfev=20000)
        res = y - logi(a, *p); w = n / n.sum(); rmse = np.sqrt(np.sum(w * res ** 2)); se = np.sqrt(np.diag(cov))
        out.append(dict(measure=meas, sex=sex, K_cm=round(p[0], 2),
                        K_at_bound='yes' if abs(p[0] - 1.10 * ymax) < 1e-3 * ymax else 'no',
                        r_per_yr=round(p[1], 3), r_se=round(se[1], 3), a50_yr=round(p[2], 2),
                        wRMSE_cm=round(rmse, 2), RMSE_pct_of_K=round(100 * rmse / p[0], 1),
                        max_abs_resid_cm=round(np.abs(res).max(), 2)))
with open('growth_fits_logistic.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
for o in out: print(o)
