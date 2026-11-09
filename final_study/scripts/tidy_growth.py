"""Tidy individual-level growth data for the overnight RMSE study.
Snyder 1977 (NIST AnthroKids individuals.csv): mm -> cm, 0.1 kg -> kg, zeros = not measured, sex 1 = male.
NHANES 1999-2000 and 2017-2020: ages 2-19, cm.
Usage: python -I tidy_growth.py <snyder_csv> <nhanes_dir> <out_dir>
"""
import sys, csv, os
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
sny, nh, out = sys.argv[1], sys.argv[2], sys.argv[3]

d = pd.read_csv(sny, encoding='latin-1', skiprows=[1])   # row 2 holds data types (and has extra fields)
num = lambda c: pd.to_numeric(d[c], errors='coerce').replace(0, float('nan'))
age = pd.to_numeric(d['AGE IN YEARS'], errors='coerce') / 1000.0       # stored in thousandths of a year
agem = pd.to_numeric(d['AGE IN MONTHS'], errors='coerce') / 12.0
print('age check (years from thousandths vs months/12), max abs diff:', float((age - agem).abs().max()))
sex = pd.to_numeric(d['SEX'], errors='coerce').map({1: 'male', 2: 'female'})
mm = lambda c: num(c) / 10.0
T = pd.DataFrame({
    'age_yr': age, 'sex': sex,
    'forearm_length': mm('RADIALE-STYLION LENGTH'),
    'forearm_circ': mm('FOREARM CIRCUMFERENCE'),
    'elbow_hand': mm('ELBOW-HAND LENGTH'),
    'upper_arm': mm('SHOULDER-ELBOW LENGTH'),
    'shank_length': mm('TIBIALE HEIGHT') - mm('SPHYRION HEIGHT'),
    'thigh_length': mm('TROCHANTERIC HEIGHT') - mm('TIBIALE HEIGHT'),
    'leg_length': mm('TROCHANTERIC HEIGHT') - mm('SPHYRION HEIGHT'),
    'calf_circ': mm('CALF CIRCUMFERENCE'),
    'thigh_circ': mm('UPPER THIGH CIRCUMFERENCE'),
    'body_mass': num('WEIGHT') / 10.0,
})
T = T[(T.age_yr >= 2) & (T.age_yr < 19) & T.sex.notna()]
T.to_csv(os.path.join(out, 'snyder_individuals_tidy.csv'), index=False)
print('Snyder individuals', len(T))
print(T.groupby('sex').count().T)
# sanity: means at ages 10-11 should match the published age-group table
print(T[(T.age_yr >= 10) & (T.age_yr < 11)].groupby('sex')[['forearm_length', 'shank_length', 'calf_circ', 'thigh_length', 'thigh_circ']].mean().round(1))

rows = []
for tag, bmx, demo, mcol in [('1999', 'bmx_1999.xpt', 'demo_1999.xpt', 'RIDAGEEX'), ('2017', 'p_bmx_2017_2020.xpt', 'p_demo_2017_2020.xpt', None)]:
    b = pd.read_sas(os.path.join(nh, bmx), format='xport')
    m = pd.read_sas(os.path.join(nh, demo), format='xport')
    x = b.merge(m, on='SEQN')
    a = x['RIDAGEYR'] + 0.5
    if mcol and mcol in x:
        a = x[mcol].where(x[mcol].notna(), x['RIDAGEYR'] * 12 + 6) / 12.0
    D = pd.DataFrame({'cycle': tag, 'age_yr': a, 'sex': x['RIAGENDR'].map({1: 'male', 2: 'female'})})
    for v, n in [('BMXCALF', 'calf_circ'), ('BMXTHICR', 'thigh_circ'), ('BMXLEG', 'upper_leg_length'),
                 ('BMXARMC', 'upper_arm_circ'), ('BMXARML', 'upper_arm_length')]:
        D[n] = x[v] if v in x else pd.NA
    rows.append(D[(D.age_yr >= 2) & (D.age_yr < 19.5)])
N = pd.concat(rows)
N.to_csv(os.path.join(out, 'nhanes_children_tidy.csv'), index=False)
print('NHANES children', len(N))
print(N.groupby('cycle').count().T)
print('youngest age with data:', {c: round(N.loc[N[c].notna(), 'age_yr'].min(), 1) for c in ['calf_circ', 'thigh_circ', 'upper_leg_length', 'upper_arm_circ']})
