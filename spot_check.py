import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import t as t_dist

# Load imputed matrix and the DE results
imputed = pd.read_csv('corrected_DE_results/tables/imputed_matrix.csv', index_col=0)
de_res = pd.read_csv('corrected_DE_results/tables/de_HFD_vs_CTRL.csv')

# Pick the top protein
top_prot = de_res.iloc[0]['Protein']
print(f'Spot-checking protein: {top_prot}')
print(f'Reported log2FC: {de_res.iloc[0]["log2FC"]:.6f}')
print(f'Reported p-value: {de_res.iloc[0]["p_value"]:.2e}')
print(f'Reported adj_p: {de_res.iloc[0]["adj_p_value"]:.6f}')

# Manually rebuild the design
wave1_prefixes = {'AC1', 'AHFD1', 'AHFD3', 'AHFD5', 'AHFD7'}
meta_records = []
for c in imputed.columns:
    name = c.replace('.raw', '')
    if 'Pio' in name: group = 'HFD_Pio'
    elif 'IP' in name: group = 'HFD_IP'
    elif 'Cr' in name: group = 'HFD_Cr'
    elif 'AHFD' in name: group = 'HFD'
    elif 'AC' in name: group = 'CTRL'
    else: group = 'Unknown'
    if name.startswith('AC'): prefix = 'AC' + name[2]
    elif name.startswith('AHFD'): prefix = 'AHFD' + name[4]
    else: prefix = 'UNK'
    wave = 'Wave1' if prefix in wave1_prefixes else 'Wave2'
    meta_records.append({'sample': c, 'group': group, 'wave': wave})
meta = pd.DataFrame(meta_records).set_index('sample')

design = pd.DataFrame(index=meta.index)
for g in sorted(meta['group'].unique()):
    design[f'group_{g}'] = (meta['group'] == g).astype(int)
design['wave2'] = (meta['wave'] == 'Wave2').astype(int)

y = imputed.loc[top_prot].values.astype(float)
X = design.values.astype(float)
model = sm.OLS(y, X).fit()

col_list = design.columns.tolist()
print(f'Design columns: {col_list}')
contrast = {'group_HFD': 1, 'group_CTRL': -1}
c_vec = np.array([contrast.get(col, 0.0) for col in col_list])
print(f'Contrast vector: {c_vec}')

estimate = float(c_vec @ model.params)
se = float(np.sqrt(c_vec @ model.cov_params() @ c_vec))
t_stat = estimate / se
df_resid = model.df_resid
p_value = float(2 * t_dist.sf(np.abs(t_stat), df_resid))

print(f'Manual log2FC: {estimate:.6f}')
print(f'Manual SE: {se:.6f}')
print(f'Manual t_stat: {t_stat:.6f}')
print(f'Manual df_resid: {df_resid}')
print(f'Manual p_value: {p_value:.2e}')

# Also manually compute group means
ctrl_vals = y[meta['group'] == 'CTRL']
hfd_vals = y[meta['group'] == 'HFD']
print(f'\nCTRL mean: {ctrl_vals.mean():.4f}, HFD mean: {hfd_vals.mean():.4f}')
print(f'Naive diff (HFD-CTRL): {hfd_vals.mean() - ctrl_vals.mean():.6f}')
print(f'Model log2FC: {estimate:.6f}')
print('(Small difference is due to wave adjustment)')
