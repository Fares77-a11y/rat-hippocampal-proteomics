#!/usr/bin/env python3
"""Re-examine sample structure with the correct interpretation:
30 columns = 30 individual rats (biological replicates).
The digit after the group prefix (1 vs 2) = Wave/batch.
R1, R2, R3 = individual rats within that wave.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

# Load data
df = pd.read_csv('diann_matrix.tsv', sep='\t')
df = df.set_index('Protein.Ids')
df = df.drop(columns=['Description'], errors='ignore')

# Build CORRECTED metadata
# Wave 1: AC1, AHFD1, AHFD3, AHFD5, AHFD7
# Wave 2: AC2, AHFD2, AHFD4, AHFD6, AHFD8
wave1_prefixes = {'AC1', 'AHFD1', 'AHFD3', 'AHFD5', 'AHFD7'}
wave2_prefixes = {'AC2', 'AHFD2', 'AHFD4', 'AHFD6', 'AHFD8'}

sample_info = []
for c in df.columns:
    name = c.replace('.raw', '')
    # Determine group
    if 'Pio' in name:
        group = 'HFD_Pio'
    elif 'IP' in name:
        group = 'HFD_IP'
    elif 'Cr' in name:
        group = 'HFD_Cr'
    elif 'AHFD' in name:
        group = 'HFD'
    elif 'AC' in name:
        group = 'CTRL'
    else:
        group = 'Unknown'

    # Determine wave from the prefix
    # Extract the prefix up to (but not including) the treatment suffix and R
    if name.startswith('AC'):
        prefix = 'AC' + name[2]  # AC1 or AC2
    elif name.startswith('AHFD'):
        prefix = 'AHFD' + name[4]  # AHFD1..AHFD8
    else:
        prefix = 'UNK'

    if prefix in wave1_prefixes:
        wave = 'Wave1'
    elif prefix in wave2_prefixes:
        wave = 'Wave2'
    else:
        wave = 'Unknown'

    sample_info.append({
        'sample': c, 'name': name, 'group': group,
        'wave': wave, 'prefix': prefix
    })

meta = pd.DataFrame(sample_info)

print('=== CORRECTED SAMPLE STRUCTURE ===')
print(meta[['name', 'group', 'wave', 'prefix']].to_string(index=False))
print()
print('=== GROUP x WAVE BALANCE ===')
print(meta.groupby(['group', 'wave']).size().unstack(fill_value=0))
print()
n_groups = meta['group'].nunique()
print(f'Total samples: {len(meta)}')
print(f'Unique groups: {n_groups}')
sizes = meta.groupby('group').size().to_dict()
print(f'Samples per group: {sizes}')
print()

# === Log2 transform ===
log2_df = np.log2(df)

# === Median normalization ===
medians = log2_df.median(axis=0)
global_median = log2_df.median().median()
norm_df = log2_df - medians + global_median

complete = norm_df.dropna()
print(f'Complete cases for PCA: {complete.shape[0]} / {norm_df.shape[0]}')

# === PCA colored by GROUP and by WAVE ===
X = complete.T.values
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=5)
pcs = pca.fit_transform(X_scaled)
ev = pca.explained_variance_ratio_
print(f'Explained variance (PC1-5): {[round(x*100,1) for x in ev]}%')

res = pd.DataFrame({
    'PC1': pcs[:, 0], 'PC2': pcs[:, 1], 'PC3': pcs[:, 2],
    'PC4': pcs[:, 3], 'PC5': pcs[:, 4],
    'Group': meta['group'].values,
    'Wave': meta['wave'].values,
    'Sample': meta['name'].values,
    'Prefix': meta['prefix'].values,
})

# --- Figure 1: PCA by Group vs by Wave ---
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
ev1, ev2 = ev[0]*100, ev[1]*100

sns.scatterplot(data=res, x='PC1', y='PC2', hue='Group', style='Wave',
                s=120, ax=axes[0], palette='Set1')
for i, row in res.iterrows():
    axes[0].text(row['PC1']+0.1, row['PC2']+0.1, row['Sample'], fontsize=6)
axes[0].set_title(f'PCA: Color=Group, Shape=Wave\nPC1={ev1:.1f}%, PC2={ev2:.1f}%')
axes[0].set_xlabel(f'PC1 ({ev1:.1f}%)')
axes[0].set_ylabel(f'PC2 ({ev2:.1f}%)')

sns.scatterplot(data=res, x='PC1', y='PC2', hue='Wave', style='Group',
                s=120, ax=axes[1], palette='Set2')
for i, row in res.iterrows():
    axes[1].text(row['PC1']+0.1, row['PC2']+0.1, row['Sample'], fontsize=6)
axes[1].set_title(f'PCA: Color=Wave, Shape=Group\nPC1={ev1:.1f}%, PC2={ev2:.1f}%')
axes[1].set_xlabel(f'PC1 ({ev1:.1f}%)')
axes[1].set_ylabel(f'PC2 ({ev2:.1f}%)')

plt.tight_layout()
plt.savefig('pca_group_wave.png', dpi=150)
print('Saved pca_group_wave.png')

# --- Figure 2: Hierarchical clustering with wave annotation ---
dist_matrix = pdist(complete.T.values, metric='correlation')
linkage_matrix = linkage(dist_matrix, method='average')

leaf_labels = []
for _, r in meta.iterrows():
    leaf_labels.append(f"{r['name']} [{r['wave']}]")

fig, ax = plt.subplots(figsize=(14, 6))
dn = dendrogram(linkage_matrix, labels=leaf_labels,
                leaf_rotation=90, leaf_font_size=8, ax=ax)
plt.title('Hierarchical Clustering (1 - Pearson r)\nLabels show [Wave]')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('dendrogram_wave.png', dpi=150)
print('Saved dendrogram_wave.png')

# --- Figure 3: Correlation heatmap sorted by group then wave ---
sort_order = meta.sort_values(['group', 'wave', 'name']).index
sorted_cols = [meta.loc[i, 'sample'] for i in sort_order]
sorted_labels = [meta.loc[i, 'name'] + ' [' + meta.loc[i, 'wave'] + ']' for i in sort_order]

corr_sorted = complete[sorted_cols].corr()
fig, ax = plt.subplots(figsize=(16, 14))
sns.heatmap(corr_sorted, annot=False, xticklabels=sorted_labels,
            yticklabels=sorted_labels,
            cmap='RdYlBu_r', vmin=0.90, vmax=1.0, ax=ax)
plt.title('Correlation Heatmap (sorted by Group, then Wave)')
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
plt.savefig('correlation_heatmap_sorted.png', dpi=150)
print('Saved correlation_heatmap_sorted.png')

# === Variance partition: Group vs Wave vs Residual ===
# Simple ANOVA-style decomposition
from scipy import stats as sp_stats

group_var_fracs = []
wave_var_fracs = []
resid_var_fracs = []

groups = meta['group'].values
waves = meta['wave'].values

for prot_id in complete.index:
    vals = complete.loc[prot_id].values
    total_var = np.var(vals, ddof=1)
    if total_var == 0:
        continue

    # Group means
    group_means = {}
    for g in np.unique(groups):
        mask = groups == g
        group_means[g] = np.mean(vals[mask])
    grand_mean = np.mean(vals)

    # SS_group
    ss_group = sum(np.sum(groups == g) * (group_means[g] - grand_mean)**2
                   for g in np.unique(groups))

    # Wave means
    wave_means = {}
    for w in np.unique(waves):
        mask = waves == w
        wave_means[w] = np.mean(vals[mask])

    # SS_wave
    ss_wave = sum(np.sum(waves == w) * (wave_means[w] - grand_mean)**2
                  for w in np.unique(waves))

    ss_total = np.sum((vals - grand_mean)**2)
    ss_resid = ss_total - ss_group - ss_wave
    if ss_resid < 0:
        ss_resid = 0  # rounding

    group_var_fracs.append(ss_group / ss_total)
    wave_var_fracs.append(ss_wave / ss_total)
    resid_var_fracs.append(ss_resid / ss_total)

print()
print('=== VARIANCE PARTITION (Group / Wave / Residual) ===')
print(f'  Group:    median = {np.median(group_var_fracs)*100:.1f}%, '
      f'mean = {np.mean(group_var_fracs)*100:.1f}%')
print(f'  Wave:     median = {np.median(wave_var_fracs)*100:.1f}%, '
      f'mean = {np.mean(wave_var_fracs)*100:.1f}%')
print(f'  Residual: median = {np.median(resid_var_fracs)*100:.1f}%, '
      f'mean = {np.mean(resid_var_fracs)*100:.1f}%')

# === Within-group correlations (now correctly: 6 bio reps per group) ===
print()
print('=== WITHIN-GROUP PAIRWISE CORRELATIONS ===')
for g in sorted(meta['group'].unique()):
    cols = meta[meta['group'] == g]['sample'].tolist()
    sub_corr = complete[cols].corr()
    mask = np.ones(sub_corr.shape, dtype=bool)
    np.fill_diagonal(mask, False)
    vals = sub_corr.values[mask]
    print(f'  {g}: mean r = {vals.mean():.4f}, min r = {vals.min():.4f}, '
          f'max r = {vals.max():.4f}')

# === Within-wave vs between-wave correlations ===
print()
print('=== WITHIN-WAVE vs BETWEEN-WAVE CORRELATIONS (within same group) ===')
for g in sorted(meta['group'].unique()):
    g_meta = meta[meta['group'] == g]
    w1_cols = g_meta[g_meta['wave'] == 'Wave1']['sample'].tolist()
    w2_cols = g_meta[g_meta['wave'] == 'Wave2']['sample'].tolist()

    # Within-wave
    within_vals = []
    for wave_cols in [w1_cols, w2_cols]:
        if len(wave_cols) > 1:
            sc = complete[wave_cols].corr()
            m = np.ones(sc.shape, dtype=bool)
            np.fill_diagonal(m, False)
            within_vals.extend(sc.values[m].tolist())

    # Between-wave
    between_vals = []
    for c1 in w1_cols:
        for c2 in w2_cols:
            between_vals.append(complete[[c1, c2]].corr().iloc[0, 1])

    if within_vals and between_vals:
        print(f'  {g}: within-wave r = {np.mean(within_vals):.4f}, '
              f'between-wave r = {np.mean(between_vals):.4f}, '
              f'delta = {np.mean(within_vals) - np.mean(between_vals):.4f}')

print()
print('Done!')
