#!/usr/bin/env python3
"""Comprehensive QC diagnostics for rat hippocampal proteomics data."""

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

# Log2 transform
log2_df = np.log2(df)

# Median normalization
medians = log2_df.median(axis=0)
global_median = log2_df.median().median()
norm_df = log2_df - medians + global_median

# Drop rows with NaN for complete-case analysis
complete = norm_df.dropna()
print(f'Complete cases: {complete.shape[0]} / {norm_df.shape[0]} proteins')

# Build sample metadata
sample_info = []
for c in complete.columns:
    name = c.replace('.raw', '')
    if 'Pio' in name:
        group = 'HFD_Pio'
        animal = name.split('Pio')[0]  # e.g. AHFD5
        rep = name.split('Pio')[1]     # e.g. R1
    elif 'IP' in name:
        group = 'HFD_IP'
        animal = name.split('IP')[0]
        rep = name.split('IP')[1]
    elif 'Cr' in name:
        group = 'HFD_Cr'
        animal = name.split('Cr')[0]
        rep = name.split('Cr')[1]
    elif 'AHFD' in name:
        group = 'HFD'
        parts = name.replace('AHFD', '')
        animal = 'AHFD' + parts[0]
        rep = parts[1:]
    elif 'AC' in name:
        group = 'CTRL'
        parts = name.replace('AC', '')
        animal = 'AC' + parts[0]
        rep = parts[1:]
    else:
        group = 'Unknown'
        animal = name
        rep = ''
    sample_info.append({
        'sample': c, 'name': name, 'group': group,
        'animal': animal, 'replicate': rep
    })

meta = pd.DataFrame(sample_info)
print()
print('=== SAMPLE STRUCTURE ===')
print(meta[['name', 'group', 'animal', 'replicate']].to_string(index=False))
print()
print('=== ANIMALS PER GROUP ===')
for g in meta['group'].unique():
    animals = meta[meta['group'] == g]['animal'].unique()
    n_animals = len(animals)
    n_samples = len(meta[meta['group'] == g])
    print(f'  {g}: {n_animals} animals, {n_samples} samples -> animals: {list(animals)}')

# === DIAGNOSTIC 1: Correlation heatmap ===
corr = complete.corr(method='pearson')
labels = []
for _, r in meta.iterrows():
    labels.append(r['name'] + ' (' + r['group'] + ')')

fig, ax = plt.subplots(figsize=(16, 14))
sns.heatmap(corr, annot=False, xticklabels=labels, yticklabels=labels,
            cmap='RdYlBu_r', vmin=0.9, vmax=1.0, ax=ax)
plt.title('Pearson Correlation Heatmap (Median-Normalized Log2)')
plt.xticks(rotation=90, fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()
plt.savefig('correlation_heatmap.png', dpi=150)
print('Saved correlation_heatmap.png')

# === DIAGNOSTIC 2: Hierarchical clustering dendrogram ===
dist_matrix = pdist(complete.T.values, metric='correlation')
linkage_matrix = linkage(dist_matrix, method='average')

fig, ax = plt.subplots(figsize=(14, 6))
dn = dendrogram(linkage_matrix,
                labels=[r['name'] for _, r in meta.iterrows()],
                leaf_rotation=90, leaf_font_size=9, ax=ax)
plt.title('Hierarchical Clustering Dendrogram (1 - Pearson correlation)')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('dendrogram.png', dpi=150)
print('Saved dendrogram.png')

# === DIAGNOSTIC 3: PCA colored by ANIMAL vs GROUP ===
X = complete.T.values
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=5)
pcs = pca.fit_transform(X_scaled)
print(f'\nExplained variance ratios (PC1-5): {pca.explained_variance_ratio_}')

res = pd.DataFrame({
    'PC1': pcs[:, 0], 'PC2': pcs[:, 1], 'PC3': pcs[:, 2],
    'PC4': pcs[:, 3], 'PC5': pcs[:, 4],
    'Group': meta['group'].values,
    'Animal': meta['animal'].values,
    'Sample': meta['name'].values
})

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# PCA colored by GROUP
ev1 = pca.explained_variance_ratio_[0] * 100
ev2 = pca.explained_variance_ratio_[1] * 100
sns.scatterplot(data=res, x='PC1', y='PC2', hue='Group', s=120,
                ax=axes[0], palette='Set1')
for i, row in res.iterrows():
    axes[0].text(row['PC1'] + 0.1, row['PC2'] + 0.1, row['Sample'], fontsize=7)
axes[0].set_title(f'PCA colored by Treatment Group\nPC1={ev1:.1f}%, PC2={ev2:.1f}%')
axes[0].set_xlabel(f'PC1 ({ev1:.1f}%)')
axes[0].set_ylabel(f'PC2 ({ev2:.1f}%)')

# PCA colored by ANIMAL
sns.scatterplot(data=res, x='PC1', y='PC2', hue='Animal', s=120,
                ax=axes[1], palette='tab10')
for i, row in res.iterrows():
    axes[1].text(row['PC1'] + 0.1, row['PC2'] + 0.1, row['Sample'], fontsize=7)
axes[1].set_title(f'PCA colored by Individual Animal\nPC1={ev1:.1f}%, PC2={ev2:.1f}%')
axes[1].set_xlabel(f'PC1 ({ev1:.1f}%)')
axes[1].set_ylabel(f'PC2 ({ev2:.1f}%)')

plt.tight_layout()
plt.savefig('pca_group_vs_animal.png', dpi=150)
print('Saved pca_group_vs_animal.png')

# === DIAGNOSTIC 4: PC3 vs PC4 ===
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
ev3 = pca.explained_variance_ratio_[2] * 100
ev4 = pca.explained_variance_ratio_[3] * 100

sns.scatterplot(data=res, x='PC3', y='PC4', hue='Group', s=120,
                ax=axes[0], palette='Set1')
for i, row in res.iterrows():
    axes[0].text(row['PC3'] + 0.1, row['PC4'] + 0.1, row['Sample'], fontsize=7)
axes[0].set_title(f'PCA colored by Treatment Group\nPC3={ev3:.1f}%, PC4={ev4:.1f}%')

sns.scatterplot(data=res, x='PC3', y='PC4', hue='Animal', s=120,
                ax=axes[1], palette='tab10')
for i, row in res.iterrows():
    axes[1].text(row['PC3'] + 0.1, row['PC4'] + 0.1, row['Sample'], fontsize=7)
axes[1].set_title(f'PCA colored by Individual Animal\nPC3={ev3:.1f}%, PC4={ev4:.1f}%')

plt.tight_layout()
plt.savefig('pca_pc3_pc4.png', dpi=150)
print('Saved pca_pc3_pc4.png')

# === DIAGNOSTIC 5: Technical replicate reproducibility ===
print()
print('=== TECHNICAL REPLICATE CORRELATIONS (within-animal) ===')
for animal in sorted(meta['animal'].unique()):
    cols = meta[meta['animal'] == animal]['sample'].tolist()
    grp = meta[meta['animal'] == animal]['group'].values[0]
    if len(cols) > 1:
        sub_corr = complete[cols].corr()
        mask = np.ones(sub_corr.shape, dtype=bool)
        np.fill_diagonal(mask, False)
        mean_corr = sub_corr.values[mask].mean()
        min_corr = sub_corr.values[mask].min()
        print(f'  {animal} ({grp}): mean r = {mean_corr:.4f}, min r = {min_corr:.4f}')

# === DIAGNOSTIC 6: Coefficient of variation per group ===
print()
print('=== INTRA-GROUP CV (median across proteins) ===')
for g in sorted(meta['group'].unique()):
    cols = meta[meta['group'] == g]['sample'].tolist()
    group_data = complete[cols]
    cv = group_data.std(axis=1) / group_data.mean(axis=1)
    print(f'  {g}: median CV = {cv.median():.4f}, mean CV = {cv.mean():.4f}')

# === DIAGNOSTIC 7: Intra-animal vs inter-animal variance ===
print()
print('=== VARIANCE DECOMPOSITION ===')
# For each protein, compute variance within animals vs between animals
within_vars = []
between_vars = []
for prot_id in complete.index:
    row = complete.loc[prot_id]
    animal_means = []
    for animal in meta['animal'].unique():
        cols = meta[meta['animal'] == animal]['sample'].tolist()
        vals = row[cols].values
        animal_means.append(np.mean(vals))
        within_vars.append(np.var(vals, ddof=0))
    between_vars.append(np.var(animal_means, ddof=0))

mean_within = np.mean(within_vars)
mean_between = np.mean(between_vars)
print(f'  Mean within-animal variance:  {mean_within:.4f}')
print(f'  Mean between-animal variance: {mean_between:.4f}')
print(f'  Ratio (between/within):       {mean_between/mean_within:.2f}')
if mean_between / mean_within > 5:
    print('  -> Between-animal variance dominates. This is expected for biological replicates.')
else:
    print('  -> Within-animal (technical) variance is relatively high.')

print()
print('Done with all diagnostics!')
