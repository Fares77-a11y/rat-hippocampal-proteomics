#!/usr/bin/env python3
"""
Corrected Differential Expression Analysis for Rat Hippocampal Proteomics
=========================================================================
- 30 rats (biological replicates), 6 per group, balanced 3-per-wave
- Cyclic loess normalization
- KNN imputation (k=5)
- Linear model: ~ 0 + group + wave
- Contrasts: HFD vs CTRL, HFD_Pio vs HFD, HFD_IP vs HFD, HFD_Cr vs HFD
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import f as f_dist
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = Path('corrected_DE_results')
OUTPUT_DIR.mkdir(exist_ok=True)
(OUTPUT_DIR / 'figures').mkdir(exist_ok=True)
(OUTPUT_DIR / 'tables').mkdir(exist_ok=True)


# ===========================================================================
# 1. LOAD DATA AND BUILD METADATA
# ===========================================================================
print('=' * 70)
print('STEP 1: Loading data and building metadata')
print('=' * 70)

df = pd.read_csv('diann_matrix.tsv', sep='\t')
df = df.set_index('Protein.Ids')
desc_col = df.pop('Description') if 'Description' in df.columns else None

# Build metadata
wave1_prefixes = {'AC1', 'AHFD1', 'AHFD3', 'AHFD5', 'AHFD7'}
wave2_prefixes = {'AC2', 'AHFD2', 'AHFD4', 'AHFD6', 'AHFD8'}

meta_records = []
for c in df.columns:
    name = c.replace('.raw', '')
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

    if name.startswith('AC'):
        prefix = 'AC' + name[2]
    elif name.startswith('AHFD'):
        prefix = 'AHFD' + name[4]
    else:
        prefix = 'UNK'

    wave = 'Wave1' if prefix in wave1_prefixes else 'Wave2'
    meta_records.append({
        'sample': c, 'name': name, 'group': group, 'wave': wave
    })

meta = pd.DataFrame(meta_records).set_index('sample')
print(f'Loaded {df.shape[0]} proteins x {df.shape[1]} samples')
print(f'Groups: {meta["group"].value_counts().to_dict()}')
print(f'Waves: {meta["wave"].value_counts().to_dict()}')


# ===========================================================================
# 2. CONTAMINANT REMOVAL
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 2: Contaminant removal')
print('=' * 70)

# Common contaminant accessions for rat proteomics
n_before = df.shape[0]
if desc_col is not None:
    contaminant_keywords = ['keratin', 'trypsin', 'serum albumin']
    contaminant_mask = desc_col.str.lower().str.contains(
        '|'.join([k.lower() for k in contaminant_keywords]), na=False
    )
    n_contam = contaminant_mask.sum()
    df = df[~contaminant_mask]
    if desc_col is not None:
        desc_col = desc_col[~contaminant_mask]
    print(f'Removed {n_contam} contaminant proteins')
else:
    print('No Description column available, skipping keyword-based contaminant removal')

print(f'Proteins remaining: {df.shape[0]}')


# ===========================================================================
# 3. LOG2 TRANSFORMATION
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 3: Log2 transformation')
print('=' * 70)

# Replace 0 with NaN before log2
df = df.replace(0, np.nan)
log2_df = np.log2(df)

missing_pct = log2_df.isna().sum().sum() / (log2_df.shape[0] * log2_df.shape[1]) * 100
print(f'Total missingness after log2: {missing_pct:.2f}%')


# ===========================================================================
# 4. CYCLIC LOESS NORMALIZATION
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 4: Cyclic loess normalization')
print('=' * 70)


# Call out to R for limma::normalizeCyclicLoess
log2_df.to_csv('temp_log2.csv')
print('Running limma::normalizeCyclicLoess in R...')
import subprocess
subprocess.run(['Rscript', 'run_limma_loess.R', 'temp_log2.csv', 'temp_norm.csv'], check=True)
norm_df = pd.read_csv('temp_norm.csv', index_col=0)
import os
os.remove('temp_log2.csv')
os.remove('temp_norm.csv')
print('Loess normalization complete.')

# Save boxplot comparison
fig, axes = plt.subplots(1, 2, figsize=(18, 6))
short_labels = [c.replace('.raw', '') for c in log2_df.columns]

bp1 = axes[0].boxplot([log2_df.iloc[:, i].dropna().values for i in range(log2_df.shape[1])],
                       tick_labels=short_labels, vert=True)
axes[0].set_title('Before Normalization (Log2)')
axes[0].tick_params(axis='x', rotation=90, labelsize=7)
axes[0].set_ylabel('Log2 Intensity')

bp2 = axes[1].boxplot([norm_df.iloc[:, i].dropna().values for i in range(norm_df.shape[1])],
                       tick_labels=short_labels, vert=True)
axes[1].set_title('After Cyclic Loess Normalization')
axes[1].tick_params(axis='x', rotation=90, labelsize=7)
axes[1].set_ylabel('Log2 Intensity')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'figures' / 'normalization_boxplots.png', dpi=150)
print('Saved normalization_boxplots.png')


# ===========================================================================
# 5. KNN IMPUTATION (k=5)
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 5: KNN imputation (k=5)')
print('=' * 70)

imputer = KNNImputer(n_neighbors=5)
imputed_values = imputer.fit_transform(norm_df.values)
imputed_df = pd.DataFrame(imputed_values, index=norm_df.index, columns=norm_df.columns)

n_imputed = norm_df.isna().sum().sum()
print(f'Imputed {n_imputed} missing values ({n_imputed / (norm_df.shape[0] * norm_df.shape[1]) * 100:.2f}%)')

# Save imputed matrix
imputed_df.to_csv(OUTPUT_DIR / 'tables' / 'imputed_matrix.csv')
print('Saved imputed_matrix.csv')


# ===========================================================================
# 6. PCA ON NORMALIZED + IMPUTED DATA
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 6: PCA on normalized + imputed data')
print('=' * 70)

X = imputed_df.T.values
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=5)
pcs = pca.fit_transform(X_scaled)
ev = pca.explained_variance_ratio_
print(f'Explained variance (PC1-5): {[round(x * 100, 1) for x in ev]}%')

pca_res = pd.DataFrame({
    'PC1': pcs[:, 0], 'PC2': pcs[:, 1], 'PC3': pcs[:, 2],
    'Group': meta['group'].values,
    'Wave': meta['wave'].values,
    'Sample': meta['name'].values,
})

fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(data=pca_res, x='PC1', y='PC2', hue='Group', style='Wave',
                s=150, ax=ax, palette='Set1')
for _, row in pca_res.iterrows():
    ax.text(row['PC1'] + 0.15, row['PC2'] + 0.15, row['Sample'], fontsize=7)
ax.set_xlabel(f'PC1 ({ev[0] * 100:.1f}%)')
ax.set_ylabel(f'PC2 ({ev[1] * 100:.1f}%)')
ax.set_title(f'PCA after Loess Normalization + KNN Imputation\n'
             f'PC1={ev[0] * 100:.1f}%, PC2={ev[1] * 100:.1f}%')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'figures' / 'pca_final.png', dpi=150)
print('Saved pca_final.png')


# ===========================================================================
# 7. DIFFERENTIAL EXPRESSION: limma eBayes
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 7: Differential expression (limma eBayes)')
print('=' * 70)

# Build design matrix: ~ 0 + group + wave (cell-means model)
design = pd.DataFrame(index=meta.index)
for g in sorted(meta['group'].unique()):
    design[f'group_{g}'] = (meta['group'] == g).astype(int)
# Wave covariate (Wave2 indicator)
design['wave2'] = (meta['wave'] == 'Wave2').astype(int)

print(f'Design matrix shape: {design.shape}')
design.to_csv('temp_design.csv')

# Define contrasts metadata
contrasts = {
    'HFD_vs_CTRL': 'HFD vs CTRL',
    'HFD_Pio_vs_HFD': 'HFD+Pio vs HFD',
    'HFD_IP_vs_HFD': 'HFD+IP vs HFD',
    'HFD_Cr_vs_HFD': 'HFD+Cr vs HFD',
}

# Run limma eBayes in R
print('Running limma::eBayes in R...')
subprocess.run(['Rscript', 'run_limma_ebayes.R', str(OUTPUT_DIR / 'tables' / 'imputed_matrix.csv'), 'temp_design.csv', str(OUTPUT_DIR / 'tables')], check=True)
os.remove('temp_design.csv')

# Load and format results
all_results = {}
for contrast_key, contrast_label in contrasts.items():
    print(f'\n  Processing contrast: {contrast_label}')
    # Read the output from R
    res = pd.read_csv(OUTPUT_DIR / 'tables' / f'limma_res_{contrast_key}.csv')
    res.rename(columns={'Unnamed: 0': 'Protein', 'logFC': 'log2FC', 'P.Value': 'p_value', 'adj.P.Val': 'adj_p_value'}, inplace=True)
    
    # Add description if available
    if desc_col is not None:
        desc_map = desc_col.to_dict()
        res['Description'] = res['Protein'].map(desc_map)

    # Significance flags (FDR < 0.05 AND |log2FC| >= 1.0)
    res['significant_FDR05'] = res['adj_p_value'] < 0.05
    res['significant_FC1'] = (res['adj_p_value'] < 0.05) & (res['log2FC'].abs() >= 1.0)
    
    n_sig_fdr = res['significant_FDR05'].sum()
    n_sig_fc = res['significant_FC1'].sum()
    print(f'    Significant (FDR<0.05): {n_sig_fdr}')
    print(f'    Significant (FDR<0.05 & |log2FC|>=1.0): {n_sig_fc}')
    
    res.to_csv(OUTPUT_DIR / 'tables' / f'de_{contrast_key}.csv', index=False)
    # Delete the raw R output to keep tables clean
    os.remove(OUTPUT_DIR / 'tables' / f'limma_res_{contrast_key}.csv')
    
    all_results[contrast_key] = res


# ===========================================================================
# 8. VOLCANO PLOTS
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 8: Volcano plots')
print('=' * 70)

fig, axes = plt.subplots(2, 2, figsize=(16, 14))
axes_flat = axes.flatten()

for idx, (contrast_key, contrast_label) in enumerate(contrasts.items()):
    ax = axes_flat[idx]
    res = all_results[contrast_key]

    # Compute -log10(adj_p)
    res['neg_log10_padj'] = -np.log10(res['adj_p_value'].clip(lower=1e-300))

    # Color: significant up/down with fold-change threshold
    colors = []
    for _, r in res.iterrows():
        if r['adj_p_value'] < 0.05 and r['log2FC'] >= 1.0:
            colors.append('#e74c3c')  # red = up
        elif r['adj_p_value'] < 0.05 and r['log2FC'] <= -1.0:
            colors.append('#3498db')  # blue = down
        else:
            colors.append('#bdc3c7')  # grey = NS (or not meeting FC threshold)

    ax.scatter(res['log2FC'], res['neg_log10_padj'],
               c=colors, s=15, alpha=0.7, edgecolors='none')
    
    # Threshold lines
    ax.axhline(-np.log10(0.05), ls='--', color='grey', lw=0.8, label='FDR=0.05')
    ax.axvline(1.0, ls=':', color='grey', lw=0.8, label='|log2FC|=1.0')
    ax.axvline(-1.0, ls=':', color='grey', lw=0.8)
    ax.axvline(0, ls='-', color='black', lw=0.5)

    # Label top significant proteins (by absolute fold change among those that are significant)
    sig = res[res['significant_FC1']].copy()
    sig['abs_fc'] = sig['log2FC'].abs()
    top_sig = sig.sort_values(['abs_fc', 'adj_p_value'], ascending=[False, True]).head(10)
    for _, r in top_sig.iterrows():
        ax.text(r['log2FC'], r['neg_log10_padj'], r['Protein'],
                fontsize=7, ha='center', va='bottom')

    n_up = ((res['adj_p_value'] < 0.05) & (res['log2FC'] >= 1.0)).sum()
    n_down = ((res['adj_p_value'] < 0.05) & (res['log2FC'] <= -1.0)).sum()
    ax.set_title(f'{contrast_label}\n'
                 f'{n_up} up, {n_down} down (FDR<0.05 & |log2FC|>=1.0)',
                 fontsize=11)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10(adj. p-value)')
    if idx == 0:
        ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'figures' / 'volcano_plots.png', dpi=150)
print('Saved volcano_plots.png')


# ===========================================================================
# 9. HEATMAP OF SIGNIFICANT PROTEINS (HFD vs CTRL)
# ===========================================================================
print('\n' + '=' * 70)
print('STEP 9: Heatmap of HFD-signature proteins')
print('=' * 70)

hfd_sig = all_results['HFD_vs_CTRL']
# Only include those that passed BOTH FDR and FC thresholds
sig_proteins = hfd_sig[hfd_sig['significant_FC1']]['Protein'].tolist()

if len(sig_proteins) > 1:
    heatmap_data = imputed_df.loc[sig_proteins]

    # Z-score per protein (row)
    heatmap_z = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0).divide(
        heatmap_data.std(axis=1), axis=0
    )

    # Sort columns by group
    col_order = meta.sort_values(['group', 'wave']).index.tolist()
    heatmap_z = heatmap_z[col_order]

    # Column colors
    group_colors_map = {
        'CTRL': '#2ecc71', 'HFD': '#e74c3c',
        'HFD_IP': '#f39c12', 'HFD_Pio': '#9b59b6', 'HFD_Cr': '#3498db'
    }
    col_colors = [group_colors_map.get(meta.loc[c, 'group'], 'grey') for c in col_order]
    col_labels = [c.replace('.raw', '') for c in col_order]
    heatmap_z.columns = col_labels

    # Clustermap
    g = sns.clustermap(heatmap_z, col_cluster=False, row_cluster=True,
                       cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                       col_colors=col_colors,
                       figsize=(14, max(8, len(sig_proteins) * 0.15)),
                       xticklabels=True, yticklabels=(len(sig_proteins) <= 80),
                       dendrogram_ratio=(0.1, 0.05))
    g.fig.suptitle(f'HFD Signature Proteins (n={len(sig_proteins)}, FDR<0.05 & |log2FC|>=1.0)\n'
                   f'Z-scored, columns ordered by group',
                   y=1.02, fontsize=12)
    g.savefig(OUTPUT_DIR / 'figures' / 'heatmap_hfd_signature.png',
              dpi=150, bbox_inches='tight')
    print(f'Saved heatmap with {len(sig_proteins)} significant proteins')
else:
    print(f'Not enough significant proteins ({len(sig_proteins)}) for clustering heatmap.')


# ===========================================================================
# 10. SUMMARY
# ===========================================================================
print('\n' + '=' * 70)
print('SUMMARY')
print('=' * 70)

summary_rows = []
for contrast_key, contrast_label in contrasts.items():
    res = all_results[contrast_key]
    n_sig_fdr = res['significant_FDR05'].sum()
    n_sig_fc = res['significant_FC1'].sum()
    n_up = ((res['adj_p_value'] < 0.05) & (res['log2FC'] >= 1.0)).sum()
    n_down = ((res['adj_p_value'] < 0.05) & (res['log2FC'] <= -1.0)).sum()
    summary_rows.append({
        'Contrast': contrast_label,
        'Sig(FDR<0.05)': n_sig_fdr,
        'Sig(FDR<0.05 & |log2FC|>=1.0)': n_sig_fc,
        'Up-regulated': n_up,
        'Down-regulated': n_down,
    })

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))
summary_df.to_csv(OUTPUT_DIR / 'tables' / 'summary.csv', index=False)

print(f'\nAll results saved to: {OUTPUT_DIR.resolve()}')
print('Done!')
