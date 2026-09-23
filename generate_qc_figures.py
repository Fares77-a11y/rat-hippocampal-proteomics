import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

OUTPUT_DIR = Path('corrected_DE_results')
FIG_DIR = OUTPUT_DIR / 'figures'

contrasts = {
    'HFD_vs_CTRL': 'HFD vs CTRL',
    'HFD_Pio_vs_HFD': 'HFD+Pio vs HFD',
    'HFD_IP_vs_HFD': 'HFD+IP vs HFD',
    'HFD_Cr_vs_HFD': 'HFD+Cr vs HFD',
}

# ===========================================================================
# 1. P-VALUE HISTOGRAMS
# ===========================================================================
print("Generating P-value histograms...")
fig, axes = plt.subplots(1, 4, figsize=(20, 4))
axes_flat = axes.flatten()

for i, (contrast_key, label) in enumerate(contrasts.items()):
    df = pd.read_csv(OUTPUT_DIR / 'tables' / f'de_{contrast_key}.csv')
    pvals = df['p_value'].dropna()
    
    ax = axes_flat[i]
    ax.hist(pvals, bins=40, color='steelblue', edgecolor='white', density=True)
    ax.axhline(1.0, color='red', ls='--', lw=1.5, label='Uniform (null)')
    
    ax.set_title(label, fontsize=12)
    ax.set_xlabel('Raw p-value', fontsize=10)
    if i == 0:
        ax.set_ylabel('Density', fontsize=10)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=9, loc='upper right')

plt.tight_layout()
fig.savefig(FIG_DIR / 'qc_pvalue_histograms.png', dpi=150)
print("  Saved qc_pvalue_histograms.png")


# ===========================================================================
# 2. MA PLOTS
# ===========================================================================
print("\nGenerating MA plots...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes_flat = axes.flatten()

for i, (contrast_key, label) in enumerate(contrasts.items()):
    df = pd.read_csv(OUTPUT_DIR / 'tables' / f'de_{contrast_key}.csv')
    ax = axes_flat[i]
    
    # limma output should have AveExpr, but if not we can just use M=log2FC
    # Let's check if AveExpr is in columns
    A = df['AveExpr'] if 'AveExpr' in df.columns else df.mean(axis=1) # Fallback if AveExpr missing, though limma topTable includes it
    M = df['log2FC']
    
    # Colors: red for significant (FDR<0.05 & |log2FC|>=1.0), grey otherwise
    is_sig = (df['adj_p_value'] < 0.05) & (df['log2FC'].abs() >= 1.0)
    colors = ['#e74c3c' if sig else '#bdc3c7' for sig in is_sig]
    alphas = [0.8 if sig else 0.4 for sig in is_sig]
    
    ax.scatter(A[~is_sig], M[~is_sig], c='#bdc3c7', s=10, alpha=0.4, edgecolors='none', label='NS')
    if is_sig.sum() > 0:
        ax.scatter(A[is_sig], M[is_sig], c='#e74c3c', s=20, alpha=0.9, edgecolors='none', label='Significant')
        
    ax.axhline(0, color='black', lw=1, ls='-')
    ax.axhline(1.0, color='grey', lw=0.8, ls=':')
    ax.axhline(-1.0, color='grey', lw=0.8, ls=':')
    
    # Fit a smooth trend line to check for intensity-dependent bias
    from statsmodels.nonparametric.smoothers_lowess import lowess
    trend = lowess(M, A, frac=0.3, return_sorted=True)
    ax.plot(trend[:, 0], trend[:, 1], color='blue', lw=2, label='Trend (loess)')
    
    ax.set_title(label, fontsize=12)
    ax.set_xlabel('Average Expression (A)', fontsize=10)
    ax.set_ylabel('Log2 Fold Change (M)', fontsize=10)
    if i == 0:
        ax.legend(fontsize=9)

plt.tight_layout()
fig.savefig(FIG_DIR / 'qc_ma_plots.png', dpi=150)
print("  Saved qc_ma_plots.png")


# ===========================================================================
# 3. CORRELATION HEATMAP (Post-normalization)
# ===========================================================================
print("\nGenerating Correlation Heatmap...")
imputed_df = pd.read_csv(OUTPUT_DIR / 'tables' / 'imputed_matrix.csv', index_col=0)

# Compute Pearson correlation between samples
corr_matrix = imputed_df.corr(method='pearson')

# Extract metadata for color coding
wave1_prefixes = {'AC1', 'AHFD1', 'AHFD3', 'AHFD5', 'AHFD7'}
meta_records = []
for c in imputed_df.columns:
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
    meta_records.append({'sample': c, 'group': group, 'wave': wave, 'name': name})

meta = pd.DataFrame(meta_records).set_index('sample')

# Set up colors
group_colors_map = {
    'CTRL': '#2ecc71', 'HFD': '#e74c3c',
    'HFD_IP': '#f39c12', 'HFD_Pio': '#9b59b6', 'HFD_Cr': '#3498db'
}
wave_colors_map = {'Wave1': '#34495e', 'Wave2': '#95a5a6'}

col_colors = pd.DataFrame({
    'Group': meta['group'].map(group_colors_map),
    'Wave': meta['wave'].map(wave_colors_map)
})

# Make labels shorter
corr_matrix.columns = meta['name']
corr_matrix.index = meta['name']
col_colors.index = meta['name']

g = sns.clustermap(corr_matrix, 
                   cmap='viridis',
                   row_colors=col_colors, 
                   col_colors=col_colors,
                   figsize=(10, 10),
                   vmin=corr_matrix.min().min(), vmax=1.0,
                   dendrogram_ratio=0.15)

g.fig.suptitle('Sample Correlation Matrix (Post-Normalization + Imputation)\nPearson R', y=1.02, fontsize=14)
g.savefig(FIG_DIR / 'qc_correlation_heatmap.png', dpi=150, bbox_inches='tight')
print("  Saved qc_correlation_heatmap.png")

print("\nDone!")
