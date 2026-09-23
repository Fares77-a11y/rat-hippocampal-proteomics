import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 1. Check p-value distributions (critical QC)
fig, axes = plt.subplots(1, 4, figsize=(20, 4))
contrasts = ['de_HFD_vs_CTRL', 'de_HFD_Pio_vs_HFD', 'de_HFD_IP_vs_HFD', 'de_HFD_Cr_vs_HFD']
labels = ['HFD vs CTRL', 'HFD+Pio vs HFD', 'HFD+IP vs HFD', 'HFD+Cr vs HFD']

for i, (cname, label) in enumerate(zip(contrasts, labels)):
    df = pd.read_csv(f'corrected_DE_results/tables/{cname}.csv')
    pvals = df['p_value'].dropna()
    axes[i].hist(pvals, bins=50, color='steelblue', edgecolor='white', density=True)
    axes[i].axhline(1.0, color='red', ls='--', lw=1, label='Uniform (null)')
    axes[i].set_title(label, fontsize=10)
    axes[i].set_xlabel('p-value')
    axes[i].set_ylabel('Density')
    axes[i].set_xlim(0, 1)
    axes[i].legend(fontsize=7)
    
    # Count how many have very small FC
    sig = df[df['adj_p_value'] < 0.05]
    tiny_fc = (sig['log2FC'].abs() < 0.5).sum()
    total_sig = len(sig)
    print(f'{label}: {total_sig} sig (FDR<0.05), of which {tiny_fc} have |log2FC| < 0.5')

plt.tight_layout()
plt.savefig('corrected_DE_results/figures/pvalue_histograms.png', dpi=150)
print('\nSaved p-value histograms')

# 2. Check if there's any protein filtering by coverage
imputed = pd.read_csv('corrected_DE_results/tables/imputed_matrix.csv', index_col=0)
print(f'\nFinal matrix: {imputed.shape[0]} proteins x {imputed.shape[1]} samples')

# 3. Check the top hit fold change
de = pd.read_csv('corrected_DE_results/tables/de_HFD_vs_CTRL.csv')
print(f'\nTop 5 HFD vs CTRL hits by p-value:')
for _, row in de.head(5).iterrows():
    fc = 2**abs(row['log2FC'])
    print(f"  {row['Protein']}: log2FC={row['log2FC']:.4f} ({fc:.2f}-fold), p={row['p_value']:.2e}, adj_p={row['adj_p_value']:.6f}")
