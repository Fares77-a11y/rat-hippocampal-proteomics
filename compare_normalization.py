import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_matrix(df, name):
    print('=' * 60)
    print(f'DATASET: {name}')
    print(f'Shape: {df.shape}')
    
    # Check for negative values
    n_neg = (df < 0).sum().sum()
    print(f'Negative values: {n_neg}')
    
    # Basic stats per sample
    medians = df.median()
    means = df.mean()
    variances = df.var()
    sums = df.sum()
    
    print(f'Range of Sample Medians: {medians.min():.4f} to {medians.max():.4f} (diff = {medians.max()-medians.min():.4g})')
    print(f'Range of Sample Means: {means.min():.4f} to {means.max():.4f} (diff = {means.max()-means.min():.4g})')
    print(f'Range of Sample Sums: {sums.min():.4e} to {sums.max():.4e} (diff = {sums.max()-sums.min():.4e})')
    
    # Check quantiles (for quantile normalization)
    q25 = df.quantile(0.25)
    q75 = df.quantile(0.75)
    print(f'Range of 25th percentiles: {q25.min():.4f} to {q25.max():.4f} (diff = {q25.max()-q25.min():.4g})')
    print(f'Range of 75th percentiles: {q75.min():.4f} to {q75.max():.4f} (diff = {q75.max()-q75.min():.4g})')
    
    return df

# 1. Hippocampus - Our Loess Normalized
print("Loading Hippocampus (Our Loess Normalized)...")
hippo_df = pd.read_csv('corrected_DE_results/tables/imputed_matrix.csv', index_col=0)
analyze_matrix(hippo_df, 'Hippocampus (Loess + KNN)')

# 2. Serum - Training Anonymized
print("\nLoading Serum Training (Normalized Data)...")
serum1_full = pd.read_excel('Rat_Serum_Training_Anonymized.xlsx', sheet_name='Normalized Data')
idx_col = 'Accession' if 'Accession' in serum1_full.columns else serum1_full.columns[0]
serum1_df = serum1_full.set_index(idx_col)
col_to_drop = [c for c in ['metabolite_id', 'Description'] if c in serum1_df.columns]
serum1_df = serum1_df.drop(columns=col_to_drop)
serum1_df = serum1_df.apply(pd.to_numeric, errors='coerce')
analyze_matrix(serum1_df, 'Serum Training (Pre-normalized)')

# 3. Serum - Pre-Diabetes HFD
print("\nLoading Serum Pre-Diabetes (Original Normalized)...")
serum2_full = pd.read_excel('Rat_Serum_Pre_Diabetes_HFD_Statistical.xlsx', sheet_name='Original Normalized')
idx_col2 = 'Accession' if 'Accession' in serum2_full.columns else serum2_full.columns[0]
serum2_df = serum2_full.set_index(idx_col2)
col_to_drop2 = [c for c in ['metabolite_id', 'Description'] if c in serum2_df.columns]
serum2_df = serum2_df.drop(columns=col_to_drop2)
serum2_df = serum2_df.apply(pd.to_numeric, errors='coerce')
analyze_matrix(serum2_df, 'Serum Pre-Diabetes (Pre-normalized)')

# Let's also check if serum1 and serum2 are actually the same dataset just renamed
print("\n" + "="*60)
print("Comparing Serum 1 and Serum 2")
print(f"Serum 1 shape: {serum1_df.shape}")
print(f"Serum 2 shape: {serum2_df.shape}")
if serum1_df.shape == serum2_df.shape:
    # check if values are close
    val1 = np.sort(serum1_df.values.flatten())
    val2 = np.sort(serum2_df.values.flatten())
    diff = np.nanmax(np.abs(val1 - val2))
    print(f"Max difference between sorted flattened values: {diff}")
    
# Plot boxplots to visualize
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
axes[0].boxplot([hippo_df[c].dropna() for c in hippo_df.columns], labels=['']*hippo_df.shape[1])
axes[0].set_title('Hippocampus (Our Loess)')
axes[1].boxplot([serum1_df[c].dropna() for c in serum1_df.columns], labels=['']*serum1_df.shape[1])
axes[1].set_title('Serum Training (Pre-normalized)')
axes[2].boxplot([serum2_df[c].dropna() for c in serum2_df.columns], labels=['']*serum2_df.shape[1])
axes[2].set_title('Serum Pre-Diabetes (Pre-normalized)')
plt.tight_layout()
plt.savefig('normalization_comparison.png')
print("\nSaved boxplot comparison to normalization_comparison.png")
