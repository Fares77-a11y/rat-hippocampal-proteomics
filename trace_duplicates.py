import pandas as pd
import numpy as np

# Check if these are identical in the RAW data too
raw = pd.read_csv('diann_matrix.tsv', sep='\t')
raw = raw.set_index('Protein.Ids')
desc = raw.pop('Description') if 'Description' in raw.columns else None

dupes = ['P48056', 'P68035', 'P62775', 'P80432', 'P06302',
         'P62138', 'P20759', 'P0C0S7', 'Q63942', 'P04041',
         'P28572', 'P0CC09', 'P13832', 'Q91XV6', 'P49186']

print("=== RAW DATA (before any processing) ===")
for p in dupes[:5]:
    if p in raw.index:
        vals = raw.loc[p].values[:6]
        print(f"  {p}: {[f'{v:.2f}' for v in vals]}")
    else:
        print(f"  {p}: NOT FOUND in raw data")

# Check if they're duplicates in raw
raw_sub = raw.loc[[p for p in dupes if p in raw.index]]
print(f"\n{len(raw_sub)} proteins found in raw data")
nunique_raw = raw_sub.drop_duplicates().shape[0]
print(f"Unique expression profiles in RAW: {nunique_raw}")

# Check the missingness pattern for these proteins
print("\nMissingness pattern in raw data:")
for p in dupes[:5]:
    if p in raw.index:
        n_zero = (raw.loc[p] == 0).sum()
        n_na = raw.loc[p].isna().sum()
        print(f"  {p}: {n_zero} zeros, {n_na} NaN")

# Also check how many total proteins have duplicate expression profiles in the imputed matrix
imputed = pd.read_csv('corrected_DE_results/tables/imputed_matrix.csv', index_col=0)
print(f"\n=== DUPLICATE CHECK ON FULL IMPUTED MATRIX ===")
print(f"Total proteins: {imputed.shape[0]}")
n_unique = imputed.drop_duplicates().shape[0]
print(f"Unique profiles: {n_unique}")
print(f"Duplicate groups: {imputed.shape[0] - n_unique} proteins share a profile with another")

# Find all duplicate groups
dup_mask = imputed.duplicated(keep=False)
dup_df = imputed[dup_mask]
print(f"\nAll proteins involved in duplicates: {dup_df.shape[0]}")

# Group them
dup_groups = dup_df.groupby(list(dup_df.columns)).apply(lambda x: list(x.index)).reset_index(drop=True)
for i, group in enumerate(dup_groups):
    print(f"  Group {i+1} ({len(group)} proteins): {group}")
