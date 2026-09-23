import pandas as pd
import numpy as np

de = pd.read_csv('corrected_DE_results/tables/de_HFD_vs_CTRL.csv')

# Check: how many proteins share the exact same log2FC and p-value?
grouped = de.groupby(['log2FC', 'p_value']).size().reset_index(name='count')
duplicates = grouped[grouped['count'] > 1].sort_values('count', ascending=False)
print("Proteins with IDENTICAL log2FC and p-value:")
print(duplicates.head(20).to_string(index=False))

# Look at the top cluster
top_fc = duplicates.iloc[0]['log2FC']
top_p = duplicates.iloc[0]['p_value']
cluster = de[(de['log2FC'] == top_fc) & (de['p_value'] == top_p)]
print(f"\nThe largest cluster ({len(cluster)} proteins) all have log2FC={top_fc:.6f}, p={top_p:.2e}")
print("Proteins in this cluster:")
for _, row in cluster.iterrows():
    desc = str(row.get('Description', ''))[:80]
    print(f"  {row['Protein']}: {desc}")

# Also check the actual expression values for these proteins
imputed = pd.read_csv('corrected_DE_results/tables/imputed_matrix.csv', index_col=0)
cluster_prots = cluster['Protein'].tolist()[:5]
print(f"\nExpression values for first 5 proteins in cluster:")
for p in cluster_prots:
    vals = imputed.loc[p].values
    print(f"  {p}: {[f'{v:.4f}' for v in vals[:6]]} ...")

# Check if they are literally identical rows
print(f"\nAre these proteins identical rows in the imputed matrix?")
mat = imputed.loc[cluster['Protein'].tolist()]
nunique = mat.drop_duplicates().shape[0]
print(f"  {len(cluster)} proteins, {nunique} unique expression profiles")
