import pandas as pd, re

df = pd.read_csv('corrected_DE_results/tables/de_HFD_vs_CTRL.csv')
print('Sample Protein IDs:')
print(df['Protein'].head(10).tolist())
print()

has_semicolon = df['Protein'].str.contains(';').sum()
print(f'Protein groups (IDs with semicolons): {has_semicolon}')
print(f'Total proteins: {len(df)}')
print()

print('Sample descriptions with gene names:')
for _, row in df.head(10).iterrows():
    d = str(row['Description'])
    gn = re.search(r'GN=(\S+)', d)
    gene = gn.group(1) if gn else 'MISSING'
    print(f"  {row['Protein']}: GN={gene}")

# Count how many have gene names
n_with_gn = sum(1 for d in df['Description'].dropna() if 'GN=' in str(d))
print(f"\nProteins with GN= in description: {n_with_gn} / {len(df)}")
