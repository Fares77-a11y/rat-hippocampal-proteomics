import pandas as pd, re

for contrast in ['de_HFD_vs_CTRL', 'de_HFD_Pio_vs_HFD', 'de_HFD_IP_vs_HFD', 'de_HFD_Cr_vs_HFD']:
    df = pd.read_csv(f'corrected_DE_results/tables/{contrast}.csv')
    sig_fc = df[df['significant_FC1'] == True]
    sig_fdr = df[df['significant_FDR05'] == True]
    print(f"=== {contrast} ===")
    print(f"  FDR<0.05 only: {len(sig_fdr)}")
    print(f"  FDR<0.05 & |log2FC|>=1: {len(sig_fc)}")
    for _, row in sig_fc.iterrows():
        gn = re.search(r'GN=(\S+)', str(row['Description']))
        gene = gn.group(1) if gn else 'NA'
        direction = 'UP' if row['log2FC'] > 0 else 'DOWN'
        print(f"    {row['Protein']} ({gene}): log2FC={row['log2FC']:.3f}, adj_p={row['adj_p_value']:.2e} [{direction}]")
    print()
