# Proteomics Differential Expression Report

## Summary
- Samples: 30
- Proteins pre-filter: 1097
- Proteins post-filter/imputation: 1097
- Contrast: `HFD_Cr vs HFD`
- Significant proteins: 36
  - Upregulated: 36
  - Downregulated: 0

## Preprocessing
- Imputation distribution: `figures/imputation_distribution.png`
- PCA (sample clustering): `figures/pca.png`

## Differential Expression Results
- Full results: `tables/de_results.csv`
- Imputed protein matrix: `tables/imputed_proteinGroups.csv`
- Volcano plot: `figures/volcano.png`

### Top 10 Upregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| P34064 | 1.645 | 1.945 |
| P00762 | 1.597 | 2.032 |
| P16446 | 1.377 | 3.812 |
| Q6DGG1 | 1.320 | 3.670 |
| P13086 | 1.285 | 3.710 |
| Q03626 | 1.284 | 4.857 |
| Q01066 | 1.247 | 1.694 |
| Q5XHZ0 | 1.183 | 3.216 |
| Q62745 | 0.914 | 2.605 |
| O35952 | 0.900 | 1.995 |

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| No downregulated proteins | | |

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`
- Provenance: `ro-crate-metadata.json`

## Disclaimer
ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
