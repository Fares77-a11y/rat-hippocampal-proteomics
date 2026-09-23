# Proteomics Differential Expression Report

## Summary
- Samples: 30
- Proteins pre-filter: 1097
- Proteins post-filter/imputation: 1097
- Contrast: `HFD_Pio vs HFD`
- Significant proteins: 16
  - Upregulated: 16
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
| P00762 | 2.080 | 3.383 |
| P04550 | 1.661 | 1.781 |
| P16446 | 1.622 | 5.809 |
| Q6DGG1 | 1.595 | 5.220 |
| P13086 | 1.543 | 5.306 |
| P34064 | 1.378 | 1.724 |
| Q03626 | 1.274 | 4.227 |
| Q5XHZ0 | 1.206 | 3.227 |
| Q62645 | 1.061 | 2.690 |
| Q62745 | 0.988 | 3.742 |

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
