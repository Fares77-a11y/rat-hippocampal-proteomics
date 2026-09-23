# Proteomics Differential Expression Report

## Summary
- Samples: 30
- Proteins pre-filter: 1097
- Proteins post-filter/imputation: 1097
- Contrast: `HFD_IP vs HFD`
- Significant proteins: 4
  - Upregulated: 4
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
| P20171 | 0.930 | 2.497 |
| Q62745 | 0.884 | 3.916 |
| Q08877 | 0.737 | 2.567 |
| Q5XIA3 | 0.657 | 2.120 |

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
