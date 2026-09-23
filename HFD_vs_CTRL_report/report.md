# Proteomics Differential Expression Report

## Summary
- Samples: 30
- Proteins pre-filter: 1097
- Proteins post-filter/imputation: 1097
- Contrast: `HFD vs CTRL`
- Significant proteins: 65
  - Upregulated: 63
  - Downregulated: 2

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
| Q62881 | 4.191 | 1.646 |
| P28570 | 1.632 | 3.483 |
| B0K020 | 1.537 | 1.675 |
| P11517 | 1.423 | 2.320 |
| P08932 | 1.298 | 3.401 |
| P62959 | 1.220 | 2.220 |
| P31647 | 1.178 | 2.207 |
| P35434 | 1.123 | 2.963 |
| P83868 | 1.123 | 2.310 |
| P40241 | 1.116 | 2.028 |

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| P51156 | -0.741 | 1.964 |
| P24268 | -0.683 | 2.311 |

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`
- Provenance: `ro-crate-metadata.json`

## Disclaimer
ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
