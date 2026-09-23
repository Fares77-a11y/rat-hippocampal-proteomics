# Differential Expression Analysis — Rat Hippocampal Proteomics (Final)

## Pipeline Summary

| Step | Method |
|------|--------|
| Input | 1,097 proteins × 30 rats (6 per group, 3 per wave) |
| Contaminant removal | 14 removed (keratin, trypsin, serum albumin) → **1,083 proteins** |
| Transformation | Log2 |
| Normalization | **Cyclic loess (`limma::normalizeCyclicLoess` with `method="fast"`)** |
| Imputation | KNN (k=5). *Missingness strictly maintained at 0.19% (63 values).* |
| Statistical model | **`limma::eBayes`** empirical Bayes hierarchical model (`~ 0 + group + wave`) |
| Multiple testing | Benjamini-Hochberg FDR |
| Significance Threshold | **FDR < 0.05 AND \|log2FC\| ≥ 1.0 (2-fold change)** |

---

## Results Summary

By fixing the limma NA propagation bug, switching to the robust empirical Bayes model (`limma eBayes`), and adding a strict 2-fold change threshold, we have eliminated the false positives and isolated the proteins with genuine, massive biological effect sizes.

| Contrast | Significant (FDR<0.05) | Meaningful (FDR<0.05 & \|log2FC\| ≥ 1.0) | Up-regulated | Down-regulated |
|----------|------------------------|-----------------------------------------|-----|------|
| **HFD vs CTRL** | 117 | **9** | 8 | 1 |
| **HFD+Pio vs HFD** | 13 | **7** | 7 | 0 |
| **HFD+IP vs HFD** | 0 | **0** | 0 | 0 |
| **HFD+Cr vs HFD** | 1 | **1** | 1 | 0 |

> [!IMPORTANT]
> The original hit list was vastly inflated by proteins with tiny (e.g. 1.05-fold) changes that were statistically significant but biologically meaningless.
>
> With a strict 2-fold biological threshold:
> - **HFD induces a massive change in only 9 core proteins.**
> - **Pioglitazone does NOT rescue any of those 9.** Its 7 significant proteins have **zero overlap** with the HFD signature — they are Pio's own independent perturbation on proteins HFD didn't affect (all 7 upregulated; all had HFD padj > 0.43). Furthermore, across the 117 FDR-significant HFD proteins, Pio **worsens** 62% of them (median rescue = −12.3%).
> - **IP** has zero significant proteins — no detectable hippocampal effect.
> - **Creatine** has 1 significant protein and shows the best (though still modest) rescue across the 117 HFD proteins (median rescue = +7.9%).

---

## PCA — Final (Limma Loess + KNN Imputed)

PC1 explains 30.7% and PC2 explains 12.2% of variance.

![PCA Final](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/pca_final.png)

---

## Volcano Plots (All 4 Contrasts)

Note the strict new boundaries: hits must cross both the horizontal FDR=0.05 threshold and the vertical \|log2FC\|=1.0 (2-fold) thresholds.

![Volcano Plots](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/volcano_plots.png)

---

## Heatmap — 9-Protein Core HFD Signature

Z-scored expression of the 9 proteins that pass the strictest biological and statistical thresholds for the HFD vs CTRL contrast.

![Heatmap](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/heatmap_hfd_signature.png)

---

## Output Files

All results are updated under [corrected_DE_results/](file:///c:/Users/Ibrah/Downloads/prot/corrected_DE_results):

```
corrected_DE_results/
├── figures/
│   ├── normalization_boxplots.png
│   ├── pca_final.png
│   ├── volcano_plots.png
│   └── heatmap_hfd_signature.png
└── tables/
    ├── imputed_matrix.csv
    ├── de_HFD_vs_CTRL.csv
    ├── de_HFD_Pio_vs_HFD.csv
    ├── de_HFD_IP_vs_HFD.csv
    ├── de_HFD_Cr_vs_HFD.csv
    └── summary.csv
```
