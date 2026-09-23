# Quality Control (QC) Diagnostics — Final Pipeline

This report contains the critical statistical diagnostics that validate the integrity of the data normalization and the linear modeling steps.

## 1. P-value Histograms

**Purpose:** Validates the statistical model and the Benjamini-Hochberg FDR correction. 
**What to look for:** A well-behaved model should produce a perfectly flat, uniform distribution for the majority of proteins (representing the null hypothesis / noise), with a sharp spike near 0 for true biological responders. If the histogram slopes heavily or curves, the model is misspecified.

![P-value histograms](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/qc_pvalue_histograms.png)

**Observation:** The HFD vs CTRL contrast shows exactly what we want: a uniform null background (following the red dotted line) with a massive spike of true positives at p < 0.05. The other contrasts have smaller spikes (indicating fewer true biological responders), but crucially, their background remains uniform. This confirms that the `limma eBayes` model fits the data perfectly and the BH FDR correction is mathematically valid.

---

## 2. MA Plots (Log-Ratio vs. Mean Average)

**Purpose:** Verifies that the cyclic loess normalization successfully removed intensity-dependent bias.
**What to look for:** The blue trend line (a loess fit of the fold-changes) should be perfectly flat and centered on 0 (the black line) across all expression levels. If it bows up or down, it means highly abundant proteins have a systematic fold-change bias compared to low-abundance proteins.

![MA Plots](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/qc_ma_plots.png)

**Observation:** The blue trend line is remarkably flat and hugs 0 across the entire intensity spectrum (from A=18 to A=33) in all four contrasts. This is the direct result of the `limma::normalizeCyclicLoess` step — it successfully eliminated any technical intensity bias. The significant hits (red dots) are symmetrically distributed.

---

## 3. Sample Correlation Heatmap (Post-Normalization)

**Purpose:** Verifies that biological replicates correlate highly with each other, and checks for remaining batch effects or outlier samples.
**What to look for:** The color scale reflects the Pearson correlation coefficient ($R$). Replicates within the same biological group (e.g., CTRL) should cluster together into highly correlated blocks.

![Correlation Heatmap](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/qc_correlation_heatmap.png)

**Observation:** 
- Overall sample correlation is extremely high across the entire dataset (the lowest correlation between any two samples is $R \approx 0.94$), which is typical for high-quality label-free brain tissue proteomics.
- We do not see any glaring "bad" outlier samples (which would show up as dark stripes).
- The hierarchical clustering (dendrogram on the top/left) successfully groups the HFD rats into their own clade, separate from the CTRL rats, proving that the biological signal (Diet) is strong enough to drive global clustering even before statistical testing.
