# Analysis Plan Review — Critical Assessment

## Overall Verdict

The pipeline structure (normalize → impute → model → correct) follows the correct standard workflow for label-free proteomics. However, this review uncovered **one critical data corruption issue** and several important gaps.

---

## 🔴 Critical Issue: 15 Proteins Collapsed to Identical Profiles

During the code review, I flagged Bug #3 (limma NA propagation) as "low severity." I was wrong — it's **critical**.

**What happened:**
1. The raw data had 63 missing values (0.20% missingness).
2. `limma::normalizeCyclicLoess(method="pairs")` propagated NAs during pairwise loess fitting, inflating missingness from 63 → 450 values.
3. KNN imputation (k=5) then filled all 450 missing values.
4. **15 proteins that originally had high individual missingness (5–9 NaN each) were all imputed to the exact same expression profile across all 30 samples.** They now share identical log2FC = −0.066, p = 9.5e-8, adj_p = 7e-6.

These 15 proteins include biologically unrelated entries (Actin, Histone H2A, Glutathione peroxidase, MAPK9, etc.) — they should not have identical expression. They are false positives inflating the hit count.

**Fix:** Either:
- Switch limma to `method="fast"` (reference-based, avoids NA propagation), or
- Impute first (on log2 data), then normalize (on complete data), or
- Filter out proteins with high missingness (e.g., require ≥50% valid values) before normalization

---

## 🟡 Missing: Fold-Change Threshold

Currently, a protein is called "significant" based solely on FDR < 0.05, with **no minimum fold-change requirement**. This is a problem:

| Contrast | Total Sig (FDR<0.05) | With \|log2FC\| < 0.5 | With \|log2FC\| ≥ 0.5 |
|----------|---------------------|----------------------|----------------------|
| HFD vs CTRL | 142 | **57 (40%)** | 85 |
| HFD+Pio vs HFD | 49 | **30 (61%)** | 19 |
| HFD+IP vs HFD | 15 | **15 (100%)** | 0 |
| HFD+Cr vs HFD | 33 | **27 (82%)** | 6 |

The top 5 "most significant" HFD vs CTRL proteins all have log2FC = −0.066, which is a **1.05-fold change** — essentially no biological change. The entire HFD+IP contrast has zero proteins with even a 1.4-fold change.

> [!WARNING]
> Without a fold-change cutoff, the hit list is dominated by statistically significant but biologically meaningless changes. Standard practice is to require **|log2FC| ≥ 0.5** (1.4-fold) or **|log2FC| ≥ 1.0** (2-fold) alongside FDR < 0.05.

---

## 🟡 Missing QC Figures

The current pipeline produces 4 figures. A publication-ready analysis needs several more:

| Figure | Status | Purpose |
|--------|--------|---------|
| Normalization boxplots | ✅ Done | Verify median alignment |
| PCA (group + wave) | ✅ Done | Sample clustering |
| Volcano plots (×4) | ✅ Done | Visualize DE results |
| Heatmap (HFD signature) | ✅ Done | Expression patterns |
| **P-value histograms** | ⚠️ Just generated | **Critical** — verify the null distribution is uniform (validates that the BH correction is appropriate) |
| **MA plots** | ❌ Missing | Verify normalization removed intensity-dependent bias |
| **Sample correlation heatmap** (post-norm) | ❌ Missing | Verify within-group samples cluster together after normalization |
| **Variance partition bar chart** | ❌ Missing | Show % variance explained by Group vs Wave vs Residual |
| **UpSet plot** | ❌ Missing | Show overlaps between the 4 contrasts |

### P-value Histograms (just generated)

![P-value histograms](C:/Users/Ibrah/.gemini/antigravity/brain/88a797f7-6694-4716-8a69-133e377d6677/figures/pvalue_histograms.png)

A well-behaved p-value distribution should show a flat (uniform) distribution with a spike near 0 for true positives. Deviations from uniformity indicate model misspecification.

---

## 🟡 Statistical Model Could Be Improved

### Current: Per-protein OLS
We fit an independent OLS model for each of 1,075 proteins. Each protein's variance is estimated from only its own 30 data points.

### Better: limma's Empirical Bayes (eBayes)
`limma` uses a hierarchical model that "borrows strength" across all proteins. For each protein, it shrinks the per-protein variance estimate toward the global trend. This has two benefits:
- **Prevents false positives** from proteins that accidentally have near-zero variance (which inflates the t-statistic)
- **Increases power** for proteins with noisy variance estimates

Since we already call R for the loess normalization, it would be straightforward to also run the full `limma::lmFit` + `limma::eBayes` pipeline in R rather than fitting OLS in Python. This is the gold standard for proteomics DE.

---

## 🟢 What Was Done Well

- ✅ Correct sample-to-group and wave mapping (verified exhaustively)
- ✅ Correct cell-means design matrix with wave covariate
- ✅ Correct contrast extraction (c'β and √(c'Σc))
- ✅ Exact scipy t-distribution p-values
- ✅ Standard BH FDR correction
- ✅ Using limma for normalization (the gold standard)

---

## Downstream Analyses Roadmap

Once the bugs above are fixed, the following analyses would complete a publication-ready story:

### Tier 1 — Essential (directly answer the biological question)

| Analysis | Purpose |
|----------|---------|
| **Rescue analysis** | For each HFD-signature protein, check if each intervention moves it back toward CTRL. Quantify % rescue per intervention. |
| **GO / Pathway enrichment** | What biological processes are disrupted by HFD in the hippocampus? (Reactome, KEGG, GO Biological Process) |
| **UpSet / Venn diagram** | Which proteins overlap between interventions? Are Pio, IP, and Cr hitting the same or different pathways? |

### Tier 2 — Strongly Recommended

| Analysis | Purpose |
|----------|---------|
| **GSEA (Gene Set Enrichment Analysis)** | Uses ALL proteins ranked by fold-change (not just significant ones). More powerful than threshold-based enrichment. |
| **STRING PPI network** | Are the HFD-signature proteins physically interacting? Do they form functional clusters? |
| **Heatmap of intervention responders** | Show how the Pio/IP/Cr-responsive proteins behave across ALL groups (not just the two being compared). |

### Tier 3 — Nice to Have

| Analysis | Purpose |
|----------|---------|
| **Kinase/TF enrichment** | Predict which upstream kinases or transcription factors drive the HFD signature. |
| **Comparison with published HFD brain proteomics** | Validate findings against literature. |
| **Dose-response directionality plot** | For each HFD-altered protein, plot its expression across CTRL → HFD → HFD+Drug as a "trajectory." |

---

## Recommended Action Plan

1. **Fix the NA propagation bug** — switch to `method="fast"` or impute before normalizing
2. **Add |log2FC| ≥ 0.5 threshold** alongside FDR < 0.05
3. **Switch to limma eBayes** for the statistical testing (do the full analysis in R)
4. **Add the missing QC figures** (p-value histograms, MA plots, correlation heatmap)
5. **Run rescue analysis + pathway enrichment** as the primary downstream analyses

Shall I implement these fixes and run the corrected pipeline?
