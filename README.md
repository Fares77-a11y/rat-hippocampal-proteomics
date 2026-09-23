# Rat Hippocampal Proteomics Analysis: HFD and Pharmacological Interventions

## Overview
This repository contains the complete analysis pipeline, raw data, and downstream interpretation for a large-scale proteomics study on rat hippocampi. The study investigates the proteomic perturbations induced by a **High-Fat Diet (HFD)** and evaluates the efficacy of three distinct pharmacological interventions:
1. **Pioglitazone (Pio)**: A PPARγ agonist (commonly used for insulin sensitization).
2. **Inorganic Phosphate (IP)**
3. **Creatine (Cr)**

## Repository Structure
```text
rat-hippocampal-proteomics/
├── data/              # Input DIA-NN intensity matrices and sample metadata
├── scripts/           # Python and R scripts for pipeline execution
│   ├── corrected_de_analysis.py   # Main DE pipeline (normalization, imputation)
│   ├── run_limma_ebayes.R         # R script for empirical Bayes statistics
│   ├── run_limma_loess.R          # Cyclic Loess normalization
│   └── downstream_analysis.R      # Functional enrichment, trajectories, & rescue
├── results/
│   ├── tables/        # Final CSVs of Differential Expression & Rescue Analysis
│   └── figures/       # High-res PNGs of all visualizations
└── reports/           # Detailed markdown reports evaluating the pipeline and biology
    ├── comprehensive_review.md    # Critical review of the biological narrative & bugs
    └── proteomics_summary.md      # High-level statistical summary
```

## The Biological Story & Interpretations

### 1. The HFD Insult: Widespread but Moderate
HFD disrupted ~11% of the detectable hippocampal proteome (117 out of 1,083 proteins at FDR < 0.05). However, only 9 "core" proteins passed a strict 2-fold change threshold (|log2FC| ≥ 1.0). 
Key biological themes altered by HFD include:
*   **Neurotransmitter Transport**: Upregulation of `Slc6a12` (GABA reuptake) and `Slc6a8` (Creatine transporter), indicating altered GABAergic signaling and energy buffering.
*   **Mitochondrial & Metabolic Stress**: Multiple metabolic enzymes (e.g., `Cisd1`, `Cox5a`, `Aldoa`) trended upward, reflecting compensatory energy production in the hippocampus.
*   **Vascular/Blood-Brain Barrier (BBB) Compromise**: Proteins like `Hbb-bs` (hemoglobin) and `Kng1` (kininogen) were elevated, pointing either to BBB leakage or neuroinflammation.

### 2. Pioglitazone: The "Paradoxical Worsening" Effect
The most striking finding of this study is that **Pioglitazone does NOT rescue the HFD-induced hippocampal damage.**
*   **No Overlap**: Pio significantly altered 7 proteins—**zero** of which overlap with the 9 core HFD-damaged proteins. Pio introduces its own independent pharmacological perturbation (e.g., upregulating mitochondrial enzymes and serine proteases) that HFD did not affect.
*   **Exacerbation**: Across the 117 HFD-dysregulated proteins, Pio actually *worsened* 62% of them, pushing their expression further away from healthy CTRL baseline levels (Median Rescue = -12.3%).

### 3. IP and Creatine Interventions
*   **Inorganic Phosphate (IP)**: Showed zero statistically significant protein changes (FDR < 0.05). It had no detectable restorative effect on the hippocampal proteome.
*   **Creatine (Cr)**: The best performer, though still modest (Median Rescue = +7.9%). Creatine partially normalized proteins like `Slc6a8` (the creatine transporter), suggesting a partial restoration of the disrupted creatine-phosphocreatine energy buffering system.

## Significant Visualizations

You can view the full suite of visualizations in the [`results/figures/`](results/figures/) directory. Highlights include:

1. **[Principal Component Analysis (PCA)](results/figures/pca_final.png)**
   Demonstrates the global variance structure after cyclic loess normalization and KNN imputation.
2. **[HFD Signature Heatmap](results/figures/heatmap_hfd_signature.png)**
   Z-scored expression clustering of the 9 core HFD signature proteins showing clean separation between healthy CTRL and HFD brains.
3. **[Intervention Heatmap](results/figures/intervention_heatmap.png)**
   Shows how the three treatments (Pio, IP, Cr) modulate the broader HFD-dysregulated proteome. Notice how the `HFD_Pio` column visually resembles the damaged HFD state more than the healthy CTRL state.
4. **[Protein Trajectories](results/figures/trajectory_plot.png)**
   Longitudinal mean intensity plots for individual core proteins across the 5 experimental groups, illustrating exact therapeutic responses (or exacerbations).
5. **[Rescue Effect Barplot](results/figures/rescue_barplot.png)**
   A macroscopic view of the median % rescue across all 117 HFD-signature proteins, emphasizing the negative effect of Pioglitazone.

## Analysis Methodology
*   **Missing Value Handling**: Conservative K-Nearest Neighbors (KNN, k=5) imputation strictly maintained to impact only 0.19% of the dataset.
*   **Normalization**: Cyclic Loess (`limma::normalizeCyclicLoess`) utilizing `method="fast"` to correct intensity-dependent biases without NA propagation.
*   **Statistical Modeling**: Cell-means parameterization with a wave covariate block (`~ 0 + group + wave`) fitted using `limma::eBayes` for empirical Bayes variance shrinkage, preventing false positives from artificially low-variance proteins.
*   **Significance**: Proteins were rigorously filtered using both Benjamini-Hochberg FDR < 0.05 and a strict 2-fold biological threshold (|log2FC| ≥ 1.0).
