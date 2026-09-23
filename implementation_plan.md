# Downstream Analysis Plan — Rat Hippocampal Proteomics

## Context

### Input Data

All input files are under `C:/Users/Ibrah/Downloads/prot/corrected_DE_results/tables/`:

| File | Contents | Key Columns |
|------|----------|-------------|
| `imputed_matrix.csv` | 1,083 × 30 log2 intensity matrix (post-loess, post-KNN) | Row index = UniProt accession, columns = sample names |
| `de_HFD_vs_CTRL.csv` | limma eBayes results for HFD vs CTRL | `Protein`, `log2FC`, `AveExpr`, `t`, `p_value`, `adj_p_value`, `B`, `Description`, `significant_FDR05`, `significant_FC1` |
| `de_HFD_Pio_vs_HFD.csv` | limma eBayes results for HFD+Pio vs HFD | Same columns |
| `de_HFD_IP_vs_HFD.csv` | limma eBayes results for HFD+IP vs HFD | Same columns |
| `de_HFD_Cr_vs_HFD.csv` | limma eBayes results for HFD+Cr vs HFD | Same columns |

### Protein IDs
- All protein IDs are UniProt accessions (e.g., `P48056`, `Q62881`).
- No protein groups (no semicolons in IDs).
- Gene names are embedded in the `Description` column as `GN=GeneName` (1,073 of 1,083 proteins have gene names).
- The organism is *Rattus norvegicus* (taxon 10116).

### Current Significant Hits (FDR < 0.05 & |log2FC| ≥ 1.0)

| Contrast | # Hits | Direction |
|----------|--------|-----------|
| HFD vs CTRL | 9 | 8 UP, 1 DOWN |
| HFD+Pio vs HFD | 7 | 7 UP |
| HFD+IP vs HFD | 0 | — |
| HFD+Cr vs HFD | 1 | 1 UP |

> [!IMPORTANT]
> For enrichment analyses (GO, KEGG, GSEA), we will use the broader **FDR < 0.05 set (117 proteins for HFD vs CTRL)**, not the strict |log2FC| ≥ 1.0 set (9 proteins), because enrichment methods require a reasonable gene list size (≥20) to have statistical power. The fold-change threshold is appropriate for identifying individual high-confidence biological hits, but enrichment analyses look for aggregate pathway-level shifts composed of many moderate-effect-size proteins.

---

## R Package Installation

The following packages need to be installed before running the analyses. Run this once:

```r
# CRAN packages
install.packages(c("UpSetR", "ggplot2", "ggrepel"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "clusterProfiler",   # GO/KEGG ORA and GSEA
    "org.Rn.eg.db",      # Rat gene annotation (UniProt ↔ Entrez ↔ Symbol)
    "enrichplot",        # Enrichment visualization (dot plots, ridge plots, etc.)
    "DOSE",              # Disease Ontology (dependency of clusterProfiler)
    "ReactomePA",        # Reactome pathway analysis
    "pathview",          # KEGG pathway diagram rendering
    "fgsea",             # Fast preranked GSEA
    "msigdbr",           # MSigDB gene sets (Hallmark, C2, C5, etc.)
    "ComplexHeatmap",    # Publication-quality heatmaps
    "STRINGdb"           # STRING PPI network queries
))
```

Save this as [install_packages.R](file:///c:/Users/Ibrah/Downloads/prot/install_packages.R) and run with `Rscript install_packages.R`.

---

## Analyses

All analyses will be implemented in a single R script: **`downstream_analysis.R`**.
All output will go to `corrected_DE_results/downstream/`.

---

### Analysis 1: ID Mapping (Prerequisite for Everything Else)

**Purpose:** Convert UniProt accessions to Entrez Gene IDs and gene symbols, which are required by all enrichment tools.

**Package:** `org.Rn.eg.db` (Bioconductor), `AnnotationDbi` (already installed)

**Algorithm:**
1. Load all 1,083 UniProt accessions from `imputed_matrix.csv`.
2. Call `AnnotationDbi::mapIds(org.Rn.eg.db, keys = uniprot_ids, keytype = "UNIPROT", column = "ENTREZID")` to get Entrez IDs.
3. Call `AnnotationDbi::mapIds(org.Rn.eg.db, keys = uniprot_ids, keytype = "UNIPROT", column = "SYMBOL")` to get gene symbols.
4. As a fallback for any unmapped proteins, parse the `GN=` field from the `Description` column.
5. Save the mapping as `id_mapping.csv` with columns: `Protein`, `EntrezID`, `Symbol`.
6. Print how many mapped successfully (expect ≥90%).

**Output:** `downstream/id_mapping.csv`

---

### Analysis 2: Rescue Analysis

**Purpose:** For each of the 117 HFD-signature proteins (FDR < 0.05), quantify whether each intervention (Pio, IP, Cr) moves the protein back toward CTRL. This directly answers: *"Does the drug rescue the HFD-induced proteomic damage?"*

**Package:** Base R, `ggplot2`, `ggrepel`

**Algorithm:**
1. Load the HFD vs CTRL results. Filter to the 117 FDR < 0.05 proteins.
2. For each protein, extract:
   - `FC_hfd` = log2FC from HFD vs CTRL (the "disease effect")
   - `FC_pio` = log2FC from HFD+Pio vs HFD (the "drug effect")
   - `FC_ip` = log2FC from HFD+IP vs HFD
   - `FC_cr` = log2FC from HFD+Cr vs HFD
3. Compute **% rescue** for each intervention:
   - `rescue_pct = -FC_drug / FC_hfd * 100`
   - If rescue_pct > 0 and ≤ 100, the drug partially rescues.
   - If rescue_pct > 100, the drug over-corrects.
   - If rescue_pct < 0, the drug worsens the HFD effect.
4. Generate a **scatter plot** for each drug: X-axis = `FC_hfd` (disease effect), Y-axis = `FC_drug` (drug effect). If the drug rescues, points should fall in the opposite quadrant from the diagonal. Add a diagonal identity line (y = -x) representing "perfect 100% rescue."
5. Generate a **bar chart** of median % rescue across all 117 proteins for each drug.
6. Save the full table as `rescue_analysis.csv`.

**Outputs:**
- `downstream/rescue_scatter_plots.png`
- `downstream/rescue_barplot.png`
- `downstream/rescue_analysis.csv`

---

### Analysis 3: GO / KEGG / Reactome Over-Representation Analysis (ORA)

**Purpose:** Identify the biological processes, molecular functions, cellular compartments, and metabolic pathways enriched in the HFD-disrupted proteome.

**Packages:** `clusterProfiler`, `org.Rn.eg.db`, `enrichplot`, `ReactomePA`, `pathview`

**Algorithm:**
1. Take the **117 FDR < 0.05 Entrez IDs** from the HFD vs CTRL contrast as the "gene list."
2. Take the **full 1,083 measured Entrez IDs** as the "universe" (background).
3. Run three separate ORA tests:
   - **GO:** `clusterProfiler::enrichGO(gene = sig_entrez, universe = all_entrez, OrgDb = org.Rn.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)`
     - `ont = "ALL"` runs BP, MF, and CC simultaneously.
   - **KEGG:** `clusterProfiler::enrichKEGG(gene = sig_entrez, universe = all_entrez, organism = "rno", pAdjustMethod = "BH", pvalueCutoff = 0.05)`
   - **Reactome:** `ReactomePA::enrichPathway(gene = sig_entrez, universe = all_entrez, organism = "rat", pAdjustMethod = "BH", pvalueCutoff = 0.05)`
4. Visualizations (using `enrichplot`):
   - **Dot plot** of the top 15 terms per ontology: `enrichplot::dotplot(go_result, showCategory=15, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")`
   - **KEGG dot plot:** `enrichplot::dotplot(kegg_result, showCategory=15)`
   - **Reactome dot plot:** `enrichplot::dotplot(reactome_result, showCategory=15)`
   - **cnetplot** (concept network): `enrichplot::cnetplot(go_result, categorySize="pvalue", foldChange=named_fc_vector)` — shows which genes belong to which GO terms.
5. For the top 3 significant KEGG pathways, render the pathway diagram with fold-change coloring: `pathview::pathview(gene.data = named_fc_vector, pathway.id = "rnoXXXXX", species = "rno")`
6. Save all enrichment results as CSV tables.

**Outputs:**
- `downstream/go_enrichment.csv`, `downstream/kegg_enrichment.csv`, `downstream/reactome_enrichment.csv`
- `downstream/go_dotplot.png`
- `downstream/kegg_dotplot.png`
- `downstream/reactome_dotplot.png`
- `downstream/go_cnetplot.png`
- `downstream/pathview_*.png` (KEGG pathway diagrams)

---

### Analysis 4: UpSet Plot (Contrast Overlaps)

**Purpose:** Visualize the overlaps between the significant protein sets from each contrast. This reveals whether the three interventions are hitting the same proteins or entirely different ones.

**Package:** `UpSetR`

**Algorithm:**
1. For each of the 4 contrasts, collect the set of FDR < 0.05 significant protein IDs. (For HFD+IP, the set is empty, so it will show a bar of 0.)
2. Build a binary membership matrix (1,083 rows × 4 columns).
3. Call `UpSetR::upset(data, sets = c("HFD_vs_CTRL", "HFD_Pio_vs_HFD", "HFD_IP_vs_HFD", "HFD_Cr_vs_HFD"), order.by = "freq", main.bar.color = "steelblue")`.
4. This will show how many proteins are unique to HFD, how many overlap between HFD and Pio, etc.

**Output:** `downstream/upset_plot.png`

---

### Analysis 5: Gene Set Enrichment Analysis (GSEA)

**Purpose:** Unlike ORA (which only uses significant proteins), GSEA uses **all 1,083 proteins ranked by their t-statistic**. This is more powerful because it detects coordinated small shifts across entire pathways that individually don't reach significance.

**Packages:** `fgsea`, `msigdbr`, `clusterProfiler`

**Algorithm:**
1. Build the ranked gene list: extract the `t` statistic from the HFD vs CTRL results, name the vector by Entrez Gene ID, sort descending.
2. Retrieve gene sets from MSigDB for rat:
   - `msigdbr::msigdbr(species = "Rattus norvegicus", category = "H")` → Hallmark gene sets
   - `msigdbr::msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = "CP:KEGG")` → KEGG curated pathways
   - `msigdbr::msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "GO:BP")` → GO Biological Process
3. Convert to the list format required by fgsea: `split(x = gene_sets$entrez_gene, f = gene_sets$gs_name)`.
4. Run GSEA: `fgsea::fgsea(pathways = genesets_list, stats = ranked_vector, minSize = 15, maxSize = 500)`.
5. Visualizations:
   - **Table plot** of the top 20 enriched/depleted pathways: `fgsea::plotGseaTable(pathways_list, ranked_vector, fgsea_res)`
   - **Enrichment curve** for the top 3 pathways: `fgsea::plotEnrichment(pathway, ranked_vector)`
   - **Dot plot** using clusterProfiler's `GSEA()` wrapper for more polished output.
6. Save full results table.

**Outputs:**
- `downstream/gsea_hallmark.csv`, `downstream/gsea_kegg.csv`, `downstream/gsea_gobp.csv`
- `downstream/gsea_hallmark_dotplot.png`
- `downstream/gsea_enrichment_curves.png`

---

### Analysis 6: STRING PPI Network

**Purpose:** Check whether the HFD-signature proteins physically interact with each other, forming functional clusters (rather than being random, unrelated hits).

**Package:** `STRINGdb`

**Algorithm:**
1. Initialize STRING: `string_db <- STRINGdb$new(version="12.0", species=10116, score_threshold=400)` (10116 = *Rattus norvegicus*, 400 = medium confidence).
2. Map the 117 FDR < 0.05 proteins: `mapped <- string_db$map(df, "Symbol", removeUnmappedRows=TRUE)`.
3. Get the interaction network: `string_db$plot_network(mapped$STRING_id)`.
4. Run functional enrichment within STRING: `enrichment <- string_db$get_enrichment(mapped$STRING_id)`.
5. Get PPI enrichment p-value to test whether the network has significantly more interactions than expected by chance.
6. Extract clusters: `clustersList <- string_db$get_clusters(mapped$STRING_id, algorithm="mcl")`.
7. Save the network image and the enrichment table.

> [!NOTE]
> The `STRINGdb` R package queries the STRING API, which requires an internet connection. If the API is unreachable, fall back to exporting the gene list and providing a URL to the STRING web interface: `https://string-db.org/cgi/network?identifiers=Slc6a8%0dNol3%0d...&species=10116`.

**Outputs:**
- `downstream/string_network.png`
- `downstream/string_enrichment.csv`
- `downstream/string_clusters.csv`

---

### Analysis 7: Intervention Response Heatmap

**Purpose:** For all proteins that are significant in *any* contrast (union of all FDR < 0.05 hits), show their expression across ALL 5 groups. This reveals whether the interventions push HFD-altered proteins back toward CTRL.

**Package:** `ComplexHeatmap`, `circlize`

**Algorithm:**
1. Collect the union of all FDR < 0.05 proteins across all 4 contrasts.
2. Extract their rows from `imputed_matrix.csv`.
3. Compute **group means** (average the 6 biological replicates per group) → a matrix of N_sig × 5 groups.
4. Z-score each row (protein) across the 5 group means.
5. Build the heatmap:
   - Rows: proteins, clustered by hierarchical clustering.
   - Columns: CTRL, HFD, HFD+IP, HFD+Pio, HFD+Cr (fixed order, NOT clustered).
   - Color: `circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))`.
   - Row annotation: a bar showing the log2FC from HFD vs CTRL, colored by direction.
   - Column annotation: group color bar.
6. Use `ComplexHeatmap::Heatmap()` for full control.

**Output:** `downstream/intervention_heatmap.png`

---

### Analysis 8: Directionality Trajectory Plot

**Purpose:** For each of the 9 core HFD-signature proteins (FDR < 0.05 & |log2FC| ≥ 1.0), plot the mean expression across the 5 experimental groups as a trajectory. This visually shows whether the interventions push the protein back toward the CTRL level.

**Package:** `ggplot2`, `ggrepel`

**Algorithm:**
1. For each of the 9 significant proteins, compute the group mean log2 intensity from `imputed_matrix.csv` (average across 6 rats per group).
2. Define the x-axis order as: CTRL → HFD → HFD+IP → HFD+Pio → HFD+Cr.
3. For each protein, plot a line connecting the 5 group means, with error bars (±SEM from the 6 biological replicates).
4. Arrange into a 3×3 faceted panel using `ggplot2::facet_wrap(~protein, scales="free_y", ncol=3)`.
5. Add a horizontal dashed line at the CTRL mean to make rescue visually obvious.

**Output:** `downstream/trajectory_plot.png`

---

## Verification Plan

### Automated Checks
After all analyses complete, the script will print a summary:
1. ID mapping success rate (expect ≥ 90%).
2. Number of significant GO / KEGG / Reactome terms (expect ≥ 5 for GO).
3. GSEA: number of significant Hallmark sets at FDR < 0.25 (the standard GSEA threshold).
4. STRING PPI enrichment p-value (expect < 0.05 if the hit list is functionally coherent).
5. Confirm all output files exist in `downstream/`.

### Manual Verification
- Check that the rescue scatter plots show Pio points trending toward the "rescue diagonal."
- Check that GO terms make biological sense for an HFD brain experiment (expect terms like "mitochondrial function," "oxidative phosphorylation," "synaptic signaling," "lipid metabolism").
- Check that the trajectory plots show Pio pulling protein levels back toward CTRL for the 9 core hits.

---

## File Structure After Completion

```
corrected_DE_results/
├── downstream/
│   ├── id_mapping.csv
│   ├── rescue_analysis.csv
│   ├── rescue_scatter_plots.png
│   ├── rescue_barplot.png
│   ├── go_enrichment.csv
│   ├── kegg_enrichment.csv
│   ├── reactome_enrichment.csv
│   ├── go_dotplot.png
│   ├── kegg_dotplot.png
│   ├── reactome_dotplot.png
│   ├── go_cnetplot.png
│   ├── pathview_*.png
│   ├── upset_plot.png
│   ├── gsea_hallmark.csv
│   ├── gsea_kegg.csv
│   ├── gsea_gobp.csv
│   ├── gsea_hallmark_dotplot.png
│   ├── gsea_enrichment_curves.png
│   ├── string_network.png
│   ├── string_enrichment.csv
│   ├── string_clusters.csv
│   ├── intervention_heatmap.png
│   └── trajectory_plot.png
└── ...
```
