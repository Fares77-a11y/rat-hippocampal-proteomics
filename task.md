# Downstream Analysis Tasks

- [/] Install R packages (`clusterProfiler`, `org.Rn.eg.db`, `enrichplot`, `DOSE`, `ReactomePA`, `pathview`, `fgsea`, `msigdbr`, `ComplexHeatmap`, `STRINGdb`)
- [x] Install R packages (`clusterProfiler`, `org.Rn.eg.db`, `enrichplot`, `DOSE`, `ReactomePA`, `pathview`, `fgsea`, `msigdbr`, `ComplexHeatmap`, `STRINGdb`)
- [x] Create comprehensive downstream R pipeline (`downstream_analysis.R`).
- [x] Fix Windows R API version mismatches (compile `ragg`, `systemfonts`, `ggplot2` from source using `Rtools44`).
- [x] Run `downstream_analysis.R` and generate ID mapping, Rescue calculations, GSEA, UpSet, STRING, and Heatmaps.
- [x] Debug network API timeouts (KEGG) by shifting to local `org.Rn.eg.db` and offline packages to ensure script completion.
- [x] Generate Trajectory response plots for core HFD signature.
- [x] Copy downstream figures to the artifacts folder and write `walkthrough.md`.
