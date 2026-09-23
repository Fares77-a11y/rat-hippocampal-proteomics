library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(pathview)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(STRINGdb)
library(UpSetR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(circlize)

dir.create("corrected_DE_results/downstream", showWarnings = FALSE)

# ===========================================================================
# 0. Load Data
# ===========================================================================
cat("Loading data...\n")
imputed_df <- read.csv("corrected_DE_results/tables/imputed_matrix.csv", row.names = 1)
de_hfd <- read.csv("corrected_DE_results/tables/de_HFD_vs_CTRL.csv")
de_pio <- read.csv("corrected_DE_results/tables/de_HFD_Pio_vs_HFD.csv")
de_ip <- read.csv("corrected_DE_results/tables/de_HFD_IP_vs_HFD.csv")
de_cr <- read.csv("corrected_DE_results/tables/de_HFD_Cr_vs_HFD.csv")

sig_hfd <- de_hfd %>% filter(adj_p_value < 0.05)
uniprot_universe <- de_hfd$Protein

# ===========================================================================
# 1. ID Mapping
# ===========================================================================
cat("Analysis 1: ID Mapping...\n")
map_df <- data.frame(Protein = uniprot_universe)
map_df$EntrezID <- mapIds(org.Rn.eg.db, keys = map_df$Protein, keytype = "UNIPROT", column = "ENTREZID", multiVals = "first")
map_df$Symbol <- mapIds(org.Rn.eg.db, keys = map_df$Protein, keytype = "UNIPROT", column = "SYMBOL", multiVals = "first")

fallback_symbols <- sapply(de_hfd$Description, function(x) {
    m <- regexpr("GN=\\S+", x)
    if (m > 0) return(sub("GN=", "", regmatches(x, m))) else return(NA)
})
names(fallback_symbols) <- de_hfd$Protein
map_df$Symbol <- ifelse(is.na(map_df$Symbol), fallback_symbols[map_df$Protein], map_df$Symbol)

write.csv(map_df, "corrected_DE_results/downstream/id_mapping.csv", row.names = FALSE)

de_hfd <- left_join(de_hfd, map_df, by = "Protein")
sig_hfd <- left_join(sig_hfd, map_df, by = "Protein")


# ===========================================================================
# 2. Rescue Analysis
# ===========================================================================
cat("\nAnalysis 2: Rescue Analysis...\n")
rescue_df <- data.frame(Protein = sig_hfd$Protein, Symbol = sig_hfd$Symbol, FC_hfd = sig_hfd$log2FC)

rescue_df$FC_pio <- de_pio$log2FC[match(rescue_df$Protein, de_pio$Protein)]
rescue_df$FC_ip <- de_ip$log2FC[match(rescue_df$Protein, de_ip$Protein)]
rescue_df$FC_cr <- de_cr$log2FC[match(rescue_df$Protein, de_cr$Protein)]

rescue_df$RescuePct_Pio <- -rescue_df$FC_pio / rescue_df$FC_hfd * 100
rescue_df$RescuePct_IP <- -rescue_df$FC_ip / rescue_df$FC_hfd * 100
rescue_df$RescuePct_Cr <- -rescue_df$FC_cr / rescue_df$FC_hfd * 100

write.csv(rescue_df, "corrected_DE_results/downstream/rescue_analysis.csv", row.names = FALSE)

p_pio <- ggplot(rescue_df, aes(x = FC_hfd, y = FC_pio)) +
    geom_point(alpha = 0.7, color = "#9b59b6") +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_text_repel(data = rescue_df %>% filter(abs(FC_pio) > 0.5), aes(label = Symbol), size = 3) +
    theme_minimal() +
    labs(title = "HFD+Pio vs HFD (Drug Effect) vs HFD vs CTRL (Disease Effect)",
         x = "Log2FC (HFD vs CTRL)", y = "Log2FC (HFD+Pio vs HFD)")

ggsave("corrected_DE_results/downstream/rescue_scatter_pio.png", p_pio, width = 6, height = 5, dpi = 150)

median_rescue <- data.frame(
    Intervention = c("Pio", "IP", "Cr"),
    Median_Rescue_Pct = c(median(rescue_df$RescuePct_Pio, na.rm=TRUE),
                          median(rescue_df$RescuePct_IP, na.rm=TRUE),
                          median(rescue_df$RescuePct_Cr, na.rm=TRUE))
)

p_bar <- ggplot(median_rescue, aes(x = Intervention, y = Median_Rescue_Pct, fill = Intervention)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Pio" = "#9b59b6", "IP" = "#f39c12", "Cr" = "#3498db")) +
    theme_minimal() +
    labs(title = "Median % Rescue of HFD Signature Proteins", y = "Median % Rescue", x = "")

ggsave("corrected_DE_results/downstream/rescue_barplot.png", p_bar, width = 5, height = 4, dpi = 150)


# ===========================================================================
# 3. GO / Reactome ORA
# ===========================================================================
cat("\nAnalysis 3: ORA...\n")
sig_entrez <- na.omit(sig_hfd$EntrezID)
all_entrez <- na.omit(de_hfd$EntrezID)

# GO Biological Process
go_res <- enrichGO(gene = sig_entrez, universe = all_entrez, OrgDb = org.Rn.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(go_res) && nrow(go_res) > 0) {
    write.csv(as.data.frame(go_res), "corrected_DE_results/downstream/go_enrichment.csv", row.names = FALSE)
    p_go <- dotplot(go_res, showCategory=15)
    ggsave("corrected_DE_results/downstream/go_dotplot.png", p_go, width=8, height=6, dpi=150)
}

# Reactome
reactome_res <- enrichPathway(gene = sig_entrez, universe = all_entrez, organism = "rat", pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(reactome_res) && nrow(reactome_res) > 0) {
    write.csv(as.data.frame(reactome_res), "corrected_DE_results/downstream/reactome_enrichment.csv", row.names = FALSE)
    p_reactome <- dotplot(reactome_res, showCategory=15)
    ggsave("corrected_DE_results/downstream/reactome_dotplot.png", p_reactome, width=8, height=6, dpi=150)
}


# ===========================================================================
# 4. UpSet Plot
# ===========================================================================
cat("\nAnalysis 4: UpSet Plot...\n")
list_sig <- list(
    HFD = de_hfd %>% filter(adj_p_value < 0.05) %>% pull(Protein),
    HFD_Pio = de_pio %>% filter(adj_p_value < 0.05) %>% pull(Protein),
    HFD_IP = de_ip %>% filter(adj_p_value < 0.05) %>% pull(Protein),
    HFD_Cr = de_cr %>% filter(adj_p_value < 0.05) %>% pull(Protein)
)

png("corrected_DE_results/downstream/upset_plot.png", width=800, height=600, res=150)
upset(fromList(list_sig), order.by = "freq", main.bar.color = "steelblue", sets.bar.color = c("#e74c3c", "#9b59b6", "#f39c12", "#3498db"))
dev.off()


# ===========================================================================
# 5. GSEA
# ===========================================================================
cat("\nAnalysis 5: GSEA...\n")
# Rank by logFC (or t-statistic)
ranked_df <- de_hfd %>% filter(!is.na(EntrezID)) %>% arrange(desc(t))
ranked_vec <- ranked_df$t
names(ranked_vec) <- ranked_df$EntrezID

# Run GSEA using GO
gsea_go <- gseGO(geneList = ranked_vec, ont = "BP", OrgDb = org.Rn.eg.db, pvalueCutoff = 1, pAdjustMethod = "BH")
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
    write.csv(as.data.frame(gsea_go), "corrected_DE_results/downstream/gsea_go_bp.csv", row.names = FALSE)
}


# ===========================================================================
# 6. STRING PPI Network
# ===========================================================================
cat("\nAnalysis 6: STRING PPI Network...\n")
string_db <- STRINGdb$new(version="12.0", species=10116, score_threshold=400, input_directory="")
mapped <- string_db$map(sig_hfd, "Protein", removeUnmappedRows=TRUE)

png("corrected_DE_results/downstream/string_network.png", width=800, height=800, res=150)
string_db$plot_network(mapped$STRING_id)
dev.off()

ppi_enrich <- string_db$get_enrichment(mapped$STRING_id)
write.csv(ppi_enrich, "corrected_DE_results/downstream/string_enrichment.csv", row.names = FALSE)

# clustersList <- string_db$get_clusters(mapped$STRING_id, algorithm="mcl")
# cluster_df <- data.frame(STRING_id = unlist(clustersList), 
#                          Cluster = rep(seq_along(clustersList), lengths(clustersList)))
# mapped_clusters <- left_join(mapped, cluster_df, by="STRING_id")
# write.csv(mapped_clusters, "corrected_DE_results/downstream/string_clusters.csv", row.names = FALSE)


# ===========================================================================
# 7. Intervention Heatmap
# ===========================================================================
cat("\nAnalysis 7: Intervention Heatmap...\n")
union_sig <- unique(c(list_sig$HFD, list_sig$HFD_Pio, list_sig$HFD_IP, list_sig$HFD_Cr))
union_sig <- intersect(union_sig, rownames(imputed_df))

if (length(union_sig) > 0) {
    # Compute group means
    meta <- data.frame(sample = colnames(imputed_df))
    meta$group <- case_when(
        grepl("Pio", meta$sample) ~ "HFD_Pio",
        grepl("IP", meta$sample) ~ "HFD_IP",
        grepl("Cr", meta$sample) ~ "HFD_Cr",
        grepl("AHFD", meta$sample) ~ "HFD",
        grepl("AC", meta$sample) ~ "CTRL",
        TRUE ~ "Unknown"
    )
    
    group_means <- sapply(c("CTRL", "HFD", "HFD_Pio", "HFD_IP", "HFD_Cr"), function(g) {
        cols <- meta$sample[meta$group == g]
        rowMeans(imputed_df[union_sig, cols, drop = FALSE])
    })
    
    # Z-score
    z_mat <- t(scale(t(group_means)))
    
    # Annotations
    fc_annot <- de_hfd$log2FC[match(union_sig, de_hfd$Protein)]
    row_ha <- rowAnnotation(log2FC_HFD = fc_annot, col = list(log2FC_HFD = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))))
    
    col_colors <- c("CTRL" = "#2ecc71", "HFD" = "#e74c3c", "HFD_Pio" = "#9b59b6", "HFD_IP" = "#f39c12", "HFD_Cr" = "#3498db")
    col_ha <- HeatmapAnnotation(Group = colnames(z_mat), col = list(Group = col_colors))
    
    png("corrected_DE_results/downstream/intervention_heatmap.png", width=800, height=800, res=150)
    ht <- Heatmap(z_mat, name = "Z-score", cluster_columns = FALSE, 
                  top_annotation = col_ha, right_annotation = row_ha,
                  row_title = sprintf("Union Significant Proteins (n=%d)", length(union_sig)),
                  show_row_names = (length(union_sig) <= 50))
    draw(ht)
    dev.off()
}


# ===========================================================================
# 8. Trajectory Plots
# ===========================================================================
cat("\nAnalysis 8: Trajectory Plots...\n")
core_hfd <- de_hfd %>% filter(significant_FC1 == "True") %>% pull(Protein)
if (length(core_hfd) > 0) {
    long_df <- imputed_df[core_hfd, , drop=FALSE]
    long_df$Protein <- rownames(long_df)
    long_df <- pivot_longer(long_df, cols = -Protein, names_to = "Sample", values_to = "Log2Int")
    
    long_df$Group <- case_when(
        grepl("Pio", long_df$Sample) ~ "HFD_Pio",
        grepl("IP", long_df$Sample) ~ "HFD_IP",
        grepl("Cr", long_df$Sample) ~ "HFD_Cr",
        grepl("AHFD", long_df$Sample) ~ "HFD",
        grepl("AC", long_df$Sample) ~ "CTRL"
    )
    long_df$Group <- factor(long_df$Group, levels = c("CTRL", "HFD", "HFD_IP", "HFD_Pio", "HFD_Cr"))
    
    summ_df <- long_df %>% group_by(Protein, Group) %>% summarize(Mean = mean(Log2Int), SE = sd(Log2Int)/sqrt(n()), .groups="drop")
    
    # Add symbols
    summ_df$Symbol <- map_df$Symbol[match(summ_df$Protein, map_df$Protein)]
    summ_df$Label <- paste0(summ_df$Protein, " (", summ_df$Symbol, ")")
    
    ctrl_means <- summ_df %>% filter(Group == "CTRL") %>% select(Label, Mean) %>% rename(CTRL_Mean = Mean)
    summ_df <- left_join(summ_df, ctrl_means, by = "Label")
    
    p_traj <- ggplot(summ_df, aes(x = Group, y = Mean, group = Label)) +
        geom_hline(aes(yintercept = CTRL_Mean), linetype = "dashed", color = "gray50") +
        geom_line(color = "black") +
        geom_point(aes(color = Group), size = 3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
        scale_color_manual(values = col_colors) +
        facet_wrap(~ Label, scales = "free_y", ncol = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Trajectory of Core HFD Signature Proteins (|log2FC| >= 1.0)", y = "Mean Log2 Intensity", x = "")
    
    ggsave("corrected_DE_results/downstream/trajectory_plot.png", p_traj, width = 10, height = max(4, ceiling(length(core_hfd)/3)*3), dpi = 150)
}

cat("\nDone!\n")
