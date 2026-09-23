library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Rn.eg.db)

cat("Loading data...\n")
imputed_df <- read.csv("corrected_DE_results/tables/imputed_matrix.csv", row.names = 1)
de_hfd <- read.csv("corrected_DE_results/tables/de_HFD_vs_CTRL.csv")

map_df <- data.frame(Protein = de_hfd$Protein)
map_df$Symbol <- mapIds(org.Rn.eg.db, keys = map_df$Protein, keytype = "UNIPROT", column = "SYMBOL", multiVals = "first")

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
    
    summ_df$Symbol <- map_df$Symbol[match(summ_df$Protein, map_df$Protein)]
    summ_df$Label <- paste0(summ_df$Protein, " (", summ_df$Symbol, ")")
    
    ctrl_means <- summ_df %>% filter(Group == "CTRL") %>% dplyr::select(Label, Mean) %>% dplyr::rename(CTRL_Mean = Mean)
    summ_df <- left_join(summ_df, ctrl_means, by = "Label")
    
    col_colors <- c("CTRL" = "#2ecc71", "HFD" = "#e74c3c", "HFD_Pio" = "#9b59b6", "HFD_IP" = "#f39c12", "HFD_Cr" = "#3498db")
    
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
