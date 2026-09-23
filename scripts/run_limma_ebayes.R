library(limma)

args = commandArgs(trailingOnly=TRUE)
imputed_file = args[1]
design_file = args[2]
out_dir = args[3]

# Read data
df = read.csv(imputed_file, row.names=1)
design_df = read.csv(design_file, row.names=1)

# Ensure design matrix row order matches expression matrix column order
design = as.matrix(design_df[colnames(df), ])

cat("Running limma::lmFit...\n")
fit = lmFit(df, design)

cat("Building contrasts...\n")
cont.matrix = makeContrasts(
  HFD_vs_CTRL = group_HFD - group_CTRL,
  HFD_Pio_vs_HFD = group_HFD_Pio - group_HFD,
  HFD_IP_vs_HFD = group_HFD_IP - group_HFD,
  HFD_Cr_vs_HFD = group_HFD_Cr - group_HFD,
  levels=design
)

cat("Running limma::contrasts.fit and eBayes...\n")
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

# Save results for each contrast
for (coef_name in colnames(cont.matrix)) {
    res = topTable(fit2, coef=coef_name, number=Inf, sort.by="none")
    out_path = file.path(out_dir, paste0("limma_res_", coef_name, ".csv"))
    write.csv(res, out_path)
    cat("Saved", out_path, "\n")
}
