options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("UpSetR", "ggplot2", "ggrepel")
for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_pkgs <- c(
    "clusterProfiler",
    "org.Rn.eg.db",
    "enrichplot",
    "DOSE",
    "ReactomePA",
    "pathview",
    "fgsea",
    "msigdbr",
    "ComplexHeatmap",
    "STRINGdb"
)

for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
}
