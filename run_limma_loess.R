library(limma)

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile = args[2]

# Read data, assuming first column is the index (Protein Accession)
df = read.csv(infile, row.names=1)
mat = as.matrix(df)

# Run limma's cyclic loess normalization using method="fast" to avoid NA propagation
cat("Running limma::normalizeCyclicLoess(method='fast')...\n")
norm_mat = normalizeCyclicLoess(mat, method="fast")

# Save to CSV
write.csv(norm_mat, outfile)
cat("Saved to", outfile, "\n")
