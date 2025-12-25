res <- read.csv("results/tables/deseq_bact_WT_24h_vs_12h.csv")

genes_condition_dep <- subset(
  res,
  padj < 0.05 &
  abs(log2FoldChange) >= 1 &
  baseMean >= 50
)

genes_condition_dep <- genes_condition_dep[
  order(genes_condition_dep$baseMean, decreasing = TRUE),
]

# create a file 
out_dir <- file.path("results", "tables", "key_genes")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# save the file
write.csv(
  genes_condition_dep,
  file.path(out_dir, "bacteria_WT_24h_vs_12h_key_genes.csv"),
  row.names = FALSE
)