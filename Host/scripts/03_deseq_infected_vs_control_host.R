# scripts/03_deseq_infected_vs_control.R
# Differential expression: infected vs control (HOST)

options(stringsAsFactors = FALSE)

# ---- Inputs ----
counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

# ---- Outputs ----
out_dir <- file.path("results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res_csv <- file.path(out_dir, "deseq_infected_vs_control_results.csv")

# ---- Load data ----
stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts_df <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_csv, check.names = FALSE)

# Ensure matching order
counts_df <- counts_df[, meta$sample_id]
stopifnot(all(colnames(counts_df) == meta$sample_id))

# ---- DESeq2 ----
suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts_df,
  colData   = meta,
  design    = ~ condition
)

# Filter low counts (standard, minimal)
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "infected", "control"))

# ---- Save results ----
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: infected vs control")
message("Results saved to: ", res_csv)
