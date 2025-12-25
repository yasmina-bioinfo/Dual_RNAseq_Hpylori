# Bacteria/scripts/07_deseq_WT_24h_vs_12h_bacteria.R
# DESeq2 (BACTERIA): WT 24h vs WT 12h

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

out_dir <- file.path("results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_bact_WT_24h_vs_12h.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# ---- Subset: WT only ----
keep <- meta$strain == "WT" & meta$time %in% c(12, 24)
meta_wt <- meta[keep, ]
meta_wt$time <- factor(meta_wt$time, levels = c(12, 24))
counts_wt <- counts[, meta_wt$sample_id, drop = FALSE]

# Checks
stopifnot(all(colnames(counts_wt) == meta_wt$sample_id))
stopifnot(length(unique(meta_wt$time)) == 2)

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts_wt,
  colData   = meta_wt,
  design    = ~ time
)

# Minimal filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Contrast: 24h vs 12h
res <- results(dds, contrast = c("time", "24", "12"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: Bacteria WT 24h vs 12h")
message("Results saved to: ", res_csv)
