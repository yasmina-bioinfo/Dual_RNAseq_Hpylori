# scripts/12_deseq_WT_24h_vs_12h_host.R
# Differential expression (HOST): WT 24h vs WT 12h

options(stringsAsFactors = FALSE)

# ---- Inputs ----
counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

# ---- Outputs ----
out_dir <- file.path("results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_WT_24h_vs_12h_host.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# ---- Subset: infected, WT only, 12h vs 24h ----
keep <- meta$condition == "infected" & meta$strain == "WT" & meta$time %in% c(12, 24)
metaWT <- meta[keep, ]
countsWT <- counts[, metaWT$sample_id]

stopifnot(all(colnames(countsWT) == metaWT$sample_id))
stopifnot(length(unique(metaWT$time)) == 2)

suppressMessages(library(DESeq2))

metaWT$time <- factor(metaWT$time, levels = c(12, 24))

dds <- DESeqDataSetFromMatrix(
  countData = countsWT,
  colData   = metaWT,
  design    = ~ time
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Contrast: 24h vs 12h
res <- results(dds, contrast = c("time", "24", "12"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: HOST WT 24h vs 12h")