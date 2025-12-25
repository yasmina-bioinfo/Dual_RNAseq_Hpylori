# Bacteria/scripts/09_deseq_KO_24h_vs_12h_bacteria.R
# DESeq2 (BACTERIA): KO 24h vs KO 12h

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

out_dir <- file.path("results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_bact_KO_24h_vs_12h.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# ---- Subset: KO only ----
keep <- meta$strain == "KO" & meta$time %in% c(12, 24)
meta_ko <- meta[keep, ]
meta_ko$time <- factor(meta_ko$time, levels = c(12, 24))

counts_ko <- counts[, meta_ko$sample_id, drop = FALSE]

# Checks
stopifnot(all(colnames(counts_ko) == meta_ko$sample_id))
stopifnot(length(unique(meta_ko$time)) == 2)

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts_ko,
  colData   = meta_ko,
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

message("✅ DESeq2 completed: Bacteria KO 24h vs 12h")
message("Results saved to: ", res_csv)
