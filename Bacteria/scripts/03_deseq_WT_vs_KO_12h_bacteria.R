# Bacteria/scripts/03_deseq_WT_vs_KO_12h_bacteria.R
# DESeq2 (BACTERIA): WT vs KO at 12h

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

out_dir <- file.path("results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_bact_WT_vs_KO_12h.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# ---- Subset: 12h only ----
keep <- meta$time == 12 & meta$strain %in% c("WT", "KO")
meta12 <- meta[keep, ]
counts12 <- counts[, meta12$sample_id, drop = FALSE]

# Checks
stopifnot(all(colnames(counts12) == meta12$sample_id))
stopifnot(length(unique(meta12$strain)) == 2)

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts12,
  colData   = meta12,
  design    = ~ strain
)

# Minimal filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Contrast: WT vs KO
res <- results(dds, contrast = c("strain", "WT", "KO"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: Bacteria WT vs KO at 12h")
message("Results saved to: ", res_csv)
