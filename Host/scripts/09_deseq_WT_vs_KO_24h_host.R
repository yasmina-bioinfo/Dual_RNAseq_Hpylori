# scripts/09_deseq_WT_vs_KO_24h_host.R
# Differential expression (HOST): WT vs KO at 24h

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

out_dir <- file.path("results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_WT_vs_KO_24h_host.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# Subset: infected at 24h, WT or KO
keep <- meta$condition == "infected" & meta$time == 24 & meta$strain %in% c("WT", "KO")
meta24 <- meta[keep, ]
counts24 <- counts[, meta24$sample_id]

# Safety checks
stopifnot(all(colnames(counts24) == meta24$sample_id))
stopifnot(all(meta24$time == 24))
stopifnot(length(unique(meta24$strain)) == 2)

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts24,
  colData   = meta24,
  design    = ~ strain
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Contrast: WT vs KO
res <- results(dds, contrast = c("strain", "WT", "KO"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: HOST WT vs KO at 24h")
message("Results saved to: ", res_csv)