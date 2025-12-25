# scripts/08_deseq_KO_12h_vs_24h_host.R
# Differential expression (HOST): KO 24h vs KO 12h

options(stringsAsFactors = FALSE)

# ---- Inputs ----
counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

# ---- Outputs ----
out_dir <- file.path("results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
res_csv <- file.path(out_dir, "deseq_KO_24h_vs_12h_host.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# ---- Subset: infected, KO only, time 12h or 24h ----
keep <- meta$condition == "infected" & meta$strain == "KO" & meta$time %in% c(12, 24)
metaKO <- meta[keep, ]
countsKO <- counts[, metaKO$sample_id]

# Safety checks
stopifnot(all(colnames(countsKO) == metaKO$sample_id))
stopifnot(length(unique(metaKO$time)) == 2)

suppressMessages(library(DESeq2))

# Set time as factor to control contrast direction
metaKO$time <- factor(metaKO$time, levels = c(12, 24))

dds <- DESeqDataSetFromMatrix(
  countData = countsKO,
  colData   = metaKO,
  design    = ~ time
)

# Minimal filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# Contrast: 24h vs 12h
res <- results(dds, contrast = c("time", "24", "12"))

# Save
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, res_csv, row.names = FALSE)

message("✅ DESeq2 completed: HOST KO 24h vs 12h")
message("Results saved to: ", res_csv)