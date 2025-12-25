# Bacteria/scripts/11_pca_bacteria_global.R
# Global PCA for bacteria (WT/KO x time) using DESeq2 VST

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# Ensure ordering
counts <- counts[, meta$sample_id, drop = FALSE]
stopifnot(all(colnames(counts) == meta$sample_id))

# Factors for plotting
meta$strain <- factor(meta$strain, levels = c("WT", "KO"))
meta$time   <- factor(meta$time, levels = c(12, 24))

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ strain + time
)

# Minimal filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Transform for PCA (no DE needed for PCA)
vsd <- vst(dds, blind = TRUE)

# PCA
mat <- t(assay(vsd))
pca <- prcomp(mat)

pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# Output figure
fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "pca_bacteria_global.png")

png(out_png, width = 900, height = 800)
plot(
  pca$x[,1], pca$x[,2],
  pch = 16,
  col = ifelse(meta$strain == "WT", "red", "blue"),
  xlab = paste0("PC1 (", round(pct[1], 1), "%)"),
  ylab = paste0("PC2 (", round(pct[2], 1), "%)"),
  main = "PCA — Bacteria RNA-seq (WT/KO x time)"
)
text(pca$x[,1], pca$x[,2], labels = paste(meta$strain, meta$time, sep = "_"), pos = 3, cex = 0.7)
legend("topright", legend = c("WT", "KO"), col = c("red", "blue"), pch = 16, bty = "n")
dev.off()

message("✅ PCA saved: ", out_png)
