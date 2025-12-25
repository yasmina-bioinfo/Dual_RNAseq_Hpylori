# scripts/05_pca_host.R
# PCA on host RNA-seq data (DESeq2 VST)

options(stringsAsFactors = FALSE)

# ---- Inputs ----
counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_csv   <- file.path("data", "metadata", "samples.csv")

# ---- Output ----
fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "PCA_host_infected_vs_control.png")

stopifnot(file.exists(counts_csv))
stopifnot(file.exists(meta_csv))

counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_csv, check.names = FALSE)

# Ensure same order
counts <- counts[, meta$sample_id]
stopifnot(all(colnames(counts) == meta$sample_id))

suppressMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

# Variance stabilizing transformation (for PCA only)
vsd <- vst(dds, blind = TRUE)

# PCA
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

# Colors by condition
cols <- ifelse(meta$condition == "infected", "red", "blue")

png(out_png, width = 900, height = 800)
plot(
  pca$x[,1], pca$x[,2],
  col = cols,
  pch = 16,
  xlab = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
  ylab = paste0("PC2 (", round(var_expl[2]*100, 1), "%)"),
  main = "PCA — Host RNA-seq (infected vs control)"
)
legend(
  "topright",
  legend = c("Infected", "Control"),
  col = c("red", "blue"),
  pch = 16,
  bty = "n"
)
dev.off()

message("✅ PCA saved: ", out_png)
