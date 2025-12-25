# scripts/07_volcano_WT_vs_KO_12h_host.R
# Volcano plot for HOST DESeq2 results: WT vs KO at 12h

options(stringsAsFactors = FALSE)

# ---- Input ----
res_csv <- file.path("results", "deseq_WT_vs_KO_12h_host.csv")
stopifnot(file.exists(res_csv))

# ---- Output ----
fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "volcano_WT_vs_KO_12h_host.png")

res <- read.csv(res_csv)

# Remove NA adjusted p-values
res <- res[!is.na(res$padj), ]

# Thresholds (keep same as before for comparability)
padj_cut <- 0.05
lfc_cut  <- 1

# Classify points
res$class <- "Not significant"
res$class[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "Up in WT"
res$class[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Up in KO"

# ---- Plot (base R) ----
png(out_png, width = 900, height = 800)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 16,
  col = ifelse(
    res$class == "Up in WT", "red",
    ifelse(res$class == "Up in KO", "blue", "grey")
  ),
  xlab = "log2 Fold Change (WT vs KO, host)",
  ylab = "-log10 adjusted p-value",
  main = "Volcano plot — Host response: WT vs KO at 12h"
)
abline(v = c(-lfc_cut, lfc_cut), lty = 2)
abline(h = -log10(padj_cut), lty = 2)
legend(
  "topright",
  legend = c("Up in WT", "Up in KO", "Not significant"),
  col = c("red", "blue", "grey"),
  pch = 16,
  bty = "n"
)
dev.off()

message("✅ Volcano plot saved: ", out_png)