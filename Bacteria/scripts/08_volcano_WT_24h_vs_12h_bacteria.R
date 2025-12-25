# Bacteria/scripts/08_volcano_WT_24h_vs_12h_bacteria.R
# Volcano plot — BACTERIA: WT 24h vs WT 12h

options(stringsAsFactors = FALSE)

res_csv <- file.path("results", "tables", "deseq_bact_WT_24h_vs_12h.csv")
stopifnot(file.exists(res_csv))

fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "volcano_bact_WT_24h_vs_12h.png")

res <- read.csv(res_csv)
res <- res[!is.na(res$padj), ]

padj_cut <- 0.05
lfc_cut  <- 1

res$class <- "Not significant"
res$class[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "Up at 24h"
res$class[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Up at 12h"

png(out_png, width = 900, height = 800)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 16,
  col = ifelse(
    res$class == "Up at 24h", "red",
    ifelse(res$class == "Up at 12h", "blue", "grey")
  ),
  xlab = "log2 Fold Change (WT 24h vs 12h, bacteria)",
  ylab = "-log10 adjusted p-value",
  main = "Volcano — Bacteria WT: 24h vs 12h"
)
abline(v = c(-lfc_cut, lfc_cut), lty = 2)
abline(h = -log10(padj_cut), lty = 2)
legend(
  "topright",
  legend = c("Up at 24h", "Up at 12h", "Not significant"),
  col = c("red", "blue", "grey"),
  pch = 16,
  bty = "n"
)
dev.off()

message("✅ Volcano saved: ", out_png)
