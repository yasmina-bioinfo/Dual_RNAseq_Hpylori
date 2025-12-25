# scripts/11_volcano_WT_vs_KO_24h_host.R
# Volcano plot — HOST: WT vs KO at 24h

options(stringsAsFactors = FALSE)

res_csv <- file.path("results", "tables", "deseq_WT_vs_KO_24h_host.csv")

fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "volcano_WT_vs_KO_24h_host.png")

res <- read.csv(res_csv)
res <- res[!is.na(res$padj), ]

padj_cut <- 0.05
lfc_cut  <- 1

res$class <- "Not significant"
res$class[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "Up in WT"
res$class[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Up in KO"

png(out_png, width = 900, height = 800)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 16,
  col = ifelse(
    res$class == "Up in WT", "red",
    ifelse(res$class == "Up in KO", "blue", "grey")
  ),
  xlab = "log2 Fold Change (WT vs KO, host, 24h)",
  ylab = "-log10 adjusted p-value",
  main = "Volcano — Host response: WT vs KO at 24h"
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

message("✅ Volcano saved: ", out_png)
