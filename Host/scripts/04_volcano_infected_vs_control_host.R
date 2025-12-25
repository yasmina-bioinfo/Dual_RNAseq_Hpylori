# scripts/04_volcano_infected_vs_control.R
# Volcano plot for DESeq2 results (infected vs control)

options(stringsAsFactors = FALSE)

# ---- Input ----
res_csv <- file.path("results", "deseq_infected_vs_control_results.csv")
stopifnot(file.exists(res_csv))

# ---- Output ----
fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
out_png <- file.path(fig_dir, "volcano_infected_vs_control.png")

res <- read.csv(res_csv)

# Remove NA adjusted p-values
res <- res[!is.na(res$padj), ]

# Thresholds
padj_cut <- 0.05
lfc_cut  <- 1

# Classify points
res$class <- "Not significant"
res$class[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "Up in infected"
res$class[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Down in infected"

# ---- Plot (base R) ----
png(out_png, width = 900, height = 800)
plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 16,
  col = ifelse(
    res$class == "Up in infected", "red",
    ifelse(res$class == "Down in infected", "blue", "grey")
  ),
  xlab = "log2 Fold Change (infected vs control)",
  ylab = "-log10 adjusted p-value",
  main = "Volcano plot — infected vs control"
)
abline(v = c(-lfc_cut, lfc_cut), lty = 2)
abline(h = -log10(padj_cut), lty = 2)
legend(
  "topright",
  legend = c("Up in infected", "Down in infected", "Not significant"),
  col = c("red", "blue", "grey"),
  pch = 16,
  bty = "n"
)
dev.off()

message("✅ Volcano plot saved: ", out_png)
