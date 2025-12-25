# Host/scripts/99_extract_key_genes_host.R
# Auto-extract key genes from ALL DESeq2 CSVs in results/tables
# Key genes = padj < 0.05 AND |log2FC| >= 1 AND baseMean >= 50

options(stringsAsFactors = FALSE)

padj_cut <- 0.05
lfc_cut  <- 1
base_cut <- 50

in_dir  <- file.path("results", "tables")
out_dir <- file.path("results", "tables", "key_genes")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# take all csv files in results/tables except the key_genes folder outputs
all_csv <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(all_csv) == 0) {
  stop("No CSV files found in: ", in_dir, "\nRun: ls results/tables")
}

extract_one <- function(path) {
  fname <- basename(path)

  res <- read.csv(path, check.names = FALSE)

  needed <- c("baseMean", "log2FoldChange", "padj")
  if (!all(needed %in% colnames(res))) {
    message("⏭️ Skipped (not a DESeq2 results table): ", fname)
    return(invisible(NULL))
  }

  # gene_id may be present; if not, try to infer
  if (!("gene_id" %in% colnames(res))) {
    # common fallback: rownames were written as first column named X
    if ("X" %in% colnames(res)) {
      res$gene_id <- res$X
    } else {
      res$gene_id <- NA
    }
  }

  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange) & !is.na(res$baseMean), ]

  key <- subset(
    res,
    padj < padj_cut &
      abs(log2FoldChange) >= lfc_cut &
      baseMean >= base_cut
  )

  key <- key[order(key$baseMean, decreasing = TRUE), ]

  out_name <- sub("\\.csv$", "", fname)
  out_path <- file.path(out_dir, paste0(out_name, "_key_genes.csv"))

  write.csv(key, out_path, row.names = FALSE)

  message("✅ Saved: ", basename(out_path), " (n=", nrow(key), ")")
}

invisible(lapply(all_csv, extract_one))

message("Done. Key genes are in: ", out_dir)
message("Thresholds: padj<", padj_cut, ", |log2FC|>=", lfc_cut, ", baseMean>=", base_cut)