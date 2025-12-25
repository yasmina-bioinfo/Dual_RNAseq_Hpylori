# bacteria/scripts/01_prepare_counts_bacteria.R
# Prepare count matrix for DESeq2 (BACTERIA) from XLSX

options(stringsAsFactors = FALSE)

library(readxl)

in_xlsx <- file.path("data", "raw", "GSE243405_Bacteria-genes.expression.xlsx")
out_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")

stopifnot(file.exists(in_xlsx))

df <- read_xlsx(in_xlsx)

# ---- Gene IDs: first column ----
gene_ids <- df[[1]]

# ---- Select count columns ----
count_cols <- grep("_count$", colnames(df), value = TRUE)
if (length(count_cols) == 0) {
  stop("No columns ending with '_count' found in bacteria file.")
}

counts_df <- df[, count_cols]

# ---- Ensure numeric & integer ----
counts_mat <- as.matrix(counts_df)
mode(counts_mat) <- "numeric"

if (any(is.na(counts_mat))) {
  stop("NA values detected in bacteria count matrix.")
}

counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- gene_ids

write.csv(
  data.frame(gene_id = rownames(counts_mat), counts_mat),
  out_csv,
  row.names = FALSE
)

message("✅ Bacteria count matrix saved: ", out_csv)
