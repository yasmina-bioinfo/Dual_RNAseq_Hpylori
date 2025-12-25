# scripts/01_prepare_counts_host.R
# Prepare count matrix for DESeq2 (HOST)
# Inputs: data/processed/host_expression.csv
# Outputs: data/processed/host_counts_matrix.csv
#          data/processed/prepare_counts_report.txt

options(stringsAsFactors = FALSE)

in_csv  <- file.path("data", "processed", "host_expression.csv")
out_csv <- file.path("data", "processed", "host_counts_matrix.csv")
rep_txt <- file.path("data", "processed", "prepare_counts_report.txt")

stopifnot(file.exists(in_csv))

df <- read.csv(in_csv, check.names = FALSE)

# ---- Identify gene ID column (assume first column with ENSG pattern) ----
gene_col_idx <- which(grepl("^ENSG", as.character(df[[1]])))
if (length(gene_col_idx) == 0) {
    stop("No ENSG gene IDs detected in the first column.")
}

gene_ids <- df[[1]]

# ---- Select count columns only ----
count_cols <- grep("_count$", colnames(df), value = TRUE)
if (length(count_cols) == 0) {
    stop("No columns ending with '_count' were found.")
}

counts_df <- df[, count_cols, drop = FALSE]

# ---- Coerce to numeric, then round to integers (DESeq2 requirement) ----
counts_mat <- as.matrix(counts_df)
mode(counts_mat) <- "numeric"

if (any(is.na(counts_mat))) {
  stop("NA values detected in count matrix.")
}

# Check if there are non-integers
diff_from_int <- max(abs(counts_mat - round(counts_mat)))
if (diff_from_int > 1e-6) {
  message("⚠️ Non-integer counts detected (max diff from integer = ", diff_from_int, "). Rounding counts for DESeq2.")
}

counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

# ---- Set rownames as gene IDs ----
rownames(counts_mat) <- gene_ids

# ---- Write outputs ----
write.csv(
    data.frame(gene_id = rownames(counts_mat), counts_mat),
    out_csv,
    row.names = FALSE
)

# ---- Minimal report ----
con <- file(rep_txt, open = "wt")
on.exit(close(con), add = TRUE)

writeLines("=== Prepare counts report (HOST) ===", con)
writeLines(paste("Date:", Sys.time()), con)
writeLines(paste("Genes:", nrow(counts_mat)), con)
writeLines(paste("Samples:", ncol(counts_mat)), con)
writeLines("Sample columns:", con)
writeLines(paste(colnames(counts_mat), collapse = ", "), con)

message("✅ Host count matrix ready:")
message(" - ", out_csv)
message("📄 Report: ", rep_txt)
