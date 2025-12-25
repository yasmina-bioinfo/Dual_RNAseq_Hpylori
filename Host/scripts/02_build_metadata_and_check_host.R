# scripts/02_build_metadata_and_check.R
# Build samples.csv from count column names and validate against count matrix

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "host_counts_matrix.csv")
meta_out   <- file.path("data", "metadata", "samples.csv")
rep_out    <- file.path("data", "metadata", "metadata_check_report.txt")

stopifnot(file.exists(counts_csv))
if (!dir.exists(file.path("data", "metadata"))) dir.create(file.path("data", "metadata"), recursive = TRUE)

# IMPORTANT: keep original column names (do not convert '-' to '.')
counts_df <- read.csv(counts_csv, check.names = FALSE)

if (!("gene_id" %in% colnames(counts_df))) {
  stop("Column 'gene_id' not found in host_counts_matrix.csv")
}

sample_cols <- setdiff(colnames(counts_df), "gene_id")
if (length(sample_cols) == 0) stop("No sample columns found in count matrix.")

parse_sample <- function(s) {
  # s can be like WT-12h-1_count OR WT.12h.1_count (if names were converted somewhere)
  base_raw <- sub("_count$", "", s)
  base <- gsub("\\.", "-", base_raw)  # normalize dots to hyphens for parsing only

  # Controls: Con-1, Con-2, Con-3 ...
  if (grepl("^Con-", base)) {
    rep_str <- sub("^Con-([0-9]+)$", "\\1", base)
    return(data.frame(
      sample_id = s,              # keep EXACT column name
      condition = "control",
      strain    = NA,
      time      = NA,
      replicate = as.integer(rep_str),
      stringsAsFactors = FALSE
    ))
  }

  # Infected: WT-12h-1, KO-24h-3 ...
  m <- regexec("^(WT|KO)-([0-9]+)h-([0-9]+)$", base)
  reg <- regmatches(base, m)[[1]]
  if (length(reg) == 0) stop("Unrecognized sample name pattern: ", s)

  data.frame(
    sample_id = s,               # keep EXACT column name
    condition = "infected",
    strain    = reg[2],
    time      = as.integer(reg[3]),
    replicate = as.integer(reg[4]),
    stringsAsFactors = FALSE
  )
}

# ---- Build metadata ----
meta_list <- lapply(sample_cols, parse_sample)
meta <- do.call(rbind, meta_list)

# ---- Safety checks ----
if (any(duplicated(meta$sample_id))) stop("Duplicate sample_id detected in metadata.")

missing_in_counts <- setdiff(meta$sample_id, sample_cols)
if (length(missing_in_counts) > 0) {
  stop("These metadata sample_id are missing in counts: ", paste(missing_in_counts, collapse = ", "))
}

missing_in_meta <- setdiff(sample_cols, meta$sample_id)
if (length(missing_in_meta) > 0) {
  stop("These count columns are missing in metadata: ", paste(missing_in_meta, collapse = ", "))
}

# ---- Write outputs ----
write.csv(meta, meta_out, row.names = FALSE)

con <- file(rep_out, open = "wt")
on.exit(close(con), add = TRUE)
writeLines("=== Metadata build & check report ===", con)
writeLines(paste("Date:", Sys.time()), con)
writeLines(paste("Samples:", nrow(meta)), con)
writeLines("\nCounts by condition:", con)
writeLines(capture.output(print(table(meta$condition))), con)
writeLines("\nCounts by strain (including NA for controls):", con)
writeLines(capture.output(print(table(meta$strain, useNA = "ifany"))), con)
writeLines("\nCounts by time (including NA for controls):", con)
writeLines(capture.output(print(table(meta$time, useNA = "ifany"))), con)

message("✅ Metadata created and validated:")
message(" - ", meta_out)
message("📄 Report: ", rep_out)
