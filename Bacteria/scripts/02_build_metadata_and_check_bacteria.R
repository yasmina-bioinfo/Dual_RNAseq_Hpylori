# Bacteria/scripts/02_build_metadata_and_check_bacteria.R
# Build bacteria metadata from sample names (supports WT-12h-1_count; ignores WT-1_count)

options(stringsAsFactors = FALSE)

counts_csv <- file.path("data", "processed", "bacteria_counts_matrix.csv")
out_dir    <- file.path("data", "metadata")
out_csv    <- file.path(out_dir, "samples.csv")

stopifnot(file.exists(counts_csv))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

counts <- read.csv(counts_csv, check.names = FALSE)
sample_ids <- colnames(counts)[-1]

parse_one <- function(x) {
  base <- sub("_count$", "", x)
  base <- gsub("\\.", "-", base)  # tolerate dot-converted names

  parts <- strsplit(base, "-", fixed = TRUE)[[1]]

  # Expect exactly: STRAIN - 12h/24h - REP
  if (length(parts) != 3) return(NULL)

  strain <- parts[1]
  time_s <- parts[2]              # "12h"
  rep_s  <- parts[3]

  time <- suppressWarnings(as.integer(sub("h$", "", time_s)))
  rep  <- suppressWarnings(as.integer(rep_s))

  if (!(strain %in% c("WT", "KO"))) return(NULL)
  if (is.na(time) || !(time %in% c(12, 24))) return(NULL)
  if (is.na(rep)) return(NULL)

  data.frame(
    sample_id = x,     # keep original column name (must match counts columns)
    strain = strain,
    time = time,
    replicate = rep,
    condition = "infected",
    stringsAsFactors = FALSE
  )
}

meta <- do.call(rbind, lapply(sample_ids, parse_one))

# Show what got ignored (e.g., WT-1_count, KO-2_count)
ignored <- setdiff(sample_ids, meta$sample_id)
if (length(ignored) > 0) {
  message("⚠️ Ignored columns (no time info):")
  print(ignored)
}

message("Samples per strain:"); print(table(meta$strain))
message("Samples per time:");   print(table(meta$time))
message("Strain x time:");      print(table(meta$strain, meta$time))

write.csv(meta, out_csv, row.names = FALSE)
message("✅ Bacteria metadata saved: ", out_csv)
