# scripts/00_import_xlsx_to_csv.R
# Import GEO processed .xlsx -> export .csv (reproducible)
# Requirements: base R + readxl

options(stringsAsFactors = FALSE)

# ---- 0) Safety checks ----
if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
    "Package 'readxl' is not installed.\n",
    "If you already have it, restart R and try again.\n",
    "If not installed, install it once with:\n",
    "install.packages('readxl')\n"
    )
}

raw_dir <- file.path("data", "raw")
out_dir <- file.path("data", "processed")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

host_xlsx <- file.path(raw_dir, "GSE243405_Host_cell-genes.expression.xlsx")
bac_xlsx  <- file.path(raw_dir, "GSE243405_Bacteria-genes.expression.xlsx")

stopifnot(file.exists(host_xlsx))
stopifnot(file.exists(bac_xlsx))

# ---- 1) Helper: read first sheet by default (safe) ----
read_first_sheet <- function(path) {
    sheets <- readxl::excel_sheets(path)
    if (length(sheets) < 1) stop("No sheets found in: ", path)
    message("File: ", basename(path), " | Sheets: ", paste(sheets, collapse = ", "))
    readxl::read_excel(path, sheet = sheets[1])
}

# ---- 2) Read ----
message("\n[1/3] Reading host xlsx...")
host_df <- read_first_sheet(host_xlsx)
message("[2/3] Reading bacteria xlsx...")
bac_df  <- read_first_sheet(bac_xlsx)

# ---- 3) Basic reporting ----
report_path <- file.path(out_dir, "import_report.txt")
con <- file(report_path, open = "wt", encoding = "UTF-8")
on.exit(close(con), add = TRUE)

writeLines("=== Import report: GSE243405 processed tables ===", con)
writeLines(paste("Date:", as.character(Sys.time())), con)

summarize_df <- function(df, name) {
    writeLines("\n---", con)
    writeLines(paste("Table:", name), con)
    writeLines(paste("Rows:", nrow(df), " | Cols:", ncol(df)), con)
    writeLines("First 12 column names:", con)
    writeLines(paste(utils::head(colnames(df), 12), collapse = " | "), con)
    writeLines(paste("NA count (total):", sum(is.na(df))), con)
}

summarize_df(host_df, "HOST")
summarize_df(bac_df,  "BACTERIA")

# ---- 4) Export to CSV ----
host_csv <- file.path(out_dir, "host_expression.csv")
bac_csv  <- file.path(out_dir, "bacteria_expression.csv")

# Ensure data.frame for write.csv
host_df <- as.data.frame(host_df)
bac_df  <- as.data.frame(bac_df)

utils::write.csv(host_df, host_csv, row.names = FALSE)
utils::write.csv(bac_df,  bac_csv,  row.names = FALSE)

message("\n✅ Export done:")
message(" - ", host_csv)
message(" - ", bac_csv)
message("📄 Report: ", report_path)
