# Open Targets gene ID <-> symbol mapping
# Loads from local cache (data/ot_genes.tsv) if available.
# Falls back to BigQuery if cache is missing and bigrquery is installed.
# To regenerate cache: Rscript scripts/cache_external_data.R

library(dplyr)

cache_path <- file.path(
  Sys.getenv("MRCOLOC_ROOT", getwd()),
  "data", "ot_genes.tsv"
)

if (file.exists(cache_path)) {
  genes <- readr::read_tsv(cache_path, show_col_types = FALSE)
  message("   Loaded ", nrow(genes), " genes from cache: ", cache_path)
} else {
  message("   Cache not found, querying BigQuery...")
  if (!requireNamespace("bigrquery", quietly = TRUE)) {
    stop("Cache file ", cache_path, " not found and bigrquery is not installed.\n",
         "Run: Rscript scripts/cache_external_data.R (requires BigQuery access)")
  }
  library(bigrquery)
  bq_project <- Sys.getenv("BQ_PROJECT_ID", "")
  bq_email <- Sys.getenv("BQ_EMAIL", "")
  if (bq_project == "" || bq_email == "") {
    stop("Set BQ_PROJECT_ID and BQ_EMAIL environment variables for BigQuery access.\n",
         "Or run: Rscript scripts/cache_external_data.R to generate the cache.")
  }
  bq_auth(email = bq_email)
  sql <- 'SELECT id, approvedSymbol FROM `open-targets-prod.platform.target`'
  tb <- bq_project_query(bq_project, sql)
  genes <- bq_table_download(tb)
  # Save cache for future runs
  readr::write_tsv(genes, cache_path)
  message("   Queried ", nrow(genes), " genes and saved cache to: ", cache_path)
}
