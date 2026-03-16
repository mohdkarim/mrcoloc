#!/usr/bin/env Rscript
# ============================================================================
# Cache External Data Dependencies
# ============================================================================
#
# This script fetches data from Google Sheets and BigQuery and saves local
# cache files so that the analysis pipeline can run without network access.
#
# Run this ONCE with Google Sheets and BigQuery access:
#   Rscript scripts/cache_external_data.R
#
# After running, the following cache files are created in data/:
#   - olink_complete.tsv           (Olink panel UniProt IDs)
#   - olink_complete_extended.tsv  (Olink panel with gene names)
#   - excluded_traits.tsv          (Traits to exclude from analysis)
#   - outcome_gwas_old.tsv         (Legacy outcome GWAS datasets)
#   - outcome_gwas_new.tsv         (New outcome GWAS datasets)
#   - ot_genes.tsv                 (Open Targets gene ID/symbol mapping)
#
# Note: The Open Targets triangulation cache (ot_triangulation.rds) is
# generated during the supplement pipeline run, not here, because it
# depends on analysis-specific gene/disease IDs.
#
# ============================================================================

cat("
================================================================================
  Caching External Data Dependencies
================================================================================
\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(googlesheets4)
})

gs4_deauth()

# --- Paths ---
project_root <- Sys.getenv("MRCOLOC_ROOT", getwd())
data_dir <- file.path(project_root, "data")

# ============================================================================
# 1. Google Sheets
# ============================================================================

message("[1/6] Caching olink_complete...")
olink <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1DBHpr_Y3pFja4tMju3ZDJV8Gv-oTLq6wEuS0HtYjGbQ",
  sheet = "olink_complete"
)
write_tsv(olink, file.path(data_dir, "olink_complete.tsv"))
cat("   Saved:", nrow(olink), "rows\n")

message("[2/6] Caching olink_complete_extended...")
olink_ext <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1DBHpr_Y3pFja4tMju3ZDJV8Gv-oTLq6wEuS0HtYjGbQ",
  sheet = "olink_complete_with_more_fields"
)
write_tsv(olink_ext, file.path(data_dir, "olink_complete_extended.tsv"))
cat("   Saved:", nrow(olink_ext), "rows\n")

message("[3/6] Caching excluded_traits...")
toremove <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1TDz8oRI5H-DMHOTs0bZgm2dcisYeuvdSCybn4NcHESw",
  sheet = "v4"
)
write_tsv(toremove, file.path(data_dir, "excluded_traits.tsv"))
cat("   Saved:", nrow(toremove), "rows (", sum(toremove$to_remove == "Y", na.rm = TRUE), "marked for removal)\n")

message("[4/6] Caching outcome_gwas_old...")
old_st2 <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1xQ9ojspUea7stOTXSO_vZbNPPg_BU_y2pyM2SigV7UI",
  sheet = "ST1 - Outcome_GWAS", skip = 1
)
write_tsv(old_st2, file.path(data_dir, "outcome_gwas_old.tsv"))
cat("   Saved:", nrow(old_st2), "rows\n")

message("[5/6] Caching outcome_gwas_new...")
new_st2 <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1aRTkoyPUrV9jOrUQAhvXl-38QYS4CE0G0aHlhry9LKg"
)
write_tsv(new_st2, file.path(data_dir, "outcome_gwas_new.tsv"))
cat("   Saved:", nrow(new_st2), "rows\n")

# ============================================================================
# 2. BigQuery: Open Targets gene table
# ============================================================================

message("[6/6] Caching Open Targets gene table from BigQuery...")

if (requireNamespace("bigrquery", quietly = TRUE)) {
  library(bigrquery)

  tryCatch({
    bq_project <- Sys.getenv("BQ_PROJECT_ID", "")
    bq_email <- Sys.getenv("BQ_EMAIL", "")
    if (bq_project == "" || bq_email == "") {
      stop("Set BQ_PROJECT_ID and BQ_EMAIL environment variables.\n",
           "Example: BQ_PROJECT_ID=your-gcp-project BQ_EMAIL=you@example.com Rscript scripts/cache_external_data.R")
    }
    bq_auth(email = bq_email)
    sql <- 'SELECT id, approvedSymbol FROM `open-targets-prod.platform.target`'
    tb <- bq_project_query(bq_project, sql)
    genes <- bq_table_download(tb)
    write_tsv(genes, file.path(data_dir, "ot_genes.tsv"))
    cat("   Saved:", nrow(genes), "genes\n")
  }, error = function(e) {
    warning("BigQuery query failed: ", e$message,
            "\n   You can still run the pipeline if data/ot_genes.tsv already exists.")
  })
} else {
  warning("bigrquery not installed. Skipping BigQuery cache.",
          "\n   You can still run the pipeline if data/ot_genes.tsv already exists.")
}

cat("
================================================================================
  Caching complete. Files saved to: ", data_dir, "
================================================================================
\n")
