#!/usr/bin/env Rscript
# ============================================================================
# Download Pan-UKB Phenotype Manifest
# ============================================================================
#
# This script downloads and processes the Pan-UKB phenotype manifest from
# a public Google Sheet. The output (panukb.rds) is used to map UKB trait
# codes to EFO terms.
#
# Source: Pan-UKB phenotype manifest
# https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8
#
# Note: This file is small (~92KB) and is included in the data/ folder.
# This script is provided for transparency and reproducibility.
#
# Usage:
#   Rscript scripts/download_panukb.R
#
# ============================================================================

cat("Downloading Pan-UKB phenotype manifest...\n\n")

# Check for required packages
if (!requireNamespace("googlesheets4", quietly = TRUE)) {
  install.packages("googlesheets4", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}

library(googlesheets4)
library(dplyr)

# Deauthorize to access public sheet without authentication
gs4_deauth()

# Download from public Google Sheet
sheet_url <- "https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8"

cat("Reading from Google Sheet...\n")
panukb <- read_sheet(sheet_url, sheet = "phenotype_manifest") %>%
  select(coding_description, filename, trait_efo_terms, trait_efos)

cat("  Rows:", nrow(panukb), "\n")
cat("  Columns:", paste(names(panukb), collapse = ", "), "\n\n")

# Determine output path
if (file.exists("scripts/download_panukb.R")) {
  output_path <- "data/panukb.rds"
} else if (file.exists("download_panukb.R")) {
  output_path <- "../data/panukb.rds"
} else {
  output_path <- "panukb.rds"
}

# Save
saveRDS(panukb, output_path)
cat("Saved to:", output_path, "\n")
cat("File size:", round(file.size(output_path) / 1024, 1), "KB\n")