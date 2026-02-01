#!/usr/bin/env Rscript
# ============================================================================
# Setup Script - Prepare Environment for mrcoloc Analysis
# ============================================================================
#
# This script:
#   1. Checks/installs required R packages
#   2. Downloads data files from Zenodo
#   3. Verifies all files are present
#   4. Runs a quick test to ensure everything works
#
# Usage:
#   Rscript scripts/setup.R
#
# ============================================================================

cat("
================================================================================
  mrcoloc Analysis Setup
================================================================================

This script will prepare your environment to reproduce the analyses from:

  'Impact of proteogenomic evidence on clinical success'
  Karim et al., Nature Genetics (2025)

================================================================================
\n")

# ============================================================================
# STEP 1: CHECK R VERSION
# ============================================================================

cat("[1/4] Checking R version...\n")

r_version <- getRversion()
if (r_version < "4.0.0") {
  stop("R version 4.0.0 or higher required. You have: ", r_version)
}
cat("  R version:", as.character(r_version), "✓\n\n")

# ============================================================================
# STEP 2: INSTALL REQUIRED PACKAGES
# ============================================================================

cat("[2/4] Checking required packages...\n")

# CRAN packages
cran_packages <- c(
  "tidyverse",
  "data.table", 
  "openxlsx",
  "DescTools",
  "ggplot2",
  "UpSetR",
  "DiagrammeR",
  "DiagrammeRsvg",
  "googlesheets4",
  "htmlwidgets",
  "Rmpfr",
  "stringr"
)

# Bioconductor packages
bioc_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "biomaRt"
)

# Check and install CRAN packages
missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("  Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Check and install Bioconductor packages
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("  Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# Verify all packages installed
all_packages <- c(cran_packages, bioc_packages)
still_missing <- all_packages[!sapply(all_packages, requireNamespace, quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}
cat("  All", length(all_packages), "packages available ✓\n\n")

# ============================================================================
# STEP 3: DOWNLOAD DATA
# ============================================================================

cat("[3/4] Downloading data files...\n\n")

# Source the download script
download_script <- "scripts/download_data.R"
if (!file.exists(download_script)) {
  download_script <- file.path(dirname(getwd()), "scripts/download_data.R")
}

if (file.exists(download_script)) {
  source(download_script)
} else {
  cat("  Warning: download_data.R not found. Skipping data download.\n")
  cat("  Run 'Rscript scripts/download_data.R' manually.\n\n")
}

# ============================================================================
# STEP 4: VERIFY SETUP
# ============================================================================

cat("\n[4/4] Verifying setup...\n")

# Check data_raw directory
data_raw <- "data_raw"
if (!dir.exists(data_raw)) {
  data_raw <- file.path(dirname(getwd()), "data_raw")
}

required_files <- c(
  "pqtl_mrcoloc_2025.rds",
  "chembl.rds",
  "ttpairs_tested.rds"
)

if (dir.exists(data_raw)) {
  present <- sum(file.exists(file.path(data_raw, required_files)))
  cat("  Core data files:", present, "/", length(required_files), "\n")
} else {
  cat("  Warning: data_raw directory not found\n")
}

# Check data directory (reference files)
data_dir <- "data"
if (!dir.exists(data_dir)) {
  data_dir <- file.path(dirname(getwd()), "data")
}

ref_files <- c("areas.tsv", "indic.tsv", "indic_topl_match.tsv")
if (dir.exists(data_dir)) {
  present <- sum(file.exists(file.path(data_dir, ref_files)))
  cat("  Reference files:", present, "/", length(ref_files), "\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  Setup Complete!\n")
cat("================================================================================\n\n")

cat("Next steps:\n")
cat("  1. If any Minikel et al. files are missing, download them manually\n")
cat("     (see data/README.md for instructions)\n\n")
cat("  2. Run the analyses:\n")
cat("     source('scripts/generate_mrcoloc_supplement.R')  # Supplementary tables\n")
cat("     source('scripts/mrcoloc_paper_2025_main_figures.R')  # Main figures\n")
cat("     source('scripts/mrcoloc_paper_2025_supp_figures.R')  # Supp figures\n\n")

cat("================================================================================\n")