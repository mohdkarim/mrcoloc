#!/usr/bin/env Rscript
# ============================================================================
# Master Setup Script
# ============================================================================
#
# This script sets up the complete analysis environment for:
#
#   "Impact of proteogenomic evidence on clinical success"
#   Karim et al., Manuscript in preparation (2025)
#
# Steps:
#   1. Install required R packages
#   2. Download data from Zenodo and GitHub
#   3. Create derived data files
#
# Usage:
#   Rscript scripts/setup.R
#
# Or from R:
#   source("scripts/setup.R")
#
# ============================================================================

cat("
================================================================================
  mrcoloc Analysis Setup
================================================================================

  Paper: Impact of proteogenomic evidence on clinical success
  Authors: Karim et al., Manuscript in preparation (2025)
  
  GitHub: https://github.com/mohdkarim/mrcoloc
  Zenodo: https://doi.org/10.5281/zenodo.18451758

================================================================================
\n")

# ============================================================================
# DETERMINE PROJECT ROOT
# ============================================================================

if (Sys.getenv("MRCOLOC_ROOT") != "") {
  project_root <- Sys.getenv("MRCOLOC_ROOT")
} else {
  if (file.exists("scripts/setup.R")) {
    project_root <- getwd()
  } else if (file.exists("setup.R")) {
    project_root <- dirname(getwd())
  } else {
    stop("Cannot determine project root. Run from project directory or set MRCOLOC_ROOT.")
  }
}

cat("Project root:", project_root, "\n\n")

# ============================================================================
# STEP 1: INSTALL REQUIRED PACKAGES
# ============================================================================

cat("--- Step 1: Checking R packages ---\n\n")

required_packages <- c(
  # Data manipulation
  "tidyverse",
  "data.table",
  "fuzzyjoin",
  "stringdist",
  
  # Statistics
  "DescTools",
  
  # Visualization
  "ggplot2",
  "patchwork",
  "UpSetR",
  "ckbplotr",
  
  # Bioconductor
  "AnnotationDbi",
  "org.Hs.eg.db",
  
  # Google Sheets (for reference data)
  "googlesheets4",
  
  # Other
  "grid",
  "forcats"
)

# Check which packages are missing
installed <- installed.packages()[, "Package"]
missing <- required_packages[!required_packages %in% installed]

if (length(missing) > 0) {
  cat("Installing missing packages:\n")
  cat("  ", paste(missing, collapse = ", "), "\n\n")
  
  # Separate CRAN and Bioconductor packages
  bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db")
  cran_missing <- missing[!missing %in% bioc_packages]
  bioc_missing <- missing[missing %in% bioc_packages]
  
  # Install CRAN packages
  if (length(cran_missing) > 0) {
    cat("Installing from CRAN...\n")
    install.packages(cran_missing, repos = "https://cloud.r-project.org")
  }
  
  # Install Bioconductor packages
  if (length(bioc_missing) > 0) {
    cat("Installing from Bioconductor...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
  }
  
  cat("\n")
} else {
  cat("All required packages are installed.\n\n")
}

# Verify installation
still_missing <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(still_missing) > 0) {
  warning("Some packages could not be installed: ", paste(still_missing, collapse = ", "))
}

# ============================================================================
# STEP 2: DOWNLOAD DATA
# ============================================================================

cat("--- Step 2: Downloading data files ---\n\n")

download_script <- file.path(project_root, "scripts", "download_data.R")

if (!file.exists(download_script)) {
  stop("download_data.R not found at: ", download_script)
}

source(download_script)

cat("\n")

# ============================================================================
# STEP 3: CREATE DERIVED DATA
# ============================================================================

cat("--- Step 3: Creating derived data files ---\n\n")

derived_script <- file.path(project_root, "scripts", "create_derived_data.R")

if (!file.exists(derived_script)) {
  stop("create_derived_data.R not found at: ", derived_script)
}

source(derived_script)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  Setup Complete!\n")
cat("================================================================================\n\n")

# Check what we have
data_raw <- file.path(project_root, "data_raw")
data_minikel <- file.path(project_root, "data", "minikel")

zenodo_files <- c(
  "ukb_ppp_mr_coloc_results.rds",
  "mr_prot_filtered_dataset_v1_v2.rds",
  "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds",
  "chembl.rds",
  "trans_genes.rds"
)

derived_files <- c(
  "pqtl_mrcoloc_2025.rds"
)

minikel_files <- c(
  "merge2.tsv.gz",
  "assoc.tsv.gz"
)

zenodo_present <- sum(file.exists(file.path(data_raw, zenodo_files)))
derived_present <- sum(file.exists(file.path(data_raw, derived_files)))
minikel_present <- sum(file.exists(file.path(data_minikel, minikel_files)))

cat("Data files status:\n")
cat("  Zenodo downloads:  ", zenodo_present, "/", length(zenodo_files), "\n")
cat("  Derived files:     ", derived_present, "/", length(derived_files), "\n")
cat("  Minikel et al.:    ", minikel_present, "/", length(minikel_files), "\n\n")

all_ready <- (zenodo_present == length(zenodo_files) && 
                derived_present == length(derived_files) &&
                minikel_present == length(minikel_files))

if (all_ready) {
  cat("All files present! You can now run the analysis:\n\n")
  cat("  # Generate main figures\n")
  cat("  Rscript scripts/mrcoloc_paper_2025_main_figures.R\n\n")
  cat("  # Generate supplementary tables\n")
  cat("  Rscript scripts/generate_mrcoloc_supplement.R\n\n")
  cat("  # Generate supplementary figures\n")
  cat("  Rscript scripts/mrcoloc_paper_2025_supp_figures.R\n\n")
} else {
  cat("WARNING: Some files are missing. Check the output above for errors.\n\n")
}

cat("================================================================================\n")