#!/usr/bin/env Rscript
# ============================================================================
# Download Required Data Files
# ============================================================================
#
# This script downloads the required data files from Zenodo to reproduce
# the analyses in:
#
#   "Impact of proteogenomic evidence on clinical success"
#   Karim et al., Nature Genetics (2025)
#
# Usage:
#   Rscript scripts/download_data.R
#
# Or from R:
#   source("scripts/download_data.R")
#
# Files are downloaded to: data_raw/
#
# ============================================================================

cat("
================================================================================
  Downloading Data Files for mrcoloc Analysis
================================================================================
\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Zenodo DOI and record ID
ZENODO_DOI <- "10.5281/zenodo.18451758"
ZENODO_RECORD <- "18451758"

# Base URL for Zenodo downloads
ZENODO_BASE_URL <- paste0("https://zenodo.org/records/", ZENODO_RECORD, "/files/")

# Files to download from Zenodo
ZENODO_FILES <- c(
  "pqtl_mrcoloc_2025.rds",
  "ukb_ppp_mr_coloc_results.rds",
  "mr_prot_filtered_dataset_v1_v2.rds",
  "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds",
  "ttpairs_tested.rds",
  "trans_genes.rds",
  "chembl.rds",
  "panukb.rds",
  "ST7_all_MR_pairs.rds"
)

# Minikel et al. files (users must download manually from Nature supplement)
MINIKEL_FILES <- c(
  "merge2.tsv.gz",
  "assoc.tsv.gz"
)

# ============================================================================
# SETUP
# ============================================================================

# Determine project root
if (Sys.getenv("PQTL_ENRICH_ROOT") != "") {
  project_root <- Sys.getenv("PQTL_ENRICH_ROOT")
} else {
  # Assume script is run from project root or scripts/ directory
  if (file.exists("scripts/download_data.R")) {
    project_root <- getwd()
  } else if (file.exists("download_data.R")) {
    project_root <- dirname(getwd())
  } else {
    stop("Cannot determine project root. Set PQTL_ENRICH_ROOT environment variable.")
  }
}

# Create data directories
data_raw_dir <- file.path(project_root, "data_raw")
data_minikel_dir <- file.path(project_root, "data", "minikel")

dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_minikel_dir, showWarnings = FALSE, recursive = TRUE)

cat("Project root:", project_root, "\n")
cat("Data directory:", data_raw_dir, "\n\n")

# ============================================================================
# DOWNLOAD FUNCTIONS
# ============================================================================
download_file <- function(url, destfile, description = NULL) {
  if (file.exists(destfile)) {
    cat("  [SKIP] Already exists:", basename(destfile), "\n")
    return(TRUE)
  }
  
  if (!is.null(description)) {
    cat("  [DOWNLOAD]", description, "\n")
  } else {
    cat("  [DOWNLOAD]", basename(destfile), "\n")
  }
  
  tryCatch({
    # Use curl for better progress and reliability
    if (Sys.which("curl") != "") {
      system2("curl", c("-L", "-o", destfile, url), stdout = FALSE)
    } else {
      download.file(url, destfile, mode = "wb", quiet = FALSE)
    }
    
    if (file.exists(destfile) && file.size(destfile) > 0) {
      cat("    -> Success:", round(file.size(destfile) / 1e6, 1), "MB\n")
      return(TRUE)
    } else {
      cat("    -> FAILED: Empty or missing file\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("    -> ERROR:", e$message, "\n")
    return(FALSE)
  })
}

# ============================================================================
# DOWNLOAD ZENODO FILES
# ============================================================================

cat("Downloading files from Zenodo (DOI:", ZENODO_DOI, ")...\n\n")

zenodo_success <- TRUE
for (filename in ZENODO_FILES) {
  url <- paste0(ZENODO_BASE_URL, filename, "?download=1")
  destfile <- file.path(data_raw_dir, filename)
  success <- download_file(url, destfile, filename)
  if (!success) zenodo_success <- FALSE
}

# ============================================================================
# CHECK MINIKEL FILES
# ============================================================================

cat("\n")
cat("Checking for Minikel et al. data files...\n\n")

minikel_missing <- character(0)
for (filename in MINIKEL_FILES) {
  filepath <- file.path(data_minikel_dir, filename)
  if (file.exists(filepath)) {
    cat("  [OK]", filename, "\n")
  } else {
    cat("  [MISSING]", filename, "\n")
    minikel_missing <- c(minikel_missing, filename)
  }
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  Download Summary\n")
cat("================================================================================\n\n")

# Check what we have
zenodo_present <- sum(file.exists(file.path(data_raw_dir, ZENODO_FILES)))
minikel_present <- sum(file.exists(file.path(data_minikel_dir, MINIKEL_FILES)))

cat("Zenodo files:", zenodo_present, "/", length(ZENODO_FILES), "\n")
cat("Minikel files:", minikel_present, "/", length(MINIKEL_FILES), "\n\n")

if (zenodo_present == length(ZENODO_FILES) && minikel_present == length(MINIKEL_FILES)) {
  cat("✓ All required files present! Ready to run analyses.\n\n")
} else {
  if (zenodo_present < length(ZENODO_FILES)) {
    cat("⚠ Some Zenodo downloads failed. Try running this script again.\n\n")
  }
  
  if (length(minikel_missing) > 0) {
    cat("⚠ Minikel et al. files must be downloaded manually:\n\n")
    cat("  1. Go to: https://doi.org/10.1038/s41586-024-07316-0\n")
    cat("  2. Download Supplementary Data files\n")
    cat("  3. Place these files in:", data_minikel_dir, "\n")
    cat("     -", paste(minikel_missing, collapse = "\n     - "), "\n\n")
  }
}

cat("================================================================================\n")