#!/usr/bin/env Rscript
# ============================================================================
# Download Required Data Files
# ============================================================================
#
# This script downloads the required data files to reproduce the analyses in:
#
#   "Impact of proteogenomic evidence on clinical success"
#   Karim et al., Nature Genetics (2025)
#
# Data sources:
#   1. Zenodo: Primary MR/coloc datasets (DOI: 10.5281/zenodo.18451758)
#   2. GitHub: Minikel et al. Nature 2024 reference data
#
# Usage:
#   Rscript scripts/download_data.R
#
# Or from R:
#   source("scripts/download_data.R")
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

# Files to download from Zenodo (5 files, ~922 MB total)
ZENODO_FILES <- c(
  "ukb_ppp_mr_coloc_results.rds",
  "mr_prot_filtered_dataset_v1_v2.rds",
  "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds",
  "chembl.rds",
  "trans_genes.rds"
)

# Minikel et al. Nature 2024 - GitHub repository
# Paper: https://doi.org/10.1038/s41586-024-07316-0
# Repo: https://github.com/ericminikel/genetic_support
MINIKEL_BASE_URL <- "https://raw.githubusercontent.com/ericminikel/genetic_support/main/data/"

# Files to download from Minikel et al. GitHub
# Core files go to data/minikel/, reference files go to data/
MINIKEL_FILES <- c(
  "merge2.tsv.gz",
  "assoc.tsv.gz",
  "areas.tsv",
  "drug_phase_summary.tsv",
  "drugs_per_indic.tsv",
  "gwas_catalog-ancestry_r2024-04-16.tsv",
  "indic.tsv",
  "indic_topl_match.tsv",
  "intogen_genes.tsv",
  "mesh_2023_topl_maps.tsv",
  "mesh_all_vocab.tsv.gz",
  "mesh_best_names.tsv.gz",
  "mesh_remappings.tsv",
  "mesh_scr_to_best.tsv.gz",
  "mesh_tree.tsv",
  "meta_acat.tsv",
  "meta_ccat.tsv",
  "meta_hcat.tsv",
  "nelson_2015_by_category.tsv",
  "nelson_2015_categories.tsv",
  "nelson_2015_category_match.tsv",
  "nelson_2015_supplementary_dataset_2.tsv",
  "nelson_2015_supplementary_dataset_3.tsv",
  "nelson_2015_table_s6.tsv",
  "nelson_2015_table_s9.tsv",
  "otg_heavy_hitter_curation.tsv",
  "pp.tsv",
  "pp_launched_alltime.tsv",
  "sim.tsv.gz",
  "universe.tsv"
)

# Gene list files (in data/gene_lists/)
GENE_LIST_FILES <- c(
  "ab_tractable.tsv",
  "catalytic_receptors.tsv",
  "enzymes.tsv",
  "ion_channels.tsv",
  "kinases.tsv",
  "nuclear_receptors.tsv",
  "rhodop_gpcr.tsv",
  "sm_tractable.tsv",
  "transporters.tsv"
)

# ============================================================================
# SETUP
# ============================================================================

# Determine project root
if (Sys.getenv("MRCOLOC_ROOT") != "") {
  project_root <- Sys.getenv("MRCOLOC_ROOT")
} else {
  if (file.exists("scripts/download_data.R")) {
    project_root <- getwd()
  } else if (file.exists("download_data.R")) {
    project_root <- dirname(getwd())
  } else {
    stop("Cannot determine project root. Run from project directory or set MRCOLOC_ROOT.")
  }
}

# Create data directories
data_raw_dir <- file.path(project_root, "data_raw")
data_dir <- file.path(project_root, "data")
data_minikel_dir <- file.path(data_dir, "minikel")
gene_lists_dir <- file.path(data_dir, "gene_lists")

dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_minikel_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gene_lists_dir, showWarnings = FALSE, recursive = TRUE)

cat("Project root:", project_root, "\n")
cat("Zenodo data directory:", data_raw_dir, "\n")
cat("Minikel data directory:", data_dir, "\n\n")

# ============================================================================
# DOWNLOAD FUNCTION
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
    if (Sys.which("curl") != "") {
      exit_code <- system2("curl", 
                           c("-L", "-f", "-o", shQuote(destfile), shQuote(url)), 
                           stdout = FALSE, stderr = FALSE)
      success <- (exit_code == 0)
    } else {
      download.file(url, destfile, mode = "wb", quiet = TRUE)
      success <- TRUE
    }
    
    if (success && file.exists(destfile) && file.size(destfile) > 0) {
      size_mb <- file.size(destfile) / 1e6
      if (size_mb >= 1) {
        cat("    -> Success:", round(size_mb, 1), "MB\n")
      } else {
        cat("    -> Success:", round(file.size(destfile) / 1e3, 1), "KB\n")
      }
      return(TRUE)
    } else {
      cat("    -> FAILED: Empty or missing file\n")
      if (file.exists(destfile)) file.remove(destfile)
      return(FALSE)
    }
  }, error = function(e) {
    cat("    -> ERROR:", e$message, "\n")
    if (file.exists(destfile)) file.remove(destfile)
    return(FALSE)
  })
}

# ============================================================================
# DOWNLOAD ZENODO FILES
# ============================================================================

cat("--- Downloading from Zenodo (DOI:", ZENODO_DOI, ") ---\n\n")

zenodo_results <- sapply(ZENODO_FILES, function(filename) {
  url <- paste0(ZENODO_BASE_URL, filename, "?download=1")
  destfile <- file.path(data_raw_dir, filename)
  download_file(url, destfile, filename)
})

# ============================================================================
# DOWNLOAD MINIKEL ET AL. FILES
# ============================================================================

cat("\n--- Downloading from Minikel et al. GitHub ---\n")
cat("    (https://github.com/ericminikel/genetic_support)\n\n")

# Files that go in data/minikel/ (core merge files)
MINIKEL_CORE_FILES <- c("merge2.tsv.gz", "assoc.tsv.gz")

minikel_results <- sapply(MINIKEL_FILES, function(filename) {
  url <- paste0(MINIKEL_BASE_URL, filename)
  # Core files go to data/minikel/, others go to data/
  if (filename %in% MINIKEL_CORE_FILES) {
    destfile <- file.path(data_minikel_dir, filename)
  } else {
    destfile <- file.path(data_dir, filename)
  }
  download_file(url, destfile, filename)
})

# ============================================================================
# DOWNLOAD GENE LIST FILES
# ============================================================================

cat("\n--- Downloading gene list files ---\n\n")

GENE_LIST_BASE_URL <- "https://raw.githubusercontent.com/ericminikel/genetic_support/main/data/gene_lists/"

gene_list_results <- sapply(GENE_LIST_FILES, function(filename) {
  url <- paste0(GENE_LIST_BASE_URL, filename)
  destfile <- file.path(gene_lists_dir, filename)
  download_file(url, destfile, filename)
})

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  Download Summary\n")
cat("================================================================================\n\n")

# Count successes
zenodo_success <- sum(zenodo_results)
minikel_success <- sum(minikel_results)
gene_list_success <- sum(gene_list_results)

cat("Zenodo files:  ", zenodo_success, "/", length(ZENODO_FILES), "\n")
cat("Minikel files: ", minikel_success, "/", length(MINIKEL_FILES), "\n")
cat("Gene lists:    ", gene_list_success, "/", length(GENE_LIST_FILES), "\n\n")

all_success <- zenodo_success == length(ZENODO_FILES) && 
               minikel_success == length(MINIKEL_FILES) &&
               gene_list_success == length(GENE_LIST_FILES)

if (all_success) {
  cat("All downloads complete!\n\n")
  cat("Next step: Run the derived data script to create pqtl_mrcoloc_2025.rds:\n")
  cat("  Rscript scripts/create_derived_data.R\n\n")
  cat("Or run the full setup:\n")
  cat("  Rscript scripts/setup.R\n\n")
} else {
  if (zenodo_success < length(ZENODO_FILES)) {
    cat("WARNING: Some Zenodo downloads failed. Try running this script again.\n")
    cat("  Failed files:\n")
    for (f in names(zenodo_results)[!zenodo_results]) {
      cat("    -", f, "\n")
    }
    cat("\n")
  }
  
  if (minikel_success < length(MINIKEL_FILES)) {
    cat("WARNING: Some Minikel et al. downloads failed.\n")
    cat("  You can manually download from:\n")
    cat("  https://github.com/ericminikel/genetic_support/tree/main/data\n\n")
  }
  
  if (gene_list_success < length(GENE_LIST_FILES)) {
    cat("WARNING: Some gene list downloads failed.\n")
    cat("  You can manually download from:\n")
    cat("  https://github.com/ericminikel/genetic_support/tree/main/data/gene_lists\n\n")
  }
}

cat("================================================================================\n")

# Return status invisibly
invisible(list(
  zenodo = zenodo_results,
  minikel = minikel_results,
  gene_lists = gene_list_results,
  success = all_success
))
