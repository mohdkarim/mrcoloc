#!/usr/bin/env Rscript
# ============================================================================
# Create Derived Data Files
# ============================================================================
#
# This script creates derived datasets from the raw Zenodo downloads:
#
#   1. pqtl_mrcoloc_2025.rds - Combined MR/coloc results (from mr_prot_unfiltered + ukb_ppp)
#
# Prerequisites:
#   Run scripts/download_data.R first to download raw data files.
#
# Usage:
#   Rscript scripts/create_derived_data.R
#
# ============================================================================

cat("
================================================================================
  Creating Derived Data Files
================================================================================
\n")

# ============================================================================
# SETUP
# ============================================================================

# Check for required packages
required_packages <- c("tidyverse", "data.table")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Determine project root
if (Sys.getenv("MRCOLOC_ROOT") != "") {
  project_root <- Sys.getenv("MRCOLOC_ROOT")
} else {
  if (file.exists("scripts/create_derived_data.R")) {
    project_root <- getwd()
  } else if (file.exists("create_derived_data.R")) {
    project_root <- dirname(getwd())
  } else {
    stop("Cannot determine project root. Run from project directory or set MRCOLOC_ROOT.")
  }
}

data_raw <- file.path(project_root, "data_raw")

cat("Project root:", project_root, "\n")
cat("Data directory:", data_raw, "\n\n")

# ============================================================================
# CHECK REQUIRED FILES
# ============================================================================

required_files <- c(
  
  "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds",
  "ukb_ppp_mr_coloc_results.rds"
)

cat("Checking required input files...\n")
missing_files <- character(0)

for (f in required_files) {
  filepath <- file.path(data_raw, f)
  if (file.exists(filepath)) {
    cat("  [OK]", f, "\n")
  } else {
    cat("  [MISSING]", f, "\n")
    missing_files <- c(missing_files, f)
  }
}

if (length(missing_files) > 0) {
  stop("\nMissing required files. Run scripts/download_data.R first.")
}

cat("\n")

# ============================================================================
# CREATE pqtl_mrcoloc_2025.rds
# ============================================================================

cat("--- Creating pqtl_mrcoloc_2025.rds ---\n\n")

output_file <- file.path(data_raw, "pqtl_mrcoloc_2025.rds")

if (file.exists(output_file)) {
  cat("  [SKIP] Already exists:", basename(output_file), "\n")
  cat("         Delete the file and re-run to regenerate.\n\n")
} else {
  
  # --- Part 1: Process older pQTL datasets ---
  cat("  Loading mr_prot_unfiltered dataset...\n")
  df <- readRDS(file.path(data_raw, 
                          "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds"))
  
  cat("    Rows:", format(nrow(df), big.mark = ","), "\n")
  
  # Calculate FDR
  cat("  Calculating FDR-adjusted p-values...\n")
  df$p_fdr <- p.adjust(df$bxy_pval, method = "fdr")
  
  # Filter to significant MR results (FDR < 5%)
  cat("  Filtering to FDR < 0.05...\n")
  dfmrcoloc_mr <- df %>% 
    filter(!is.na(exp_out_gsmr) & !is.na(ensid) & p_fdr < 0.05)
  
  cat("    Rows after filtering:", format(nrow(dfmrcoloc_mr), big.mark = ","), "\n")
  
  # Create key and select columns
  dfmrcoloc_mr$key <- with(dfmrcoloc_mr, paste0(hgnc_protein, "_", outcome))
  
  pqtl1 <- dfmrcoloc_mr %>%
    select(key, Data, nsnp, cis_trans_mr, bxy, bxy_pval, coloc_cis, coloc_h4_cis, 
           varid_left, coloc_trans, coloc_h4_trans, SNP_trans_coloc) %>%
    rename_with(~str_replace_all(., c(
      "nsnp" = "nsnps",
      "varid_left" = "snp_ciscoloc",
      "SNP_trans_coloc" = "snp_transcoloc"
    ))) %>%
    distinct()
  
  # Standardize cis_trans_mr values
  pqtl1$cis_trans_mr <- with(pqtl1, 
                             ifelse(cis_trans_mr == "cis", "Cis",
                                    ifelse(cis_trans_mr == "trans", "Trans",
                                           ifelse(cis_trans_mr == "mixed", "Mixed", cis_trans_mr))))
  
  cat("    pqtl1 rows:", format(nrow(pqtl1), big.mark = ","), "\n")
  
  # Clean up memory
  rm(df, dfmrcoloc_mr)
  gc(verbose = FALSE)
  
  # --- Part 2: Process UKB-PPP data ---
  cat("  Loading ukb_ppp_mr_coloc_results...\n")
  res <- readRDS(file.path(data_raw, "ukb_ppp_mr_coloc_results.rds"))
  
  cat("    Rows:", format(nrow(res), big.mark = ","), "\n")
  
  res <- res %>%
    mutate(across(where(is.character), ~na_if(.x, ""))) %>%
    mutate(
      outcome = coalesce(accession, trait),
      hgnc_protein = str_extract(prot, "^[^_]+"),
      Data = "UKBPPP_2023"
    ) %>%
    rename_with(~str_replace_all(., c(
      "IVs" = "cis_trans_mr",
      "pp4_cis" = "coloc_h4_cis",
      "pp4_trans" = "coloc_h4_trans"
    ))) %>%
    mutate(
      coloc_cis = if_else(!is.na(pp1_cis), "Yes", "No"),
      coloc_trans = if_else(!is.na(pp1_trans), "Yes", "No"),
      key = paste0(hgnc_protein, "_", outcome)
    )
  
  pqtl2 <- res %>%
    select(key, Data, nsnps, cis_trans_mr, bxy, bxy_pval, coloc_cis, coloc_h4_cis, 
           snp_ciscoloc, coloc_trans, coloc_h4_trans, snp_transcoloc) %>%
    distinct()
  
  cat("    pqtl2 rows:", format(nrow(pqtl2), big.mark = ","), "\n")
  
  # Clean up memory
  rm(res)
  gc(verbose = FALSE)
  
  # --- Merge datasets ---
  cat("  Merging pqtl1 and pqtl2...\n")
  pqtl <- rbind(pqtl1, pqtl2)
  
  cat("    Combined rows:", format(nrow(pqtl), big.mark = ","), "\n")
  
  # Save
  cat("  Saving pqtl_mrcoloc_2025.rds...\n")
  saveRDS(pqtl, output_file)
  
  file_size <- file.size(output_file) / 1e6
  cat("    -> Success:", round(file_size, 1), "MB\n\n")
  
  # Clean up
  rm(pqtl1, pqtl2, pqtl)
  gc(verbose = FALSE)
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("================================================================================\n")
cat("  Derived Data Summary\n")
cat("================================================================================\n\n")

# List all files in data_raw
cat("Files in data_raw/:\n\n")

files <- list.files(data_raw, pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  size <- file.size(f) / 1e6
  cat(sprintf("  %-60s %8.1f MB\n", basename(f), size))
}

cat("\n")
cat("Derived data creation complete!\n\n")
cat("You can now run the analysis scripts:\n")
cat("  Rscript scripts/mrcoloc_paper_2025_main_figures.R\n")
cat("  Rscript scripts/generate_mrcoloc_supplement.R\n\n")
cat("================================================================================\n")