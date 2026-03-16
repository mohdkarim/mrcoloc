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

t0 <- Sys.time()
t_step <- Sys.time()

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

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
cat("--- Creating pqtl_mrcoloc_2025.rds ---\n\n")
t_step <- Sys.time()

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
# CREATE pgenes.rds (background gene set)
# ============================================================================

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
cat("--- Creating pgenes.rds ---\n\n")
t_step <- Sys.time()

pgenes_file <- file.path(data_raw, "pgenes.rds")

if (file.exists(pgenes_file)) {
  cat("  [SKIP] Already exists:", basename(pgenes_file), "\n\n")
} else {
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })

  # Olink genes
  cat("  Loading Olink panel genes...\n")
  olink <- read_tsv(file.path(project_root, "data", "olink_complete.tsv"), show_col_types = FALSE)
  olink2 <- olink %>%
    separate_rows(`Uniprot ID`, sep = ",") %>%
    mutate(hgnc_protein = AnnotationDbi::mapIds(org.Hs.eg.db, keys = `Uniprot ID`,
           column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")) %>%
    filter(!is.na(hgnc_protein))
  olink_genes <- unique(as.character(olink2$hgnc_protein))

  # Other genes from unfiltered dataset
  cat("  Extracting gene list from unfiltered dataset...\n")
  df_unfiltered <- readRDS(file.path(data_raw,
    "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds"))
  othergenes <- unique(df_unfiltered$hgnc_protein[!is.na(df_unfiltered$hgnc_protein)])
  rm(df_unfiltered); gc(verbose = FALSE)

  pgenes <- unique(c(olink_genes, othergenes))
  saveRDS(pgenes, pgenes_file)
  cat("    -> Success:", length(pgenes), "genes\n\n")
  rm(olink, olink2, olink_genes, othergenes, pgenes)
  gc(verbose = FALSE)
}

# ============================================================================
# CREATE merge3_pqtl.rds (merged therapeutic index + pQTL data)
# ============================================================================

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
cat("--- Creating merge3_pqtl.rds ---\n\n")
t_step <- Sys.time()

merge3_file <- file.path(data_raw, "merge3_pqtl.rds")

if (file.exists(merge3_file)) {
  cat("  [SKIP] Already exists:", basename(merge3_file), "\n\n")
} else {
  mr_pval_threshold <- 0.05 / 47e6

  cat("  Loading pqtl_mrcoloc_2025.rds and filtering by Bonferroni...\n")
  pqtl2 <- readRDS(file.path(data_raw, "pqtl_mrcoloc_2025.rds")) %>%
    filter(bxy_pval <= mr_pval_threshold)

  cat("  Loading merge2.tsv.gz...\n")
  minikel_dir <- file.path(project_root, "data", "minikel")
  merge2 <- read_tsv(file.path(minikel_dir, "merge2.tsv.gz"), show_col_types = FALSE) %>%
    mutate(
      otg_study = if_else(assoc_source == "OTG",
        str_remove(original_link, "https://genetics.opentargets.org/study/"), NA_character_),
      otg_study = str_remove(otg_study, "FINNGEN_R6_"),
      key = paste0(gene, "_", otg_study)
    )

  cat("  Merging pQTL data with therapeutic index...\n")
  pqtl_cols <- c("nsnps", "cis_trans_mr", "bxy", "bxy_pval",
                 "coloc_cis", "coloc_h4_cis", "snp_ciscoloc",
                 "coloc_trans", "coloc_h4_trans", "snp_transcoloc")

  merge2_with_pqtl <- merge2 %>% left_join(pqtl2, by = "key")
  pqtl_rows <- merge2_with_pqtl %>% filter(!is.na(cis_trans_mr)) %>% mutate(original_link = "pqtl")
  merge2_cleaned <- merge2_with_pqtl %>% mutate(across(all_of(pqtl_cols), ~ if_else(!is.na(cis_trans_mr), NA, .)))
  merge3_pqtl <- bind_rows(merge2_cleaned, pqtl_rows)

  # Platform annotations
  merge3_pqtl <- merge3_pqtl %>%
    mutate(platform = case_when(
      Data %in% c("UKBPPP_2023", "SCALLOP_2020", "HILLARY_2019", "FOLKERSEN_2017") ~ "Olink",
      Data %in% c("SUN_2018", "SUHRE_2017", "PIETZNER_2020") ~ "Somascan",
      Data == "OLLI_2017" ~ "Other",
      TRUE ~ NA_character_
    ))

  # Therapeutic area
  area <- fread(file.path(project_root, "data", "areas.tsv"))
  topl <- fread(file.path(project_root, "data", "indic_topl_match.tsv"))
  ta <- merge(area, topl, by = "topl")
  merge3_pqtl$therapeutic_area <- ta$area[match(merge3_pqtl$indication_mesh_id, ta$indication_mesh_id)]

  saveRDS(merge3_pqtl, merge3_file)
  file_size <- file.size(merge3_file) / 1e6
  cat("    -> Success:", format(nrow(merge3_pqtl), big.mark = ","), "rows,", round(file_size, 1), "MB\n\n")
  rm(pqtl2, merge2, merge2_with_pqtl, pqtl_rows, merge2_cleaned, merge3_pqtl, area, topl, ta)
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

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("\nTotal time: ", round(difftime(Sys.time(), t0, units="mins"), 1), " minutes")