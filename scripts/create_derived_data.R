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

# Auto-install missing CRAN packages
cran_pkgs <- c("tidyverse", "data.table")
missing <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if (length(missing) > 0) {
  cat("Installing missing CRAN packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# Bioconductor packages (needed for pgenes construction)
bioc_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[,"Package"]]
if (length(missing_bioc) > 0) {
  cat("Installing missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
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

# ---------------------------------------------------------------------------
# [R2.10 FIX] Revised for the Nature Medicine revision, in response to
# Reviewer 2 comment 10. Two defects in the previous version:
#
#  (1) Olink panel accessions were mapped to symbols by a single route
#      (org.Hs.eg.db) and unmapped accessions were then silently dropped by
#      filter(!is.na(...)). 40 of 2,960 accessions fail that mapping, and two of
#      them - LPA and PSCA - are drug targets assayed by UKB-PPP with
#      Bonferroni-significant pQTL MR associations. They were therefore missing
#      from the very background they belong in, which inflated the headline
#      relative success estimate. Accessions are now resolved through BOTH
#      available routes (the local data/olink_complete_extended.tsv UniProt ->
#      Gene name table, and org.Hs.eg.db) and the union taken; 22 symbols
#      resolve only via the local table and 103 only via org.Hs.eg.db, so
#      neither route alone is sufficient. Accessions unresolvable by both are
#      now REPORTED rather than silently discarded.
#
#  (2) No per-platform measured-protein sets existed, so the per-platform rows
#      of Figure 1a fell back to `platform == p` in merge3_pqtl. That column is
#      only populated on Bonferroni-significant rows, so those backgrounds were
#      conditioned on the very evidence under test. We now also emit
#      pgenes_platform.rds, which decomposes pgenes by platform, and the figure
#      script uses it for the platform backgrounds.
#
# A protein tested for pQTL MR is additionally unioned in as a completeness
# backstop, restricted to symbols that appear as targets in the therapeutic
# index - only those can affect any target-indication estimate, and the
# restriction discards parsing junk by construction (sub("_.*$", "", key)
# truncates compound protein names such as CKMT1A_CKMT1B).
#
# LIMITATION: the unfiltered MR dataset does not contain UKBPPP_2023, which is
# why the Olink panel manifest is needed to represent UKB-PPP. For the seven
# other studies the protein lists derive from proteins that had a testable
# instrument, so assayed proteins with no detected pQTL are under-represented
# for the SomaScan studies. Closing that would require their assay manifests,
# which are not available here. See peer_review/background_bug_and_fix.txt.
# ---------------------------------------------------------------------------

pgenes_file   <- file.path(data_raw, "pgenes.rds")
platform_file <- file.path(data_raw, "pgenes_platform.rds")
rebuild_pgenes <- isTRUE(as.logical(Sys.getenv("MRCOLOC_REBUILD_PGENES", "FALSE")))

if (file.exists(pgenes_file) && file.exists(platform_file) && !rebuild_pgenes) {
  cat("  [SKIP] Already exists:", basename(pgenes_file), "+", basename(platform_file), "\n\n")
} else {
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })

  # Sentinel strings that must never enter a gene set. "NULL" was present in the
  # pre-revision pgenes.rds, from as.character() on an unmapped accession.
  BAD_SYMBOLS <- c("NA", "NULL", "", "NaN")

  plat_of <- function(x) case_when(
    x %in% c("UKBPPP_2023", "SCALLOP_2020", "HILLARY_2019", "FOLKERSEN_2017") ~ "Olink",
    x %in% c("SUN_2018", "SUHRE_2017", "PIETZNER_2020")                       ~ "Somascan",
    x == "OLLI_2017"                                                          ~ "Other",
    TRUE ~ NA_character_
  )

  # --- 1. Olink assay manifest, resolved by both available routes ------------
  cat("  Resolving Olink panel accessions to gene symbols (two routes)...\n")
  acc <- read_tsv(file.path(project_root, "data", "olink_complete.tsv"),
                  show_col_types = FALSE) %>%
    separate_rows(`Uniprot ID`, sep = ",") %>%
    mutate(up = trimws(`Uniprot ID`)) %>%
    filter(!is.na(up), up != "") %>%
    distinct(up)

  local_map <- read_tsv(file.path(project_root, "data", "olink_complete_extended.tsv"),
                        show_col_types = FALSE) %>%
    transmute(up = `UniProt ID`, sym = `Gene name`) %>%
    filter(!is.na(up), !is.na(sym)) %>%
    distinct()

  res <- acc %>%
    left_join(local_map, by = "up") %>%
    group_by(up) %>%
    summarise(sym_local = dplyr::first(na.omit(sym)), .groups = "drop")

  res$sym_orgdb <- tryCatch(
    as.character(suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = res$up, column = "SYMBOL",
      keytype = "UNIPROT", multiVals = "first"))),
    error = function(e) { cat("    [warn] org.Hs.eg.db UNIPROT mapping unavailable\n"); NA_character_ }
  )

  olink_panel <- setdiff(unique(na.omit(c(res$sym_local, res$sym_orgdb))), BAD_SYMBOLS)
  cat(sprintf("    accessions %d -> symbols %d | local-only %d, orgdb-only %d, unresolved %d\n",
              nrow(res), length(olink_panel),
              sum(!is.na(res$sym_local) & is.na(res$sym_orgdb)),
              sum(is.na(res$sym_local) & !is.na(res$sym_orgdb)),
              sum(is.na(res$sym_local) & is.na(res$sym_orgdb))))

  # --- 2. Per-study protein lists from the unfiltered MR dataset -------------
  cat("  Extracting per-study protein lists (large file, ~3 min)...\n")
  df_unfiltered <- readRDS(file.path(data_raw,
    "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds"))
  study_prot <- df_unfiltered %>%
    transmute(prot = hgnc_protein, study = as.character(Data)) %>%
    filter(!is.na(prot), prot != "") %>%
    mutate(platform = plat_of(study)) %>%
    distinct()
  rm(df_unfiltered); gc(verbose = FALSE)

  # --- 3. Completeness backstop: proteins tested for pQTL MR ----------------
  # Vocabulary is merge2's gene column (the therapeutic index). Verified
  # identical to unique(merge3_pqtl$gene) - 2,517 symbols both ways - and used
  # here because merge3_pqtl.rds is created later in this script.
  cat("  Extracting tested proteins from pqtl_mrcoloc_2025.rds...\n")
  ti_genes <- read_tsv(file.path(project_root, "data", "minikel", "merge2.tsv.gz"),
                       show_col_types = FALSE) %>%
    filter(!is.na(gene), gene != "") %>% distinct(gene) %>% pull(gene)

  pq <- readRDS(file.path(data_raw, "pqtl_mrcoloc_2025.rds"))
  tested <- tibble(prot = sub("_.*$", "", pq$key), study = as.character(pq$Data)) %>%
    filter(!is.na(prot), prot != "") %>%
    mutate(platform = plat_of(study)) %>%
    distinct()
  rm(pq); gc(verbose = FALSE)

  # --- 4. Assemble ----------------------------------------------------------
  pset <- function(p) setdiff(unique(c(
    study_prot$prot[study_prot$platform == p],
    intersect(tested$prot[tested$platform == p], ti_genes))), BAD_SYMBOLS)

  pgenes_platform <- list(
    Olink    = setdiff(unique(c(olink_panel, pset("Olink"))), BAD_SYMBOLS),
    Somascan = pset("Somascan"),
    Other    = pset("Other")
  )
  pgenes <- unique(unlist(pgenes_platform, use.names = FALSE))

  stopifnot(setequal(pgenes, unique(unlist(pgenes_platform))))

  saveRDS(pgenes, pgenes_file)
  saveRDS(pgenes_platform, platform_file)
  cat(sprintf("    -> Success: %d genes (Olink %d, Somascan %d, Other %d)\n\n",
              length(pgenes), length(pgenes_platform$Olink),
              length(pgenes_platform$Somascan), length(pgenes_platform$Other)))
  rm(acc, local_map, res, olink_panel, study_prot, tested, ti_genes,
     pgenes, pgenes_platform)
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