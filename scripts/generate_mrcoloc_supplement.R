#!/usr/bin/env Rscript
# ============================================================================
# Generate Supplementary Tables for Proteogenomics Paper
# Consolidated script to produce mrcoloc_supplement.xlsx
# ============================================================================
#
# This script generates all supplementary tables matching the structure of
# mrcoloc_supplement_v3.xlsx:
#
#   Key           - Definitions and column descriptions
#   ST1           - Outcome GWAS datasets
#   ST2           - Proteomic GWAS datasets
#   ST3           - Excluded traits
#   ST4           - Figure 1a enrichment data
#   ST5           - Figure 1b UpSet plot data
#   ST6           - Figure 1c gene family enrichment data
#   ST7           - Chi-square distribution (pQTL across TAs)
#   ST8           - TA Enrichment Universe 1
#   ST9           - TA Enrichment Universe 2
#   ST10          - Spearman Summary
#   ST11          - Per-TA RS Universe 1
#   ST12          - Per-TA RS Universe 2
#   ST13          - Breslow-Day Summary
#   ST14          - LOO Universe 1
#   ST15          - LOO Universe 2
#   ST16          - All MR target-trait pairs
#   ST17          - pQTL success TI pairs
#
# Usage:
#   Rscript generate_mrcoloc_supplement.R
#
# Output:
#   mrcoloc_supplement.xlsx
#
# Author: Mohd Anisul Karim <mohd@variantbio.com>
# ============================================================================

cat("
================================================================================
  Generating Supplementary Tables (mrcoloc_supplement.xlsx)
================================================================================
\n")

# ============================================================================
# SECTION 0: SETUP
# ============================================================================

message("[0/12] Loading packages and configuration...")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(tidyverse)
  library(data.table)
  library(Rmpfr)
  library(stringr)
  library(openxlsx)
  library(DescTools)
  library(googlesheets4)
})

gs4_deauth()

# --- Paths ---
project_root <- Sys.getenv("PQTL_ENRICH_ROOT", 
                           "/home/mohd/mohd-sandbox/pQTL_enrichment/mrcoloc_paper2025")
data_raw     <- file.path(project_root, "data_raw")
data_dir     <- file.path(project_root, "genetic_support-main", "data")
output_dir   <- file.path(project_root, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Thresholds ---
BONFERRONI_P     <- 0.05 / 47e6
SIMILARITY_THR   <- 0.8
MIN_L2G_SHARE    <- 0.5
COLOC_H4_THR     <- 0.8

# --- Source helpers ---
source(file.path(project_root, "R/genes_tbl.R"))

# ============================================================================
# SECTION 1: HELPER FUNCTIONS
# ============================================================================

message("[1/12] Setting up helper functions...")

# Clean dataframe for export (remove list columns)
clean_for_export <- function(df) {
  df %>% ungroup() %>% distinct() %>% select(where(~ !is.list(.x)))
}

# Compute mantissa/exponent for very small p-values
compute_mantissa_exponent <- function(df, pval_cols) {
  for (pval_col in pval_cols) {
    beta_col <- gsub("_pval|_p", "", pval_col)
    se_col   <- gsub("_pval|_p", "_se", pval_col)
    
    pvals <- df[[pval_col]]
    betas <- df[[beta_col]]
    ses   <- df[[se_col]]
    
    mantissa <- exponent <- numeric(length(pvals))
    
    nz <- pvals > 0 & !is.na(pvals)
    if (any(nz)) {
      exponent[nz] <- floor(log10(pvals[nz]))
      mantissa[nz] <- pvals[nz] / (10^exponent[nz])
    }
    
    z <- pvals == 0 & !is.na(pvals)
    if (any(z)) {
      z_scores <- abs(betas[z] / ses[z])
      pval_mpfr <- 2 * pnorm(mpfr(z_scores, precBits = 200), lower.tail = FALSE)
      str_vals <- strsplit(format(pval_mpfr, scientific = TRUE), "e")
      mantissa[z] <- as.numeric(sapply(str_vals, `[[`, 1))
      exponent[z] <- as.integer(sapply(str_vals, `[[`, 2))
    }
    
    df[[paste0(pval_col, "_mantissa")]] <- mantissa
    df[[paste0(pval_col, "_exponent")]] <- exponent
  }
  df
}

# Harmonize column types between two dataframes
harmonize_column_types <- function(df1, df2, verbose = FALSE) {
  for (col in intersect(names(df1), names(df2))) {
    t1 <- class(df1[[col]])[1]; t2 <- class(df2[[col]])[1]
    if (t1 != t2) {
      if (verbose) message("Type mismatch: ", col, " (", t1, " vs ", t2, ")")
      if (t1 == "list" && all(lengths(df1[[col]]) <= 1)) df1[[col]] <- unlist(df1[[col]])
      if (t2 == "list" && all(lengths(df2[[col]]) <= 1)) df2[[col]] <- unlist(df2[[col]])
      if (class(df1[[col]])[1] %in% c("numeric","integer","double") ||
          class(df2[[col]])[1] %in% c("numeric","integer","double")) {
        df1[[col]] <- suppressWarnings(as.numeric(df1[[col]]))
        df2[[col]] <- suppressWarnings(as.numeric(df2[[col]]))
      } else {
        df1[[col]] <- as.character(df1[[col]])
        df2[[col]] <- as.character(df2[[col]])
      }
    }
  }
  list(df1 = df1, df2 = df2)
}

# Relative risk with Katz CI
rr_katz <- function(x1, n1, x2, n2, add_half = FALSE) {
  if (anyNA(c(x1, n1, x2, n2)) || any(c(n1, n2) <= 0) || 
      x1 < 0 || x2 < 0 || x1 > n1 || x2 > n2) {
    return(tibble(rr = NA_real_, lwr = NA_real_, upr = NA_real_))
  }
  if (add_half && (x1 == 0 || x2 == 0)) {
    x1 <- x1 + 0.5; x2 <- x2 + 0.5; n1 <- n1 + 1; n2 <- n2 + 1
  }
  if (x2 == 0) return(tibble(rr = Inf, lwr = NA_real_, upr = NA_real_))
  if (x1 == 0) return(tibble(rr = 0, lwr = NA_real_, upr = NA_real_))
  out <- as.data.frame(DescTools::BinomRatioCI(x1 = x1, n1 = n1, x2 = x2, n2 = n2, method = "katz"))
  tibble(rr = out$est, lwr = out$lwr.ci, upr = out$upr.ci)
}

# Build 3D array for Breslow-Day test
build_3d_array <- function(ta_table) {
  K <- nrow(ta_table)
  arr <- array(0, dim = c(2, 2, K),
               dimnames = list(support = c("pqtl", "no_pqtl"),
                               outcome = c("success", "failure"),
                               ta = ta_table$therapeutic_area))
  for (i in seq_len(K)) {
    arr["pqtl", "success", i]    <- ta_table$pqtl_success[i]
    arr["pqtl", "failure", i]    <- ta_table$pqtl_total[i] - ta_table$pqtl_success[i]
    arr["no_pqtl", "success", i] <- ta_table$nopqtl_success[i]
    arr["no_pqtl", "failure", i] <- ta_table$nopqtl_total[i] - ta_table$nopqtl_success[i]
  }
  arr
}

# Breslow-Day test for heterogeneity
breslow_day_test <- function(array_3d) {
  K <- dim(array_3d)[3]
  if (K < 2) return(list(statistic = NA, df = NA, p.value = NA))
  mh_result <- tryCatch(mantelhaen.test(array_3d, correct = FALSE), error = function(e) NULL)
  if (is.null(mh_result)) return(list(statistic = NA, df = NA, p.value = NA))
  mh_or <- mh_result$estimate
  
  chi_sq <- 0; valid_strata <- 0
  for (i in seq_len(K)) {
    a <- array_3d[1,1,i]; b <- array_3d[1,2,i]; c <- array_3d[2,1,i]; d <- array_3d[2,2,i]
    n <- a + b + c + d
    if (n == 0) next
    n1 <- a + b; m1 <- a + c
    coef_a <- 1 - mh_or; coef_b <- mh_or * (n1 + m1) + (n - n1 - m1); coef_c <- -mh_or * n1 * m1
    if (abs(coef_a) < 1e-10) { a_exp <- -coef_c / coef_b
    } else {
      disc <- coef_b^2 - 4 * coef_a * coef_c
      if (disc < 0) next
      a_exp <- (-coef_b + sqrt(disc)) / (2 * coef_a)
    }
    b_exp <- n1 - a_exp; c_exp <- m1 - a_exp; d_exp <- n - n1 - m1 + a_exp
    if (any(c(a_exp, b_exp, c_exp, d_exp) <= 0)) next
    var_a <- 1 / (1/a_exp + 1/b_exp + 1/c_exp + 1/d_exp)
    if (!is.na(var_a) && var_a > 0) { chi_sq <- chi_sq + (a - a_exp)^2 / var_a; valid_strata <- valid_strata + 1 }
  }
  df <- valid_strata - 1
  if (df < 1) return(list(statistic = NA, df = NA, p.value = NA))
  list(statistic = chi_sq, df = df, p.value = pchisq(chi_sq, df, lower.tail = FALSE))
}

# Helper: Get best T-I pair per target-indication
get_ti_best <- function(df, support_tis, indic_df) {
  df %>%
    filter(!is.na(gene), gene != "", !is.na(indication_mesh_id), !is.na(ccat)) %>%
    left_join(indic_df %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
    filter(genetic_insight != "none") %>%
    mutate(
      highest_phase = case_when(
        !is.na(succ_3_a) ~ 3,
        !is.na(succ_2_3) ~ 2,
        !is.na(succ_1_2) ~ 1,
        !is.na(succ_p_1) ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
    group_by(ti_uid) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gensup = ti_uid %in% support_tis$ti_uid)
}

# Helper: Compute relative success (RS)
compute_rs <- function(ti_best, pgenes_filter = NULL) {
  if (!is.null(pgenes_filter)) {
    ti_best <- ti_best %>% filter(gene %in% pgenes_filter)
  }
  
  launched <- ti_best %>% filter(!is.na(succ_3_a)) %>% mutate(success = succ_3_a)
  phase1 <- ti_best %>% filter(!is.na(succ_1_2))
  
  succ_gs <- sum(launched$gensup & launched$success, na.rm = TRUE)
  succ_nogs <- sum(!launched$gensup & launched$success, na.rm = TRUE)
  total_gs <- sum(phase1$gensup, na.rm = TRUE)
  total_nogs <- sum(!phase1$gensup, na.rm = TRUE)
  
  if (any(c(total_gs, total_nogs) == 0)) {
    return(tibble(est = NA, lwr.ci = NA, upr.ci = NA, n = ""))
  }
  
  out <- as.data.frame(BinomRatioCI(succ_gs, total_gs, succ_nogs, total_nogs, method = "katz"))
  out$n <- paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
  as_tibble(out)
}

cat("  Helper functions loaded.\n")
# ============================================================================
# SECTION 2: LOAD AND PREPARE DATA
# ============================================================================

message("[2/12] Loading and preparing data...")

# Load pQTL data
pqtl2 <- readRDS(file.path(data_raw, "pqtl_mrcoloc_2025.rds")) %>%
  filter(bxy_pval <= BONFERRONI_P)

# Load merge2 (genetic support data)
merge2 <- read_tsv(file.path(data_dir, "merge2.tsv.gz"), show_col_types = FALSE) %>%
  mutate(
    otg_study = if_else(assoc_source == "OTG",
                        str_remove(original_link, "https://genetics.opentargets.org/study/"),
                        NA_character_),
    otg_study = str_remove(otg_study, "FINNGEN_R6_"),
    key = paste0(gene, "_", otg_study)
  )

# Load indications
indic <- read_tsv(file.path(data_dir, "indic.tsv"), show_col_types = FALSE)

# Merge pQTL data with genetic support data
pqtl_cols <- c("nsnps", "cis_trans_mr", "bxy", "bxy_pval", "coloc_cis", "coloc_h4_cis",
               "snp_ciscoloc", "coloc_trans", "coloc_h4_trans", "snp_transcoloc")

merge2_with_pqtl <- merge2 %>% left_join(pqtl2, by = "key")

# Create separate pQTL rows
pqtl_rows <- merge2_with_pqtl %>%
  filter(!is.na(cis_trans_mr)) %>%
  mutate(original_link = "pqtl")

# Clean merged data
merge2_clean <- merge2_with_pqtl %>%
  mutate(across(all_of(pqtl_cols), ~ if_else(!is.na(cis_trans_mr), NA, .)))

# Combine
merge3_pqtl <- bind_rows(merge2_clean, pqtl_rows)

# Add therapeutic area annotation
area <- fread(file.path(data_dir, "areas.tsv"))
topl <- fread(file.path(data_dir, "indic_topl_match.tsv"))
ta_map <- merge(area, topl, by = "topl")
merge3_pqtl$therapeutic_area <- ta_map$area[match(merge3_pqtl$indication_mesh_id, ta_map$indication_mesh_id)]

# Load background genes (pgenes)
olink <- suppressMessages(read_sheet(
  "https://docs.google.com/spreadsheets/d/1DBHpr_Y3pFja4tMju3ZDJV8Gv-oTLq6wEuS0HtYjGbQ",
  sheet = "olink_complete"))
olink_genes <- olink %>%
  separate_rows(`Uniprot ID`, sep = ",") %>%
  mutate(hgnc_protein = mapIds(org.Hs.eg.db, keys = `Uniprot ID`, column = "SYMBOL",
                               keytype = "UNIPROT", multiVals = "first")) %>%
  filter(!is.na(hgnc_protein)) %>% pull(hgnc_protein) %>% unique()
otherttpairs <- readRDS(file.path(data_raw, "ttpairs_tested.rds"))
pgenes <- unique(c(olink_genes, unique(gsub("_.*", "", otherttpairs))))

cat(sprintf("  Background genes: %d\n", length(pgenes)))

# ============================================================================
# SECTION 3: DEFINE SUPPORT SETS
# ============================================================================

message("[3/12] Defining support sets...")

# pQTL support
df_pqtl_support <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
         gene %in% pgenes) %>%
  distinct(ti_uid)

# Any genetic support
df_any_genetic_support <- bind_rows(
  merge3_pqtl %>% filter(assoc_source == "OTG", comb_norm >= SIMILARITY_THR,
                         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE) %>% distinct(ti_uid),
  merge3_pqtl %>% filter(assoc_source == "OMIM", comb_norm >= SIMILARITY_THR) %>% distinct(ti_uid),
  merge3_pqtl %>% filter(assoc_source == "PICCOLO", comb_norm >= SIMILARITY_THR,
                         !is.na(pic_h4), pic_h4 >= 0.9) %>% distinct(ti_uid),
  merge3_pqtl %>% filter(assoc_source == "Genebass", comb_norm >= SIMILARITY_THR) %>% distinct(ti_uid)
) %>% distinct(ti_uid)

cat(sprintf("  pQTL-supported TIs: %d\n", nrow(df_pqtl_support)))
cat(sprintf("  Any genetic support TIs: %d\n", nrow(df_any_genetic_support)))

# ============================================================================
# SECTION 4: BUILD TI_BEST_BASE
# ============================================================================

message("[4/12] Building ti_best_base...")

ti_best_base <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "", !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat), gene %in% pgenes, !is.na(succ_1_2), !is.na(therapeutic_area)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2, TRUE ~ 1)) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup()

cat(sprintf("  ti_best_base: %d TIs\n", nrow(ti_best_base)))
# ============================================================================
# SECTION 5: THERAPEUTIC AREA HETEROGENEITY ANALYSIS (ST7-ST15, ST17)
# ============================================================================

message("[5/12] Running therapeutic area heterogeneity analysis...")

run_full_analysis <- function(ti_best, pqtl_uids, universe_label) {
  
  df <- ti_best %>%
    mutate(pqtl = ti_uid %in% pqtl_uids,
           launched = !is.na(succ_3_a) & succ_3_a)
  
  # Per-TA table
  ta_all <- df %>%
    group_by(therapeutic_area) %>%
    summarise(pqtl_total = sum(pqtl), nopqtl_total = sum(!pqtl),
              pqtl_success = sum(pqtl & launched), nopqtl_success = sum(!pqtl & launched),
              .groups = "drop") %>%
    mutate(pqtl_failure = pqtl_total - pqtl_success,
           nopqtl_failure = nopqtl_total - nopqtl_success)
  
  # Filter valid strata
  ta_valid <- ta_all %>%
    filter(pqtl_total > 0, nopqtl_total > 0,
           (pqtl_success + nopqtl_success) > 0,
           (pqtl_failure + nopqtl_failure) > 0)
  
  # Add RS and baseline metrics
  ta_table <- ta_valid %>%
    mutate(baseline_sr = nopqtl_success / nopqtl_total,
           pqtl_prevalence = pqtl_total / (pqtl_total + nopqtl_total)) %>%
    rowwise() %>%
    mutate(rs_tbl = list(rr_katz(pqtl_success, pqtl_total, nopqtl_success, nopqtl_total, add_half = TRUE))) %>%
    unnest(rs_tbl) %>% ungroup() %>%
    mutate(rs_str = sprintf("%.2f (%.2f-%.2f)", rr, lwr, upr),
           sparse = pqtl_total < 5) %>%
    arrange(desc(rr))
  
  # Overall RS
  overall <- ta_valid %>%
    summarise(pqtl_success = sum(pqtl_success), pqtl_total = sum(pqtl_total),
              nopqtl_success = sum(nopqtl_success), nopqtl_total = sum(nopqtl_total))
  overall_rs <- rr_katz(overall$pqtl_success, overall$pqtl_total,
                        overall$nopqtl_success, overall$nopqtl_total)
  
  # Breslow-Day test
  bd <- breslow_day_test(build_3d_array(ta_valid))
  
  # Spearman correlation
  spearman <- cor.test(ta_table$baseline_sr, ta_table$pqtl_prevalence, method = "spearman")
  
  # Leave-one-out
  loo <- map_dfr(ta_valid$therapeutic_area, function(drop_ta) {
    ta_sub <- ta_valid %>% filter(therapeutic_area != drop_ta)
    pooled <- ta_sub %>%
      summarise(pqtl_success = sum(pqtl_success), pqtl_total = sum(pqtl_total),
                nopqtl_success = sum(nopqtl_success), nopqtl_total = sum(nopqtl_total))
    rs_tbl <- rr_katz(pooled$pqtl_success, pooled$pqtl_total,
                      pooled$nopqtl_success, pooled$nopqtl_total)
    tibble(dropped_ta = drop_ta,
           pqtl = sprintf("%d/%d", pooled$pqtl_success, pooled$pqtl_total),
           nopqtl = sprintf("%d/%d", pooled$nopqtl_success, pooled$nopqtl_total),
           rs = rs_tbl$rr, lwr = rs_tbl$lwr, upr = rs_tbl$upr,
           rs_str = sprintf("%.2f (%.2f-%.2f)", rs_tbl$rr, rs_tbl$lwr, rs_tbl$upr))
  })
  
  loo_full <- bind_rows(
    tibble(dropped_ta = "(none - overall)",
           pqtl = sprintf("%d/%d", overall$pqtl_success, overall$pqtl_total),
           nopqtl = sprintf("%d/%d", overall$nopqtl_success, overall$nopqtl_total),
           rs = overall_rs$rr, lwr = overall_rs$lwr, upr = overall_rs$upr,
           rs_str = sprintf("%.2f (%.2f-%.2f)", overall_rs$rr, overall_rs$lwr, overall_rs$upr)),
    loo
  )
  
  # Chi-square test
  ta_with_expected <- ta_table %>%
    mutate(
      ta_total = pqtl_total + nopqtl_total,
      expected_pqtl = ta_total * sum(pqtl_total) / sum(ta_total),
      obs_exp_ratio = pqtl_total / expected_pqtl,
      enrichment = case_when(
        obs_exp_ratio > 1.5 ~ "enriched",
        obs_exp_ratio < 0.67 ~ "depleted",
        TRUE ~ "as expected"
      )
    )
  
  chisq_result <- tryCatch(
    chisq.test(x = ta_table$pqtl_total, 
               p = (ta_table$pqtl_total + ta_table$nopqtl_total) / sum(ta_table$pqtl_total + ta_table$nopqtl_total)),
    error = function(e) list(statistic = NA, parameter = NA, p.value = NA)
  )
  
  list(
    universe = universe_label,
    ta_table = ta_table %>%
      select(therapeutic_area, pqtl_total, pqtl_success, nopqtl_total, nopqtl_success, 
             rr, lwr, upr, rs_str, baseline_sr, pqtl_prevalence, sparse) %>%
      arrange(desc(pqtl_total)),
    ta_distribution = ta_with_expected %>%
      select(therapeutic_area, ta_total, pqtl_total, expected_pqtl, obs_exp_ratio, enrichment) %>%
      arrange(desc(obs_exp_ratio)),
    chisq = tibble(universe = universe_label, 
                   chisq = as.numeric(chisq_result$statistic), 
                   df = as.numeric(chisq_result$parameter), 
                   p = chisq_result$p.value),
    bd = tibble(universe = universe_label, bd_chisq = bd$statistic, bd_df = bd$df, bd_p = bd$p.value),
    spearman = tibble(universe = universe_label, rho = spearman$estimate, p = spearman$p.value),
    loo = loo_full
  )
}

# Run both universes
ti_best_u1 <- ti_best_base
ti_best_u2 <- ti_best_base %>% filter(ti_uid %in% df_any_genetic_support$ti_uid)

results_u1 <- run_full_analysis(ti_best_u1, df_pqtl_support$ti_uid,
                                "UNIVERSE 1: All Phase-I-entered TIs (measured proteins)")
results_u2 <- run_full_analysis(ti_best_u2, df_pqtl_support$ti_uid,
                                "UNIVERSE 2: Phase-I-entered + ANY genetic evidence")

# Extract successful pQTL+ TIs (ST17)
ST17 <- ti_best_base %>%
  filter(ti_uid %in% df_pqtl_support$ti_uid,
         !is.na(succ_3_a) & succ_3_a) %>%
  select(gene, indication_mesh_id, indication_mesh_term, therapeutic_area, ti_uid) %>%
  arrange(gene, indication_mesh_id)

cat(sprintf("  Successful pQTL+ TIs (ST17): %d\n", nrow(ST17)))

# Assign tables
ST7  <- bind_rows(results_u1$chisq, results_u2$chisq)
ST8  <- results_u1$ta_distribution
ST9  <- results_u2$ta_distribution
ST10 <- bind_rows(results_u1$spearman, results_u2$spearman)
ST11 <- results_u1$ta_table
ST12 <- results_u2$ta_table
ST13 <- bind_rows(results_u1$bd, results_u2$bd)
ST14 <- results_u1$loo
ST15 <- results_u2$loo

cat("  ST7-ST15, ST17 complete.\n")
# ============================================================================
# SECTION 6: BUILD ST16 - ALL MR TARGET-TRAIT PAIRS
# ============================================================================

message("[6/12] Building ST16 (All MR target-trait pairs)...")

# Load new UKBPPP MR-coloc results
new <- readRDS(file.path(data_raw, "ukb_ppp_mr_coloc_results.rds")) %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  mutate(
    outcome = coalesce(accession, trait),
    hgnc_protein = str_extract(prot, "^[^_]+"),
    outcome_trait_efo = ifelse(grepl("^GCST", accession), gsub("^.*?-.*?-(.*)$", "\\1", trait), NA),
    Data = "UKBPPP_2023"
  ) %>%
  rename(
    cis_trans_mr = IVs, outcome_trait = pheno,
    coloc_h1_cis = pp1_cis, coloc_h2_cis = pp2_cis, coloc_h3_cis = pp3_cis, 
    coloc_h4_cis = pp4_cis, coloc_h4_h3_cis = pp43_cis,
    coloc_h1_trans = pp1_trans, coloc_h2_trans = pp2_trans, coloc_h3_trans = pp3_trans,
    coloc_h4_trans = pp4_trans, coloc_h4_h3_trans = pp43_trans
  ) %>%
  select(hgnc_protein, outcome, outcome_trait, outcome_trait_efo, Data, nsnps, 
         cis_trans_mr, starts_with("bxy"), starts_with("coloc_"), starts_with("snp_"), 
         Samples, Cases) %>%
  distinct()

new$pav_cismr <- NA
new$ensid <- genes$id[match(new$hgnc_protein, genes$approvedSymbol)]

olink_full <- read_sheet("https://docs.google.com/spreadsheets/d/1DBHpr_Y3pFja4tMju3ZDJV8Gv-oTLq6wEuS0HtYjGbQ", 
                         sheet = "olink_complete_with_more_fields")

olink2 <- olink_full %>%
  separate_rows(`UniProt ID`, sep = ",") %>%
  mutate(hgnc = unlist(mapIds(org.Hs.eg.db, `UniProt ID`, "SYMBOL", "UNIPROT", multiVals="first"))) %>%
  filter(!is.na(hgnc))

missgenes <- new %>% filter(is.na(ensid)) %>% distinct(hgnc_protein) %>%
  mutate(hgnc2 = as.character(olink2$hgnc[match(hgnc_protein, olink_full$`Gene name`)]),
         hgnc3 = protein_map[hgnc_protein],
         final = coalesce(hgnc3, hgnc2))

lookup <- setNames(missgenes$final, missgenes$hgnc_protein)
idx <- match(new$hgnc_protein, names(lookup))
repl <- lookup[idx]
new$hgnc_protein[!is.na(repl) & repl != new$hgnc_protein] <- repl[!is.na(repl) & repl != new$hgnc_protein]
new$ensid <- genes$id[match(new$hgnc_protein, genes$approvedSymbol)]

new <- compute_mantissa_exponent(new, c("bxy_pval", "bxy_pval_egger", "bxy_pval_median"))

# Load legacy MR-coloc results
df1 <- readRDS(file.path(data_raw, "mr_prot_filtered_dataset_v1_v2.rds"))
df2 <- readRDS(file.path(data_raw, "mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds"))

old <- df1 %>%
  left_join(df2 %>% select(exp_out_gsmr_coloc, coloc_h1_trans:coloc_h4_h3_trans, SNP_trans_coloc), 
            by = "exp_out_gsmr_coloc")

# Recover missing gene symbols/ensids for legacy data
missgenes_old <- old %>% 
  filter(is.na(ensid) | is.na(hgnc_protein)) %>% 
  select(protein_trait, ensid, hgnc_protein) %>% 
  distinct()

missgenes_old <- missgenes_old %>%
  mutate(
    hgnc_protein = coalesce(
      hgnc_protein,
      genes$approvedSymbol[match(protein_trait, genes$approvedSymbol)],
      genes$approvedSymbol[match(ensid, genes$id)]
    ),
    ensid = coalesce(
      ensid,
      genes$id[match(protein_trait, genes$approvedSymbol)]
    )
  )

old <- old %>%
  left_join(missgenes_old[, c("protein_trait", "hgnc_protein", "ensid")],
            by = "protein_trait", suffix = c("", ".mg")) %>%
  mutate(
    hgnc_protein = coalesce(hgnc_protein, hgnc_protein.mg,
                            missgenes_old$hgnc_protein[match(ensid, missgenes_old$ensid)]),
    ensid = coalesce(ensid, ensid.mg)
  ) %>%
  select(-ends_with(".mg")) %>%
  filter(!is.na(hgnc_protein))

old <- old %>%
  rename(nsnps=nsnp, Cases=n_cases, Samples=n_initial,
         bxy_egger=mr_egger, bxy_se_egger=mr_egger_se, bxy_pval_egger=mr_egger_p,
         bxy_median=mr_wm, bxy_se_median=mr_wm_se, bxy_pval_median=mr_wm_p,
         bxy_pval_egger_mantissa=mr_egger_p_mantissa, bxy_pval_egger_exponent=mr_egger_p_exponent,
         bxy_pval_median_mantissa=mr_wm_p_mantissa, bxy_pval_median_exponent=mr_wm_p_exponent,
         snp_ciscoloc=varid_left, snp_transcoloc=SNP_trans_coloc,
         coloc_h1_cis=coloc_h1, coloc_h2_cis=coloc_h2, coloc_h3_cis=coloc_h3, 
         coloc_h4_cis=coloc_h4, coloc_h4_h3_cis=coloc_h4_h3) %>%
  select(Data, hgnc_protein, ensid, outcome, outcome_trait, outcome_trait_efo, cis_trans_mr, nsnps,
         starts_with("bxy"), -bxy_lci, -bxy_uci, starts_with("coloc_"), -coloc_n_vars,
         starts_with("snp_"), Samples, Cases, pav_cismr) %>%
  mutate(cis_trans_mr = case_when(cis_trans_mr=="cis"~"Cis", cis_trans_mr=="trans"~"Trans", 
                                  cis_trans_mr=="mixed"~"Mixed", TRUE~cis_trans_mr))

# Harmonize and combine datasets
harm <- harmonize_column_types(old, new)
comb <- bind_rows(harm$df1, harm$df2) %>%
  filter(bxy_pval <= BONFERRONI_P) %>%
  arrange(bxy_pval_exponent, bxy_pval_mantissa)

# Add trait keys
panukb <- readRDS(file.path(data_raw, "panukb.rds"))
panukb[panukb == "NA"] <- NA
panukb$trait_code <- sub("\\.tsv\\.bgz$", "", panukb$filename)

pos <- match(comb$outcome, panukb$trait_code)
comb$outcome_trait2 <- ifelse(is.na(panukb$trait_efo_terms[pos]) | panukb$trait_efo_terms[pos]=="NA",
                              NA, panukb$trait_efo_terms[pos])
comb$outcome_trait_efo2 <- ifelse(is.na(panukb$trait_efos[pos]) | panukb$trait_efos[pos]=="NA",
                                  NA, panukb$trait_efos[pos])
comb$outcome_trait2 <- ifelse(!is.na(comb$outcome_trait_efo2) & is.na(comb$outcome_trait2),
                              panukb$coding_description[pos], comb$outcome_trait2)

mesh <- fread(file.path(data_dir, "assoc.tsv.gz"))
mesh$otg_study <- ""
mesh$otg_study[mesh$source=="OTG"] <- gsub("https://genetics.opentargets.org/study/", "",
                                           mesh$original_link[mesh$source=="OTG"])
mesh$otg_study <- gsub("FINNGEN_R6_", "", mesh$otg_study)

pos <- match(comb$outcome, mesh$otg_study)
comb$mesh_id <- mesh$mesh_id[pos]
comb$mesh_term <- mesh$mesh_term[pos]

comb$outcome_trait[grepl("pheno 48 / pheno 49", comb$outcome_trait)] <- "waist-hip ratio"
comb$trait_key_term <- coalesce(comb$mesh_term, comb$outcome_trait)
comb$trait_key_term <- ifelse(!is.na(comb$outcome_trait2) & grepl("categorical-", comb$outcome),
                              comb$outcome_trait2, comb$trait_key_term)
comb$trait_key <- coalesce(comb$mesh_id, comb$outcome_trait_efo, comb$outcome_trait_efo2, comb$trait_key_term)
comb$ttpair <- paste(comb$hgnc_protein, comb$trait_key, sep = "_")

# Remove irrelevant traits
toremove <- read_sheet("https://docs.google.com/spreadsheets/d/1TDz8oRI5H-DMHOTs0bZgm2dcisYeuvdSCybn4NcHESw", 
                       sheet = "v4") %>% filter(to_remove == "Y")
comb <- comb[!comb$trait_key_term %in% toremove$trait_key_term, ]

cat(sprintf("  Unique target-trait pairs: %d\n", n_distinct(comb$ttpair)))

# Trans-coloc genes
tbl <- readRDS(file.path(data_raw, "trans_genes.rds"))
comb$trans_coloc_genes <- tbl$nearest_gene_symbols[match(comb$snp_transcoloc, tbl$variantId)]

# Replication & multi-method support
comb_summary <- comb %>%
  group_by(ttpair) %>%
  summarise(
    replication = n_distinct(Data) > 1 | n_distinct(outcome) > 1,
    cis_mr = any(cis_trans_mr=="Cis"), trans_mr = any(cis_trans_mr=="Trans"), mixed_mr = any(cis_trans_mr=="Mixed"),
    cis_coloc = any(coloc_h4_cis >= 0.8, na.rm=TRUE), trans_coloc = any(coloc_h4_trans >= 0.8, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(mr_coloc_types = pmap_chr(list(cis_mr, trans_mr, mixed_mr, cis_coloc, trans_coloc), ~{
    labs <- c("cis-MR","trans-MR","mixed-MR","cis-coloc","trans-coloc")[c(...)]
    if (length(labs)==0) NA_character_ else paste(labs, collapse=", ")
  }), multi_method_count = ifelse(is.na(mr_coloc_types), 0L, str_count(mr_coloc_types, ",")+1L)) %>%
  select(ttpair, replication, mr_coloc_types, multi_method_count)

comb2 <- comb %>% left_join(comb_summary, by = "ttpair")

# Triangulation
efo <- comb2 %>% distinct(outcome_trait_efo, outcome_trait_efo2, ensid, ttpair)
split_ids <- function(x) { if(is.na(x)) return(NA_character_); str_trim(unlist(str_split(str_remove_all(str_remove_all(x, '^c\\(|\\)$'), '"'), ","))) }
efo_clean <- efo %>%
  mutate(all_ids = map2(outcome_trait_efo, outcome_trait_efo2, ~{ ids <- c(split_ids(.x), split_ids(.y)); ids[!is.na(ids) & ids != ""] })) %>%
  unnest(all_ids) %>% rename(efo_id = all_ids) %>%
  mutate(key = paste0(ensid, "_", efo_id))

# Extract unique ensid and efo values for triangulate_ot.R
ensid <- unique(efo_clean$ensid)
efo <- unique(efo_clean$efo_id)

source(file.path(project_root, "R/triangulate_ot.R"))

efo_clean$triangulation <- tr$concatenatedDatatypeIds[match(efo_clean$key, tr$key)]
efo_clean$harmonic_genetic_score <- tr$harmonic_genetic_score[match(efo_clean$key, tr$key)]

comb2$triangulation <- efo_clean$triangulation[match(comb2$ttpair, efo_clean$ttpair)]
comb2$harmonic_genetic_score <- efo_clean$harmonic_genetic_score[match(comb2$ttpair, efo_clean$ttpair)]

# Drug target matching
merge3 <- merge2 %>% filter(comb_norm >= 0.8) %>% select(indication_mesh_term, ccat, key, gene, l2g_share) %>% distinct()

comb2$drug_target_pp <- ifelse(comb2$hgnc_protein %in% unique(merge3$gene), "yes", "no")
comb2$ppkey <- paste0(comb2$hgnc_protein, "_", comb2$outcome)
pos <- match(comb2$ppkey, merge3$key)
comb2$drug_target_pp_indication_match <- merge3$indication_mesh_term[pos]
comb2$drug_target_pp_highest_phase_indication_match <- merge3$ccat[pos]
comb2$l2g_share <- merge3$l2g_share[pos]

chembl <- readRDS(file.path(data_raw, "chembl.rds"))
chembl_ttpair <- chembl %>% group_by(targetId, diseaseFromSourceMappedId) %>%
  summarise(chembl_highest_phase = max(clinicalPhase), diseaseFromSource = first(diseaseFromSource), .groups="drop") %>%
  mutate(key = paste0(targetId, "_", diseaseFromSourceMappedId))

comb2$drug_target_chembl <- ifelse(comb2$ensid %in% unique(chembl$targetId), "yes", "no")
pos <- match(efo_clean$key, chembl_ttpair$key)
efo_clean$drug_target_chembl_indication_match <- chembl_ttpair$diseaseFromSource[pos]
efo_clean$drug_target_chembl_highest_phase_indication_match <- chembl_ttpair$chembl_highest_phase[pos]
comb2$drug_target_chembl_indication_match <- efo_clean$drug_target_chembl_indication_match[match(comb2$ttpair, efo_clean$ttpair)]
comb2$drug_target_chembl_highest_phase_indication_match <- efo_clean$drug_target_chembl_highest_phase_indication_match[match(comb2$ttpair, efo_clean$ttpair)]

comb2$positive_control <- ifelse(!is.na(comb2$drug_target_pp_indication_match) | !is.na(comb2$drug_target_chembl_indication_match), "yes", "no")
comb2$repositioning_opportunity <- ifelse((comb2$drug_target_pp=="yes" | comb2$drug_target_chembl=="yes") & comb2$positive_control=="no", "yes", "no")

# HLA annotation
hla_chr <- "6"
hla_start <- 25000000
hla_end <- 34000000

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = unique(comb2$hgnc_protein),
  mart = ensembl
)

gene_coords$hla <- with(gene_coords, 
                        ifelse(chromosome_name == hla_chr & 
                                 start_position >= hla_start & 
                                 end_position <= hla_end, 
                               "yes", "no"))

hla_lookup <- gene_coords %>%
  select(hgnc_symbol, hla) %>%
  group_by(hgnc_symbol) %>%
  summarise(hla = ifelse(any(hla == "yes"), "yes", "no"), .groups = "drop")

comb2 <- comb2 %>%
  left_join(hla_lookup, by = c("hgnc_protein" = "hgnc_symbol")) %>%
  mutate(hla = ifelse(is.na(hla), "no", hla))

# --- TOP-LD Functions ---
to_topld_id <- function(x) {
  x <- str_trim(x)
  if (str_detect(x, "^chr[0-9XYM]+:[0-9]+:[ACGT]+:[ACGT]+$")) return(x)
  parts <- str_split(x, "[_:]", simplify = TRUE)
  if (ncol(parts) != 4) return(NA_character_)
  paste0("chr", parts[1], ":", parts[2], ":", parts[3], ":", parts[4])
}

run_topld_chunks <- function(snp_ids, topld_bin, topld_dir, chunk_size=200, 
                             ld_threshold=0.6, population="EUR", maf=0.01, rerun=FALSE) {
  n_chunks <- ceiling(length(snp_ids) / chunk_size)
  ld_list <- vector("list", n_chunks)
  old_wd <- getwd(); setwd(topld_dir); on.exit(setwd(old_wd))
  
  for (i in seq_len(n_chunks)) {
    output_fn <- paste0("outputLD_chunk_", i, ".txt")
    if (file.exists(output_fn) && !rerun) {
      tryCatch({ ld_list[[i]] <- fread(output_fn, data.table=FALSE) }, error=function(e) NULL)
    } else {
      idx_start <- (i-1)*chunk_size + 1; idx_end <- min(i*chunk_size, length(snp_ids))
      writeLines(snp_ids[idx_start:idx_end], paste0("input_chunk_", i, ".txt"))
      system2(topld_bin, c("-thres", ld_threshold, "-pop", population, "-maf", maf,
                           "-inFile", paste0("input_chunk_", i, ".txt"),
                           "-outputLD", output_fn, "-outputInfo", paste0("outputInfo_chunk_", i, ".txt")))
      if (file.exists(output_fn)) ld_list[[i]] <- fread(output_fn, data.table=FALSE)
    }
  }
  bind_rows(ld_list[!sapply(ld_list, is.null)])
}

extract_pav_summary <- function(ld_data) {
  if (is.null(ld_data) || nrow(ld_data) == 0) return(NULL)
  pav_terms <- c("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained",
                 "frameshift_variant","stop_lost","start_lost","transcript_amplification",
                 "feature_elongation","feature_truncation","inframe_insertion","inframe_deletion",
                 "missense_variant","protein_altering_variant")
  
  ld_data %>%
    mutate(index_snp = paste0("chr", CHROM, ":", `POS1(gb38)`, ":", REF1, ":", ALT1),
           gene2 = `marker2-GeneName`, cons2 = `marker2-LocationRelativeToGene`) %>%
    separate_rows(gene2, cons2, sep = "\\|") %>%
    separate_rows(cons2, sep = ",") %>%
    mutate(gene2 = str_trim(gene2), cons2 = str_trim(cons2)) %>%
    filter(!is.na(gene2), gene2 != "", !is.na(cons2), cons2 != "") %>%
    distinct(index_snp, gene = gene2, so_term = cons2) %>%
    group_by(index_snp, gene) %>%
    summarise(pav_flag = any(so_term %in% pav_terms),
              pav_terms = paste(sort(unique(so_term[so_term %in% pav_terms])), collapse=";"), 
              .groups="drop") %>%
    mutate(pav_terms = na_if(pav_terms, ""))
}

# Run PAV annotation
topld_dir <- file.path(project_root, "topld_api")
topld_bin <- file.path(topld_dir, "topld_api")

cis_snps <- comb %>% filter(cis_trans_mr=="Cis", !is.na(snp_ciscoloc), snp_ciscoloc!="", !grepl("mapp", snp_ciscoloc)) %>%
  distinct(snp_ciscoloc) %>% pull()

cis_snp_map <- tibble(index_snp_internal = cis_snps) %>%
  mutate(index_snp_topld = map_chr(index_snp_internal, to_topld_id)) %>%
  filter(!is.na(index_snp_topld))

if (nrow(cis_snp_map) > 0) {
  ld_all <- run_topld_chunks(unique(cis_snp_map$index_snp_topld), topld_bin, topld_dir)
  cis_pav_summary <- extract_pav_summary(ld_all)
  
  if (!is.null(cis_pav_summary)) {
    comb2 <- comb2 %>%
      left_join(cis_snp_map, by = c("snp_ciscoloc" = "index_snp_internal")) %>%
      left_join(cis_pav_summary, by = c("index_snp_topld" = "index_snp", "hgnc_protein" = "gene")) %>%
      mutate(pav_cismr = case_when(cis_trans_mr=="Cis" & !is.na(pav_flag) & pav_flag ~ "yes",
                                   cis_trans_mr=="Cis" ~ "no", TRUE ~ NA_character_))
  }
}

# Ensure columns exist even if no PAV data
if (!"pav_terms" %in% names(comb2)) comb2$pav_terms <- NA_character_
if (!"pav_cismr" %in% names(comb2)) comb2$pav_cismr <- NA_character_

# Select best row per target-trait pair
ST16 <- comb2 %>%
  mutate(coloc_best = pmax(coloc_h4_cis, coloc_h4_trans, na.rm=TRUE),
         has_p = !is.na(bxy_pval_exponent) & !is.na(bxy_pval_mantissa)) %>%
  group_by(ttpair) %>%
  mutate(
    snp_ciscoloc_all = toString(unique(na.omit(ifelse(coloc_h4_cis >= 0.8, snp_ciscoloc, NA)))),
    snp_transcoloc_all = toString(unique(na.omit(ifelse(coloc_h4_trans >= 0.8, snp_transcoloc, NA)))),
    trans_coloc_genes_all = toString(unique(na.omit(ifelse(!is.na(trans_coloc_genes), paste0(snp_transcoloc, " (", trans_coloc_genes, ")"), NA)))),
    data_sources = toString(unique(Data)), outcomes = toString(unique(outcome)),
    pav_cismr_any = any(pav_cismr == "yes", na.rm=TRUE), pav_terms_all = toString(unique(na.omit(pav_terms)))
  ) %>%
  arrange(desc(has_p), bxy_pval_exponent, bxy_pval_mantissa, desc(coalesce(coloc_best, -Inf))) %>%
  slice_head(n = 1) %>% ungroup() %>%
  mutate(pav_cismr = case_when(pav_cismr_any ~ "yes", !is.na(pav_cismr) ~ "no", TRUE ~ NA_character_)) %>%
  dplyr::select(hgnc_protein, ttpair, data_sources, outcomes, outcome_trait, trait_key_term,
                starts_with("bxy"), snp_ciscoloc_all, snp_transcoloc_all, trans_coloc_genes_all,
                pav_cismr, pav_terms_all, l2g_share, hla, harmonic_genetic_score, replication, triangulation, multi_method_count, mr_coloc_types,
                positive_control, repositioning_opportunity)

cat(sprintf("  ST16 (All MR pairs): %d rows\n", nrow(ST16)))
# ============================================================================
# SECTION 7: BUILD ST1-ST3 (Simple tables)
# ============================================================================

message("[7/12] Building ST1-ST3...")

# --- ST1: Outcome GWAS ---
old_st2 <- read_sheet("https://docs.google.com/spreadsheets/d/1xQ9ojspUea7stOTXSO_vZbNPPg_BU_y2pyM2SigV7UI", 
                      sheet = "ST1 - Outcome_GWAS", skip = 1) %>%
  transmute(outcome_datasets = as.character(outcome_datasets), outcome_trait = as.character(outcome_trait),
            outcome_trait_efo = as.character(outcome_trait_efo), outcome_data_source = as.character(outcome_data_source))

new_st2 <- read_sheet("https://docs.google.com/spreadsheets/d/1aRTkoyPUrV9jOrUQAhvXl-38QYS4CE0G0aHlhry9LKg") %>%
  filter(!is.na(outcome_datasets)) %>%
  mutate(outcome_data_source = case_when(str_detect(outcome_datasets, "^GCST") ~ "GWAS Catalog",
                                         str_detect(outcome_datasets, "^FINNGEN") ~ "FinnGen", TRUE ~ "pan-UK Biobank"))

ST1 <- bind_rows(old_st2, new_st2 %>% select(any_of(names(old_st2)))) %>%
  mutate(outcome_datasets = str_replace(outcome_datasets, "^FINNGEN_R[0-9]+_", "")) %>%
  group_by(outcome_datasets, outcome_data_source) %>%
  summarise(across(c(outcome_trait, outcome_trait_efo), ~first(na.omit(.))), .groups = "drop")

# --- ST2: Proteomic GWAS ---
ST2 <- tribble(
  ~proteomic_dataset, ~platform, ~sample_size, ~`no. of proteins measured`,
  "OLLI_2017",        "BioRad",        8293,    41,
  "FOLKERSEN_2017",   "Olink",         3394,    83,
  "HILLARY_2019",     "Olink",          750,    92,
  "SCALLOP_2020",     "Olink",        30931,    90,
  "SUN_2018",         "Somalogic",     3300,  3622,
  "SUHRE_2017",       "Somalogic",     1000,  1124,
  "PIETZNER_2020",    "Somalogic",    10708,   186,
  "UKBPPP_2023",      "Olink",        54219,  2923
)

# --- ST3: Excluded Traits ---
ST3 <- toremove %>% select(trait_key_term)

cat(sprintf("  ST1 (Outcome GWAS): %d rows\n", nrow(ST1)))
cat(sprintf("  ST2 (Proteomic GWAS): %d rows\n", nrow(ST2)))
cat(sprintf("  ST3 (Excluded Traits): %d rows\n", nrow(ST3)))

# ============================================================================
# SECTION 8: BUILD ST4-ST6 (Enrichment tables)
# ============================================================================

message("[8/12] Computing ST4-ST6 (enrichment tables)...")

# --- ST4: Enrichment by source ---
# pQTL support
pqtl_support_st4 <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE
  ) %>%
  distinct(ti_uid)

ti_best_pqtl <- get_ti_best(merge3_pqtl, pqtl_support_st4, indic)
rs_pqtl <- compute_rs(ti_best_pqtl, pgenes) %>%
  mutate(source = "pQTL", label = "By genetic evidence source")

# Other genetic sources
sources <- c("OMIM" = "omim", "All OTG" = "otg", "Genebass" = "genebass", 
             "FinnGen" = "finngen", "GWAS catalog" = "gwascat")

rs_other <- map_dfr(names(sources), function(src_name) {
  src <- sources[src_name]
  
  support <- merge3_pqtl %>%
    filter(!is.na(gene), !is.na(indication_mesh_id), !is.na(ccat)) %>%
    left_join(indic %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
    filter(genetic_insight != "none") %>%
    filter(case_when(
      src == "omim" ~ assoc_source == "OMIM" & comb_norm >= SIMILARITY_THR,
      src == "otg" ~ assoc_source == "OTG" & comb_norm >= SIMILARITY_THR & 
        !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE,
      src == "genebass" ~ assoc_source == "Genebass" & comb_norm >= SIMILARITY_THR,
      src == "finngen" ~ assoc_source == "OTG" & grepl("FINNGEN", original_link) & 
        comb_norm >= SIMILARITY_THR & !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE,
      src == "gwascat" ~ assoc_source == "OTG" & grepl("GCST", original_link) & 
        comb_norm >= SIMILARITY_THR & !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE
    )) %>%
    distinct(ti_uid)
  
  ti_best <- get_ti_best(merge3_pqtl, support, indic)
  compute_rs(ti_best) %>% mutate(source = src_name, label = "By genetic evidence source")
})

# L2G thresholds only
source(file.path(project_root, "R/pipeline_best.R"))
source(file.path(project_root, "R/advancement_rr.R"))

rs_l2g <- map_dfr(c(0.25, 0.50, 0.75), function(thr) {
  pb <- pipeline_best(merge2, phase = "combined", basis = "ti", associations = "OTG",
                      share_mode = "L2G", min_share = thr, verbose = FALSE)
  rr <- advancement_rr(pb) %>% filter(phase == "I-Launch")
  tibble(
    est = rr$rs_mean, lwr.ci = rr$rs_l, upr.ci = rr$rs_u,
    source = paste0("L2G share: >= ", thr),
    n = paste0("(", rr$x_yes, "/", rr$n_yes, ")/(", rr$x_no, "/", rr$n_no, ")"),
    label = "By L2G share only"
  )
})

# pQTL + L2G thresholds
rs_pqtl_l2g <- map_dfr(c(0.25, 0.50, 0.75), function(thr) {
  support <- merge3_pqtl %>%
    filter(
      grepl("pqtl", original_link, ignore.case = TRUE),
      comb_norm >= SIMILARITY_THR,
      !is.na(l2g_share), l2g_share >= thr,
      gene %in% pgenes
    ) %>%
    distinct(ti_uid)
  
  ti_best <- get_ti_best(merge3_pqtl %>% filter(gene %in% pgenes), support, indic)
  compute_rs(ti_best) %>%
    mutate(source = paste0("L2G share: >= ", thr), label = "By pQTL + L2G share")
})

# MR types
rs_mrtype <- map_dfr(c("Cis", "Trans", "Mixed"), function(mrtype) {
  support <- merge3_pqtl %>%
    filter(
      grepl("pqtl", original_link, ignore.case = TRUE),
      comb_norm >= SIMILARITY_THR,
      !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
      cis_trans_mr == mrtype,
      gene %in% pgenes
    ) %>%
    distinct(ti_uid)
  
  ti_best <- get_ti_best(merge3_pqtl %>% filter(gene %in% pgenes), support, indic)
  compute_rs(ti_best) %>%
    mutate(source = paste0(mrtype, "-MR"), label = "By pQTL MR-coloc type")
})

# MR + coloc
support_coloc <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
    (coloc_h4_cis >= COLOC_H4_THR | coloc_h4_trans >= COLOC_H4_THR),
    gene %in% pgenes
  ) %>%
  distinct(ti_uid)

rs_mrcoloc <- compute_rs(get_ti_best(merge3_pqtl %>% filter(gene %in% pgenes), support_coloc, indic)) %>%
  mutate(source = "MR+coloc", label = "By pQTL MR-coloc type")

# Combine all for ST4
ST4 <- bind_rows(rs_l2g, rs_pqtl_l2g, rs_pqtl, rs_other, rs_mrtype, rs_mrcoloc) %>%
  mutate(subgroup = row_number()) %>%
  select(panel_group = label, source_label = source, rs_estimate = est,
         rs_lwr_95ci = lwr.ci, rs_upr_95ci = upr.ci, count_string = n, row_order = subgroup)

cat(sprintf("  ST4 (Figure1a): %d rows\n", nrow(ST4)))

# --- ST5: UpSet plot data ---
pqtl_launched <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
    succ_3_a == TRUE
  ) %>%
  distinct(ti_uid, gene, indication_mesh_id, indication_mesh_term, 
           cis_trans_mr, coloc_h4_cis, coloc_h4_trans)

ST5 <- pqtl_launched %>%
  mutate(target_indication_id = paste0(gene, "-", indication_mesh_id)) %>%
  group_by(target_indication_id, gene, indication_mesh_id, indication_mesh_term) %>%
  summarise(
    cis_mr = as.integer(any(cis_trans_mr == "Cis")),
    trans_mr = as.integer(any(cis_trans_mr == "Trans")),
    mixed_mr = as.integer(any(cis_trans_mr == "Mixed")),
    coloc_cis = as.integer(any(coloc_h4_cis >= COLOC_H4_THR, na.rm = TRUE)),
    coloc_trans = as.integer(any(coloc_h4_trans >= COLOC_H4_THR, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(n_methods_supported = cis_mr + trans_mr + mixed_mr + coloc_cis + coloc_trans) %>%
  arrange(desc(n_methods_supported), gene)

cat(sprintf("  ST5 (Figure1b): %d rows\n", nrow(ST5)))

# --- ST6: Gene family enrichment ---
gene_list_dir <- file.path(project_root, "genetic_support-main/data/gene_lists")
gene_families <- list.files(gene_list_dir, pattern = "\\.tsv$", full.names = TRUE) %>%
  map_dfr(~ read_tsv(.x, col_names = "gene", show_col_types = FALSE) %>%
            mutate(family = tools::file_path_sans_ext(basename(.x)))) %>%
  distinct(gene, family)

ti_best_fam <- merge3_pqtl %>%
  filter(!is.na(gene), !is.na(indication_mesh_id), !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(
    highest_phase = case_when(
      succ_3_a ~ 3,
      succ_2_3 ~ 2,
      succ_1_2 ~ 1,
      succ_p_1 ~ 0
    ),
    pqtl_support = grepl("pqtl", original_link, ignore.case = TRUE) & 
      comb_norm >= SIMILARITY_THR & 
      !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE
  ) %>%
  arrange(ti_uid, desc(pqtl_support), desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup() %>%
  inner_join(gene_families, by = "gene", relationship = "many-to-many")

# L2G only enrichment
fam_l2g <- ti_best_fam %>%
  mutate(supported = !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE) %>%
  group_by(family, supported) %>%
  summarise(
    x = sum(!is.na(succ_3_a) & succ_3_a),
    n = sum(!is.na(succ_1_2)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = supported, values_from = c(x, n), values_fill = 0) %>%
  rowwise() %>%
  mutate(
    rs = list(as.data.frame(BinomRatioCI(x_TRUE, n_TRUE, x_FALSE, n_FALSE, method = "katz"))),
    est = rs$est, lwr.ci = rs$lwr.ci, upr.ci = rs$upr.ci
  ) %>%
  ungroup() %>%
  mutate(
    analysis = "L2G only",
    n_str = paste0("(", x_TRUE, "/", n_TRUE, ") / (", x_FALSE, "/", n_FALSE, ")")
  )

# pQTL enrichment
fam_pqtl <- ti_best_fam %>%
  filter(gene %in% pgenes) %>%
  mutate(supported = grepl("pqtl", original_link, ignore.case = TRUE)) %>%
  group_by(family, supported) %>%
  summarise(
    x = sum(!is.na(succ_3_a) & succ_3_a),
    n = sum(!is.na(succ_1_2)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = supported, values_from = c(x, n), values_fill = 0) %>%
  rowwise() %>%
  mutate(
    rs = list(as.data.frame(BinomRatioCI(x_TRUE, n_TRUE, x_FALSE, n_FALSE, method = "katz"))),
    est = rs$est, lwr.ci = rs$lwr.ci, upr.ci = rs$upr.ci
  ) %>%
  ungroup() %>%
  mutate(
    analysis = "L2G + pQTL",
    n_str = paste0("(", x_TRUE, "/", n_TRUE, ") / (", x_FALSE, "/", n_FALSE, ")")
  )

ST6 <- bind_rows(fam_l2g, fam_pqtl) %>%
  select(analysis_type = analysis, family, rs_estimate = est,
         rs_lwr_95ci = lwr.ci, rs_upr_95ci = upr.ci, count_string = n_str) %>%
  arrange(analysis_type, desc(rs_estimate))

cat(sprintf("  ST6 (Figure1c): %d rows\n", nrow(ST6)))
# ============================================================================
# SECTION 9: BUILD KEY SHEET (Definitions + Column Descriptions)
# ============================================================================

message("[9/12] Building Key sheet...")

# Create Key sheet with definitions, universes, formulas, and column definitions
Key <- tribble(
  ~Section, ~Item, ~Description,
  
  # DEFINITIONS
  "DEFINITIONS", "pQTL+ TIs require ALL of:", NA,
  NA, "- MR p-value", "< 0.05/47M (Bonferroni-corrected)",
  NA, "- Indication-trait similarity", ">= 0.8 (semantic similarity between indication and GWAS trait)",
  NA, "- L2G share", ">= 0.5 (locus-to-gene score from Open Targets Genetics)",
  NA, "- Gene platform", "Gene measured on Olink or SomaScan proteomics platform",
  NA, "- Clinical stage", "Reached at least Phase I",
  NA, "- For MR+coloc", "Additionally requires colocalization H4 >= 0.8",
  NA, NA, NA,
  NA, "pQTL- TIs:", NA,
  NA, "- Definition", "All TIs in the universe that do not meet pQTL+ criteria",
  NA, NA, NA,
  
  # UNIVERSES
  "UNIVERSES", "Universe 1", "All Phase I-entered TIs where target is a measured protein (includes TIs with/without traditional genetic evidence)",
  NA, "Universe 2", "Subset of Universe 1 with existing genetic evidence (OTG L2G >= 0.5, OMIM, PICCOLO H4 >= 0.9, or Genebass)",
  NA, "- Note", "pQTL+ TIs are identical across universes; only the pQTL- comparator group changes",
  NA, NA, NA,
  
  # FORMULAS
  "FORMULAS", "Relative Success (RS)", "RS = (pqtl_success / pqtl_total) / (nopqtl_success / nopqtl_total)",
  NA, "Baseline Success Rate", "baseline_sr = nopqtl_success / nopqtl_total",
  NA, "pQTL Prevalence", "pqtl_prevalence = pqtl_total / (pqtl_total + nopqtl_total)",
  NA, "Spearman Correlation", "Rank correlation between baseline_sr and pqtl_prevalence across TAs; tests if pQTL is enriched in 'easy' TAs",
  NA, "Breslow-Day Test", "Tests homogeneity of odds ratios across TA strata; non-significant p indicates consistent effect across TAs",
  NA, "Chi-square Test", "Tests if pQTL+ TIs are non-randomly distributed across TAs; compares observed vs expected under proportional distribution",
  NA, "Observed/Expected Ratio", "obs_exp_ratio = pqtl_total / expected_pqtl, where expected_pqtl = ta_total * sum(pqtl_total) / sum(ta_total)",
  NA, NA, NA,
  
  # COLUMN DEFINITIONS BY SHEET
  "COLUMN DEFINITIONS", "--- ST1: Outcome_GWAS ---", NA,
  NA, "outcome_datasets", "Study identifier for outcome GWAS",
  NA, "outcome_data_source", "Source database (GWAS Catalog, FinnGen, pan-UK Biobank)",
  NA, "outcome_trait", "Description of the outcome trait",
  NA, "outcome_trait_efo", "Experimental Factor Ontology (EFO) ID(s) for the trait",
  NA, NA, NA,
  
  NA, "--- ST2: Proteomic_GWAS ---", NA,
  NA, "proteomic_dataset", "Name of the proteomic GWAS dataset",
  NA, "platform", "Proteomics platform used (Olink, Somalogic, BioRad)",
  NA, "sample_size", "Number of samples in the study",
  NA, "no. of proteins measured", "Number of proteins measured",
  NA, NA, NA,
  
  NA, "--- ST3: Excluded_Traits ---", NA,
  NA, "trait_key_term", "Trait excluded from analysis (non-disease or QC trait)",
  NA, NA, NA,
  
  NA, "--- ST4: Figure1a ---", NA,
  NA, "panel_group", "Grouping of enrichment analysis (by genetic evidence source, L2G share, etc.)",
  NA, "source_label", "Label for the genetic evidence source or analysis type",
  NA, "rs_estimate", "Relative success estimate (ratio of success rates)",
  NA, "rs_lwr_95ci", "Lower bound of 95% confidence interval",
  NA, "rs_upr_95ci", "Upper bound of 95% confidence interval",
  NA, "count_string", "(supported_successes/supported_total)/(unsupported_successes/unsupported_total)",
  NA, "row_order", "Order for display in figure",
  NA, NA, NA,
  
  NA, "--- ST5: Figure1b ---", NA,
  NA, "target_indication_id", "Unique identifier (gene-indication_mesh_id)",
  NA, "gene", "HGNC gene symbol",
  NA, "indication_mesh_id", "MeSH identifier for the indication",
  NA, "indication_mesh_term", "MeSH term for the indication",
  NA, "cis_mr", "Binary: cis-MR evidence (1=yes, 0=no)",
  NA, "trans_mr", "Binary: trans-MR evidence",
  NA, "mixed_mr", "Binary: mixed-MR evidence",
  NA, "coloc_cis", "Binary: cis-colocalization evidence (H4>=0.8)",
  NA, "coloc_trans", "Binary: trans-colocalization evidence (H4>=0.8)",
  NA, "n_methods_supported", "Total number of methods supporting this pair",
  NA, NA, NA,
  
  NA, "--- ST6: Figure1c ---", NA,
  NA, "analysis_type", "Type of analysis (L2G only vs L2G + pQTL)",
  NA, "family", "Gene family name",
  NA, "rs_estimate", "Relative success estimate for this gene family",
  NA, "rs_lwr_95ci", "Lower bound of 95% CI",
  NA, "rs_upr_95ci", "Upper bound of 95% CI",
  NA, "count_string", "(supported_successes/supported_total)/(unsupported_successes/unsupported_total)",
  NA, NA, NA,
  
  NA, "--- ST7: ChiSq_Distribution ---", NA,
  NA, "universe", "Universe label",
  NA, "chisq", "Chi-square statistic",
  NA, "df", "Degrees of freedom",
  NA, "p", "P-value; significant indicates non-random pQTL distribution across TAs",
  NA, NA, NA,
  
  NA, "--- ST8/ST9: TA_Enrichment ---", NA,
  NA, "therapeutic_area", "Therapeutic area name",
  NA, "ta_total", "Total TIs in this TA (pqtl_total + nopqtl_total)",
  NA, "pqtl_total", "Observed number of pQTL+ TIs",
  NA, "expected_pqtl", "Expected pQTL+ TIs if distributed proportionally",
  NA, "obs_exp_ratio", "Observed/Expected ratio; >1 = enriched, <1 = depleted",
  NA, "enrichment", "Label: 'enriched' (>1.5), 'depleted' (<0.67), or 'as expected'",
  NA, NA, NA,
  
  NA, "--- ST10: Spearman_Summary ---", NA,
  NA, "universe", "Universe label",
  NA, "rho", "Spearman rank correlation coefficient",
  NA, "p", "P-value; non-significant indicates no selection bias toward easy TAs",
  NA, NA, NA,
  
  NA, "--- ST11/ST12: PerTA_RS ---", NA,
  NA, "therapeutic_area", "Therapeutic area name",
  NA, "pqtl_total", "Number of pQTL+ TIs in this TA",
  NA, "pqtl_success", "Number of pQTL+ TIs that launched",
  NA, "nopqtl_total", "Number of pQTL- TIs in this TA",
  NA, "nopqtl_success", "Number of pQTL- TIs that launched",
  NA, "rr", "Relative Success point estimate (Katz method)",
  NA, "lwr", "Lower 95% CI",
  NA, "upr", "Upper 95% CI",
  NA, "rs_str", "RS formatted as 'estimate (lower-upper)'",
  NA, "baseline_sr", "Baseline success rate = nopqtl_success / nopqtl_total",
  NA, "pqtl_prevalence", "Proportion with pQTL support = pqtl_total / (pqtl_total + nopqtl_total)",
  NA, "sparse", "TRUE if pqtl_total < 5 (limited information)",
  NA, NA, NA,
  
  NA, "--- ST13: BreslowDay_Summary ---", NA,
  NA, "universe", "Universe label",
  NA, "bd_chisq", "Breslow-Day chi-square statistic",
  NA, "bd_df", "Degrees of freedom (number of valid TAs - 1)",
  NA, "bd_p", "P-value; non-significant indicates homogeneous effect across TAs",
  NA, NA, NA,
  
  NA, "--- ST14/ST15: LOO ---", NA,
  NA, "dropped_ta", "TA dropped for sensitivity analysis ('none - overall' = no TA dropped)",
  NA, "pqtl", "pQTL+ successes/total after dropping TA",
  NA, "nopqtl", "pQTL- successes/total after dropping TA",
  NA, "rs", "RS point estimate after dropping TA",
  NA, "lwr", "Lower 95% CI",
  NA, "upr", "Upper 95% CI",
  NA, "rs_str", "RS formatted as 'estimate (lower-upper)'",
  NA, NA, NA,
  
  NA, "--- ST16: All_MR_pairs ---", NA,
  NA, "hgnc_protein", "HGNC gene symbol for the protein",
  NA, "ttpair", "Unique target-trait pair identifier",
  NA, "data_sources", "Comma-separated list of data sources",
  NA, "outcomes", "Comma-separated list of outcome identifiers",
  NA, "outcome_trait", "Description of the outcome trait",
  NA, "trait_key_term", "Standardized trait key term",
  NA, "bxy", "MR effect estimate (beta)",
  NA, "bxy_se", "Standard error of MR effect",
  NA, "bxy_pval", "P-value for MR effect",
  NA, "bxy_pval_mantissa", "Mantissa of p-value (for very small p-values)",
  NA, "bxy_pval_exponent", "Exponent of p-value",
  NA, "bxy_egger", "MR-Egger effect estimate",
  NA, "bxy_se_egger", "SE of MR-Egger estimate",
  NA, "bxy_pval_egger", "P-value for MR-Egger test",
  NA, "bxy_median", "Weighted median MR effect estimate",
  NA, "bxy_se_median", "SE of weighted median estimate",
  NA, "bxy_pval_median", "P-value for weighted median test",
  NA, "snp_ciscoloc_all", "SNPs with cis-colocalization (H4>=0.8)",
  NA, "snp_transcoloc_all", "SNPs with trans-colocalization (H4>=0.8)",
  NA, "trans_coloc_genes_all", "Genes near trans-colocalization SNPs",
  NA, "pav_cismr", "Protein-altering variant annotation (yes/no)",
  NA, "pav_terms_all", "Sequence ontology terms for PAVs",
  NA, "l2g_share", "Locus-to-gene score share",
  NA, "hla", "Gene in HLA region (yes/no)",
  NA, "harmonic_genetic_score", "Open Targets harmonic sum of genetic evidence",
  NA, "replication", "Finding replicated across datasets (TRUE/FALSE)",
  NA, "triangulation", "Open Targets data type IDs supporting this association",
  NA, "multi_method_count", "Number of MR/coloc methods supporting this pair",
  NA, "mr_coloc_types", "Comma-separated list of MR/coloc types",
  NA, "positive_control", "Known drug target for this indication (yes/no)",
  NA, "repositioning_opportunity", "Drug target for other indications (yes/no)",
  NA, NA, NA,
  
  NA, "--- ST17: pqtl_success_ti_pairs ---", NA,
  NA, "gene", "HGNC gene symbol",
  NA, "indication_mesh_id", "MeSH identifier for the indication",
  NA, "indication_mesh_term", "MeSH term for the indication",
  NA, "therapeutic_area", "Therapeutic area",
  NA, "ti_uid", "Unique target-indication identifier"
)

message("  Key sheet built with ", nrow(Key), " rows")
# ============================================================================
# SECTION 11: CREATE EXCEL WORKBOOK
# ============================================================================

message("[11/12] Creating Excel workbook...")

# Prepare all tables
all_tables <- list(
  "Key" = Key,
  "ST1 - Outcome_GWAS" = ST1,
  "ST2 - Proteomic_GWAS" = ST2,
  "ST3 - Excluded_Traits" = ST3,
  "ST4 - Figure1a" = ST4,
  "ST5 - Figure1b" = ST5,
  "ST6 - Figure1c" = ST6,
  "ST7 - ChiSq_Distribution" = ST7,
  "ST8 - TA_Enrichment_Universe1" = ST8,
  "ST9 - TA_Enrichment_Universe2" = ST9,
  "ST10 - Spearman_Summary" = ST10,
  "ST11 - PerTA_RS_Universe1" = ST11,
  "ST12 - PerTA_RS_Universe2" = ST12,
  "ST13 - BreslowDay_Summary" = ST13,
  "ST14 - LOO_Universe1" = ST14,
  "ST15 - LOO_Universe2" = ST15,
  "ST16 - All_MR_pairs" = ST16,
  "ST17 - pqtl_success_ti_pairs" = ST17
)

# Sheet titles
sheet_titles <- list(
  "Key" = "Key: Definitions, Formulas, and Column Descriptions",
  "ST1 - Outcome_GWAS" = "Supplementary Table 1: Outcome GWAS Datasets",
  "ST2 - Proteomic_GWAS" = "Supplementary Table 2: Proteomic GWAS Datasets",
  "ST3 - Excluded_Traits" = "Supplementary Table 3: Excluded Traits",
  "ST4 - Figure1a" = "Supplementary Table 4: Enrichment by Genetic Evidence Source (Figure 1a)",
  "ST5 - Figure1b" = "Supplementary Table 5: UpSet Plot Data - MR/Coloc Overlap (Figure 1b)",
  "ST6 - Figure1c" = "Supplementary Table 6: Gene Family Enrichment (Figure 1c)",
  "ST7 - ChiSq_Distribution" = "Supplementary Table 7: Chi-Square Test - pQTL Distribution Across TAs",
  "ST8 - TA_Enrichment_Universe1" = "Supplementary Table 8: Therapeutic Area Enrichment (Universe 1)",
  "ST9 - TA_Enrichment_Universe2" = "Supplementary Table 9: Therapeutic Area Enrichment (Universe 2)",
  "ST10 - Spearman_Summary" = "Supplementary Table 10: Spearman Correlation Summary",
  "ST11 - PerTA_RS_Universe1" = "Supplementary Table 11: Per-TA Relative Success (Universe 1)",
  "ST12 - PerTA_RS_Universe2" = "Supplementary Table 12: Per-TA Relative Success (Universe 2)",
  "ST13 - BreslowDay_Summary" = "Supplementary Table 13: Breslow-Day Heterogeneity Test Summary",
  "ST14 - LOO_Universe1" = "Supplementary Table 14: Leave-One-Out Sensitivity (Universe 1)",
  "ST15 - LOO_Universe2" = "Supplementary Table 15: Leave-One-Out Sensitivity (Universe 2)",
  "ST16 - All_MR_pairs" = "Supplementary Table 16: All Bonferroni-Significant MR Target-Trait Pairs",
  "ST17 - pqtl_success_ti_pairs" = "Supplementary Table 17: Successful pQTL-Supported TI Pairs"
)

# Clean all tables
all_tables_clean <- map(all_tables, clean_for_export)

# Create workbook
wb <- createWorkbook()

# Define styles
title_style <- createStyle(
  fontSize = 12,
  textDecoration = "bold",
  fgFill = "#4472C4",
  fontColour = "#FFFFFF",
  halign = "left",
  valign = "center"
)

header_style <- createStyle(
  fontSize = 11,
  textDecoration = "bold",
  fgFill = "#D9E2F3",
  border = "TopBottomLeftRight",
  borderColour = "#8EA9DB"
)

section_style <- createStyle(
  fontSize = 11,
  textDecoration = "bold",
  fgFill = "#E2EFDA"
)

# Add each sheet
for (sheet_name in names(all_tables_clean)) {
  
  df <- all_tables_clean[[sheet_name]]
  title <- sheet_titles[[sheet_name]]
  
  addWorksheet(wb, sheet_name)
  
  # Row 1: Title (merged cell)
  writeData(wb, sheet_name, title, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:max(ncol(df), 3), rows = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
  setRowHeights(wb, sheet_name, rows = 1, heights = 25)
  
  if (sheet_name == "Key") {
    # Special handling for Key sheet - no header row, just data
    writeData(wb, sheet_name, df, startRow = 3, startCol = 1, headerStyle = header_style)
    
    # Style section headers (DEFINITIONS, UNIVERSES, FORMULAS, COLUMN DEFINITIONS)
    section_rows <- which(!is.na(df$Section)) + 3  # +3 for title row + empty row + header
    for (r in section_rows) {
      addStyle(wb, sheet_name, section_style, rows = r, cols = 1:3, gridExpand = TRUE)
    }
    
    # Set column widths
    setColWidths(wb, sheet_name, cols = 1, widths = 25)
    setColWidths(wb, sheet_name, cols = 2, widths = 35)
    setColWidths(wb, sheet_name, cols = 3, widths = 100)
    
  } else {
    # Standard sheets: Header row + data
    writeData(wb, sheet_name, df, startRow = 3, startCol = 1, headerStyle = header_style)
    
    # Apply header style
    addStyle(wb, sheet_name, header_style, rows = 3, cols = 1:ncol(df), gridExpand = TRUE)
    
    # Auto-width columns (with max width cap)
    setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
    
    # Freeze panes (title + header)
    freezePane(wb, sheet_name, firstActiveRow = 4, firstActiveCol = 1)
  }
  
  cat(sprintf("  Added sheet: %s (%d rows)\n", sheet_name, nrow(df)))
}

# ============================================================================
# SECTION 12: SAVE WORKBOOK
# ============================================================================

message("[12/12] Saving workbook...")

output_file <- file.path(output_dir, "mrcoloc_supplement.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(sprintf("\n  Saved: %s\n", output_file))

# ============================================================================
# SUMMARY
# ============================================================================

cat("
================================================================================
  Supplementary Tables Complete!
================================================================================

Output file: ", output_file, "

Sheet summary:
")

for (sheet_name in names(all_tables_clean)) {
  df <- all_tables_clean[[sheet_name]]
  cat(sprintf("  %-30s %6d rows x %2d cols\n", sheet_name, nrow(df), ncol(df)))
}

cat("
Key statistics:
  ST1  - Outcome GWAS datasets:        ", nrow(ST1), "
  ST2  - Proteomic GWAS datasets:      ", nrow(ST2), "
  ST3  - Excluded traits:              ", nrow(ST3), "
  ST4  - Enrichment by source:         ", nrow(ST4), "
  ST5  - UpSet plot data:              ", nrow(ST5), "
  ST6  - Gene family enrichment:       ", nrow(ST6), "
  ST7  - Chi-square distribution:      ", nrow(ST7), "
  ST8  - TA Enrichment (U1):           ", nrow(ST8), "
  ST9  - TA Enrichment (U2):           ", nrow(ST9), "
  ST10 - Spearman summary:             ", nrow(ST10), "
  ST11 - Per-TA RS (U1):               ", nrow(ST11), "
  ST12 - Per-TA RS (U2):               ", nrow(ST12), "
  ST13 - Breslow-Day summary:          ", nrow(ST13), "
  ST14 - LOO (U1):                     ", nrow(ST14), "
  ST15 - LOO (U2):                     ", nrow(ST15), "
  ST16 - All MR pairs:                 ", nrow(ST16), "
  ST17 - pQTL success TI pairs:        ", nrow(ST17), "

================================================================================
")