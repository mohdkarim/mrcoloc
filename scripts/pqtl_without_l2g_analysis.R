#!/usr/bin/env Rscript
# ============================================================================
# Analysis: pQTL-supported pairs stratified by L2G status
# ============================================================================
# Addresses Reviewer 1 comment 5: "what is the relative success of
# pQTL-supported pairs that lack L2G support?"
#
# Three analyses:
# 1. pQTL support WITH L2G >= 0.5 (current definition)
# 2. pQTL support WITHOUT L2G >= 0.5 (L2G < 0.5 or missing)
# 3. pQTL support regardless of L2G (all pQTL MR-significant)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DescTools)
})

project_root <- Sys.getenv("MRCOLOC_ROOT", getwd())
data_dir     <- file.path(project_root, "data")
minikel_dir  <- file.path(project_root, "data", "minikel")

# Load data
merge3_pqtl <- readRDS(file.path(project_root, "data_raw/merge3_pqtl.rds"))
indic       <- read_tsv(file.path(data_dir, "indic.tsv"), show_col_types = FALSE)
pgenes      <- readRDS(file.path(project_root, "data_raw/pgenes.rds"))

SIMILARITY_THR <- 0.8
MIN_L2G_SHARE  <- 0.5

cat("\n=== pQTL without L2G Analysis ===\n\n")
cat("Background gene set (pgenes):", length(pgenes), "genes\n")
cat("Total rows in merge3_pqtl:", nrow(merge3_pqtl), "\n\n")

# --- Helper: compute RS for a given set of supported TIs ---
compute_rs_from_support <- function(support_tis, universe_df, label) {
  ti_best <- universe_df %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat)) %>%
    left_join(indic %>% select(indication_mesh_id, genetic_insight),
              by = "indication_mesh_id") %>%
    filter(genetic_insight != "none") %>%
    mutate(highest_phase = case_when(
      !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
      !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_
    )) %>%
    arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
    group_by(ti_uid) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gensup = ti_uid %in% support_tis$ti_uid)

  long     <- ti_best %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_best %>% filter(!is.na(succ_1_2))

  succ_gs    <- sum(long$gensup & long$success, na.rm = TRUE)
  total_gs   <- sum(baseline$gensup, na.rm = TRUE)
  succ_nogs  <- sum(!long$gensup & long$success, na.rm = TRUE)
  total_nogs <- sum(!baseline$gensup, na.rm = TRUE)

  if (total_gs == 0 || total_nogs == 0 || succ_gs > total_gs) {
    cat(sprintf("  %s: CANNOT COMPUTE (supported: %d/%d, unsupported: %d/%d)\n",
                label, succ_gs, total_gs, succ_nogs, total_nogs))
    return(tibble(label = label, est = NA, lwr.ci = NA, upr.ci = NA,
                  succ_gs = succ_gs, total_gs = total_gs,
                  succ_nogs = succ_nogs, total_nogs = total_nogs,
                  rate_supported = NA, rate_unsupported = NA))
  }

  out <- as.data.frame(BinomRatioCI(
    x1 = succ_gs, n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  ))

  rate_sup   <- succ_gs / total_gs
  rate_unsup <- succ_nogs / total_nogs

  cat(sprintf("  %s:\n", label))
  cat(sprintf("    RS = %.2f (95%% CI: %.2f - %.2f)\n", out$est, out$lwr.ci, out$upr.ci))
  cat(sprintf("    Supported:   %d/%d = %.1f%%\n", succ_gs, total_gs, rate_sup * 100))
  cat(sprintf("    Unsupported: %d/%d = %.1f%%\n", succ_nogs, total_nogs, rate_unsup * 100))

  tibble(label = label, est = out$est, lwr.ci = out$lwr.ci, upr.ci = out$upr.ci,
         succ_gs = succ_gs, total_gs = total_gs,
         succ_nogs = succ_nogs, total_nogs = total_nogs,
         rate_supported = rate_sup, rate_unsupported = rate_unsup)
}

# Restrict universe to measured proteins
universe <- merge3_pqtl %>% filter(gene %in% pgenes)

# ============================================================================
# Analysis 1: Current definition - pQTL + L2G >= 0.5
# ============================================================================
cat("--- Analysis 1: pQTL with L2G >= 0.5 (current paper definition) ---\n")

support_with_l2g <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
    gene %in% pgenes
  ) %>%
  distinct(ti_uid)

cat("  Supported TI pairs:", nrow(support_with_l2g), "\n")
r1 <- compute_rs_from_support(support_with_l2g, universe, "pQTL + L2G >= 0.5")

# ============================================================================
# Analysis 2: pQTL support WITHOUT L2G (L2G < 0.5 or missing)
# ============================================================================
cat("\n--- Analysis 2: pQTL without L2G (L2G < 0.5 or missing) ---\n")

support_no_l2g <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    (is.na(l2g_share) | l2g_share < MIN_L2G_SHARE),
    gene %in% pgenes
  ) %>%
  distinct(ti_uid)

cat("  Supported TI pairs:", nrow(support_no_l2g), "\n")
r2 <- compute_rs_from_support(support_no_l2g, universe, "pQTL without L2G (< 0.5)")

# ============================================================================
# Analysis 3: pQTL support regardless of L2G
# ============================================================================
cat("\n--- Analysis 3: pQTL regardless of L2G (any l2g_share) ---\n")

support_any_l2g <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    gene %in% pgenes
  ) %>%
  distinct(ti_uid)

cat("  Supported TI pairs:", nrow(support_any_l2g), "\n")
r3 <- compute_rs_from_support(support_any_l2g, universe, "pQTL any L2G")

# ============================================================================
# Analysis 4: More granular L2G stratification
# ============================================================================
cat("\n--- Analysis 4: pQTL support by L2G strata ---\n")

results_strata <- map_dfr(list(
  list(label = "pQTL + L2G >= 0.75", filter_fn = function(df) df %>% filter(!is.na(l2g_share), l2g_share >= 0.75)),
  list(label = "pQTL + L2G 0.5-0.75", filter_fn = function(df) df %>% filter(!is.na(l2g_share), l2g_share >= 0.5, l2g_share < 0.75)),
  list(label = "pQTL + L2G 0.25-0.5", filter_fn = function(df) df %>% filter(!is.na(l2g_share), l2g_share >= 0.25, l2g_share < 0.5)),
  list(label = "pQTL + L2G < 0.25", filter_fn = function(df) df %>% filter(!is.na(l2g_share), l2g_share < 0.25)),
  list(label = "pQTL + L2G missing", filter_fn = function(df) df %>% filter(is.na(l2g_share)))
), function(spec) {
  support <- merge3_pqtl %>%
    filter(
      grepl("pqtl", original_link, ignore.case = TRUE),
      comb_norm >= SIMILARITY_THR,
      gene %in% pgenes
    ) %>%
    spec$filter_fn() %>%
    distinct(ti_uid)

  cat("  ", spec$label, "- Supported TI pairs:", nrow(support), "\n")
  compute_rs_from_support(support, universe, spec$label)
})

# ============================================================================
# Analysis 5: What fraction of the 23 launched pQTL pairs have L2G < 0.5?
# ============================================================================
cat("\n--- Analysis 5: L2G distribution among launched pQTL pairs ---\n")

launched_pqtl <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= SIMILARITY_THR,
    gene %in% pgenes,
    succ_3_a == TRUE
  ) %>%
  distinct(ti_uid, gene, indication_mesh_term, l2g_share) %>%
  arrange(l2g_share)

cat("  Total launched pQTL TI pairs (with any l2g):", nrow(launched_pqtl), "\n")
cat("  With L2G >= 0.5:", sum(!is.na(launched_pqtl$l2g_share) & launched_pqtl$l2g_share >= 0.5), "\n")
cat("  With L2G < 0.5:", sum(!is.na(launched_pqtl$l2g_share) & launched_pqtl$l2g_share < 0.5), "\n")
cat("  With L2G missing:", sum(is.na(launched_pqtl$l2g_share)), "\n")
cat("\n  Launched pairs with low/missing L2G:\n")

low_l2g_launched <- launched_pqtl %>%
  filter(is.na(l2g_share) | l2g_share < 0.5)

if (nrow(low_l2g_launched) > 0) {
  print(low_l2g_launched, n = Inf)
} else {
  cat("  (none)\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n\n=== SUMMARY TABLE ===\n")
all_results <- bind_rows(r1, r2, r3, results_strata) %>%
  filter(!is.na(est))

cat(sprintf("\n%-30s %6s %10s %20s %20s\n",
            "Label", "RS", "95% CI", "Supported", "Unsupported"))
cat(strrep("-", 95), "\n")
for (i in seq_len(nrow(all_results))) {
  r <- all_results[i, ]
  cat(sprintf("%-30s %6.2f (%5.2f-%5.2f) %8d/%-8d %8d/%-8d\n",
              r$label, r$est, r$lwr.ci, r$upr.ci,
              r$succ_gs, r$total_gs, r$succ_nogs, r$total_nogs))
}
cat("\n")
