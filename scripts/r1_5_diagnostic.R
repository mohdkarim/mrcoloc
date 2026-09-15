#!/usr/bin/env Rscript
# ============================================================================
# R1.5 diagnostic: is "lacks L2G support" being defined correctly?
# ============================================================================
# Earlier analysis filtered ROWS by l2g_share band then took distinct ti_uid.
# Because l2g_share is row-level (per gene-trait-study triplet), a TI pair with
# rows at both 0.15 and 0.90 would land in BOTH the "<0.5" and ">=0.5" strata.
# The correct TI-level definition of "lacks L2G support" is max(l2g_share) < thr.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DescTools)
})

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic       <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes      <- readRDS("data_raw/pgenes.rds")

SIMILARITY_THR <- 0.8

cat("\n================================================================\n")
cat("  DIAGNOSTIC 1: row-level vs TI-level l2g_share\n")
cat("================================================================\n\n")

# pQTL support WITHOUT any L2G condition (the correct primitive)
pqtl_rows <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         gene %in% pgenes)

# TI-level max l2g_share among pQTL rows
ti_l2g <- pqtl_rows %>%
  group_by(ti_uid, gene) %>%
  summarise(
    max_l2g   = suppressWarnings(max(l2g_share, na.rm = TRUE)),
    min_l2g   = suppressWarnings(min(l2g_share, na.rm = TRUE)),
    n_rows    = n(),
    n_na_l2g  = sum(is.na(l2g_share)),
    .groups = "drop"
  ) %>%
  mutate(max_l2g = ifelse(is.infinite(max_l2g), NA_real_, max_l2g),
         min_l2g = ifelse(is.infinite(min_l2g), NA_real_, min_l2g))

cat("The genes I previously reported as 'launched with L2G < 0.5':\n\n")
check <- ti_l2g %>%
  filter(gene %in% c("PCSK9", "APOB", "ANGPTL3", "F2", "CSF3", "VWF", "IL4R", "KLK3")) %>%
  arrange(gene, ti_uid)
print(as.data.frame(check), row.names = FALSE)

cat("\n>> If max_l2g is high, these pairs DO have L2G support and were\n")
cat("   wrongly counted in the 'lacks L2G' stratum by the row-level filter.\n")

cat("\n================================================================\n")
cat("  DIAGNOSTIC 2: is l2g_share ever missing for pQTL pairs?\n")
cat("================================================================\n\n")

cat(sprintf("pQTL-supported TI pairs (MR + MeSH only): %d\n", nrow(ti_l2g)))
cat(sprintf("  with max_l2g non-missing: %d\n", sum(!is.na(ti_l2g$max_l2g))))
cat(sprintf("  with max_l2g MISSING:     %d\n", sum(is.na(ti_l2g$max_l2g))))
cat(sprintf("  rows with NA l2g_share:   %d / %d\n",
            sum(pqtl_rows$l2g_share %>% is.na()), nrow(pqtl_rows)))

cat("\nIs TNF present at all in the analysis dataset?\n")
tnf <- merge3_pqtl %>% filter(gene == "TNF") %>%
  summarise(rows = n(),
            n_pqtl_rows = sum(grepl("pqtl", original_link, ignore.case = TRUE)),
            n_na_l2g = sum(is.na(l2g_share)),
            max_l2g = suppressWarnings(max(l2g_share, na.rm = TRUE)))
print(as.data.frame(tnf), row.names = FALSE)

cat("\nIs SOST present?\n")
sost <- merge3_pqtl %>% filter(gene == "SOST") %>%
  summarise(rows = n(),
            n_pqtl_rows = sum(grepl("pqtl", original_link, ignore.case = TRUE)),
            n_na_l2g = sum(is.na(l2g_share)),
            max_l2g = suppressWarnings(max(l2g_share, na.rm = TRUE)))
print(as.data.frame(sost), row.names = FALSE)

cat("\n================================================================\n")
cat("  DIAGNOSTIC 3: distribution of TI-level max L2G among pQTL pairs\n")
cat("================================================================\n\n")

brk <- ti_l2g %>%
  mutate(band = case_when(
    is.na(max_l2g)   ~ "missing",
    max_l2g <  0.25  ~ "< 0.25",
    max_l2g <  0.50  ~ "0.25-0.50",
    max_l2g <  0.75  ~ "0.50-0.75",
    TRUE             ~ ">= 0.75"
  )) %>%
  count(band)
print(as.data.frame(brk), row.names = FALSE)

cat("\nCumulative (the reviewer's requested direction):\n")
for (thr in c(0.25, 0.50, 0.75)) {
  n_below <- sum(!is.na(ti_l2g$max_l2g) & ti_l2g$max_l2g < thr)
  cat(sprintf("  pQTL-supported pairs with max L2G < %.2f: %d\n", thr, n_below))
}

cat("\n================================================================\n")
cat("  Saving TI-level L2G table for downstream use\n")
cat("================================================================\n")
saveRDS(ti_l2g, "data_raw/ti_l2g_pqtl_maxshare.rds")
cat("  -> data_raw/ti_l2g_pqtl_maxshare.rds\n\n")
