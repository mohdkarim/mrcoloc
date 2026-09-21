#!/usr/bin/env Rscript
# ============================================================================
# R2.1: does universal pQTL coverage erode prognostic value?
# ============================================================================
# Reviewer 2's premise: "if virtually all the genome has a pQTL, then almost by
# definition that does not then offer any prognostic value".
#
# Counter: pQTL existence is necessary but nowhere near sufficient. Support
# additionally requires (i) a Bonferroni-significant MR association between the
# protein and a disease trait, (ii) MeSH similarity >= 0.8 between that trait
# and the drug's indication, and (iii) L2G share >= 0.5. This script quantifies
# how much each step discriminates.
# ============================================================================

suppressPackageStartupMessages({library(tidyverse); library(openxlsx)})

merge3 <- readRDS("data_raw/merge3_pqtl.rds")
indic  <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes <- readRDS("data_raw/pgenes.rds")
BONF   <- 0.05 / 4.7e7

# ---- scale of the underlying MR resource -----------------------------------
st16 <- read.xlsx("output/mrcoloc_supplement.xlsx",
                  sheet = "ST16 - All_MR_pairs", startRow = 3)
cat("\n=== scale of the pQTL MR resource ===\n")
cat(sprintf("  proteins measured on the platforms (background)   : %d\n", length(pgenes)))
cat(sprintf("  Bonferroni-significant MR target-trait pairs (ST16): %d\n", nrow(st16)))
cat(sprintf("  distinct proteins with >=1 significant MR assoc    : %d (%.0f%% of measured)\n",
            n_distinct(st16$hgnc_protein),
            100 * n_distinct(st16$hgnc_protein) / length(pgenes)))

# ---- filtering cascade over Phase I T-I pairs -------------------------------
base <- merge3 %>%
  filter(gene %in% pgenes, !is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "", !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none")

ti_phase1 <- base %>%
  mutate(hp = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
                        !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
  arrange(ti_uid, desc(hp), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup() %>%
  filter(!is.na(succ_1_2)) %>% pull(ti_uid)

pq <- base %>% filter(grepl("pqtl", original_link, ignore.case = TRUE))

s1 <- ti_phase1
s2 <- pq %>% filter(ti_uid %in% s1) %>% distinct(ti_uid) %>% pull(ti_uid)
s3 <- pq %>% filter(ti_uid %in% s1, !is.na(bxy_pval), bxy_pval <= BONF) %>%
      distinct(ti_uid) %>% pull(ti_uid)
s4 <- pq %>% filter(ti_uid %in% s1, !is.na(bxy_pval), bxy_pval <= BONF,
                    comb_norm >= 0.8) %>% distinct(ti_uid) %>% pull(ti_uid)
s5 <- pq %>% filter(ti_uid %in% s1, !is.na(bxy_pval), bxy_pval <= BONF,
                    comb_norm >= 0.8, !is.na(l2g_share), l2g_share >= 0.5) %>%
      distinct(ti_uid) %>% pull(ti_uid)

# NOTE: s2 ("any pQTL MR association present") and s3 ("Bonferroni-significant")
# are the SAME set here, because merge3 only ever retains Bonferroni-significant
# associations - the significance filter is applied upstream, not by this script.
# The s2 row was therefore removed from the published table: it retained 100% of
# the previous row, which reads as a filter that did nothing and invites the
# question "what did that step do?". The two are collapsed into one step whose
# label states the threshold. s2 is still computed and checked below so the fact
# is recorded in the run log. Do not reinstate the row without re-checking this.
stopifnot(setequal(s2, s3))

casc <- tibble(
  step = c("Phase I T-I pairs, target measured on a proteomics platform",
           "  + Bonferroni-significant pQTL MR association for the target (p < 1.06e-9)",
           "  + MR trait MeSH-matched to the indication (>= 0.8)",
           "  + L2G share >= 0.5  [= pQTL-supported, as published]"),
  n = c(length(s1), length(s3), length(s4), length(s5))) %>%
  mutate(pct_of_universe = sprintf("%.2f%%", 100 * n / length(s1)),
         retained_from_prev = c(NA, sprintf("%.1f%%", 100 * n[-1] / lag(n)[-1])))

cat(sprintf("\n  [check] pairs with any pQTL association = %d; with a Bonferroni-significant\n", length(s2)))
cat(sprintf("          one = %d. Identical, because merge3 is already significance-filtered,\n", length(s3)))
cat("          so these two steps are reported as one row in ST20.\n")

cat("\n=== filtering cascade ===\n\n")
print(as.data.frame(casc), row.names = FALSE, right = FALSE)

cat("\n=== the discriminating step ===\n")
cat(sprintf("  Of %d Phase I pairs whose target IS measured (i.e. pQTL data available),\n", length(s1)))
cat(sprintf("  %d (%.2f%%) meet the support definition. %d (%.1f%%) do not.\n",
            length(s5), 100*length(s5)/length(s1),
            length(s1)-length(s5), 100*(length(s1)-length(s5))/length(s1)))
cat(sprintf("\n  Even among pairs with a Bonferroni-significant pQTL MR association\n"))
cat(sprintf("  somewhere (%d pairs), only %d (%.1f%%) survive disease matching and L2G.\n",
            length(s3), length(s5), 100*length(s5)/length(s3)))

write_tsv(casc, "output/r2_1_saturation_cascade.tsv")
cat("\n-> output/r2_1_saturation_cascade.tsv\n\n")
