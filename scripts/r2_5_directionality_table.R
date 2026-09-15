#!/usr/bin/env Rscript
# ============================================================================
# R2.5: build the directional-concordance curation table
# ============================================================================
# Reviewer 2 comment 5 asks whether the direction of the pQTL effect is
# consistent with the therapeutic perturbation, and specifically about CIS
# effects ("I'm talking for all pQTLs, cis-effects"), noting trans is already
# discussed at L141-151.
#
# bxy is expressed per SD HIGHER genetically predicted plasma protein (L208-209).
#   bxy > 0  => higher protein associated with HIGHER risk of the outcome
#   bxy < 0  => higher protein associated with LOWER risk
#
# Concordance (to be evaluated on CIS effects only):
#   drug INHIBITS target  -> concordant if bxy > 0 (drug lowers a harmful protein)
#   drug ACTIVATES/replaces -> concordant if bxy < 0 (drug raises a protective protein)
#
# Trans effects are reported separately and NOT scored, because plasma levels of
# a trans-regulated protein need not track target activity in the therapeutic
# direction (see CSF3/CSF3R, L141-151).
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic       <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes      <- readRDS("data_raw/pgenes.rds")

STHR <- 0.8
L2GT <- 0.5
BONF <- 0.05 / 4.7e7      # ~1.06e-9, as stated at L229

base <- merge3_pqtl %>%
  filter(gene %in% pgenes,
         !is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "", !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none")

# pQTL-supported TI pairs (paper definition)
pqtl_ids <- base %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= STHR, !is.na(l2g_share), l2g_share >= L2GT) %>%
  distinct(ti_uid) %>% pull(ti_uid)

# Phase-I entrants among them
ti_phase1 <- base %>%
  mutate(hp = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
                        !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
  arrange(ti_uid, desc(hp), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup() %>%
  filter(ti_uid %in% pqtl_ids, !is.na(succ_1_2)) %>%
  transmute(ti_uid, gene, indication_mesh_id, indication_mesh_term,
            therapeutic_area,
            launched = !is.na(succ_3_a) & succ_3_a,
            highest_phase = hp)

cat(sprintf("\npQTL-supported Phase I pairs: %d (launched: %d)\n",
            nrow(ti_phase1), sum(ti_phase1$launched)))

# --- direction of effect per pair, split by cis / trans ---------------------
mr_rows <- base %>%
  filter(ti_uid %in% ti_phase1$ti_uid,
         grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= STHR,
         !is.na(bxy), !is.na(bxy_pval), bxy_pval <= BONF)

dir_summary <- mr_rows %>%
  mutate(kind = case_when(cis_trans_mr == "Cis"   ~ "cis",
                          cis_trans_mr == "Trans" ~ "trans",
                          cis_trans_mr == "Mixed" ~ "mixed",
                          TRUE ~ "other")) %>%
  group_by(ti_uid, kind) %>%
  summarise(
    n_assoc   = n(),
    n_pos     = sum(bxy > 0),
    n_neg     = sum(bxy < 0),
    bxy_med   = median(bxy),
    sign_str  = case_when(all(bxy > 0) ~ "+", all(bxy < 0) ~ "-", TRUE ~ "mixed sign"),
    .groups   = "drop"
  )

wide <- dir_summary %>%
  filter(kind %in% c("cis", "trans", "mixed")) %>%
  pivot_wider(id_cols = ti_uid, names_from = kind,
              values_from = c(sign_str, bxy_med, n_assoc))

tbl <- ti_phase1 %>%
  left_join(wide, by = "ti_uid") %>%
  transmute(
    gene, indication = indication_mesh_term, therapeutic_area,
    launched = if_else(launched, "launched", "did not launch"),
    highest_phase,
    cis_sign    = sign_str_cis,
    cis_bxy     = round(bxy_med_cis, 3),
    n_cis       = n_assoc_cis,
    trans_sign  = sign_str_trans,
    trans_bxy   = round(bxy_med_trans, 3),
    n_trans     = n_assoc_trans,
    mixed_sign  = sign_str_mixed,
    # --- TO BE CURATED ---
    drug_moa               = NA_character_,   # inhibitor | activator/replacement | other
    example_drug           = NA_character_,
    expected_bxy_sign      = NA_character_,   # + if inhibitor, - if activator
    cis_concordant         = NA_character_,   # yes | no | n/a
    moa_source             = NA_character_,
    notes                  = NA_character_
  ) %>%
  arrange(desc(launched == "launched"), gene, indication)

cat("\n--- unique targets requiring MoA curation ---\n")
ut <- tbl %>% distinct(gene) %>% arrange(gene) %>% pull(gene)
cat(sprintf("  %d targets: %s\n", length(ut), paste(ut, collapse = ", ")))

cat("\n--- cis direction availability ---\n")
cat(sprintf("  pairs with a cis effect:   %d\n", sum(!is.na(tbl$cis_sign))))
cat(sprintf("  pairs with only trans:     %d\n",
            sum(is.na(tbl$cis_sign) & !is.na(tbl$trans_sign))))
cat(sprintf("  pairs with mixed-only:     %d\n",
            sum(is.na(tbl$cis_sign) & is.na(tbl$trans_sign) & !is.na(tbl$mixed_sign))))
cat(sprintf("  cis sign inconsistent within pair: %d\n",
            sum(tbl$cis_sign == "mixed sign", na.rm = TRUE)))

cat("\n--- table ---\n")
print(as.data.frame(tbl %>% select(gene, indication, launched, cis_sign, cis_bxy,
                                  n_cis, trans_sign, trans_bxy, n_trans)),
      row.names = FALSE)

write_tsv(tbl, "output/r2_5_directionality_curation.tsv", na = "")
cat("\n  -> output/r2_5_directionality_curation.tsv\n")
cat("     (drug_moa / example_drug / expected_bxy_sign / cis_concordant left blank for curation)\n\n")
