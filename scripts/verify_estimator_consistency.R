#!/usr/bin/env Rscript
# ============================================================================
# Verify: do the L2G rows and pQTL rows of Figure 1a use the same estimator?
# ============================================================================
# L2G rows  -> pipeline_best() + advancement_rr()
# pQTL rows -> compute_rr() (direct BinomRatioCI on a single 2x2)
#
# Checks:
#   1. Does advancement_rr actually use a product-across-phases estimator,
#      or does always_katz short-circuit it to a direct ratio?
#   2. Reproduce ST4's L2G >= 0.5 row (127/484)/(1392/12538) to confirm the
#      pipeline is understood.
#   3. Compute the pQTL row BOTH ways on identical data and compare.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DescTools)
  library(binom)
  library(glue)
})

project_root <- getwd()
source(file.path(project_root, "R/pipeline_best.R"))
source(file.path(project_root, "R/advancement_rr.R"))

cat("\n================================================================\n")
cat("  1. Is the Wald/product branch reachable?\n")
cat("================================================================\n\n")
cat(sprintf("  always_katz = %s\n", always_katz))
cat(sprintf("  Wald branch condition is `!always_katz & x>=3 & y>=3` -> %s\n",
            (!always_katz) & TRUE))
cat("  => with always_katz=TRUE the product/Wald branch is UNREACHABLE;\n")
cat("     binom_ratio_atomic always falls through to BinomRatioCI (Katz),\n")
cat("     which uses (x_yes,n_yes,x_no,n_no) directly and IGNORES mean_yes/mean_no.\n")

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic  <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes <- readRDS("data_raw/pgenes.rds")
merge2 <- read_tsv("data/minikel/merge2.tsv.gz", show_col_types = FALSE) %>%
  mutate(
    otg_study = if_else(assoc_source == "OTG",
      str_remove(original_link, "https://genetics.opentargets.org/study/"), NA_character_),
    otg_study = str_remove(otg_study, "FINNGEN_R6_"),
    key = paste0(gene, "_", otg_study))

cat("\n================================================================\n")
cat("  2. Reproduce ST4 L2G >= 0.5 row via pipeline_best + advancement_rr\n")
cat("================================================================\n\n")
pb  <- pipeline_best(merge2, phase = "combined", basis = "ti",
                     associations = c("OTG"), share_mode = "L2G",
                     min_share = 0.5, verbose = FALSE)
rr  <- advancement_rr(pb)
il  <- rr %>% filter(phase == "I-Launch")
cat(sprintf("  advancement_rr I-Launch: (%d/%d)/(%d/%d)  RS=%.2f (%.2f-%.2f)\n",
            il$x_yes, il$n_yes, il$x_no, il$n_no, il$rs_mean, il$rs_l, il$rs_u))
cat("  ST4 published:           (127/484)/(1392/12538)  RS=2.36 (2.02-2.77)\n")
cat(sprintf("  MATCH: %s\n", identical(c(il$x_yes, il$n_yes, il$x_no, il$n_no),
                                       c(127L, 484L, 1392L, 12538L))))

cat("\n  Full per-phase output (shows what mean_yes/mean_no would have given):\n")
print(as.data.frame(rr %>% select(phase, x_yes, n_yes, x_no, n_no,
                                 mean_yes, mean_no, rs_mean, rs_l, rs_u)),
      row.names = FALSE, digits = 4)

# what the product estimator WOULD give for the L2G row
pr <- rr %>% filter(phase %in% c("I", "II", "III"))
prod_yes <- prod(pr$mean_yes); prod_no <- prod(pr$mean_no)
cat(sprintf("\n  If the product estimator were used: P(S)_yes=%.4f P(S)_no=%.4f -> RS=%.2f\n",
            prod_yes, prod_no, prod_yes / prod_no))
cat(sprintf("  Direct Katz ratio actually reported:                        RS=%.2f\n",
            il$rs_mean))

cat("\n================================================================\n")
cat("  3. pQTL row computed BOTH ways on identical data\n")
cat("================================================================\n\n")

# --- build the pQTL analysis frame exactly as the main figures script does ---
similarity_threshold <- 0.8
min_l2g_share        <- 0.5

df_pqtl_support <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= similarity_threshold,
         !is.na(l2g_share), l2g_share >= min_l2g_share) %>%
  distinct(ti_uid)

ti_best_all <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "", !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup() %>%
  mutate(gensup = ti_uid %in% df_pqtl_support$ti_uid) %>%
  filter(gene %in% pgenes)          # measured-protein universe

# --- (a) compute_rr style: single 2x2, launched vs Phase-I entrants ---------
long_a <- ti_best_all %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
base_a <- ti_best_all %>% filter(!is.na(succ_1_2))
xa1 <- sum(long_a$gensup & long_a$success, na.rm = TRUE)
na1 <- sum(base_a$gensup, na.rm = TRUE)
xa2 <- sum(!long_a$gensup & long_a$success, na.rm = TRUE)
na2 <- sum(!base_a$gensup, na.rm = TRUE)
ra  <- as.data.frame(BinomRatioCI(xa1, na1, xa2, na2, method = "katz"))
cat(sprintf("  (a) compute_rr (main figures): (%d/%d)/(%d/%d)  RS=%.2f (%.2f-%.2f)\n",
            xa1, na1, xa2, na2, ra$est, ra$lwr.ci, ra$upr.ci))

# --- (b) advancement_rr style on the SAME frame -----------------------------
pb_pqtl <- ti_best_all %>%
  mutate(target_status = if_else(gensup, "genetically supported target", "other"),
         similarity = comb_norm)
rr_b <- advancement_rr(pb_pqtl)
il_b <- rr_b %>% filter(phase == "I-Launch")
cat(sprintf("  (b) advancement_rr on same data: (%d/%d)/(%d/%d)  RS=%.2f (%.2f-%.2f)\n",
            il_b$x_yes, il_b$n_yes, il_b$x_no, il_b$n_no,
            il_b$rs_mean, il_b$rs_l, il_b$rs_u))

cat("\n  Per-phase detail for the pQTL set:\n")
print(as.data.frame(rr_b %>% select(phase, x_yes, n_yes, x_no, n_no,
                                    mean_yes, mean_no, rs_mean, rs_l, rs_u)),
      row.names = FALSE, digits = 4)

pr_b <- rr_b %>% filter(phase %in% c("I", "II", "III"))
py <- prod(pr_b$mean_yes); pn <- prod(pr_b$mean_no)
cat(sprintf("\n  Product estimator would give: P(S)_yes=%.4f P(S)_no=%.4f -> RS=%.2f\n",
            py, pn, py / pn))

cat("\n================================================================\n")
cat("  VERDICT\n")
cat("================================================================\n\n")
cat(sprintf("  L2G row  (advancement_rr) : RS = %.2f\n", il$rs_mean))
cat(sprintf("  pQTL row (compute_rr)     : RS = %.2f\n", ra$est))
cat(sprintf("  pQTL row (advancement_rr) : RS = %.2f\n", il_b$rs_mean))
cat(sprintf("\n  Do the two code paths agree on the pQTL row? %s\n",
            ifelse(abs(ra$est - il_b$rs_mean) < 0.01, "YES - same estimator",
                   sprintf("NO - differ by %.2f", abs(ra$est - il_b$rs_mean)))))
cat(sprintf("  Denominator definitions: compute_rr n_yes=%d vs advancement_rr n_yes=%d\n",
            na1, il_b$n_yes))
