#!/usr/bin/env Rscript
# ============================================================================
# R1.5: RS of pQTL-supported pairs that LACK L2G support
# ============================================================================
# Reviewer 1 asked for L2G thresholds in the "less than" direction:
#   "< 0.25, L2G < 0.50 and L2G < 0.75 - what is the relative success of
#    pQTL-supported pairs that lack L2G support?"
#
# CRITICAL: l2g_share is row-level (gene-trait-study triplet). "Lacks L2G
# support" must be defined at TI level as max(l2g_share) < threshold, i.e. NO
# study assigns this gene a high L2G share for a matched trait. Filtering rows
# and then de-duplicating ti_uid double-counts pairs across strata.
#
# Analyses:
#   A. Literal request  - RS of pQTL+ & maxL2G<thr vs standard background
#   B. Matched/conditional - within the L2G-negative stratum only, pQTL+ vs pQTL-
#   C. Disjoint bands   - dose-response without double-counting
#   D. 2x2 factorial    - (pQTL +/-) x (L2G +/-), all vs double-negative
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DescTools)
})

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic       <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes      <- readRDS("data_raw/pgenes.rds")

SIMILARITY_THR <- 0.8
MIN_L2G_SHARE  <- 0.5

# --- Universe: measured proteins, genetic-insight indications, best row per TI
universe <- merge3_pqtl %>%
  filter(gene %in% pgenes,
         !is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none")

# TI-level max L2G share across ALL rows for that pair (any evidence source).
# All-NA -> NA, treated as L2G-negative below.
ti_maxl2g <- universe %>%
  group_by(ti_uid) %>%
  summarise(max_l2g = suppressWarnings(max(l2g_share, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(max_l2g = ifelse(is.infinite(max_l2g), NA_real_, max_l2g))

# pQTL support = MR significant + MeSH matched (NO L2G condition)
pqtl_tis <- universe %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR) %>%
  distinct(ti_uid) %>% pull(ti_uid)

# Best row per TI for phase outcomes
ti_best <- universe %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_
  )) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(ti_maxl2g, by = "ti_uid") %>%
  mutate(has_pqtl = ti_uid %in% pqtl_tis)

cat(sprintf("\nUniverse: %d TI pairs | Phase I: %d | pQTL-supported: %d\n",
            nrow(ti_best),
            sum(!is.na(ti_best$succ_1_2)),
            sum(ti_best$has_pqtl)))

# --- RS helper: supported vs unsupported within a given universe -------------
rs_of <- function(dat, sup_flag, label, min_n = 5) {
  d  <- dat %>% mutate(.sup = sup_flag)
  ph <- d %>% filter(!is.na(succ_1_2))
  lg <- d %>% filter(!is.na(succ_3_a))

  x1 <- sum(lg$.sup  & lg$succ_3_a, na.rm = TRUE); n1 <- sum(ph$.sup,  na.rm = TRUE)
  x2 <- sum(!lg$.sup & lg$succ_3_a, na.rm = TRUE); n2 <- sum(!ph$.sup, na.rm = TRUE)

  if (n1 < min_n || n2 < min_n || x1 == 0 || x2 == 0) {
    cat(sprintf("  %-46s  UNDERPOWERED  (%d/%d) vs (%d/%d)\n", label, x1, n1, x2, n2))
    return(tibble(label, est = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_,
                  x1, n1, x2, n2, flag = "underpowered"))
  }
  rr <- as.data.frame(BinomRatioCI(x1, n1, x2, n2, method = "katz"))
  cat(sprintf("  %-46s  RS=%5.2f (%4.2f-%5.2f)   (%d/%d)=%4.1f%% vs (%d/%d)=%4.1f%%\n",
              label, rr$est, rr$lwr.ci, rr$upr.ci,
              x1, n1, 100 * x1 / n1, x2, n2, 100 * x2 / n2))
  tibble(label, est = rr$est, lwr.ci = rr$lwr.ci, upr.ci = rr$upr.ci,
         x1, n1, x2, n2, flag = "ok")
}

# L2G-negative indicator at a threshold (NA max_l2g counts as negative)
l2g_neg <- function(v, thr) is.na(v) | v < thr

# ============================================================================
cat("\n================================================================\n")
cat("  A. LITERAL REQUEST\n")
cat("     pQTL-supported AND max L2G < thr, vs standard background\n")
cat("     (background = all pairs not in that supported set)\n")
cat("================================================================\n\n")

A <- map_dfr(c(0.25, 0.50, 0.75), function(thr) {
  rs_of(ti_best,
        ti_best$has_pqtl & l2g_neg(ti_best$max_l2g, thr),
        sprintf("pQTL+ & max L2G < %.2f", thr))
})

cat("\n  Reference (current paper definition):\n")
A_ref <- rs_of(ti_best,
               ti_best$has_pqtl & !is.na(ti_best$max_l2g) &
                 ti_best$max_l2g >= MIN_L2G_SHARE,
               "pQTL+ & max L2G >= 0.50")

# ============================================================================
cat("\n================================================================\n")
cat("  B. MATCHED / CONDITIONAL  (the stronger test)\n")
cat("     Restrict universe to pairs L2G alone would MISS (max L2G < thr),\n")
cat("     then compare pQTL+ vs pQTL- WITHIN that stratum.\n")
cat("     Answers: among pairs L2G misses, does pQTL still pick winners?\n")
cat("================================================================\n\n")

B <- map_dfr(c(0.25, 0.50, 0.75), function(thr) {
  sub <- ti_best %>% filter(l2g_neg(max_l2g, thr))
  cat(sprintf("  [stratum max L2G < %.2f: %d pairs, %d in Phase I]\n",
              thr, nrow(sub), sum(!is.na(sub$succ_1_2))))
  rs_of(sub, sub$has_pqtl, sprintf("pQTL+ vs pQTL- | max L2G < %.2f", thr))
})

# ============================================================================
cat("\n================================================================\n")
cat("  C. DISJOINT BANDS (dose-response, no double-counting)\n")
cat("================================================================\n\n")

bands <- list(
  list(lab = "pQTL+ & L2G missing",     f = function(v) is.na(v)),
  list(lab = "pQTL+ & L2G < 0.25",      f = function(v) !is.na(v) & v <  0.25),
  list(lab = "pQTL+ & L2G 0.25-0.50",   f = function(v) !is.na(v) & v >= 0.25 & v < 0.50),
  list(lab = "pQTL+ & L2G 0.50-0.75",   f = function(v) !is.na(v) & v >= 0.50 & v < 0.75),
  list(lab = "pQTL+ & L2G >= 0.75",     f = function(v) !is.na(v) & v >= 0.75)
)
C <- map_dfr(bands, function(b)
  rs_of(ti_best, ti_best$has_pqtl & b$f(ti_best$max_l2g), b$lab))

# ============================================================================
cat("\n================================================================\n")
cat("  D. 2x2 FACTORIAL at L2G 0.50\n")
cat("     All cells vs the double-negative reference cell\n")
cat("================================================================\n\n")

ti_best <- ti_best %>%
  mutate(l2g_pos = !is.na(max_l2g) & max_l2g >= MIN_L2G_SHARE,
         cell = case_when(
           has_pqtl  &  l2g_pos ~ "pQTL+ / L2G+",
           !has_pqtl &  l2g_pos ~ "pQTL- / L2G+",
           has_pqtl  & !l2g_pos ~ "pQTL+ / L2G-",
           TRUE                 ~ "pQTL- / L2G-"
         ))

ph <- ti_best %>% filter(!is.na(succ_1_2))
lg <- ti_best %>% filter(!is.na(succ_3_a))

cell_tab <- tibble(cell = c("pQTL+ / L2G+", "pQTL+ / L2G-",
                            "pQTL- / L2G+", "pQTL- / L2G-")) %>%
  rowwise() %>%
  mutate(
    n_phase1 = sum(ph$cell == cell),
    n_launch = sum(lg$cell == cell & lg$succ_3_a, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(rate = ifelse(n_phase1 > 0, n_launch / n_phase1, NA_real_))

print(as.data.frame(cell_tab %>% mutate(rate = sprintf("%.1f%%", 100 * rate))),
      row.names = FALSE)

ref <- cell_tab %>% filter(cell == "pQTL- / L2G-")
cat(sprintf("\n  Reference cell: pQTL- / L2G-  = %d/%d = %.1f%%\n\n",
            ref$n_launch, ref$n_phase1, 100 * ref$rate))

D <- cell_tab %>% filter(cell != "pQTL- / L2G-") %>%
  rowwise() %>%
  mutate(res = list({
    if (n_phase1 < 5 || n_launch == 0) {
      cat(sprintf("  %-46s  UNDERPOWERED  (%d/%d)\n", cell, n_launch, n_phase1))
      tibble(est = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_)
    } else {
      rr <- as.data.frame(BinomRatioCI(n_launch, n_phase1,
                                       ref$n_launch, ref$n_phase1, method = "katz"))
      cat(sprintf("  %-46s  RS=%5.2f (%4.2f-%5.2f)   (%d/%d)=%4.1f%%\n",
                  paste(cell, "vs ref"), rr$est, rr$lwr.ci, rr$upr.ci,
                  n_launch, n_phase1, 100 * n_launch / n_phase1))
      tibble(est = rr$est, lwr.ci = rr$lwr.ci, upr.ci = rr$upr.ci)
    }
  })) %>%
  unnest(res) %>%
  ungroup()

# ============================================================================
cat("\n================================================================\n")
cat("  WITHDRAWN NUMBER CHECK\n")
cat("================================================================\n\n")
cat("  Previously reported (row-level filter, INVALID): RS = 2.60 (1.60-4.22)\n")
a50 <- A %>% filter(grepl("0.50", label))
if (!is.na(a50$est)) {
  cat(sprintf("  Corrected TI-level equivalent:                  RS = %.2f (%.2f-%.2f)\n",
              a50$est, a50$lwr.ci, a50$upr.ci))
} else {
  cat("  Corrected TI-level equivalent:                  UNDERPOWERED\n")
}
cat(sprintf("  (supported pairs: %d launched / %d Phase I)\n", a50$x1, a50$n1))

saveRDS(list(A = A, A_ref = A_ref, B = B, C = C, D = D, cells = cell_tab),
        "output/r1_5_results.rds")
cat("\n  -> output/r1_5_results.rds\n")

# ============================================================================
# [ST22] Supplementary table for Reviewer 1 comment 5: relative success of
# pQTL-supported pairs stratified by their TI-level maximum L2G share, plus the
# literal below-threshold comparisons the reviewer asked for.
#
# The reviewer asked for L2G < 0.25, < 0.50 and < 0.75. Those strata contain too
# few launched pQTL-supported pairs to estimate: hence the underpowered flags.
# NOTE: max L2G share is taken at TI level. Filtering rows on l2g_share and then
# taking distinct(ti_uid) yields non-disjoint strata - see r1_5_diagnostic.R.
#
# Columns follow ST4's convention for direct inclusion in the supplement.
# ============================================================================
st22 <- bind_rows(
  C %>% transmute(panel_group = "By TI-level maximum L2G share band",
                  source_label = label,
                  count_string = sprintf("(%d/%d)/(%d/%d)", x1, n1, x2, n2),
                  rs_estimate = est, rs_lwr_95ci = lwr.ci, rs_upr_95ci = upr.ci,
                  note = flag),
  A %>% transmute(panel_group = "Below-threshold, as requested by the reviewer",
                  source_label = label,
                  count_string = sprintf("(%d/%d)/(%d/%d)", x1, n1, x2, n2),
                  rs_estimate = est, rs_lwr_95ci = lwr.ci, rs_upr_95ci = upr.ci,
                  note = flag),
  B %>% transmute(panel_group = "Below-threshold, matched strata",
                  source_label = label,
                  count_string = sprintf("(%d/%d)/(%d/%d)", x1, n1, x2, n2),
                  rs_estimate = est, rs_lwr_95ci = lwr.ci, rs_upr_95ci = upr.ci,
                  note = flag)
)
write_tsv(st22, "output/ST22_l2g_band_stratification.tsv")
cat("  -> output/ST22_l2g_band_stratification.tsv\n\n")
