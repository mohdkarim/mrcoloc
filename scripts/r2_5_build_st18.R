#!/usr/bin/env Rscript
# ============================================================================
# Build ST18: directional sensitivity of the pQTL relative success estimate
# ============================================================================
# Same column format as ST4 (Figure 1a data) so the two can be read together:
#   panel_group | source_label | rs_estimate | rs_lwr_95ci | rs_upr_95ci | count_string
#
# Scenario B is the headline conservative estimate: directionally non-aligned
# launched pairs are removed from the numerator while the denominator is held at
# all pQTL-supported Phase I pairs. Direction was curated only among launched
# pairs, so discordant failures among the supported-but-not-launched pairs remain
# in the denominator - the estimate is therefore a lower bound.
#
# All counts are derived at run time; none are hardcoded. After the background
# correction (see peer_review/background_bug_and_fix.txt) the denominator is 49
# rather than 46 and scenario B gives 3.28 (2.22-4.84) rather than 3.49.
# ============================================================================

suppressPackageStartupMessages({library(tidyverse); library(DescTools); library(openxlsx)})

# [R2.10 FIX] The background counts and the supported denominator were previously
# hardcoded (X_BG <- 757; N_BG <- 7154, with 23/46 and 17/46 written into the
# rows below), so this script never read pgenes and would have kept building ST18
# from the pre-correction numbers. All four are now derived, and the aligned
# count comes from the curated alignment table rather than a literal.
merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic       <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes      <- readRDS("data_raw/pgenes.rds")

supp_ti <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= 0.8, !is.na(l2g_share), l2g_share >= 0.5) %>%
  distinct(ti_uid) %>% pull(ti_uid)

ti_best <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "", !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat), gene %in% pgenes) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup() %>%
  mutate(gensup = ti_uid %in% supp_ti)

.long <- ti_best %>% filter(!is.na(succ_3_a))
.base <- ti_best %>% filter(!is.na(succ_1_2))
X_BG  <- sum(!.long$gensup & .long$succ_3_a, na.rm = TRUE)
N_BG  <- sum(!.base$gensup, na.rm = TRUE)
N_SUP <- sum(.base$gensup, na.rm = TRUE)
X_SUP <- sum(.long$gensup & .long$succ_3_a, na.rm = TRUE)
rm(merge3_pqtl); gc(verbose = FALSE)

aln         <- read_tsv("output/r2_5_alignment_table.tsv", show_col_types = FALSE)
N_ALIGNED   <- sum(aln$alignment == "aligned")
N_NONALIGN  <- sum(aln$alignment != "aligned")

cat(sprintf("  derived: background %d/%d | supported %d/%d | aligned %d of %d launched\n",
            X_BG, N_BG, X_SUP, N_SUP, N_ALIGNED, nrow(aln)))

L2G  <- c(127, 484, 1392, 12538)   # ST4 reference row; full universe, unaffected

row_of <- function(panel, label, x, n, x2 = X_BG, n2 = N_BG) {
  rr <- as.data.frame(BinomRatioCI(x, n, x2, n2, method = "katz"))
  tibble(panel_group = panel, source_label = label,
         rs_estimate = rr$est, rs_lwr_95ci = rr$lwr.ci, rs_upr_95ci = rr$upr.ci,
         count_string = sprintf("(%d/%d)/(%d/%d)", x, n, x2, n2))
}

ST18 <- bind_rows(
  row_of("Directional sensitivity", "pQTL: all supported pairs (as published)", X_SUP, N_SUP),
  row_of("Directional sensitivity", "pQTL: directionally aligned only",         N_ALIGNED, N_SUP),
  row_of("Directional sensitivity", "pQTL: aligned, aligned denominator",       N_ALIGNED, N_SUP - N_NONALIGN),
  row_of("Reference (from ST4)",    "L2G share: >= 0.5",
         L2G[1], L2G[2], L2G[3], L2G[4])
)

print(as.data.frame(ST18 %>% mutate(across(starts_with("rs_"), ~round(., 3)))),
      row.names = FALSE, right = FALSE)

write_tsv(ST18, "output/ST18_directional_sensitivity.tsv")

wb <- createWorkbook()
sh <- "ST18 - Directional_Sensitivity"
addWorksheet(wb, sh)
title_style  <- createStyle(fontSize = 12, textDecoration = "bold")
header_style <- createStyle(textDecoration = "bold", border = "bottom")
writeData(wb, sh,
  "Supplementary Table 18: Directional Sensitivity of pQTL Relative Success",
  startRow = 1, startCol = 1)
addStyle(wb, sh, title_style, rows = 1, cols = 1)
writeData(wb, sh, ST18, startRow = 3, startCol = 1, headerStyle = header_style)
freezePane(wb, sh, firstActiveRow = 4, firstActiveCol = 1)
setColWidths(wb, sh, cols = 1:ncol(ST18), widths = "auto")
saveWorkbook(wb, "output/ST18_directional_sensitivity.xlsx", overwrite = TRUE)

cat("\n-> output/ST18_directional_sensitivity.tsv\n")
cat("-> output/ST18_directional_sensitivity.xlsx\n\n")
