#!/usr/bin/env Rscript
# ============================================================================
# R1.5: RS of pQTL support IRRESPECTIVE of L2G share
# ============================================================================
# Throughout the manuscript and supplement, pQTL+ requires L2G share >= 0.5
# (see Key sheet: "pQTL+ TIs require ALL of: ... L2G share >= 0.5"). This
# script reports the estimate when that requirement is dropped, i.e. pQTL
# support = Bonferroni-significant MR + MeSH similarity >= 0.8 only.
#
# Reported in both universes used in the supplement:
#   Universe 1 - all Phase I TIs where the target is a measured protein
#   Universe 2 - subset of Universe 1 with existing genetic evidence
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
PICCOLO_H4     <- 0.9

base <- merge3_pqtl %>%
  filter(gene %in% pgenes,
         !is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none")

# --- Two pQTL support definitions -------------------------------------------
# (a) paper: the pQTL association itself must carry l2g_share >= 0.5
pqtl_paper <- base %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE) %>%
  distinct(ti_uid) %>% pull(ti_uid)

# (b) irrespective of L2G: MR + MeSH match only
pqtl_anyl2g <- base %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR) %>%
  distinct(ti_uid) %>% pull(ti_uid)

# --- Universe 2 membership: any existing genetic evidence -------------------
gen_ev <- base %>%
  filter(comb_norm >= SIMILARITY_THR) %>%
  filter(
    (assoc_source == "OTG"      & !is.na(l2g_share) & l2g_share >= MIN_L2G_SHARE) |
    (assoc_source == "OMIM") |
    (assoc_source == "PICCOLO"  & !is.na(pic_h4)    & pic_h4    >= PICCOLO_H4) |
    (assoc_source == "Genebass")
  ) %>%
  distinct(ti_uid) %>% pull(ti_uid)

# --- Best row per TI --------------------------------------------------------
ti_best <- base %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_
  )) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup()

rs_of <- function(dat, sup_ids, label) {
  d  <- dat %>% mutate(.sup = ti_uid %in% sup_ids)
  ph <- d %>% filter(!is.na(succ_1_2))
  lg <- d %>% filter(!is.na(succ_3_a))

  x1 <- sum(lg$.sup  & lg$succ_3_a, na.rm = TRUE); n1 <- sum(ph$.sup,  na.rm = TRUE)
  x2 <- sum(!lg$.sup & lg$succ_3_a, na.rm = TRUE); n2 <- sum(!ph$.sup, na.rm = TRUE)

  rr <- as.data.frame(BinomRatioCI(x1, n1, x2, n2, method = "katz"))
  cat(sprintf("  %-42s RS=%5.2f (%4.2f-%5.2f)  (%d/%d)=%4.1f%% vs (%d/%d)=%4.1f%%\n",
              label, rr$est, rr$lwr.ci, rr$upr.ci,
              x1, n1, 100 * x1 / n1, x2, n2, 100 * x2 / n2))
  tibble(label, est = rr$est, lwr.ci = rr$lwr.ci, upr.ci = rr$upr.ci,
         x_sup = x1, n_sup = n1, x_unsup = x2, n_unsup = n2)
}

cat("\n================================================================\n")
cat("  UNIVERSE 1 (measured-protein background)\n")
cat("================================================================\n\n")
u1 <- bind_rows(
  rs_of(ti_best, pqtl_paper,  "pQTL+ requiring L2G >= 0.5 (paper)"),
  rs_of(ti_best, pqtl_anyl2g, "pQTL+ irrespective of L2G share")
)

cat("\n================================================================\n")
cat("  UNIVERSE 2 (restricted to TIs with existing genetic evidence)\n")
cat("================================================================\n\n")
ti_u2 <- ti_best %>% filter(ti_uid %in% gen_ev)
cat(sprintf("  [Universe 2: %d TI pairs, %d at Phase I]\n\n",
            nrow(ti_u2), sum(!is.na(ti_u2$succ_1_2))))
u2 <- bind_rows(
  rs_of(ti_u2, pqtl_paper,  "pQTL+ requiring L2G >= 0.5 (paper)"),
  rs_of(ti_u2, pqtl_anyl2g, "pQTL+ irrespective of L2G share")
)

cat("\n================================================================\n")
cat("  COVERAGE vs PRECISION trade-off (Universe 1)\n")
cat("================================================================\n\n")
cat(sprintf("  pQTL+ TI pairs, L2G >= 0.5 required : %4d\n", length(pqtl_paper)))
cat(sprintf("  pQTL+ TI pairs, L2G unrestricted    : %4d  (+%.0f%%)\n",
            length(pqtl_anyl2g),
            100 * (length(pqtl_anyl2g) - length(pqtl_paper)) / length(pqtl_paper)))
cat(sprintf("  all L2G>=0.5 pairs are a subset      : %s\n",
            all(pqtl_paper %in% pqtl_anyl2g)))
cat(sprintf("\n  Phase I pairs : %d -> %d\n", u1$n_sup[1], u1$n_sup[2]))
cat(sprintf("  Launched pairs: %d -> %d\n", u1$x_sup[1], u1$x_sup[2]))
cat(sprintf("  RS            : %.2f -> %.2f\n", u1$est[1], u1$est[2]))

cat("\n================================================================\n")
cat("  Why so few pQTL+ pairs survive a LOW L2G filter\n")
cat("================================================================\n\n")
ti_maxl2g <- base %>%
  group_by(ti_uid) %>%
  summarise(max_l2g = suppressWarnings(max(l2g_share, na.rm = TRUE)), .groups = "drop") %>%
  mutate(max_l2g = ifelse(is.infinite(max_l2g), NA_real_, max_l2g))

tab <- ti_best %>%
  filter(ti_uid %in% pqtl_anyl2g) %>%
  left_join(ti_maxl2g, by = "ti_uid") %>%
  mutate(band = case_when(
    is.na(max_l2g)  ~ "L2G missing",
    max_l2g < 0.25  ~ "< 0.25",
    max_l2g < 0.50  ~ "0.25-0.50",
    max_l2g < 0.75  ~ "0.50-0.75",
    TRUE            ~ ">= 0.75"
  )) %>%
  group_by(band) %>%
  summarise(
    pqtl_pairs   = n(),
    reached_ph1  = sum(!is.na(succ_1_2)),
    launched     = sum(!is.na(succ_3_a) & succ_3_a, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(tab), row.names = FALSE)
cat(sprintf("\n  TOTAL: %d pQTL+ pairs | %d reached Phase I | %d launched\n",
            sum(tab$pqtl_pairs), sum(tab$reached_ph1), sum(tab$launched)))

saveRDS(list(u1 = u1, u2 = u2, bands = tab), "output/r1_5_irrespective_l2g.rds")
cat("\n  -> output/r1_5_irrespective_l2g.rds\n\n")
