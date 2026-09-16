#!/usr/bin/env Rscript
# ============================================================================
# R1.6: Is pQTL support independent of OMIM and Genebass evidence?
# ============================================================================
# Reviewer asked to "test pQTLs in combination with the other genetic evidence
# sources, mainly Genebass and OMIM, to show that pQTLs are independent".
#
# Problems with the earlier draft analysis:
#   - "pQTL excluding OMIM" removed OMIM pairs from the SUPPORTED set but left
#     them in the BACKGROUND, which is enriched -> biased, not an independence
#     test.
#   - Union analyses (pQTL OR OMIM) are combined-evidence estimates; they say
#     nothing about independence.
#   - No factorial cross-classification, which is what "in combination" means.
#   - Figure 1a computes pQTL in the measured-protein universe but OMIM/Genebass
#     in the full universe -> denominators are not matched.
#
# This script does it properly:
#   1. Matched-universe marginal RS for each source (fixes denominators)
#   2. 2x2 factorial: pQTL x OMIM and pQTL x Genebass, all vs double-negative
#   3. Stratified RS of pQTL within each source's +/- strata
#   4. Breslow-Day homogeneity test = formal test of no interaction
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DescTools)
})

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic       <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes      <- readRDS("data_raw/pgenes.rds")

STHR  <- 0.8
L2GT  <- 0.5

# ---- Full universe (no pgenes restriction) and matched universe ------------
full_base <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "", !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none")

matched_base <- full_base %>% filter(gene %in% pgenes)

# ---- Support sets ----------------------------------------------------------
pqtl_ids <- full_base %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= STHR, !is.na(l2g_share), l2g_share >= L2GT) %>%
  distinct(ti_uid) %>% pull(ti_uid)

omim_ids <- full_base %>%
  filter(assoc_source == "OMIM", comb_norm >= STHR) %>%
  distinct(ti_uid) %>% pull(ti_uid)

genebass_ids <- full_base %>%
  filter(assoc_source == "Genebass", comb_norm >= STHR) %>%
  distinct(ti_uid) %>% pull(ti_uid)

best_of <- function(df) {
  df %>%
    mutate(hp = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
                          !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
    arrange(ti_uid, desc(hp), desc(comb_norm)) %>%
    group_by(ti_uid) %>% slice(1) %>% ungroup()
}

ti_full    <- best_of(full_base)
ti_matched <- best_of(matched_base)

# ---- Cohort framing: among Phase I entrants, did the pair launch? ----------
cohort <- function(ti) {
  ti %>% filter(!is.na(succ_1_2)) %>%
    mutate(launched = !is.na(succ_3_a) & succ_3_a)
}

rs2 <- function(dat, sup, lab, min_cell = 5) {
  d <- dat %>% mutate(.s = sup)
  x1 <- sum(d$.s & d$launched);  n1 <- sum(d$.s)
  x2 <- sum(!d$.s & d$launched); n2 <- sum(!d$.s)
  if (n1 < min_cell || x1 == 0 || x2 == 0) {
    cat(sprintf("  %-44s UNDERPOWERED (%d/%d) vs (%d/%d)\n", lab, x1, n1, x2, n2))
    return(tibble(label = lab, est = NA, lwr.ci = NA, upr.ci = NA, x1, n1, x2, n2))
  }
  rr <- as.data.frame(BinomRatioCI(x1, n1, x2, n2, method = "katz"))
  cat(sprintf("  %-44s RS=%5.2f (%4.2f-%5.2f)  (%d/%d)=%5.1f%% vs (%d/%d)=%4.1f%%\n",
              lab, rr$est, rr$lwr.ci, rr$upr.ci, x1, n1, 100*x1/n1, x2, n2, 100*x2/n2))
  tibble(label = lab, est = rr$est, lwr.ci = rr$lwr.ci, upr.ci = rr$upr.ci, x1, n1, x2, n2)
}

# ============================================================================
cat("\n================================================================\n")
cat("  1. DENOMINATOR CHECK: Figure 1a mixes universes\n")
cat("================================================================\n\n")
cat("  As published in ST4, the 'by genetic evidence source' panel uses:\n")
cat("    pQTL     background  7,154  (measured proteins)\n")
cat("    OMIM     background 12,830  (full universe)\n")
cat("    Genebass background 12,976  (full universe)\n\n")
cat("  Recomputed in the FULL universe (should reproduce ST4):\n")
cf <- cohort(ti_full)
invisible(rs2(cf, cf$ti_uid %in% pqtl_ids,     "pQTL      (full universe)"))
invisible(rs2(cf, cf$ti_uid %in% omim_ids,     "OMIM      (full universe)"))
invisible(rs2(cf, cf$ti_uid %in% genebass_ids, "Genebass  (full universe)"))

cat("\n  Recomputed in the MATCHED universe (measured proteins, like pQTL):\n")
cm <- cohort(ti_matched)
marg <- bind_rows(
  rs2(cm, cm$ti_uid %in% pqtl_ids,     "pQTL      (matched universe)"),
  rs2(cm, cm$ti_uid %in% omim_ids,     "OMIM      (matched universe)"),
  rs2(cm, cm$ti_uid %in% genebass_ids, "Genebass  (matched universe)")
)

# ============================================================================
cat("\n================================================================\n")
cat("  2. OVERLAP between evidence sources (matched universe, Phase I)\n")
cat("================================================================\n\n")
cm <- cm %>% mutate(
  pq = ti_uid %in% pqtl_ids,
  om = ti_uid %in% omim_ids,
  gb = ti_uid %in% genebass_ids
)
cat(sprintf("  pQTL+ pairs at Phase I: %d\n", sum(cm$pq)))
cat(sprintf("    also OMIM+     : %d (%.0f%%)\n", sum(cm$pq & cm$om),
            100*sum(cm$pq & cm$om)/sum(cm$pq)))
cat(sprintf("    also Genebass+ : %d (%.0f%%)\n", sum(cm$pq & cm$gb),
            100*sum(cm$pq & cm$gb)/sum(cm$pq)))
cat(sprintf("    neither        : %d (%.0f%%)\n", sum(cm$pq & !cm$om & !cm$gb),
            100*sum(cm$pq & !cm$om & !cm$gb)/sum(cm$pq)))

# ============================================================================
factorial_2x2 <- function(dat, other_flag, other_name) {
  cat(sprintf("\n================================================================\n"))
  cat(sprintf("  3. 2x2 FACTORIAL: pQTL x %s\n", other_name))
  cat(sprintf("================================================================\n\n"))
  d <- dat %>% mutate(.o = other_flag,
    cell = case_when( pq &  .o ~ sprintf("pQTL+ / %s+", other_name),
                      pq & !.o ~ sprintf("pQTL+ / %s-", other_name),
                     !pq &  .o ~ sprintf("pQTL- / %s+", other_name),
                     TRUE      ~ sprintf("pQTL- / %s-", other_name)))
  tab <- d %>% group_by(cell) %>%
    summarise(n = n(), launched = sum(launched), .groups = "drop") %>%
    mutate(rate = sprintf("%.1f%%", 100*launched/n))
  print(as.data.frame(tab), row.names = FALSE)

  refcell <- sprintf("pQTL- / %s-", other_name)
  ref <- tab %>% filter(cell == refcell)
  cat(sprintf("\n  reference: %s = %d/%d = %.1f%%\n\n",
              refcell, ref$launched, ref$n, 100*ref$launched/ref$n))
  # [ST19] cell-vs-reference estimates are now collected, not just printed
  cellrs <- bind_rows(lapply(setdiff(tab$cell, refcell), function(cc) {
    r <- tab %>% filter(cell == cc)
    if (r$n < 5 || r$launched == 0) {
      cat(sprintf("  %-30s vs ref   UNDERPOWERED (%d/%d)\n", cc, r$launched, r$n))
      return(tibble(cell = cc, est = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_))
    }
    rr <- as.data.frame(BinomRatioCI(r$launched, r$n, ref$launched, ref$n, method = "katz"))
    cat(sprintf("  %-30s vs ref   RS=%5.2f (%4.2f-%5.2f)\n",
                cc, rr$est, rr$lwr.ci, rr$upr.ci))
    tibble(cell = cc, est = rr$est, lwr.ci = rr$lwr.ci, upr.ci = rr$upr.ci)
  }))

  # --- stratified pQTL effect within each stratum of the other source
  # [ST19] results collected as well as printed
  cat(sprintf("\n  Stratified: pQTL+ vs pQTL- within each %s stratum\n\n", other_name))
  strat <- bind_rows(lapply(c(TRUE, FALSE), function(lev) {
    sub <- d %>% filter(.o == lev)
    rs2(sub, sub$pq, sprintf("pQTL effect | %s %s", other_name,
                             ifelse(lev, "POSITIVE", "NEGATIVE")))
  }))

  # --- Breslow-Day test of homogeneity of the pQTL OR across strata
  arr <- array(NA_integer_, dim = c(2,2,2),
               dimnames = list(pQTL = c("yes","no"),
                               outcome = c("launch","no"),
                               stratum = c("other_pos","other_neg")))
  for (i in seq_along(c(TRUE, FALSE))) {
    lev <- c(TRUE, FALSE)[i]
    sub <- d %>% filter(.o == lev)
    arr["yes","launch",i] <- sum(sub$pq & sub$launched)
    arr["yes","no",i]     <- sum(sub$pq & !sub$launched)
    arr["no","launch",i]  <- sum(!sub$pq & sub$launched)
    arr["no","no",i]      <- sum(!sub$pq & !sub$launched)
  }
  cat("\n  2x2x2 table:\n"); print(arr)
  bd <- try(BreslowDayTest(arr), silent = TRUE)
  bd_p <- NA_real_
  if (!inherits(bd, "try-error")) {
    bd_p <- as.numeric(bd$p.value)
    cat(sprintf("\n  Breslow-Day homogeneity of pQTL OR across %s strata: p = %.3f\n",
                other_name, bd_p))
    if (bd_p < 0.05) {
      cat("  (p < 0.05 => SIGNIFICANT heterogeneity: the pQTL effect is NOT constant\n")
      cat("   across ", other_name, " strata. Do NOT describe the sources as independent;\n", sep = "")
      cat("   use complementary-coverage wording. See CLAUDE.md error 2.)\n")
    } else {
      cat("  (p > 0.05 => no evidence of interaction, though note the test is\n")
      cat("   underpowered with this few pairs per cell)\n")
    }
  } else {
    cat("\n  Breslow-Day test not computable (empty cells)\n")
  }
  # [ST19] return everything needed to build the supplementary table
  invisible(list(other = other_name, tab = tab, ref = ref,
                 cellrs = cellrs, strat = strat, bd_p = bd_p))
}

t_om <- factorial_2x2(cm, cm$om, "OMIM")
t_gb <- factorial_2x2(cm, cm$gb, "Genebass")

# ============================================================================
cat("\n================================================================\n")
cat("  4. pQTL with NEITHER OMIM NOR Genebass support\n")
cat("================================================================\n\n")
sub <- cm %>% filter(!om & !gb)
cat(sprintf("  [stratum: %d Phase I pairs]\n", nrow(sub)))
neither <- rs2(sub, sub$pq, "pQTL effect | no OMIM and no Genebass")

saveRDS(list(marginal = marg, omim = t_om, genebass = t_gb, neither = neither),
        "output/r1_6_independence.rds")
cat("\n  -> output/r1_6_independence.rds\n")

# ============================================================================
# [ST19] Supplementary table: pQTL x OMIM and pQTL x Genebass factorial,
# stratified estimates and Breslow-Day homogeneity tests.
# Column names follow ST4's convention so the sheet can be dropped straight
# into generate_mrcoloc_supplement.R.
# ============================================================================
row_ <- function(panel, label, x, n, est = NA, lwr = NA, upr = NA, note = "") {
  tibble(panel_group = panel, source_label = label,
         count_string = if (is.na(x)) "" else sprintf("(%d/%d)", x, n),
         rate = if (is.na(x)) "" else sprintf("%.1f%%", 100 * x / n),
         rs_estimate = est, rs_lwr_95ci = lwr, rs_upr_95ci = upr, note = note)
}

st19 <- bind_rows(lapply(list(t_om, t_gb), function(f) {
  panel <- sprintf("pQTL x %s", f$other)
  cells <- f$tab %>% left_join(f$cellrs, by = "cell")
  bind_rows(
    bind_rows(lapply(seq_len(nrow(cells)), function(i) {
      r <- cells[i, ]
      row_(panel, r$cell, r$launched, r$n, r$est, r$lwr.ci, r$upr.ci,
           if (identical(r$cell, f$ref$cell)) "reference cell" else "vs reference cell")
    })),
    bind_rows(lapply(seq_len(nrow(f$strat)), function(i) {
      s <- f$strat[i, ]
      row_(sprintf("%s (stratified)", panel), s$label, s$x1, s$n1,
           s$est, s$lwr.ci, s$upr.ci, sprintf("comparator (%d/%d)", s$x2, s$n2))
    })),
    row_(sprintf("%s (stratified)", panel),
         sprintf("Breslow-Day homogeneity of the pQTL odds ratio across %s strata", f$other),
         NA, NA, NA, NA, NA, sprintf("p = %.3f", f$bd_p))
  )
}))
st19 <- bind_rows(st19,
  row_("Neither source present", neither$label, neither$x1, neither$n1,
       neither$est, neither$lwr.ci, neither$upr.ci,
       sprintf("comparator (%d/%d)", neither$x2, neither$n2)))

write_tsv(st19, "output/ST19_pqtl_omim_genebass_factorial.tsv")
cat("  -> output/ST19_pqtl_omim_genebass_factorial.tsv\n\n")
