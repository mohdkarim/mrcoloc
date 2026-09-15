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
  for (cc in setdiff(tab$cell, refcell)) {
    r <- tab %>% filter(cell == cc)
    if (r$n < 5 || r$launched == 0) {
      cat(sprintf("  %-30s vs ref   UNDERPOWERED (%d/%d)\n", cc, r$launched, r$n)); next
    }
    rr <- as.data.frame(BinomRatioCI(r$launched, r$n, ref$launched, ref$n, method = "katz"))
    cat(sprintf("  %-30s vs ref   RS=%5.2f (%4.2f-%5.2f)\n",
                cc, rr$est, rr$lwr.ci, rr$upr.ci))
  }

  # --- stratified pQTL effect within each stratum of the other source
  cat(sprintf("\n  Stratified: pQTL+ vs pQTL- within each %s stratum\n\n", other_name))
  for (lev in c(TRUE, FALSE)) {
    sub <- d %>% filter(.o == lev)
    rs2(sub, sub$pq, sprintf("pQTL effect | %s %s", other_name,
                             ifelse(lev, "POSITIVE", "NEGATIVE")))
  }

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
  if (!inherits(bd, "try-error")) {
    cat(sprintf("\n  Breslow-Day homogeneity of pQTL OR across %s strata: p = %.3f\n",
                other_name, bd$p.value))
    cat("  (p > 0.05 => no evidence of interaction => pQTL effect consistent\n")
    cat("   across strata, i.e. acts independently of ", other_name, ")\n", sep = "")
  } else {
    cat("\n  Breslow-Day test not computable (empty cells)\n")
  }
  invisible(tab)
}

t_om <- factorial_2x2(cm, cm$om, "OMIM")
t_gb <- factorial_2x2(cm, cm$gb, "Genebass")

# ============================================================================
cat("\n================================================================\n")
cat("  4. pQTL with NEITHER OMIM NOR Genebass support\n")
cat("================================================================\n\n")
sub <- cm %>% filter(!om & !gb)
cat(sprintf("  [stratum: %d Phase I pairs]\n", nrow(sub)))
invisible(rs2(sub, sub$pq, "pQTL effect | no OMIM and no Genebass"))

saveRDS(list(marginal = marg, omim = t_om, genebass = t_gb),
        "output/r1_6_independence.rds")
cat("\n  -> output/r1_6_independence.rds\n\n")
