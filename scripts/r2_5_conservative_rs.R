#!/usr/bin/env Rscript
# ============================================================================
# R2.5: conservative RS excluding directionally non-aligned launched pairs
# ============================================================================
# Rationale: the reviewer's position is that a pQTL association pointing the
# opposite way to the therapeutic mechanism should not "count" as support. We
# therefore remove non-aligned pairs from the NUMERATOR (successes) while
# holding the denominator at all pQTL-supported Phase I pairs.
#
# This is deliberately conservative. Direction was curated only among the
# launched pairs; the supported-but-not-launched pairs were not curated, so any
# discordant failures among them are NOT removed. The analysis therefore
# penalises successes only and gives no offsetting credit, producing a lower
# bound on RS.
#
# Variants reported (counts derived at run time, not hardcoded):
#   A. published    all launched supported pairs / all supported Phase I pairs
#   B. conservative aligned launched only / all supported Phase I pairs
#   C. consistent   aligned launched only / supported Phase I minus non-aligned
#                   (shown for completeness; less conservative than B)
#
# Post-correction values are A 4.44 (3.27-6.02), B 3.28 (2.22-4.84),
# C 3.74 (2.57-5.44) on 17 aligned of 23 launched, denominator 49.
# ============================================================================

suppressPackageStartupMessages({library(tidyverse); library(DescTools)})

aln <- read_tsv("output/r2_5_alignment_table.tsv", show_col_types = FALSE)

n_launched   <- nrow(aln)
non_aligned  <- aln %>% filter(alignment != "aligned")
n_aligned    <- sum(aln$alignment == "aligned")

cat(sprintf("\nlaunched pQTL-supported pairs: %d\n", n_launched))
cat(sprintf("non-aligned: %d (%s)\n", nrow(non_aligned),
            paste(non_aligned$Target, collapse=", ")))
cat(sprintf("aligned: %d\n\n", n_aligned))

# [R2.10 FIX] These three counts were previously hardcoded as
#   X_BG <- 757; N_BG <- 7154; N_SUP <- 46
# which meant this script did not read pgenes at all and would have kept
# emitting the pre-correction values (RS 3.49) after the background was fixed,
# silently desynchronising R2.5 from every other number in the response. They
# are now derived from the data, exactly as the figure and supplement scripts do.
cat("  Deriving background counts from the current pgenes...\n")
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

X_BG  <- sum(!.long$gensup & .long$succ_3_a, na.rm = TRUE)   # launched, unsupported
N_BG  <- sum(!.base$gensup, na.rm = TRUE)                    # Phase I, unsupported
N_SUP <- sum(.base$gensup, na.rm = TRUE)                      # Phase I, supported
X_SUP <- sum(.long$gensup & .long$succ_3_a, na.rm = TRUE)     # launched, supported

cat(sprintf("    background %d/%d | supported %d/%d\n", X_BG, N_BG, X_SUP, N_SUP))
rm(merge3_pqtl); gc(verbose = FALSE)

rs <- function(x, n, label, note="") {
  rr <- as.data.frame(BinomRatioCI(x, n, X_BG, N_BG, method="katz"))
  cat(sprintf("  %-42s RS = %.2f (95%% CI %.2f-%.2f)   (%d/%d = %.1f%%)  %s\n",
              label, rr$est, rr$lwr.ci, rr$upr.ci, x, n, 100*x/n, note))
  tibble(scenario=label, x=x, n=n, rs=rr$est, lwr=rr$lwr.ci, upr=rr$upr.ci)
}

cat(sprintf("=== RS against the measured-protein background (%d/%d) ===\n\n", X_BG, N_BG))
res <- bind_rows(
  rs(X_SUP, N_SUP, "A. published"),
  rs(n_aligned, N_SUP, "B. conservative (drop non-aligned from numerator)"),
  rs(n_aligned, N_SUP - nrow(non_aligned),
     "C. consistent (drop from both)", "<- less conservative than B")
)

cat("\n=== benchmarks for comparison ===\n")
# L2G rows use the full T-I universe and are unaffected by the pgenes correction.
cat("  L2G share >= 0.5 alone (ST4)              RS = 2.36 (95% CI 2.02-2.77)\n")
cat("  Minikel et al headline, any genetic support RS = 2.6\n")

cat("\n=== interpretation ===\n")
a <- res %>% filter(str_starts(scenario,"A"))
b <- res %>% filter(str_starts(scenario,"B"))
cat(sprintf("  Removing every directionally non-aligned success lowers RS from %.2f to %.2f\n", a$rs, b$rs))
cat(sprintf("  (95%% CI %.2f-%.2f). The lower bound excludes 1 and the CIs overlap the published\n", b$lwr, b$upr))
cat("  estimate. NOTE: do not claim the conservative estimate significantly exceeds L2G\n")
cat("  alone - its lower bound overlaps the L2G upper bound of 2.77.\n")

write_tsv(res, "output/r2_5_conservative_rs.tsv")
cat("\n-> output/r2_5_conservative_rs.tsv\n\n")
