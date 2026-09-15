#!/usr/bin/env Rscript
# ============================================================================
# R2.2: ascertainment check - do proteomics platform proteins have an inflated
#       baseline success rate?
# ============================================================================
# Reviewer 2, comment 2: "What is the relative success rate of proteins that were
# assayed by the platforms but did NOT carry a pQTL, but nonetheless were
# successful in advancing to phase 1?"
#
# This script was written for the Nature Medicine revision because the numbers
# quoted in our response existed only as assertions in the response letter - no
# script or saved output computed them, so they could not be regenerated after
# the background correction and were not covered by the code-availability
# statement.
#
# Definitions:
#   measured, no pQTL  = Phase I pairs whose target is in pgenes, without support
#   unmeasured         = the full T-I universe minus the measured set, taken from
#                        the L2G >= 0.5 row of Figure 1a (Phase I 13,022;
#                        launched 1,519 = 127 + 1,392). Those L2G counts use the
#                        full universe and are unaffected by the pgenes fix.
#
# Also reports the per-phase decomposition, which is held in reserve rather than
# used in the response: Preclinical -> Phase I entry is equivalent for supported
# and unsupported measured targets, and the separation appears at Phase II -> III.
#
# Outputs: output/r2_2_ascertainment.tsv
# ============================================================================

suppressPackageStartupMessages({library(tidyverse); library(DescTools)})

project_root <- Sys.getenv("MRCOLOC_ROOT", getwd())
out_dir      <- file.path(project_root, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

merge3_pqtl <- readRDS(file.path(project_root, "data_raw/merge3_pqtl.rds"))
indic       <- read_tsv(file.path(project_root, "data/indic.tsv"), show_col_types = FALSE)
pgenes      <- readRDS(file.path(project_root, "data_raw/pgenes.rds"))

SIM <- 0.8; L2G_THR <- 0.5
FULL_PH1 <- 484 + 12538      # full-universe Phase I entrants (ST4, L2G >= 0.5 row)
FULL_LAU <- 127 + 1392       # full-universe launches

cat("\n================================================================\n")
cat("  R2.2: ascertainment check\n")
cat("================================================================\n\n")

supp_ti <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIM, !is.na(l2g_share), l2g_share >= L2G_THR) %>%
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

long <- ti_best %>% filter(!is.na(succ_3_a))
base <- ti_best %>% filter(!is.na(succ_1_2))

sup_lau  <- sum(long$gensup & long$succ_3_a, na.rm = TRUE)
sup_ph1  <- sum(base$gensup, na.rm = TRUE)
uns_lau  <- sum(!long$gensup & long$succ_3_a, na.rm = TRUE)
uns_ph1  <- sum(!base$gensup, na.rm = TRUE)

meas_ph1 <- sup_ph1 + uns_ph1
meas_lau <- sup_lau + uns_lau
unm_ph1  <- FULL_PH1 - meas_ph1
unm_lau  <- FULL_LAU - meas_lau

show_rs <- function(label, x1, n1, x2, n2) {
  rr <- as.data.frame(BinomRatioCI(x1, n1, x2, n2, method = "katz"))
  cat(sprintf("  %-46s RS = %.2f (95%% CI %.2f-%.2f)\n", label, rr$est, rr$lwr.ci, rr$upr.ci))
  tibble(comparison = label, x1 = x1, n1 = n1, x2 = x2, n2 = n2,
         rs = rr$est, lwr = rr$lwr.ci, upr = rr$upr.ci)
}

cat("--- absolute launch rates from Phase I ---\n")
cat(sprintf("  measured + pQTL support : %4d/%-5d = %.1f%%\n", sup_lau, sup_ph1, 100*sup_lau/sup_ph1))
cat(sprintf("  measured, no pQTL       : %4d/%-5d = %.1f%%\n", uns_lau, uns_ph1, 100*uns_lau/uns_ph1))
cat(sprintf("  not measured            : %4d/%-5d = %.1f%%\n\n", unm_lau, unm_ph1, 100*unm_lau/unm_ph1))

cat("--- the reviewer's question: is the comparator inflated? ---\n")
res <- bind_rows(
  show_rs("measured-no-pQTL vs unmeasured (KEY TEST)", uns_lau, uns_ph1, unm_lau, unm_ph1),
  show_rs("pQTL vs measured-no-pQTL (as published)",   sup_lau, sup_ph1, uns_lau, uns_ph1),
  show_rs("pQTL vs all non-pQTL, full universe",       sup_lau, sup_ph1,
          FULL_LAU - sup_lau, FULL_PH1 - sup_ph1),
  show_rs("pQTL vs unmeasured only",                   sup_lau, sup_ph1, unm_lau, unm_ph1)
)

key <- res %>% slice(1)
cat("\n--- interpretation ---\n")
if (key$upr < 1) {
  cat("  Platform proteins WITHOUT pQTL support are LESS likely to launch than\n")
  cat("  unmeasured targets, so panel membership does not inflate the comparator;\n")
  cat("  if anything it makes the pQTL estimate conservative.\n")
} else {
  cat("  WARNING: the key test no longer excludes 1. The response wording in R2.2\n")
  cat("  ('does not inflate our comparator') must be revisited.\n")
}

# --- held in reserve: per-phase decomposition -------------------------------
cat("\n--- per-phase transitions (reserve material, not used in the response) ---\n")
# succ_* columns are logical: TRUE = advanced past that transition. The rate is
# therefore mean(col) among pairs with a non-missing value for it.
phase_cols <- c("Preclinical->I" = "succ_p_1", "I->II" = "succ_1_2",
                "II->III" = "succ_2_3", "III->Launch" = "succ_3_a")
for (nm in names(phase_cols)) {
  cc <- phase_cols[[nm]]
  d  <- ti_best %>% filter(!is.na(.data[[cc]]))
  a  <- mean(d[[cc]][d$gensup],  na.rm = TRUE)
  b  <- mean(d[[cc]][!d$gensup], na.rm = TRUE)
  cat(sprintf("  %-15s supported %.1f%% (n=%d) vs unsupported %.1f%% (n=%d)  ratio %.2f\n",
              nm, 100*a, sum(d$gensup), 100*b, sum(!d$gensup), a/b))
}

write_tsv(res, file.path(out_dir, "r2_2_ascertainment.tsv"))
cat(sprintf("\n-> output/r2_2_ascertainment.tsv\n\n"))
