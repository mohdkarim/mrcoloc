#!/usr/bin/env Rscript
# ============================================================================
# Supplementary Figures: Proteogenomic Evidence and Drug Development Success
# ============================================================================
#
# This script generates supplementary figures for the 2025 proteogenomics paper:
#   - Figure S2: Enrichment by therapeutic area
#   - Figure S3: Enrichment by clinical phase transition
#
# Usage:
#   Rscript mrcoloc_paper_2025_supp_figures.R
#
# Prerequisites:
#   Run mrcoloc_paper_2025_figure1.R first to generate required objects,
#   or this script will load/compute them independently.
#
# ============================================================================

# --- 0. Setup and Configuration ---------------------------------------------

message("
================================================================================
  Supplementary Figures: Proteogenomic Evidence and Drug Development Success
================================================================================
")

message("[0/5] Loading packages...")

# Auto-install missing packages
cran_pkgs <- c("tidyverse", "data.table", "DescTools", "ckbplotr", "fuzzyjoin", "scales",
               "binom", "epitools", "janitor", "glue", "lawstat", "weights", "openxlsx", "optparse", "MASS")
missing <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if (length(missing) > 0) {
  message("   Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DescTools)
  library(ckbplotr)
  library(grid)
  library(fuzzyjoin)
  library(scales)
})

# --- Configuration -----------------------------------------------------------

project_root <- Sys.getenv(
  "MRCOLOC_ROOT",
  getwd()
)

fig_dir       <- file.path(project_root, "figures")
data_dir      <- file.path(project_root, "data")
minikel_dir   <- file.path(project_root, "data", "minikel")
r_dir         <- project_root

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

# Analysis thresholds
mr_pval_threshold    <- 0.05 / 47e6
similarity_threshold <- 0.8
min_l2g_share        <- 0.5

message("   Project root: ", project_root)
message("   Output dir:   ", fig_dir)

t0 <- Sys.time()
t_step <- Sys.time()

# --- 1. Load Data (if not already in environment) ----------------------------

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[1/5] Loading data...")
t_step <- Sys.time()

# Load pre-computed data (created by create_derived_data.R)
merge3_pqtl <- readRDS(file.path(project_root, "data_raw/merge3_pqtl.rds"))

# Also load merge2 for full Minikel universe comparison (Figure S4)
merge2 <- read_tsv(file.path(minikel_dir, "merge2.tsv.gz"), show_col_types = FALSE) %>%
  mutate(
    otg_study = if_else(assoc_source == "OTG",
      str_remove(original_link, "https://genetics.opentargets.org/study/"), NA_character_),
    otg_study = str_remove(otg_study, "FINNGEN_R6_"),
    key = paste0(gene, "_", otg_study)
  )

# Indications
indic <- read_tsv(file.path(data_dir, "indic.tsv"), show_col_types = FALSE)

# Background gene set
pgenes <- readRDS(file.path(project_root, "data_raw/pgenes.rds"))

message("   Data loaded successfully")
message("   Background gene set (pgenes): ", length(pgenes), " genes")


# --- 2. Define pQTL-supported TIs (restricted to measured proteins) ----------

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[2/5] Defining pQTL-supported T-I pairs...")
t_step <- Sys.time()

df_pqtl_support <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= similarity_threshold,
    !is.na(l2g_share), l2g_share >= min_l2g_share,
    gene %in% pgenes
  ) %>%
  distinct(ti_uid)

message("   pQTL-supported T-I pairs: ", nrow(df_pqtl_support))


# --- 3. Figure S2: Enrichment by Therapeutic Area ----------------------------

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[3/5] Generating Figure S2 (enrichment by therapeutic area)...")
t_step <- Sys.time()

# Get all therapeutic areas (restricted to measured proteins)
tas <- merge3_pqtl %>%
  filter(!is.na(therapeutic_area), gene %in% pgenes) %>%
  distinct(therapeutic_area) %>%
  pull(therapeutic_area)

message("   Therapeutic areas found: ", length(tas))

# Compute enrichment per therapeutic area
enrichment_by_ta <- map_dfr(tas, function(ta) {
  
  # Supported TIs in this TA (restricted to measured proteins)
  df_support <- merge3_pqtl %>%
    filter(
      grepl("pqtl", original_link, ignore.case = TRUE),
      comb_norm >= similarity_threshold,
      !is.na(l2g_share), l2g_share >= min_l2g_share,
      therapeutic_area == ta,
      gene %in% pgenes
    ) %>%
    distinct(ti_uid)
  
  # Best TI row per indication (restricted to measured proteins)
  ti_best <- merge3_pqtl %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(indication_mesh_id), indication_mesh_id != "",
      !is.na(therapeutic_area), therapeutic_area == ta,
      gene %in% pgenes
    ) %>%
    left_join(indic %>% select(indication_mesh_id, genetic_insight),
              by = "indication_mesh_id") %>%
    filter(genetic_insight != "none") %>%
    arrange(ti_uid, desc(comb_norm)) %>%
    group_by(ti_uid) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gensup = ti_uid %in% df_support$ti_uid)
  
  # Extract success info
  long <- ti_best %>%
    filter(!is.na(succ_3_a)) %>%
    rename(success = succ_3_a)
  
  baseline <- ti_best %>%
    filter(!is.na(succ_1_2))
  
  # Counts
  succ_gs <- sum(long$gensup & long$success, na.rm = TRUE)
  succ_nogs <- sum(!long$gensup & long$success, na.rm = TRUE)
  total_gs <- sum(baseline$gensup, na.rm = TRUE)
  total_nogs <- sum(!baseline$gensup, na.rm = TRUE)
  
  # Guard against bad input
  if (any(c(total_gs, total_nogs) == 0) || succ_gs > total_gs || succ_nogs > total_nogs) {
    tibble(
      source = ta,
      est = NA_real_,
      lwr.ci = NA_real_,
      upr.ci = NA_real_,
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    )
  } else {
    out <- as.data.frame(BinomRatioCI(
      x1 = succ_gs, n1 = total_gs,
      x2 = succ_nogs, n2 = total_nogs,
      method = "katz"
    ))
    out$source <- ta
    out$n <- paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    as_tibble(out)
  }
})

# Prepare for forest plot
# Exclude TAs with 0 supported successes (RS = 0) or undefined RS
enrichment_by_ta <- enrichment_by_ta %>%
  filter(!is.na(est), est > 0) %>%
  arrange(desc(est)) %>%
  mutate(
    subgroup = row_number(),
    label = "By therapeutic area"
  )

row_labels_ta <- enrichment_by_ta %>%
  transmute(
    subgroup = subgroup,
    label = str_to_title(source)
  )

figS2 <- forest_plot(
  enrichment_by_ta,
  col.key           = "subgroup",
  col.lci           = "lwr.ci",
  col.uci           = "upr.ci",
  exponentiate      = FALSE,
  nullval           = 1,
  col.left          = "n",
  stroke            = 0.8,
  shape             = 16,
  pointsize         = 3,
  estcolumn         = TRUE,
  base_size         = 11,
  xlim              = c(0.5, 25),
  xticks            = c(1, 5, 10, 15, 20, 25),
  col.right.heading = "RS (95% CI)",
  col.left.heading  = "Supported/\nUnsupported",
  xlab              = "Phase I-Launch Relative Success (RS)",
  panel.headings    = NULL,
  row.labels        = row_labels_ta
)

ggsave(file.path(fig_dir, "figS2_enrichment_by_ta.pdf"),
       figS2$plot, width = 9, height = 6, device = "pdf")

message("   Saved: figS2_enrichment_by_ta.pdf")
message("   Note: Excluded 6 TAs with no pQTL-supported T-I pairs (endocrine, immune,")
message("         infection, ophthalmology, other, psychiatry) and 3 TAs with no")
message("         supported successes (neurology, oncology, signs/symptoms)")


# --- 4. Figure S3: Enrichment by Phase Transition ----------------------------

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[4/5] Generating Figure S3 (enrichment by phase transition)...")
t_step <- Sys.time()

# Phase mapping - include cumulative transitions
phase_map <- tibble(
  phase   = c("Preclinical > I", "I > II", "II > III", "III > Launch", "Preclinical > Launch", "I > Launch"),
  varname = c("succ_p_1", "succ_1_2", "succ_2_3", "succ_3_a", "succ_p_launch", "succ_1_launch")
)

# Build ti_best_all (restricted to measured proteins)
ti_best_all <- merge3_pqtl %>%
  filter(
    !is.na(gene), gene != "",
    !is.na(indication_mesh_id), indication_mesh_id != "",
    !is.na(ccat),
    gene %in% pgenes
  ) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3,
    !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1,
    !is.na(succ_p_1) ~ 0,
    TRUE ~ NA_real_
  )) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(gensup = ti_uid %in% df_pqtl_support$ti_uid)

message("   T-I pairs in background: ", nrow(ti_best_all))

# Compute enrichment by phase
enrichment_by_phase <- map_dfr(1:nrow(phase_map), function(i) {
  this_phase <- phase_map[i, ]
  
  if (this_phase$varname == "succ_p_launch") {
    # Special case: Preclinical to Launch
    long <- ti_best_all %>%
      filter(!is.na(succ_p_1)) %>%
      mutate(success = succ_3_a == TRUE) %>%
      select(ti_uid, gensup, success)
  } else if (this_phase$varname == "succ_1_launch") {
    # Special case: Phase I to Launch
    long <- ti_best_all %>%
      filter(!is.na(succ_1_2)) %>%
      mutate(success = succ_3_a == TRUE) %>%
      select(ti_uid, gensup, success)
  } else {
    # Standard phase transitions
    long <- ti_best_all %>%
      select(ti_uid, gensup, all_of(this_phase$varname)) %>%
      rename(success = all_of(this_phase$varname)) %>%
      filter(!is.na(success))
  }
  
  # Calculate totals
  total_gs <- sum(long$gensup, na.rm = TRUE)
  total_nogs <- sum(!long$gensup, na.rm = TRUE)
  succ_gs <- sum(long$gensup & long$success, na.rm = TRUE)
  succ_nogs <- sum(!long$gensup & long$success, na.rm = TRUE)
  
  if (any(c(total_gs, total_nogs) == 0) || succ_gs > total_gs || succ_nogs > total_nogs) {
    tibble(
      source = this_phase$phase,
      est = NA_real_,
      lwr.ci = NA_real_,
      upr.ci = NA_real_,
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    )
  } else {
    out <- as.data.frame(BinomRatioCI(
      x1 = succ_gs, n1 = total_gs,
      x2 = succ_nogs, n2 = total_nogs,
      method = "katz"
    ))
    out$source <- this_phase$phase
    out$n <- paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    as_tibble(out)
  }
})

# Prepare for forest plot
enrichment_by_phase <- enrichment_by_phase %>%
  filter(!is.na(est)) %>%
  mutate(
    source = factor(source, levels = phase_map$phase),
    subgroup = as.numeric(source)
  )

row_labels_phase <- enrichment_by_phase %>%
  transmute(
    subgroup = as.character(subgroup),
    label = as.character(source)
  )

figS3 <- forest_plot(
  enrichment_by_phase,
  col.key           = "subgroup",
  col.lci           = "lwr.ci",
  col.uci           = "upr.ci",
  exponentiate      = FALSE,
  nullval           = 1,
  col.left          = "n",
  stroke            = 0.8,
  shape             = 16,
  pointsize         = 3,
  estcolumn         = TRUE,
  base_size         = 11,
  xlim              = c(0.5, 7),
  xticks            = c(1, 2, 3, 4, 5, 6, 7),
  col.right.heading = "RS (95% CI)",
  col.left.heading  = "Supported/\nUnsupported",
  xlab              = "Relative Success (RS)",
  panel.headings    = NULL,
  row.labels        = row_labels_phase
)

ggsave(file.path(fig_dir, "figS3_enrichment_by_phase.pdf"),
       figS3$plot, width = 9, height = 5, device = "pdf")

message("   Saved: figS3_enrichment_by_phase.pdf")


# --- 5. Summary --------------------------------------------------------------

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[5/5] Complete!")
message("
================================================================================
  Supplementary Figures generation complete
================================================================================

Output files in: ", fig_dir, "

Figures:
  - figS2_enrichment_by_ta.pdf       Forest plot by therapeutic area (8 TAs)
  - figS3_enrichment_by_phase.pdf    Forest plot by phase transition (6 phases)

Notes:
  - Figure S2 excludes 9 TAs: 6 with no pQTL-supported T-I pairs,
    3 with no supported successes (RS=0)
  - Background restricted to ", length(pgenes), " measured proteins
  - pQTL-supported T-I pairs: ", nrow(df_pqtl_support), "
")
# ==============================================================================
# Figure S4: pQTL Enrichment with Full Minikel Background
# ==============================================================================
# Sensitivity analysis comparing pQTL enrichment using:
# 1. Measured proteins background (main analysis)
# 2. Full Minikel T-I pairs background (no pgenes restriction)
# ==============================================================================

if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")
message("[S4] Generating Figure S4 (pQTL enrichment - background comparison)...")
t_step <- Sys.time()

# pQTL-supported TI pairs (same as main analysis)
df_pqtl_support <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= 0.8,
    !is.na(l2g_share), l2g_share >= 0.5
  ) %>%
  distinct(ti_uid)

# --- Analysis 1: Measured proteins background (pgenes) ---
ti_best_pgenes <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat), gene %in% pgenes) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_
  )) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(gensup = ti_uid %in% df_pqtl_support$ti_uid)

long_pg <- ti_best_pgenes %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
base_pg <- ti_best_pgenes %>% filter(!is.na(succ_1_2))

rs_pgenes <- as.data.frame(BinomRatioCI(
  x1 = sum(long_pg$gensup & long_pg$success), n1 = sum(base_pg$gensup),
  x2 = sum(!long_pg$gensup & long_pg$success), n2 = sum(!base_pg$gensup),
  method = "katz"
)) %>% mutate(
  source = "pQTL (measured proteins)",
  n = paste0("(", sum(long_pg$gensup & long_pg$success), "/", sum(base_pg$gensup), ")/(",
             sum(!long_pg$gensup & long_pg$success), "/", sum(!base_pg$gensup), ")")
)

# --- Analysis 2: Full Minikel background (no pgenes restriction) ---
ti_best_full <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(highest_phase = case_when(
    !is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
    !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_
  )) %>%
  arrange(ti_uid, desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(gensup = ti_uid %in% df_pqtl_support$ti_uid)

long_full <- ti_best_full %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
base_full <- ti_best_full %>% filter(!is.na(succ_1_2))

rs_full <- as.data.frame(BinomRatioCI(
  x1 = sum(long_full$gensup & long_full$success), n1 = sum(base_full$gensup),
  x2 = sum(!long_full$gensup & long_full$success), n2 = sum(!base_full$gensup),
  method = "katz"
)) %>% mutate(
  source = "pQTL (all Minikel T-I pairs)",
  n = paste0("(", sum(long_full$gensup & long_full$success), "/", sum(base_full$gensup), ")/(",
             sum(!long_full$gensup & long_full$success), "/", sum(!base_full$gensup), ")")
)

# --- Combine and plot ---
figS4_data <- bind_rows(rs_pgenes, rs_full) %>%
  mutate(colour = "blue", subgroup = row_number())

row_labels_s4 <- data.frame(
  subgroup = 1:2,
  label = c("Measured proteins background", "Full Minikel background")
)

figS4 <- forest_plot(
  figS4_data,
  col.key = "subgroup", col.lci = "lwr.ci", col.uci = "upr.ci",
  exponentiate = FALSE, colour = "colour", nullval = 1,
  col.left = "n", stroke = 1, shape = 16, estcolumn = TRUE,
  base_size = 14, xlim = c(1, 9), xticks = 1:8,
  col.right.heading = "RS, 95% CI",
  col.left.heading = "(A[G]/S[G])/(A![G]/S![G])",
  xlab = "Phase I-Launch Relative Success (RS), 95% CI",
  panel.headings = "pQTL Enrichment by Background Universe",
  row.labels = row_labels_s4
)

ggsave(file.path(fig_dir, "figS4_background_comparison.pdf"), 
       plot = figS4$plot, device = "pdf", width = 8, height = 3)
ggsave(file.path(fig_dir, "figS4_background_comparison.png"), 
       plot = figS4$plot, device = "png", width = 8, height = 3, dpi = 300)

message("   Saved: figS4_background_comparison.pdf/.png")
message("   Measured proteins: RS = ", round(rs_pgenes$est, 2), 
        " (95% CI: ", round(rs_pgenes$lwr.ci, 2), "-", round(rs_pgenes$upr.ci, 2), ")")
message("   Full Minikel:      RS = ", round(rs_full$est, 2),
        " (95% CI: ", round(rs_full$lwr.ci, 2), "-", round(rs_full$upr.ci, 2), ")")
if (exists("t_step")) message("   done (", round(difftime(Sys.time(), t_step, units="secs")), "s)")

message("\nTotal time: ", round(difftime(Sys.time(), t0, units="mins"), 1), " minutes")
