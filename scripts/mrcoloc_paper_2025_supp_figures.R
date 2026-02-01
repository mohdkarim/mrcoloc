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

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(data.table)
  library(DescTools)
  library(googlesheets4)
  library(ckbplotr)
  library(grid)
  library(fuzzyjoin)
  library(scales)
})

gs4_deauth()

# --- Configuration -----------------------------------------------------------

project_root <- Sys.getenv(
  "PQTL_ENRICH_ROOT",
  "/home/mohd/mohd-sandbox/pQTL_enrichment/mrcoloc_paper2025"
)

fig_dir       <- file.path(project_root, "figures")
data_dir      <- file.path(project_root, "genetic_support-main", "data")
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

# --- 1. Load Data (if not already in environment) ----------------------------

message("[1/5] Loading data...")

# Check if merge3_pqtl exists, if not load it
if (!exists("merge3_pqtl")) {
  
  # pQTL MR-coloc results
  pqtl2 <- readRDS(file.path(project_root, "data_raw/pqtl_mrcoloc_2025.rds")) %>%
    filter(bxy_pval <= mr_pval_threshold)
  
  # Merge2 with annotations
  merge2 <- read_tsv(file.path(data_dir, "merge2.tsv.gz"),
                     show_col_types = FALSE) %>%
    mutate(
      otg_study = if_else(
        assoc_source == "OTG",
        str_remove(original_link, "https://genetics.opentargets.org/study/"),
        NA_character_
      ),
      otg_study = str_remove(otg_study, "FINNGEN_R6_"),
      key       = paste0(gene, "_", otg_study)
    )
  
  pqtl_cols <- c("nsnps", "cis_trans_mr", "bxy", "bxy_pval",
                 "coloc_cis", "coloc_h4_cis", "snp_ciscoloc",
                 "coloc_trans", "coloc_h4_trans", "snp_transcoloc")
  
  merge2_with_pqtl <- merge2 %>%
    left_join(pqtl2, by = "key")
  
  pqtl_rows <- merge2_with_pqtl %>%
    filter(!is.na(cis_trans_mr)) %>%
    mutate(original_link = "pqtl")
  
  merge2_cleaned <- merge2_with_pqtl %>%
    mutate(across(all_of(pqtl_cols),
                  ~ if_else(!is.na(cis_trans_mr), NA, .)))
  
  merge3_pqtl <- bind_rows(merge2_cleaned, pqtl_rows)
  
  # Platform and therapeutic area annotations
  merge3_pqtl <- merge3_pqtl %>%
    mutate(
      platform = case_when(
        Data %in% c("UKBPPP_2023", "SCALLOP_2020",
                    "HILLARY_2019", "FOLKERSEN_2017") ~ "Olink",
        Data %in% c("SUN_2018", "SUHRE_2017", "PIETZNER_2020") ~ "Somascan",
        Data == "OLLI_2017" ~ "Other",
        TRUE ~ NA_character_
      )
    )
  
  area <- fread(file.path(data_dir, "areas.tsv"))
  topl <- fread(file.path(data_dir, "indic_topl_match.tsv"))
  ta   <- merge(area, topl, by = "topl")
  pos  <- match(merge3_pqtl$indication_mesh_id, ta$indication_mesh_id)
  merge3_pqtl$therapeutic_area <- ta$area[pos]
}

# Load indications
indic <- read_tsv(file.path(data_dir, "indic.tsv"), show_col_types = FALSE)

# Load proteomics background genes if not present
if (!exists("pgenes")) {
  olink <- read_sheet(
    "https://docs.google.com/spreadsheets/d/1DBHpr_Y3pFja4tMju3ZDJV8Gv-oTLq6wEuS0HtYjGbQ",
    sheet = "olink_complete"
  )
  
  olink2 <- olink %>%
    separate_rows(`Uniprot ID`, sep = ",") %>%
    mutate(
      hgnc_protein = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys      = `Uniprot ID`,
        column    = "SYMBOL",
        keytype   = "UNIPROT",
        multiVals = "first"
      )
    ) %>%
    filter(!is.na(hgnc_protein))
  
  olink_genes <- unique(as.character(olink2$hgnc_protein))
  otherttpairs <- readRDS(file.path(project_root, "data_raw/ttpairs_tested.rds"))
  othergenes   <- unique(gsub("_.*", "", otherttpairs))
  pgenes       <- unique(c(olink_genes, othergenes))
}

message("   Data loaded successfully")
message("   Background gene set (pgenes): ", length(pgenes), " genes")


# --- 2. Define pQTL-supported TIs (restricted to measured proteins) ----------

message("[2/5] Defining pQTL-supported T-I pairs...")

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

message("[3/5] Generating Figure S2 (enrichment by therapeutic area)...")

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

message("[4/5] Generating Figure S3 (enrichment by phase transition)...")

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