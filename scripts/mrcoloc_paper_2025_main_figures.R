#!/usr/bin/env Rscript
# ============================================================================
# Figure 1: Proteogenomic Evidence and Drug Development Success
# ============================================================================
#
# This script generates Figure 1 for the 2025 proteogenomics paper:
#   - Figure 1a: Forest plot of relative success by genetic evidence source
#   - Figure 1b: UpSet plot of MR/coloc overlap
#   - Figure 1c: Forest plot of enrichment by gene family (L2G vs L2G+pQTL)
#
# Usage:
#   Rscript mrcoloc_paper_2025_figure1.R
#
# Environment:
#   Set PQTL_ENRICH_ROOT to your project directory, or it defaults to:
#   /home/mohd/mohd-sandbox/pQTL_enrichment/mrcoloc_paper2025
#
# Outputs (in figures/):
#   - fig1a_forest_main.pdf/.png
#   - fig1b_upset.pdf/.png
#   - fig1c_forest_gene_annot.pdf/.png
#
# ============================================================================

# --- 0. Setup and Configuration ---------------------------------------------

message("
================================================================================
  Figure 1: Proteogenomic Evidence and Drug Development Success
================================================================================
")

message("[0/7] Loading packages...")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(data.table)
  library(DescTools)
  library(UpSetR)
  library(googlesheets4)
  library(ckbplotr)
  library(grid)
  library(magick)
})

gs4_deauth()

# --- Configuration -----------------------------------------------------------

project_root <- Sys.getenv(
  "PQTL_ENRICH_ROOT",
  "/home/mohd/mohd-sandbox/pQTL_enrichment/mrcoloc_paper2025"
)

# Directory structure
fig_dir       <- file.path(project_root, "figures")
data_dir      <- file.path(project_root, "genetic_support-main", "data")
gene_list_dir <- file.path(data_dir, "gene_lists")
r_dir         <- project_root

# Create output directory if needed
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

# Analysis thresholds
mr_pval_threshold    <- 0.05 / 47e6
similarity_threshold <- 0.8
min_l2g_share        <- 0.5
coloc_pp4_threshold  <- 0.8

l2g_thresholds <- c(0.25, 0.50, 0.75)

# Phase mapping
phase_map <- tibble(
  phase   = factor(c("Preclinical", "I", "II", "III"),
                   ordered = TRUE,
                   levels  = c("Preclinical", "I", "II", "III")),
  phorder = 0:3,
  varname = c("succ_p_1", "succ_1_2", "succ_2_3", "succ_3_a")
)

# Genetic evidence sources
sources <- c(
  "OMIM (Mendelian disease)", "All OTG", "PICCOLO",
  "Genebass (gene burden)", "FinnGen", "Neale UKBB", "GWAS catalog"
)

message("   Project root: ", project_root)
message("   Output dir:   ", fig_dir)

# --- 1. Helper Functions -----------------------------------------------------

message("[1/7] Defining helper functions...")

#' Compute relative risk with confidence interval
compute_rr <- function(succ_gs, total_gs, succ_nogs, total_nogs,
                       source, label = NA_character_) {
  if (any(c(total_gs, total_nogs) == 0) ||
      succ_gs > total_gs || succ_nogs > total_nogs) {
    return(tibble(
      est    = NA_real_,
      lwr.ci = NA_real_,
      upr.ci = NA_real_,
      source = source,
      n      = paste0("(", succ_gs, "/", total_gs, ")/(", 
                      succ_nogs, "/", total_nogs, ")"),
      label  = label
    ))
  }
  
  out <- as.data.frame(BinomRatioCI(
    x1 = succ_gs,  n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  ))
  out$source <- source
  out$n      <- paste0("(", succ_gs, "/", total_gs, ")/(",
                       succ_nogs, "/", total_nogs, ")")
  out$label  <- label
  as_tibble(out)
}

#' Select best target-indication row per pair
select_best_ti_row <- function(df, indic_tbl) {
  df %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat)) %>%
    left_join(indic_tbl %>% select(indication_mesh_id, genetic_insight),
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
    ungroup()
}

#' Save plot as PDF and PNG
save_figure <- function(plot_obj = NULL, filename, width, height,
                        is_base_graphics = FALSE, base_plot_fn = NULL) {
  pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename, ".png"))
  
  if (is_base_graphics && !is.null(base_plot_fn)) {
    pdf(pdf_path, width = width, height = height, onefile = FALSE)
    base_plot_fn()
    dev.off()
  } else if (!is.null(plot_obj)) {
    ggsave(pdf_path, plot = plot_obj, device = "pdf",
           width = width, height = height, units = "in")
  }
  
  # Convert PDF to PNG using magick
  img <- image_read_pdf(pdf_path, density = 300)
  image_write(img, path = png_path, format = "png")
  
  message("   Saved: ", basename(pdf_path))
  message("   Saved: ", basename(png_path))
  
  invisible(pdf_path)
}

# --- 2. Load and Prepare Data ------------------------------------------------

message("[2/7] Loading data...")

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

# Merge pQTL data
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

# Platform annotations
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

# Therapeutic area
area <- fread(file.path(data_dir, "areas.tsv"))
topl <- fread(file.path(data_dir, "indic_topl_match.tsv"))
ta   <- merge(area, topl, by = "topl")
pos  <- match(merge3_pqtl$indication_mesh_id, ta$indication_mesh_id)
merge3_pqtl$therapeutic_area <- ta$area[pos]

# Indications
indic <- read_tsv(file.path(data_dir, "indic.tsv"), show_col_types = FALSE)

# Olink panel genes
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

# Other tested genes
otherttpairs <- readRDS(file.path(project_root, "data_raw/ttpairs_tested.rds"))
othergenes   <- unique(gsub("_.*", "", otherttpairs))
pgenes       <- unique(c(olink_genes, othergenes))

message("   Loaded ", nrow(merge3_pqtl), " T-I associations")
message("   ", length(pgenes), " genes in proteomics background")

# --- 3. Compute Enrichment Metrics -------------------------------------------

message("[3/7] Computing enrichment metrics...")

# Define pQTL-supported TI pairs
df_pqtl_support <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= similarity_threshold,
    !is.na(l2g_share), l2g_share >= min_l2g_share
  ) %>%
  distinct(ti_uid)

# Best TI per pair
ti_best_all <- merge3_pqtl %>%
  as.data.frame() %>%
  select_best_ti_row(indic_tbl = indic) %>%
  mutate(gensup = ti_uid %in% df_pqtl_support$ti_uid)

# --- 3a. pQTL enrichment within Olink/Soma genes ---

get_pgene_enrichment <- function(ti_best, pgenes) {
  ti_pgene <- ti_best %>% filter(gene %in% pgenes)
  
  long <- ti_pgene %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_pgene %>% filter(!is.na(succ_1_2))
  
  compute_rr(
    succ_gs    = sum(long$gensup & long$success, na.rm = TRUE),
    total_gs   = sum(baseline$gensup, na.rm = TRUE),
    succ_nogs  = sum(!long$gensup & long$success, na.rm = TRUE),
    total_nogs = sum(!baseline$gensup, na.rm = TRUE),
    source     = "pQTL",
    label      = "By genetic evidence source"
  )
}

pqtl_enrichment <- get_pgene_enrichment(ti_best_all, pgenes)

# --- 3b. Enrichment by genetic evidence source (Minikel-style) ---

get_enrichment_minikel_style <- function(src) {
  supported_tis <- merge3_pqtl %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat)) %>%
    left_join(indic %>% select(indication_mesh_id, genetic_insight),
              by = "indication_mesh_id") %>%
    filter(genetic_insight != "none") %>%
    filter(case_when(
      src == "OMIM (Mendelian disease)" ~ 
        assoc_source == "OMIM" & comb_norm >= similarity_threshold,
      src == "All OTG" ~ 
        assoc_source == "OTG" & comb_norm >= similarity_threshold & 
        !is.na(l2g_share) & l2g_share >= min_l2g_share,
      src == "PICCOLO" ~ 
        assoc_source == "PICCOLO" & comb_norm >= similarity_threshold & 
        !is.na(pic_h4) & pic_h4 >= 0.9,
      src == "Genebass (gene burden)" ~ 
        assoc_source == "Genebass" & comb_norm >= similarity_threshold,
      src == "FinnGen" ~ 
        assoc_source == "OTG" & grepl("FINNGEN", original_link, ignore.case = TRUE) &
        comb_norm >= similarity_threshold & !is.na(l2g_share) & l2g_share >= min_l2g_share,
      src == "Neale UKBB" ~ 
        assoc_source == "OTG" & grepl("NEALE", original_link, ignore.case = TRUE) &
        comb_norm >= similarity_threshold & !is.na(l2g_share) & l2g_share >= min_l2g_share,
      src == "GWAS catalog" ~ 
        assoc_source == "OTG" & grepl("GCST", original_link, ignore.case = TRUE) &
        comb_norm >= similarity_threshold & !is.na(l2g_share) & l2g_share >= min_l2g_share,
      TRUE ~ FALSE
    )) %>%
    distinct(ti_uid) %>%
    mutate(gensup = TRUE)
  
  ti_best <- merge3_pqtl %>%
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
    mutate(gensup = ti_uid %in% supported_tis$ti_uid)
  
  long     <- ti_best %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_best %>% filter(!is.na(succ_1_2))
  
  succ_gs    <- sum(long$gensup & long$success)
  succ_nogs  <- sum(!long$gensup & long$success)
  total_gs   <- sum(baseline$gensup)
  total_nogs <- sum(!baseline$gensup)
  
  as.data.frame(BinomRatioCI(
    x1 = succ_gs, n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  )) %>%
    mutate(
      source = src,
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    )
}

enrichment_overall <- map_dfr(sources, get_enrichment_minikel_style)

# --- 3c. L2G-only enrichment ---

source(file.path(r_dir, "R/pipeline_best.R"))
source(file.path(r_dir, "R/advancement_rr.R"))

enrichment_l2g_minikel <- map_dfr(l2g_thresholds, function(thr) {
  pb <- pipeline_best(merge2, phase = "combined", basis = "ti",
                      associations = c("OTG"), share_mode = "L2G",
                      min_share = thr, verbose = FALSE)
  rr  <- advancement_rr(pb)
  row <- rr %>% filter(phase == "I-Launch")
  
  tibble(
    est    = row$rs_mean,
    lwr.ci = row$rs_l,
    upr.ci = row$rs_u,
    source = paste0("L2G share: >= ", thr),
    n      = paste0("(", row$x_yes, "/", row$n_yes, ")/(",
                    row$x_no, "/", row$n_no, ")"),
    label  = "By L2G share only"
  )
})

# --- 3d. pQTL + L2G enrichment ---

enrichment_pqtl_l2g <- map_dfr(l2g_thresholds, function(threshold) {
  pqtl <- merge3_pqtl %>% filter(gene %in% pgenes)
  
  df <- pqtl %>%
    left_join(indic %>% select(indication_mesh_id, genetic_insight),
              by = "indication_mesh_id") %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat), genetic_insight != "none") %>%
    mutate(gs_dynamic = if_else(
      grepl("pqtl", original_link, ignore.case = TRUE) &
        comb_norm >= similarity_threshold &
        !is.na(l2g_share) & l2g_share >= threshold,
      "yes", "no"
    ))
  
  ti_best <- df %>%
    group_by(ti_uid) %>%
    arrange(desc(gs_dynamic == "yes"), desc(comb_norm)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gensup = gs_dynamic == "yes")
  
  long   <- ti_best %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  phase1 <- ti_best %>% filter(!is.na(succ_1_2))
  
  compute_rr(
    succ_gs    = sum(long$gensup & long$success, na.rm = TRUE),
    total_gs   = sum(phase1$gensup, na.rm = TRUE),
    succ_nogs  = sum(!long$gensup & long$success, na.rm = TRUE),
    total_nogs = sum(!phase1$gensup, na.rm = TRUE),
    source     = paste0("L2G share: >= ", threshold),
    label      = "By pQTL + L2G share"
  )
})

# --- 3e. Enrichment by MR type ---

enrichment_mrtype_df <- map_dfr(c("Cis", "Trans", "Mixed"), function(mrtype) {
  df_mrtype_support <- merge3_pqtl %>%
    filter(grepl("pqtl", original_link, ignore.case = TRUE),
           comb_norm >= similarity_threshold,
           !is.na(l2g_share), l2g_share >= min_l2g_share,
           cis_trans_mr == mrtype, gene %in% pgenes) %>%
    distinct(ti_uid)
  
  ti_best_mrtype <- merge3_pqtl %>%
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
    mutate(gensup = ti_uid %in% df_mrtype_support$ti_uid)
  
  long     <- ti_best_mrtype %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_best_mrtype %>% filter(!is.na(succ_1_2))
  
  succ_gs    <- sum(long$gensup & long$success)
  succ_nogs  <- sum(!long$gensup & long$success)
  total_gs   <- sum(baseline$gensup)
  total_nogs <- sum(!baseline$gensup)
  
  as.data.frame(BinomRatioCI(
    x1 = succ_gs, n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  )) %>%
    mutate(
      source = paste0(mrtype, "-MR"),
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    )
})

# --- 3f. MR + coloc enrichment ---

enrichment_mrcoloc_df <- {
  df_coloc_support <- merge3_pqtl %>%
    filter(grepl("pqtl", original_link, ignore.case = TRUE),
           comb_norm >= similarity_threshold,
           !is.na(l2g_share), l2g_share >= min_l2g_share,
           (coloc_h4_cis >= 0.8 | coloc_h4_trans >= 0.8),
           gene %in% pgenes) %>%
    distinct(ti_uid)
  
  ti_best_coloc <- merge3_pqtl %>%
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
    mutate(gensup = ti_uid %in% df_coloc_support$ti_uid)
  
  long     <- ti_best_coloc %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_best_coloc %>% filter(!is.na(succ_1_2))
  
  succ_gs    <- sum(long$gensup & long$success)
  succ_nogs  <- sum(!long$gensup & long$success)
  total_gs   <- sum(baseline$gensup)
  total_nogs <- sum(!baseline$gensup)
  
  as.data.frame(BinomRatioCI(
    x1 = succ_gs, n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  )) %>%
    mutate(
      source = "MR+coloc",
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")")
    )
}

enrichment_mrcoloc_combined_df <- bind_rows(enrichment_mrtype_df, enrichment_mrcoloc_df) %>%
  mutate(label = "By pQTL MR-coloc type")

# --- 3g. Enrichment by proteomics platform ---

enrichment_by_platform <- map_dfr(c("Somascan", "Olink"), function(p) {
  df_platform_support <- merge3_pqtl %>%
    filter(grepl("pqtl", original_link, ignore.case = TRUE),
           comb_norm >= similarity_threshold,
           !is.na(l2g_share), l2g_share >= min_l2g_share,
           platform == p) %>%
    distinct(ti_uid)
  
  ti_best_platform <- merge3_pqtl %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat), platform == p) %>%
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
    mutate(gensup = ti_uid %in% df_platform_support$ti_uid)
  
  long     <- ti_best_platform %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
  baseline <- ti_best_platform %>% filter(!is.na(succ_1_2))
  
  succ_gs    <- sum(long$gensup & long$success, na.rm = TRUE)
  succ_nogs  <- sum(!long$gensup & long$success, na.rm = TRUE)
  total_gs   <- sum(baseline$gensup, na.rm = TRUE)
  total_nogs <- sum(!baseline$gensup, na.rm = TRUE)
  
  if (any(c(total_gs, total_nogs) == 0) || succ_gs > total_gs || succ_nogs > total_nogs) {
    return(tibble(
      est = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_,
      source = paste0("pQTL: ", p),
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")"),
      label = "By protein platform"
    ))
  }
  
  as.data.frame(BinomRatioCI(
    x1 = succ_gs, n1 = total_gs,
    x2 = succ_nogs, n2 = total_nogs,
    method = "katz"
  )) %>%
    mutate(
      source = paste0("pQTL: ", p),
      n = paste0("(", succ_gs, "/", total_gs, ")/(", succ_nogs, "/", total_nogs, ")"),
      label = "By protein platform"
    )
})

message("   Computed enrichment for ", length(sources) + 1, " genetic evidence sources")
message("   Computed enrichment for ", length(l2g_thresholds), " L2G thresholds")

# --- 4. Figure 1a: Main Forest Plot ------------------------------------------

message("[4/7] Generating Figure 1a (main forest plot)...")

# Combine all enrichment results
en1 <- bind_rows(enrichment_overall, pqtl_enrichment) %>%
  mutate(label = "By genetic evidence source") %>%
  arrange(desc(est))

en2 <- bind_rows(enrichment_l2g_minikel, enrichment_pqtl_l2g)
en3 <- enrichment_mrcoloc_combined_df
en4 <- enrichment_by_platform

en <- bind_rows(en1, en2, en3, en4) %>%
  mutate(
    colour   = if_else(grepl("pqtl", source, ignore.case = TRUE) |
                         grepl("pqtl", label, ignore.case = TRUE), "blue", "black"),
    subgroup = row_number()
  )

row_labels <- en %>%
  select(subgroup, label, source) %>%
  rename(group = label, label = source)

fig1a <- forest_plot(
  en,
  col.key           = "subgroup",
  col.lci           = "lwr.ci",
  col.uci           = "upr.ci",
  exponentiate      = FALSE,
  colour            = "colour",
  nullval           = 1,
  col.left          = "n",
  stroke            = 1,
  shape             = 16,
  estcolumn         = TRUE,
  base_size         = 14,
  xlim              = c(1, 9),
  xticks            = 1:8,
  col.right.heading = "RS, 95% CI",
  col.left.heading  = "(A[G]/S[G])/\n(A![G]/S![G])",
  xlab              = "Phase I-Launch Relative Success (RS), 95% CI",
  panel.headings    = "Relative Success of pQTL- and GWAS-supported T-I Pairs",
  row.labels        = row_labels
)

save_figure(fig1a$plot, "fig1a_forest_main", width = 8, height = 8)

# --- 5. Figure 1b: UpSet Plot ------------------------------------------------

message("[5/7] Generating Figure 1b (UpSet plot)...")

pqtl_supported <- merge3_pqtl %>%
  filter(
    grepl("pqtl", original_link, ignore.case = TRUE),
    comb_norm >= similarity_threshold,
    !is.na(l2g_share), l2g_share >= min_l2g_share,
    succ_3_a == TRUE
  ) %>%
  distinct(
    ti_uid, gene, indication_mesh_term, assoc_mesh_term,
    comb_norm, l2g_share, bxy, bxy_pval, nsnps,
    cis_trans_mr, coloc_h4_cis, snp_ciscoloc, coloc_h4_trans, snp_transcoloc,
    succ_p_1, succ_1_2, succ_2_3, succ_3_a
  ) %>%
  mutate(
    cis_mr      = cis_trans_mr == "Cis",
    trans_mr    = cis_trans_mr == "Trans",
    mixed_mr    = cis_trans_mr == "Mixed",
    coloc_cis   = !is.na(coloc_h4_cis) & coloc_h4_cis >= 0.8,
    coloc_trans = !is.na(coloc_h4_trans) & coloc_h4_trans >= 0.8
  )

sets_success <- list(
  cis_mr      = pqtl_supported %>% filter(cis_mr) %>% pull(ti_uid) %>% unique(),
  trans_mr    = pqtl_supported %>% filter(trans_mr) %>% pull(ti_uid) %>% unique(),
  mixed_mr    = pqtl_supported %>% filter(mixed_mr) %>% pull(ti_uid) %>% unique(),
  coloc_cis   = pqtl_supported %>% filter(coloc_cis) %>% pull(ti_uid) %>% unique(),
  coloc_trans = pqtl_supported %>% filter(coloc_trans) %>% pull(ti_uid) %>% unique()
)

sets_success_df <- fromList(sets_success)

# Save PDF
pdf(file.path(fig_dir, "fig1b_upset.pdf"), width = 10, height = 9)
upset(sets_success_df, order.by = "freq", text.scale = 2)
grid.text(
  "Overlapping Genetic Support from MR and \nColoc Methods for Phase 1 to Launch T-I Pairs",
  x = 0.7, y = 0.97, just = "top",
  gp = gpar(fontsize = 15, fontface = "bold")
)
dev.off()

# Save PNG
png(file.path(fig_dir, "fig1b_upset.png"), width = 10, height = 9, units = "in", res = 300)
upset(sets_success_df, order.by = "freq", text.scale = 2)
grid.text(
  "Overlapping Genetic Support from MR and \nColoc Methods for Phase 1 to Launch T-I Pairs",
  x = 0.7, y = 0.97, just = "top",
  gp = gpar(fontsize = 15, fontface = "bold")
)
dev.off()

message("   Saved: fig1b_upset.pdf")
message("   Saved: fig1b_upset.png")

# --- 6. Figure 1c: Gene Family Forest Plot -----------------------------------

message("[6/7] Generating Figure 1c (gene family forest plot)...")

# Load gene families
gene_list_files <- list.files(gene_list_dir, pattern = "\\.tsv$", full.names = TRUE)

if (length(gene_list_files) == 0) {
  stop("No .tsv files found in gene_list_dir: ", gene_list_dir,
       "\nPlease check that the path exists and contains gene family files.")
}

message("   Found ", length(gene_list_files), " gene family files")

gene_families <- map_dfr(gene_list_files, ~ {
  family_name <- tools::file_path_sans_ext(basename(.x))
  read_tsv(.x, col_names = "gene", show_col_types = FALSE) %>%
    mutate(family = family_name)
}) %>%
  distinct(gene, family)

# Best T-I with pQTL support flag
ti_best_pqtl <- merge3_pqtl %>%
  filter(!is.na(gene), gene != "",
         !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat)) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight),
            by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(
    highest_phase = case_when(
      succ_3_a ~ 3, succ_2_3 ~ 2, succ_1_2 ~ 1, succ_p_1 ~ 0, TRUE ~ NA_real_
    ),
    pqtl_support = grepl("pqtl", original_link, ignore.case = TRUE) &
      comb_norm >= 0.8 & !is.na(l2g_share) & l2g_share >= 0.5
  ) %>%
  arrange(ti_uid, desc(pqtl_support), desc(highest_phase), desc(comb_norm)) %>%
  group_by(ti_uid) %>%
  slice(1) %>%
  ungroup()

ti_best_annotated <- ti_best_pqtl %>%
  inner_join(gene_families, by = "gene", relationship = "many-to-many")

# L2G-based enrichment by family
ti_family_support <- ti_best_annotated %>%
  filter(!is.na(family)) %>%
  mutate(genetic_support = ifelse(!is.na(l2g_share) & l2g_share >= 0.5, 
                                  "supported", "unsupported")) %>%
  group_by(family, ti_uid, genetic_support) %>%
  summarise(
    phase1 = any(!is.na(succ_1_2)),
    launch = any(!is.na(succ_3_a) & succ_3_a),
    .groups = "drop"
  )

family_enrichment <- ti_family_support %>%
  group_by(family, genetic_support) %>%
  summarise(x_yes = sum(launch), n_yes = sum(phase1), .groups = "drop") %>%
  pivot_wider(names_from = genetic_support, values_from = c(x_yes, n_yes), values_fill = 0) %>%
  rowwise() %>%
  mutate(
    rs_tbl = data.frame(BinomRatioCI(
      x1 = x_yes_supported, n1 = n_yes_supported,
      x2 = x_yes_unsupported, n2 = n_yes_unsupported,
      method = "katz"
    )),
    est = rs_tbl$est, lwr.ci = rs_tbl$lwr.ci, upr.ci = rs_tbl$upr.ci
  ) %>%
  ungroup() %>%
  select(family, x_yes_supported, n_yes_supported, 
         x_yes_unsupported, n_yes_unsupported, est, lwr.ci, upr.ci)

# pQTL-based enrichment by family
ti_best_annotated_pqtl <- ti_best_annotated %>%
  mutate(pqtl_supported = grepl("pqtl", original_link, ignore.case = TRUE)) %>%
  filter(gene %in% pgenes)

family_ti_pqtl <- ti_best_annotated_pqtl %>%
  filter(!is.na(family)) %>%
  group_by(family, ti_uid, pqtl_supported) %>%
  summarise(
    phase1 = any(!is.na(succ_1_2)),
    launch = any(!is.na(succ_3_a) & succ_3_a),
    .groups = "drop"
  )

family_summary <- family_ti_pqtl %>%
  group_by(family, pqtl_supported) %>%
  summarise(n_yes = sum(phase1), x_yes = sum(launch), .groups = "drop") %>%
  pivot_wider(names_from = pqtl_supported, values_from = c(n_yes, x_yes),
              names_prefix = "pqtl_", values_fill = 0)

enrichment_pqtl_by_family <- family_summary %>%
  rowwise() %>%
  mutate(
    rs_tbl = data.frame(BinomRatioCI(
      x1 = x_yes_pqtl_TRUE, n1 = n_yes_pqtl_TRUE,
      x2 = x_yes_pqtl_FALSE, n2 = n_yes_pqtl_FALSE,
      method = "katz"
    )),
    est = rs_tbl$est, lwr.ci = rs_tbl$lwr.ci, upr.ci = rs_tbl$upr.ci
  ) %>%
  ungroup() %>%
  select(family, est, lwr.ci, upr.ci,
         x_yes_pqtl_TRUE, n_yes_pqtl_TRUE, x_yes_pqtl_FALSE, n_yes_pqtl_FALSE)

# Prepare for plotting
l2g_df <- family_enrichment %>%
  mutate(source = "L2G", family = fct_reorder(family, est))

pqtl_df <- enrichment_pqtl_by_family %>%
  rename(
    x_yes_supported = x_yes_pqtl_TRUE, n_yes_supported = n_yes_pqtl_TRUE,
    x_yes_unsupported = x_yes_pqtl_FALSE, n_yes_unsupported = n_yes_pqtl_FALSE
  ) %>%
  mutate(source = "L2G + pQTL", family = fct_reorder(as.factor(family), est))

common_families <- intersect(l2g_df$family, pqtl_df$family)
l2g_df  <- l2g_df %>% filter(family %in% common_families)
pqtl_df <- pqtl_df %>% filter(family %in% common_families)

order_levels <- l2g_df %>% arrange(est) %>% pull(family) %>% as.character()
l2g_df$family  <- factor(l2g_df$family, levels = rev(order_levels))
pqtl_df$family <- factor(pqtl_df$family, levels = rev(order_levels))

l2g_clean <- l2g_df %>%
  mutate(
    analysis = "L2G only",
    subgroup = family,
    n = paste0("(", x_yes_supported, "/", n_yes_supported, ") / (",
               x_yes_unsupported, "/", n_yes_unsupported, ")")
  )

pqtl_clean <- pqtl_df %>%
  mutate(
    analysis = "L2G + pQTL",
    subgroup = family,
    n = paste0("(", x_yes_supported, "/", n_yes_supported, ") / (",
               x_yes_unsupported, "/", n_yes_unsupported, ")")
  )

forest_df <- bind_rows(l2g_clean, pqtl_clean) %>%
  mutate(analysis = factor(analysis, levels = c("L2G only", "L2G + pQTL"))) %>%
  select(analysis, subgroup, est, lwr.ci, upr.ci, n)

row_labels_1c <- data.frame(
  subgroup = unique(forest_df$subgroup),
  label    = gsub("_", " ", unique(forest_df$subgroup))
)

forest_df_clean <- forest_df %>%
  mutate(
    est    = ifelse(analysis == "L2G + pQTL" & (is.na(est) | lwr.ci == 0), NA, est),
    lwr.ci = ifelse(analysis == "L2G + pQTL" & lwr.ci == 0, NA, lwr.ci),
    upr.ci = ifelse(analysis == "L2G + pQTL" & (is.na(upr.ci) | lwr.ci == 0), NA, upr.ci)
  )

ordered_subgroups <- forest_df_clean %>%
  filter(analysis == "L2G + pQTL") %>%
  mutate(est_order = ifelse(is.na(est), -Inf, est)) %>%
  arrange(desc(est_order)) %>%
  pull(subgroup)

forest_df_clean <- forest_df_clean %>%
  mutate(subgroup = factor(subgroup, levels = ordered_subgroups))

row_labels_1c <- row_labels_1c %>%
  mutate(subgroup = factor(subgroup, levels = ordered_subgroups)) %>%
  arrange(subgroup)

forest_df_clean$colour <- with(forest_df_clean, 
                               ifelse(grepl("pqtl", analysis, ignore.case = TRUE), 
                                      "blue", "black"))

fig1c <- forest_plot(
  split(forest_df_clean, ~analysis),
  col.key           = "subgroup",
  col.lci           = "lwr.ci",
  col.uci           = "upr.ci",
  colour            = "colour",
  exponentiate      = FALSE,
  logscale          = TRUE,
  nullval           = 1,
  col.left          = "n",
  stroke            = 0.8,
  shape             = 16,
  estcolumn         = TRUE,
  xlim              = c(0.4, 10),
  xticks            = c(1, 2, 4, 8),
  col.right.heading = "RS, 95% CI",
  col.left.heading  = "(A[G]/S[G])/\n(A![G]/S![G])",
  xlab              = "Phase I-Launch Relative Success (RS), 95% CI",
  row.labels        = row_labels_1c
)

save_figure(fig1c$plot, "fig1c_forest_gene_annot", width = 10, height = 4)

# --- 7. Summary --------------------------------------------------------------

message("[7/7] Complete!")
message("
================================================================================
  Figure 1 generation complete
================================================================================

Output files in: ", fig_dir, "

Individual panels (PDF + PNG at 300 DPI):
  - fig1a_forest_main.pdf/.png       Main forest plot
  - fig1b_upset.pdf/.png             UpSet plot (MR/coloc overlap)
  - fig1c_forest_gene_annot.pdf/.png Gene family forest plot

Combine panels manually in Google Slides, PowerPoint, or Illustrator
for the final Figure 1.

")