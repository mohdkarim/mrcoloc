suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DescTools)
})

merge3_pqtl <- readRDS("data_raw/merge3_pqtl.rds")
indic <- read_tsv("data/indic.tsv", show_col_types = FALSE)
pgenes <- readRDS("data_raw/pgenes.rds")
# [R2.10 FIX] per-platform ASSAY sets; see create_derived_data.R
pgenes_platform <- readRDS("data_raw/pgenes_platform.rds")

SIMILARITY_THR <- 0.8
MIN_L2G_SHARE  <- 0.5

cat("\n================================================================\n")
cat("  REVIEWER 2 COMMENT 10: Platform counts discrepancy\n")
cat("================================================================\n\n")

cat("Platforms in data:\n")
print(table(merge3_pqtl$platform, useNA = "ifany"))

cat("\n\n--- Reproducing the per-platform RS analysis ---\n\n")

plat_res <- list()

for (p in c("Somascan", "Olink")) {
  df_platform_support <- merge3_pqtl %>%
    filter(grepl("pqtl", original_link, ignore.case = TRUE),
           comb_norm >= SIMILARITY_THR,
           !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
           platform == p) %>%
    distinct(ti_uid)

  # [R2.10 FIX] background = proteins measured on this platform, was `platform == p`
  ti_best_platform <- merge3_pqtl %>%
    filter(!is.na(gene), gene != "",
           !is.na(indication_mesh_id), indication_mesh_id != "",
           !is.na(ccat), gene %in% pgenes_platform[[p]]) %>%
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
  total_gs   <- sum(baseline$gensup, na.rm = TRUE)
  succ_nogs  <- sum(!long$gensup & long$success, na.rm = TRUE)
  total_nogs <- sum(!baseline$gensup, na.rm = TRUE)

  cat(sprintf("Platform: %s\n", p))
  cat(sprintf("  Supported TIs:   %d\n", nrow(df_platform_support)))
  cat(sprintf("  RS string:       (%d/%d)/(%d/%d)\n", succ_gs, total_gs, succ_nogs, total_nogs))
  cat(sprintf("  Universe (Phase I TIs): %d\n", nrow(baseline)))
  cat(sprintf("  Unique genes in universe: %d\n\n", length(unique(baseline$gene))))

  plat_res[[p]] <- tibble(
    platform = p,
    x_gs     = succ_gs,   n_gs   = total_gs,
    x_nogs   = succ_nogs, n_nogs = total_nogs,
    universe = nrow(baseline)
  )
}

plat_res <- bind_rows(plat_res)

cat("\n--- The main pQTL row (all platforms combined, pgenes universe) ---\n\n")

pqtl_support_all <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE) %>%
  distinct(ti_uid)

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
  mutate(gensup = ti_uid %in% pqtl_support_all$ti_uid)

long_all     <- ti_best_pgenes %>% filter(!is.na(succ_3_a)) %>% rename(success = succ_3_a)
baseline_all <- ti_best_pgenes %>% filter(!is.na(succ_1_2))

x_all    <- sum(long_all$gensup & long_all$success)
n_all    <- sum(baseline_all$gensup)
x_no_all <- sum(!long_all$gensup & long_all$success)
n_no_all <- sum(!baseline_all$gensup)

cat(sprintf("  All platforms (pgenes background):\n"))
cat(sprintf("  RS string: (%d/%d)/(%d/%d)\n", x_all, n_all, x_no_all, n_no_all))
cat(sprintf("  Universe (Phase I TIs): %d\n\n", nrow(baseline_all)))

cat("\n================================================================\n")
cat("  WHY THE NUMBERS DIFFER\n")
cat("================================================================\n\n")

cat("1) DIFFERENT UNIVERSES per platform:\n")
cat("   The per-platform row restricts BOTH numerator AND denominator to rows\n")
cat("   where platform == X. This means the background (unsupported) is only\n")
cat("   TIs for genes assayed on THAT specific platform.\n\n")

soma_genes <- merge3_pqtl %>% filter(platform == "Somascan") %>% distinct(gene) %>% pull()
olink_genes <- merge3_pqtl %>% filter(platform == "Olink") %>% distinct(gene) %>% pull()
cat(sprintf("   Somascan genes in data:  %d\n", length(soma_genes)))
cat(sprintf("   Olink genes in data:     %d\n", length(olink_genes)))
cat(sprintf("   Overlap (both):          %d\n", sum(soma_genes %in% olink_genes)))
cat(sprintf("   Union:                   %d\n", length(union(soma_genes, olink_genes))))
cat(sprintf("   pgenes (full set):       %d\n\n", length(pgenes)))

cat("2) PLATFORM OVERLAP in pQTL-supported TIs:\n")

soma_support <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
         platform == "Somascan") %>%
  distinct(ti_uid) %>% pull()

olink_support <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
         platform == "Olink") %>%
  distinct(ti_uid) %>% pull()

cat(sprintf("   pQTL-supported TIs (Somascan): %d\n", length(soma_support)))
cat(sprintf("   pQTL-supported TIs (Olink):    %d\n", length(olink_support)))
cat(sprintf("   Overlap (both platforms):      %d\n", sum(soma_support %in% olink_support)))
cat(sprintf("   Union:                         %d\n", length(union(soma_support, olink_support))))
cat(sprintf("   Main pQTL row (all):           %d\n\n", nrow(pqtl_support_all)))

cat("3) Launched pQTL-supported TIs by platform:\n\n")
launched_pqtl <- merge3_pqtl %>%
  filter(grepl("pqtl", original_link, ignore.case = TRUE),
         comb_norm >= SIMILARITY_THR,
         !is.na(l2g_share), l2g_share >= MIN_L2G_SHARE,
         succ_3_a == TRUE) %>%
  distinct(ti_uid, platform)

platform_summary <- launched_pqtl %>%
  group_by(ti_uid) %>%
  summarise(platforms = paste(sort(unique(platform)), collapse = "+"), .groups = "drop")

print(table(platform_summary$platforms))
cat(sprintf("\n   Unique launched TIs: %d\n", nrow(platform_summary)))

cat("\n4) RECONCILIATION WITH THE REVIEWER'S ARITHMETIC:\n")
cat("   Naively summing the two per-platform rows of Figure 1a reproduces the\n")
cat("   figures quoted in the comment. Shown against the headline pQTL row:\n\n")

cat(sprintf("   Summed platform rows:  (%d/%d)/(%d/%d)   [reviewer read 26/58 and 494/4102]\n",
    sum(plat_res$x_gs), sum(plat_res$n_gs), sum(plat_res$x_nogs), sum(plat_res$n_nogs)))
cat(sprintf("   Headline pQTL row:     (%d/%d)/(%d/%d)\n", x_all, n_all, x_no_all, n_no_all))
cat(sprintf("   Launched pairs double-counted by the sum: %d (%d summed - %d unique)\n",
    sum(plat_res$x_gs) - x_all, sum(plat_res$x_gs), x_all))
cat(sprintf("   Phase I universes: %s; headline %d (union of all 8 datasets)\n",
    paste(sprintf("%s %d", plat_res$platform, plat_res$universe), collapse = ", "),
    nrow(baseline_all)))
# [R2.10 FIX] This explanation was previously the wrong way round. Before the fix
# the platform backgrounds were conditioned on having a significant association,
# so their sum was SMALLER than the headline denominator. Now that they are
# assay-based they OVERLAP instead, and their sum EXCEEDS the headline.
cat(sprintf("   Summed platform universes %d vs headline %d; the difference is the\n",
    sum(plat_res$universe), nrow(baseline_all)))
cat(sprintf("   %d Phase I pairs whose target is measured on both platforms.\n",
    sum(plat_res$universe) - nrow(baseline_all)))
cat("   Neither numerators nor denominators sum, for the same reason: overlap.\n")

cat("\n================================================================\n")
cat("  ANSWER TO REVIEWER\n")
cat("================================================================\n\n")
cat("The discrepancy arises because:\n\n")
cat("a) Per-platform denominators differ from the main row because each\n")
cat("   platform analysis restricts the background to genes measured on\n")
cat("   THAT platform only. Somascan and Olink measure different protein\n")
cat("   sets, so the denominators are smaller (platform-specific subsets).\n\n")
cat("b) Per-platform numerators sum to > 23 because some TI pairs have\n")
cat("   pQTL support from BOTH Somascan and Olink studies (double-counted\n")
cat("   when summing across platforms). The main pQTL row counts each TI\n")
cat("   only once regardless of how many platforms provide evidence.\n\n")
cat("c) The main row uses pgenes (union of ALL platform genes) as\n")
cat("   background, which is larger than any single platform.\n")
