#!/usr/bin/env Rscript
# ============================================================================
# R2.6: cluster-robust confidence intervals for relative success
# ============================================================================
# Addresses Reviewer 2 comment 6: "Are the confidence intervals appropriately
# calibrated given that many TIs hit the same target?"
#
# Method: relative success is a risk ratio, so we fit a log-link Poisson GLM
#   launched ~ pQTL_support
# and take a cluster-robust (sandwich) variance clustered on TARGET GENE. This is
# the modified-Poisson-with-robust-variance estimator for a risk ratio
# (Zou 2004), with the clustered extension (Zou & Donner 2013). The published
# Katz intervals are retained as primary for comparability with Minikel et al;
# these are a sensitivity analysis.
#
# Clustering is applied to ALL observations in both arms (753 target clusters),
# not only to the 12 targets contributing launched supported pairs.
#
# Two independent implementations are computed and cross-checked:
#   (a) survey::svyglm with ids = ~gene   (linearisation)
#   (b) a hand-coded CR1 sandwich          (no package dependency)
# plus (c) a cluster bootstrap over target genes, as a third check.
#
# Outputs: output/ST21_cluster_robust_ci.{tsv,xlsx}
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
})

project_root <- Sys.getenv("MRCOLOC_ROOT", getwd())
out_dir      <- file.path(project_root, "output")
pgenes_file  <- file.path(project_root, "data_raw", "pgenes.rds")

m     <- readRDS(file.path(project_root, "data_raw/merge3_pqtl.rds"))
indic <- read_tsv(file.path(project_root, "data/indic.tsv"), show_col_types = FALSE)
pg    <- readRDS(pgenes_file)

SIM <- 0.8; L2G <- 0.5

supp_set <- function(min_l2g = L2G) {
  m %>% filter(grepl("pqtl", original_link, ignore.case = TRUE),
               comb_norm >= SIM, !is.na(l2g_share), l2g_share >= min_l2g) %>%
    distinct(ti_uid) %>% pull()
}

ti_best <- m %>%
  filter(!is.na(gene), gene != "", !is.na(indication_mesh_id), indication_mesh_id != "",
         !is.na(ccat), gene %in% pg) %>%
  left_join(indic %>% select(indication_mesh_id, genetic_insight), by = "indication_mesh_id") %>%
  filter(genetic_insight != "none") %>%
  mutate(hp = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2,
                        !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0, TRUE ~ NA_real_)) %>%
  arrange(ti_uid, desc(hp), desc(comb_norm)) %>%
  group_by(ti_uid) %>% slice(1) %>% ungroup()

# --- interval estimators -----------------------------------------------------

katz <- function(x1, n1, x2, n2) {
  r <- (x1 / n1) / (x2 / n2)
  s <- sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2)
  c(rs = r, lwr = r * exp(-1.96 * s), upr = r * exp(1.96 * s), se_log = s)
}

# hand-coded CR1 cluster-robust sandwich for a log-link Poisson GLM
cr1_sandwich <- function(fit, cluster) {
  X <- model.matrix(fit)
  u <- residuals(fit, type = "response")          # y - mu  (score weight for log link)
  bread <- solve(crossprod(X, X * fit$fitted.values))
  uX <- X * u
  meat <- rowsum(uX, group = cluster) %>% { crossprod(.) }
  G <- length(unique(cluster)); N <- nrow(X); K <- ncol(X)
  adj <- (G / (G - 1)) * ((N - 1) / (N - K))       # CR1 small-sample correction
  bread %*% meat %*% bread * adj
}

cluster_ci <- function(dat) {
  d <- dat %>% transmute(y = as.integer(success), gs = as.integer(gs), gene)
  fit <- glm(y ~ gs, family = poisson(link = "log"), data = d)
  b   <- coef(fit)["gs"]

  # (b) hand-coded CR1
  V   <- cr1_sandwich(fit, d$gene)
  se_cr1 <- sqrt(V["gs", "gs"])

  # (a) survey linearisation
  des <- svydesign(ids = ~gene, weights = ~1, data = d)
  sfit <- svyglm(y ~ gs, design = des, family = quasipoisson(link = "log"))
  se_svy <- summary(sfit)$coefficients["gs", "Std. Error"]

  # (c) cluster bootstrap over target genes.
  # Index rows by cluster ONCE - a per-replicate linear scan over 753 genes x
  # 7,215 rows is ~10^10 operations and does not finish.
  set.seed(1)
  idx_by_gene <- split(seq_len(nrow(d)), d$gene)
  genes <- names(idx_by_gene)
  yv <- d$y; gv <- d$gs
  bs <- replicate(2000, {
    g   <- sample.int(length(genes), length(genes), replace = TRUE)
    ii  <- unlist(idx_by_gene[g], use.names = FALSE)
    y1  <- yv[ii][gv[ii] == 1]; y0 <- yv[ii][gv[ii] == 0]
    if (!length(y1) || !length(y0)) return(NA_real_)
    p1 <- mean(y1); p0 <- mean(y0)
    if (p1 == 0 || p0 == 0) return(NA_real_)
    log(p1 / p0)
  })
  bs <- bs[is.finite(bs)]

  list(rs = exp(b),
       cr1 = exp(b + c(-1.96, 1.96) * se_cr1), se_cr1 = se_cr1,
       svy = exp(b + c(-1.96, 1.96) * se_svy), se_svy = se_svy,
       boot = exp(quantile(bs, c(0.025, 0.975), na.rm = TRUE)),
       n_boot_ok = length(bs),
       n_clusters = length(genes))
}

# --- rows to report ----------------------------------------------------------

rows <- list(
  list(label = "pQTL (+ L2G >= 0.5)  [headline]", min_l2g = 0.50),
  list(label = "pQTL + L2G >= 0.25",              min_l2g = 0.25),
  list(label = "pQTL + L2G >= 0.75",              min_l2g = 0.75)
)

res <- map_dfr(rows, function(r) {
  s  <- supp_set(r$min_l2g)
  b  <- ti_best %>% mutate(gs = ti_uid %in% s)
  L  <- b %>% filter(!is.na(succ_3_a)) %>% mutate(success = succ_3_a)
  B  <- b %>% filter(!is.na(succ_1_2))
  x1 <- sum(L$gs & L$success); n1 <- sum(B$gs)
  x2 <- sum(!L$gs & L$success); n2 <- sum(!B$gs)
  k  <- katz(x1, n1, x2, n2)

  # cluster analysis is on the Phase I cohort with launch as the outcome
  dat <- B %>% mutate(success = ti_uid %in% L$ti_uid[L$success])
  cl  <- cluster_ci(dat)

  tibble(
    row = r$label,
    counts = sprintf("(%d/%d)/(%d/%d)", x1, n1, x2, n2),
    rs = k["rs"],
    katz = sprintf("%.2f-%.2f", k["lwr"], k["upr"]),
    cluster_cr1 = sprintf("%.2f-%.2f", cl$cr1[1], cl$cr1[2]),
    cluster_survey = sprintf("%.2f-%.2f", cl$svy[1], cl$svy[2]),
    cluster_boot = sprintf("%.2f-%.2f", cl$boot[1], cl$boot[2]),
    se_ratio = cl$se_cr1 / k["se_log"],
    n_targets = cl$n_clusters,
    boot_ok = cl$n_boot_ok
  )
})

cat("\n=== Relative success: Katz vs cluster-robust intervals (clustered on target gene) ===\n\n")
print(as.data.frame(res %>% mutate(rs = round(rs, 2), se_ratio = round(se_ratio, 2))), row.names = FALSE)
cat("\nse_ratio = cluster-robust SE(log RS) / Katz SE(log RS). >1 means clustering widens the interval.\n")
cat("Three estimators shown as a cross-check; they should agree closely.\n")

write_tsv(res, file.path(out_dir, "r2_6_cluster_robust_ci.tsv"))
cat("\n-> output/r2_6_cluster_robust_ci.tsv\n")

# --- ST21-format table (provisional number) ---------------------------------
# Column names follow ST4's convention (panel_group, source_label, rs_estimate,
# rs_lwr_95ci, rs_upr_95ci, count_string) so it can be dropped into
# generate_mrcoloc_supplement.R unchanged when the fix is ported.

st21 <- res %>%
  mutate(
    panel_group   = "Cluster-robust sensitivity",
    source_label  = row,
    rs_estimate   = round(rs, 3),
    count_string  = counts,
    katz_95ci             = katz,
    cluster_robust_95ci   = cluster_cr1,
    cluster_bootstrap_95ci = cluster_boot,
    se_inflation_vs_katz  = round(se_ratio, 2),
    n_target_clusters     = n_targets
  ) %>%
  select(panel_group, source_label, count_string, rs_estimate, katz_95ci,
         cluster_robust_95ci, cluster_bootstrap_95ci, se_inflation_vs_katz,
         n_target_clusters)

write_tsv(st21, file.path(out_dir, "ST21_cluster_robust_ci.tsv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  sheet <- "ST21 - Cluster_robust_CI"
  addWorksheet(wb, sheet)
  title <- paste("Supplementary Table 21: Relative success with confidence intervals",
                 "accounting for repeated targets (clustered on target gene)")
  writeData(wb, sheet, title, startRow = 1, startCol = 1)
  mergeCells(wb, sheet, cols = 1:ncol(st21), rows = 1)
  addStyle(wb, sheet, createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  writeData(wb, sheet, st21, startRow = 3, startCol = 1,
            headerStyle = createStyle(textDecoration = "bold", border = "bottom"))
  setColWidths(wb, sheet, cols = 1:ncol(st21), widths = "auto")
  freezePane(wb, sheet, firstActiveRow = 4)
  note <- c(
    "Primary intervals in Figure 1a and Supplementary Table 4 are Katz intervals, following Minikel et al, so that",
    "estimates remain comparable with the L2G-based estimates reproduced in the same figure.",
    "Cluster-robust intervals use a CR1 sandwich variance on a log-link Poisson model of launch on pQTL support,",
    "clustered on target gene across all target-indication pairs in both arms (Zou 2004; Zou and Donner 2013).",
    "Cluster bootstrap resamples target genes with replacement, 2000 replicates, percentile interval.",
    "Cluster-robust intervals must not be compared against the unadjusted Katz intervals of the L2G rows."
  )
  writeData(wb, sheet, data.frame(Notes = note), startRow = nrow(st21) + 5, startCol = 1)
  saveWorkbook(wb, file.path(out_dir, "ST21_cluster_robust_ci.xlsx"), overwrite = TRUE)
  cat("-> out/ST21_cluster_robust_ci.tsv + .xlsx\n")
}
