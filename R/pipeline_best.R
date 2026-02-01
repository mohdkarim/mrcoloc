overall_start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...')

options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(binom))
suppressMessages(library(glue))
suppressMessages(library(lawstat))
suppressMessages(library(weights))
suppressMessages(library(epitools))
suppressMessages(library(DescTools))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))
suppressMessages(library(MASS)); summarize=dplyr::summarize; select=dplyr::select; rename=dplyr::rename; slice=dplyr::slice
if(interactive()) {
  setwd('/home/mohd/mohd-sandbox/pQTL_enrichment/mrcoloc_paper2025/genetic_support-main')
}


option_list = list(
  make_option(c("-o", "--oto"), action="store_true", default=FALSE, 
              help="limit to drugs with one target only (oto mode) [default %default]")
); 
opt = parse_args(OptionParser(option_list=option_list))

pp = read_tsv('data/pp.tsv', col_types=cols())
drug_phase_summary = read_tsv('data/drug_phase_summary.tsv', col_types=cols())
assoc = read_tsv('data/assoc.tsv.gz', col_types=cols())
indic = read_tsv("data/indic.tsv", col_types=cols())
indic_topl_match = read_tsv('data/indic_topl_match.tsv', col_types=cols())
universe = read_tsv('data/universe.tsv', col_types=cols())
meta_hcat = read_tsv('data/meta_hcat.tsv', col_types=cols())
meta_acat = read_tsv('data/meta_acat.tsv', col_types=cols())
meta_ccat = read_tsv('data/meta_ccat.tsv', col_types=cols())
mesh_best_names = read_tsv('data/mesh_best_names.tsv.gz', col_types=cols())
sim = read_tsv('data/sim.tsv.gz', col_types=cols())

if (opt$oto) {
  pp$hcat = pp$oto_hcat
  pp$acat = pp$oto_acat
  pp$ccat = pp$oto_ccat
  pp$hcatnum = meta_hcat$num[match(pp$hcat, meta_hcat$cat)]
  pp$acatnum = meta_acat$num[match(pp$acat, meta_hcat$cat)]
  pp$ccatnum = meta_ccat$num[match(pp$ccat, meta_hcat$cat)]
  merge2$hcat = pp$oto_hcat[match(merge2$ti_uid, pp$ti_uid)]
  merge2$acat = pp$oto_acat[match(merge2$ti_uid, pp$ti_uid)]
  merge2$ccat = pp$oto_ccat[match(merge2$ti_uid, pp$ti_uid)]
  pp = pp[!is.na(pp$ccat),]
  merge2 = merge2[!is.na(merge2$ccat),]
}

# constants

active_clinical = tibble(cat=c('Phase I','Phase II','Phase III'))

####
# Functions
####

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x < 0, '-', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

abs_or = function(odds_ratio) {
  abs_odds_ratio = odds_ratio
  flip_indices = odds_ratio < 1 & !is.na(odds_ratio)
  abs_odds_ratio[flip_indices] = 1/odds_ratio[flip_indices]
  return (abs_odds_ratio)
}


pipeline_best = function(merged_table,
                         basis='ti',
                         phase='combined',
                         require_insight=TRUE,
                         share_mode='L2G', # other option is V2G
                         min_share=0.5,
                         max_share=1.0,
                         worst_rank=Inf, # set to Inf if you want to include all
                         min_h4 = 0.9,
                         include_missing=FALSE,
                         associations=c('OMIM','GWAS'),
                         otg_subcat=c(''),
                         genebass_subcat=NULL,
                         mendelian_mechanism='',
                         min_year=2005,
                         max_year=2022,
                         firstyear=F,
                         minusomim=F,
                         lacking=NULL, # association sources required to be lacked by the T-I
                         andalso=NULL, # association sources required to *also* endorse the T-I
                         minusothersubcat=F,
                         mingenecount=0,
                         maxgenecount=Inf,
                         mapping_basis='all',
                         min_beta=0,
                         max_beta=Inf,
                         min_or=1,
                         max_or=Inf,
                         min_maf=0,
                         max_maf=1,
                         threshold=0.8,
                         network_list=NA,
                         verbose=T) {
  
  start_time = Sys.time()
  
  mtable = merged_table
  if (verbose) {
    cat(file=stderr(),'Starting row count: ',nrow(mtable),'\n')
    flush.console()
  }
  
  # add & select unique ID
  if (basis %in% c('di_mesh','drug-indication')) {
    mtable$di_uid = paste0(mtable$drugid,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$di_uid
  } else if (basis %in% c('ti','target-indication')) {
    mtable$ti_uid = paste0(mtable$gene,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$ti_uid
  } else if (basis=='drug') {
    mtable$uid = mtable$drugid
  }
  
  # assign highest level of advancement depending on phase specified
  if (phase == 'active') {
    meta = meta_acat
    mtable$cat = mtable$acat
  } else if (phase == 'historical') {
    meta = meta_hcat
    mtable$cat = mtable$hcat
  } else if (phase == 'combined') {
    meta = meta_ccat
    mtable$cat = mtable$ccat
  }
  
  mtable$catnum = meta$num[match(mtable$cat, meta$cat)]
  # remove "Other"
  mtable = mtable[mtable$cat != 'Other' ,]
  # map L2G share
  mtable$assoc_share = mtable$l2g_share
  mtable$assoc_rank = mtable$l2g_rank
  
  if (verbose) {
    cat(file=stderr(),'Selecting user-specified filters...')
    flush.console()
  }
  
  # by default, require non-missing target & indication
  if (!include_missing) {
    mtable = mtable[mtable$gene != '' & mtable$indication_mesh_id != '' & !is.na(mtable$gene) & !is.na(mtable$indication_mesh_id),]
  }
  # genetic insight requirement
  if (require_insight) {
    mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id[indic$genetic_insight != 'none'],]
  } else {
    # otherwise simply require the indication be present in the indic table
    mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id,]
  }
  
  # remove omim-supported associations if desired. only works in T-I mode
  # note that order of operations is important - this must come before associations filter
  if (minusomim) {
    # look for first year in which a target-*indication* pair was genetically supported
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% c('OMIM')) %>%
      select(ti_uid) -> omim_supported_ti
    # retain the null rows (i.e. no association) or those where hte T-I is not in OMIM
    # what gets removed? e.g. OTG associations that were already established by OMIM
    mtable %>%
      filter(is.na(assoc_source) | !(ti_uid %in% omim_supported_ti$ti_uid)) -> mtable
  }
  
  # remove any association sources required to be lacked
  if (!is.null(lacking)) {
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% lacking) %>%
      filter(!(assoc_source %in% 'OTG' & l2g_share >= min_share)) %>%
      select(ti_uid) -> lackable_supported_ti
    mtable %>%
      filter(is.na(assoc_source) | !(ti_uid %in% lackable_supported_ti$ti_uid)) -> mtable
  }
  
  if (!is.null(andalso)) {
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% andalso) %>%
      filter(!(assoc_source %in% 'OTG' & l2g_share < min_share)) %>%
      select(ti_uid) -> andalso_supported_ti
    mtable %>%
      filter(is.na(assoc_source) | (ti_uid %in% andalso_supported_ti$ti_uid)) -> mtable
  }
  
  # user-specified sources of genetic associations
  # allow user to specify either grouping terms like "GWAS", or specific sources
  source_map = tibble(source=c("OTG", "PICCOLO", "Genebass", "OMIM", "intOGen"),
                      source_name=c('GWAS','GWAS','GWAS','OMIM','Somatic'))
  if (!(identical(associations, c('OMIM','GWAS','Somatic'))) ) {
    mtable %>%
      left_join(source_map, by=c('assoc_source'='source')) %>%
      filter(is.na(source_name) | source_name %in% associations | assoc_source %in% associations) -> mtable
    
  }
  
  # further filter of subtype of OTG association
  if (otg_subcat != '') {
    # first subset to just OTG
    mtable %>%
      filter(is.na(assoc_source) | assoc_source %in% 'OTG') -> mtable
    # now pick the types
    mtable %>%
      mutate(gwas_source = case_when(grepl('GCST', original_link) ~ 'GWAS Catalog',
                                     grepl('FINNGEN', original_link) ~ 'FinnGen',
                                     grepl('NEALE',original_link) ~ 'Neale UKBB',
                                     grepl('pqtl',original_link) ~ 'pQTL',
                                     TRUE ~ 'Other')) %>%
      filter(is.na(assoc_source) | gwas_source %in% otg_subcat) -> mtable
  }
  
  # further filter of annotation & test in Genebass
  if (!is.null(genebass_subcat)) {
    grepstring = paste(genebass_subcat, collapse='|')
    # mtable %>%
    #   filter(!is.na(assoc_source) & assoc_source %in% 'Genebass') %>%
    #   filter(grepl(grepstring, assoc_info)) -> genebass_hits
    mtable %>%
      filter(is.na(assoc_source) | !(assoc_source %in% 'Genebass') | grepl(grepstring, assoc_info)) -> mtable
  }
  
  # apply user-specified filter of OMIM disease mechanism
  if (mendelian_mechanism != '') {
    mtable %>%
      filter(!(mtable$assoc_source %in% 'OMIM') | is.na(mtable$assoc_info) | grepl(mendelian_mechanism,mtable$assoc_info)) -> mtable
  }
  
  # apply user-specified OTG gene mapping share & rank minimum/maximum
  if (share_mode == 'V2G') {
    mtable$assoc_share = mtable$v2g_share
    mtable$assoc_rank = mtable$v2g_rank
    assoc$assoc_share = assoc$v2g_share # needed in assoc table too for genecount section below
  } else if (share_mode == 'L2G') {
    mtable$assoc_share = mtable$l2g_share
    mtable$assoc_rank = mtable$l2g_rank
    assoc$assoc_share = assoc$l2g_share
  }
  
  # worst rank
  if (worst_rank < Inf) {
    mtable %>%
      filter(!(assoc_source %in% 'OTG') | mtable$assoc_rank <= worst_rank) -> mtable
  }
  
  # note that among OTG associations, throw out any with NA share, as these would be zeroes (does not occur in Dec 2021 dataset anyway)
  # and note that with L2G a significant number of associations actually have 100% share, so only delete > max_share and not >= max_share
  mtable %>%
    filter(!(assoc_source %in% 'OTG') | (!is.na(assoc_share) & assoc_share >= min_share & assoc_share <= max_share)) -> mtable
  
  # apply user-specified H4 minimum / maximum
  if (min_h4 > .9) {
    mtable %>%
      filter(!(assoc_source %in% 'PICCOLO') | (!is.na(mtable$pic_h4) & mtable$pic_h4 >= min_h4)) -> mtable
  }
  
  # apply user-specified genecount minimum/maximum
  if (mingenecount > 0 | maxgenecount < Inf) {
    assoc %>%
      filter(source %in% associations) %>%
      filter(source!='OTG' | (assoc_share >= min_share & assoc_share <= max_share)) %>%
      group_by(mesh_id) %>%
      summarize(.groups='keep', n_genes=length(unique(gene))) -> gene_counts
    mtable$gene_count = gene_counts$n_genes[match(mtable$assoc_mesh_id, gene_counts$mesh_id)]
    mtable %>%
      filter(is.na(gene_count) | gene_count >= mingenecount & gene_count <= maxgenecount) -> mtable
  }
  
  
  # apply "first year" criteria if applicable
  if (firstyear) {
    # look for first year in which a target-*indication* pair was genetically supported
    # only works in T-I mode
    mtable %>%
      filter(comb_norm >= threshold) %>%
      group_by(ti_uid) %>%
      summarize(.groups='keep', min_assoc_year=min(assoc_year)) -> ti_first_sup
    mtable$min_assoc_year = ti_first_sup$min_assoc_year[match(mtable$ti_uid, ti_first_sup$ti_uid)]
    # keep entries with no assoc year (the null rows), or where assoc year = the min assoc year, i.e. this is
    # the first report of this genetic association (or tied for first)
    mtable %>%
      filter(is.na(assoc_year) | assoc_year == min_assoc_year) -> mtable
  }
  
  # avoid -Inf values in comparisons by hard coding in case of all missing values:
  if (sum(!is.na(mtable$assoc_year)) == 0) {
    mtable_intrinsic_max_year = 2021
    mtable_intrinsic_min_year = 2000
  } else {
    mtable_intrinsic_max_year = max(mtable$assoc_year, na.rm=T)
    mtable_intrinsic_min_year = min(mtable$assoc_year, na.rm=T)
  }
  
  # apply user-specified filter of association years - for OTG only
  if (min_year > mtable_intrinsic_min_year | max_year < mtable_intrinsic_max_year) {
    mtable %>%
      filter(is.na(assoc_source) | assoc_source != 'OTG' | assoc_year >= min_year & assoc_year <= max_year) -> mtable
  }
  
  # join back in beta
  mtable$abs_beta = abs(assoc$beta[match(mtable$arow, assoc$arow)])
  if (min_beta > 0 | max_beta < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_beta) & mtable$abs_beta >= min_beta & mtable$abs_beta <= max_beta)) -> mtable
  }
  
  # same as beta but for OR
  mtable$abs_or = abs_or(assoc$odds_ratio[match(mtable$arow, assoc$arow)])  
  if (min_or > 1 | max_or < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    # leave null rows (with no association source) but delete all those mapped to OTG that do not have or, or have or outside range
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_or) & mtable$abs_or >= min_or & mtable$abs_or < max_or)) -> mtable
  }
  
  
  # lead SNP maf
  if (min_maf > 0 | max_maf < 1) {
    mtable$lead_maf = pmin(mtable$af_gnomad_nfe, 1-mtable$af_gnomad_nfe)
    mtable$lead_maf[!is.na(mtable$lead_maf) & mtable$lead_maf < 0] = NA
    # lead_maf >= min_maf & lead_maf < max_maf
    # >= and < gets you "[, )" logic
    # also remove those that are GWAS where lead_maf is NA - likely in non-European populations so af_gnomad_nfe is not relevant
    mtable %>% 
      filter(is.na(assoc_source) | !(assoc_source %in% c('OTG','PICCOLO','Genebass')) | (!is.na(lead_maf) & (lead_maf >= min_maf & lead_maf < max_maf))) -> mtable 
  }
  
  
  
  
  if (verbose) {
    cat(file=stderr(),'Selecting highest phase reached and best genetic similarity...')
    flush.console()
  }
  
  
  
  suppressWarnings(mtable %>% group_by(uid) %>% summarize(maxsim = max(comb_norm, na.rm=T), maxcat=max(catnum, na.rm=T)) -> step1)
  
  if (verbose) {
    cat(file=stderr(),nrow(step1),'rows remain.\n')
    flush.console()
  }
  
  if (verbose) {
    cat(file=stderr(),'Joining back in program details...')
    flush.console()
  }
  
  # add a filter first - only slightly reduces row count
  mtable %>%
    filter(uid %in% step1$uid & comb_norm %in% unique(step1$maxsim) & catnum %in% unique(step1$maxcat)) -> mtable
  
  # use tidy to left join
  step1 %>%
    left_join(mtable, by = c("uid" = "uid", "maxsim" = "comb_norm", "maxcat" = "catnum")) %>%
    rename(similarity=maxsim, catnum=maxcat) -> step2
  

  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows found after join.\n')
    flush.console()
  }
  if (verbose) {
    cat(file=stderr(),'Removing duplicates resulting from ties...')
    flush.console()
  }
  
  # prioritize rows with known outcome at most advanced phase, then de-dup
  step2 %>%
    mutate(highest_phase_with_known_outcome = case_when(!is.na(succ_3_a) ~ 3,
                                                        !is.na(succ_2_3) ~ 2,
                                                        !is.na(succ_1_2) ~ 1,
                                                        !is.na(succ_p_1) ~ 0)) %>%
    arrange(uid, desc(highest_phase_with_known_outcome)) %>%
    group_by(uid) %>%
    slice(1) %>%
    ungroup() -> step2 # row count should drop back to ~ that of step1
  
  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows remain.\n')
    flush.console()
  }
  
  # annotate in additional info:
  step2$areas = indic$areas[match(step2$indication_mesh_id, indic$indication_mesh_id)]
  step2$genetic_insight = replace_na(indic$genetic_insight[match(step2$indication_mesh_id, indic$indication_mesh_id)],'none')
  step2$target_status = ''
  step2$target_status[step2$similarity >= threshold] = 'genetically supported target'
  step2$target_status[step2$similarity <  threshold] = 'unsupported target'
  step2$target_status[require_insight & step2$genetic_insight == 'none'] = 'indication lacks genetic insight'
  step2$target_status[is.na(step2$gene) | step2$gene == ''] = 'no target annotated'
  step2$target_status[is.na(step2$indication_mesh_id) | step2$indication_mesh_id == ''] = 'no indication annotated'
  
  
  if (verbose) {
    cat(file=stderr(),paste0('Using sim threshold ',threshold,', "genetically supported target" rows: ',sum(step2$target_status=='genetically supported target'),'....\n'))
    time_elapsed = (Sys.time() - start_time)
    cat(file=stderr(),'pipeline_best completed in',round(time_elapsed,1),units(time_elapsed),'.\n')
    flush.console()
  }
  
  return (step2)
}