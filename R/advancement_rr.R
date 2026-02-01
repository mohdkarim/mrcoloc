advancement_rr = function(best, alpha = 0.05, threshold = NA) {
  
  phase_map = tibble(phase=factor(c('Preclinical','I','II','III','I-Launch'),ordered=T,levels=c('Preclinical','I','II','III','I-Launch')), 
                     phorder = 0:4,
                     varname=c('succ_p_1','succ_1_2','succ_2_3','succ_3_a','succ_1_a'))
  
  # determine whether operating on a pipeline_best output that required genetic insight
  require_insight = 'indication lacks genetic insight' %in% best$target_status
  if (is.na(threshold)) {
    best$gensup = best$target_status=='genetically supported target'
  } else {
    best$gensup = !is.na(best$similarity) & best$similarity >= threshold
  } 
  
  best %>%
    filter(genetic_insight != 'none' | (!require_insight)) %>%
    select(ti_uid, gensup, succ_p_1:succ_3_a) %>%
    pivot_longer(succ_p_1:succ_3_a) %>%
    inner_join(phase_map, by=c('name'='varname')) %>%
    filter(!is.na(value)) %>%
    rename(success=value) %>%
    mutate(gs = ifelse(gensup,'yes','no')) %>%
    select(ti_uid, gs, phase, success) -> long
  
  denoms_structure = tibble(gs=c('no','yes'))
  
  long %>%
    mutate(gs = as.character(gs)) %>%
    filter(phase != 'Preclinical') %>%
    group_by(gs) %>%
    summarize(.groups='keep',
              denom = length(unique(ti_uid))) %>%
    ungroup() %>%
    right_join(denoms_structure, by='gs') %>%
    mutate(denom = replace_na(denom, 0)) -> denoms
  
  rr_long_structure = crossing(gs=c('yes','no'),phase=c('Preclinical','I','II','III'))
  
  # Wilson CI for single phases
  long %>%
    mutate(gs = as.character(gs)) %>%
    group_by(gs, phase) %>%
    summarize(.groups='keep',
              x = sum(success),
              n = sum(!is.na(success))) %>%
    ungroup() %>%
    right_join(rr_long_structure, by=c('gs','phase')) %>%
    mutate(x=replace_na(x, 0),
           n=replace_na(n, 0)) %>%
    mutate(binom = binom.confint(x, n, 1-alpha, method='wilson')[,c('mean','lower','upper')]) %>%
    mutate(mean = binom$mean, l=binom$lower, u=binom$upper) %>%
    select(gs, phase, x, n, mean, l, u) -> rs_long
  
  # Wald CI for product of P(S) across I-Launch
  rs_long %>%
    filter(phase != 'Preclinical' & phase != 'I-Launch') %>%
    group_by(gs) %>%
    summarize(.groups='keep',
              x = x[phase=='III'],
              m = prod(mean),
              l = prod(mean) - qnorm(1 - (alpha)/2) * sqrt(prod(mean * (1 - mean) / n + mean^2) - prod(mean)^2), 
              u = prod(mean) + qnorm(1 - (alpha)/2) * sqrt(prod(mean * (1 - mean) / n + mean^2) - prod(mean)^2)) %>%
    rename(mean=m) %>%
    ungroup() %>%
    left_join(denoms, by='gs') %>%
    mutate(n = denom) %>%
    mutate(phase='I-Launch') %>%
    select(gs, phase, x, n, mean, l, u) -> ilrows
  
  rs_long %>%
    bind_rows(ilrows) %>%
    pivot_wider(id_cols=phase, names_from = gs, names_sep='_', values_from=c(x,n,mean,l,u)) %>%
    inner_join(select(phase_map, phase, phorder), by='phase') %>%
    arrange(phorder) %>%
    select(-phorder) %>%
    ungroup() %>%
    mutate(
      binratio = suppressWarnings(binom_ratio(x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha = 0.05))) %>%
    mutate(
      rs_mean = binratio$rs_mean,
      rs_l = binratio$rs_l,
      rs_u = binratio$rs_u,
      fraction = glue("{x_yes}/{n_yes}") 
    ) %>%
    select(phase, x_yes, n_yes, x_no, n_no, mean_yes, l_yes, u_yes, mean_no, l_no, u_no, rs_mean, rs_l, rs_u, fraction) -> rr
  
  return (rr)
}

binom_ratio = function(x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha =  0.05) {
  setNames(as_tibble(t(mapply(binom_ratio_atomic, x_yes, n_yes, x_no, n_no, mean_yes, mean_no, alpha))), c('rs_mean','rs_l','rs_u'))
}

binom_ratio_atomic = function(x_yes, n_yes, x_no, n_no, mean_yes=NA, mean_no=NA, alpha =  0.05) {
  
  # Wald by default
  if (!always_katz & x_yes >= 3 & x_no >=3) {
    mean = mean_yes / mean_no
    #mean = (x_yes/n_yes) / (x_no/n_no)
    lower = mean - qnorm(1 - alpha/2) * sqrt(1/(n_yes + n_no) * mean^2 * (1/mean_yes + 1/mean_no))
    upper = mean + qnorm(1 - alpha/2) * sqrt(1/(n_yes + n_no) * mean^2 * (1/mean_yes + 1/mean_no))
  } else { # Katz for small N
    # for the I-Launch row, mean_yes or mean_no may be NA because a phase had no data
    # important to leave that NA - don't use the filled-in total denominator
    if (is.na(mean_yes) | is.na(mean_no)) {
      mean = lower = upper = as.numeric(NA)
    } else {
      binom_obj = BinomRatioCI(x_yes, n_yes, x_no, n_no, conf=1-alpha)
      mean = binom_obj[1,'est']
      lower = binom_obj[1,'lwr.ci']
      upper = binom_obj[1,'upr.ci']
    }
  }
  return (cbind(mean=mean,lower=lower,upper=upper))
}

always_katz = T