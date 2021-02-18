equiv_linear_inter <- function(meta_samples, full_samples, meta_nocov_samples,
                              dose_range = c(-1, 1), 
                              betas = c(-.6, -.2, .5, .9, .2, .5, -.3, .4),
                              balanced = c(TRUE, FALSE)){
  
  possible_doses <- seq(dose_range[1], dose_range[2], by = .2)
  n <- length(possible_doses)

  all_samps_meta <- meta_samples %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows()
  
  vars_meta <- all_samps_meta[,1:8]
  
  vals_meta <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_meta)),
                          B2 = rep(vars_meta[,2], n),
                          B3 = rep(vars_meta[,3], n),
                          B4 = rep(vars_meta[,4], n),
                          B7 = rep(vars_meta[,7], n),
                          B8 = rep(vars_meta[,8], n))
  
  vals_meta <- vals_meta %>%
    mutate(d1_meta = ((B3)*d0 - B2)/(B4))
  
  
  vals2_meta <- vals_meta %>% 
    group_by(d0) %>% 
    summarize(lower_meta= quantile(d1_meta, probs = .025), 
              upper_meta = quantile(d1_meta, probs = 0.975),
              median_meta = median(d1_meta)) %>%
    mutate(vals_true = ((betas[3])*d0 - betas[2])/(betas[4]))

  
  all_samps_meta_nocov <- meta_nocov_samples %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows()
  
  vars_meta_nocov <- all_samps_meta_nocov[,1:4]
  
  vals_meta_nocov <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_meta_nocov)),
                          B2 = rep(vars_meta_nocov[,2], n),
                          B3 = rep(vars_meta_nocov[,3], n),
                          B4 = rep(vars_meta_nocov[,4], n))
  
  vals_meta_nocov <- vals_meta_nocov %>%
    mutate(d1_meta_nocov = ((B3)*d0 - B2)/(B4))
  
  
  vals2_meta_nocov <- vals_meta_nocov %>% 
    group_by(d0) %>% 
    summarize(lower_meta_nocov = quantile(d1_meta_nocov, probs = .025), 
              upper_meta_nocov = quantile(d1_meta_nocov, probs = 0.975),
              median_meta_nocov = median(d1_meta_nocov))
  
  all_samps_full <- full_samples %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows()
  
  vars_full <- all_samps_full[,1:8]
  
  vals_full <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_full)),
                          B2 = rep(vars_full[,2], n),
                          B3 = rep(vars_full[,3], n),
                          B4 = rep(vars_full[,4], n),
                          B7 = rep(vars_full[,7], n),
                          B8 = rep(vars_full[,8], n))
  
  vals_full <- vals_full %>%
    mutate(d1_full_female = ((B3+B7)*d0 - B2)/(B4+B8),
           d1_full_male = ((B3)*d0 - B2)/(B4))
  
  
  vals2_full <- vals_full %>% 
    group_by(d0) %>% 
    summarize(lower_full_female = quantile(d1_full_female, probs = .025), 
              upper_full_female = quantile(d1_full_female, probs = 0.975),
              median_full_female = median(d1_full_female),
              lower_full_male = quantile(d1_full_male, probs = .025), 
              upper_full_male = quantile(d1_full_male, probs = 0.975),
              median_full_male = median(d1_full_male))
  
  vals2_all <- full_join(vals2_meta, vals2_full, by = "d0") %>%
    full_join(vals2_meta_nocov, by = "d0")
  
  
  return(vals2_all)
}
