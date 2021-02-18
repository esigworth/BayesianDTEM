full_run_inter <- function(balanced_trials = c(TRUE, FALSE),
                     tau2 = 0,
                     doses = doses,
                     nA = 75, nB = 75,
                     nsubjA = 100, nsubjB = 100,
                     minsubj = 50, maxsubj = 200,
                     iters = 10000, thin = 1, burnin = 5000){
  
  dat <- sim_linear_int(balanced_trials = balanced_trials, 
                        tau2 = tau2, doses = doses, nA = nA, nB = nB,
                        nsubjA = nsubjA, nsubjB = nsubjB,
                        minsubj = minsubj, maxsubj = maxsubj)
  

  mod_meta <-"
    model{
    
    #Likelihood
      for( i in 1:n)
        {
          logit_outcome[i]~dnorm(mu[i], 1/(logit_var[i]))
          mu[i]<-b[1]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+
          b[5]*cov1[i] + b[6]*cov2[i] + u[trial[i]]
        }
    
    for (j in 1:nid)
      {
        u[j]~dnorm(0, prec_u)
      }
    #priors 
    for(j in 1:6) { b[j]~dnorm(0, 1E-6)}
    
    prec_u <- pow(sigma, 2)
    sigma ~ dgamma(0.001, 0.001)
    
    
    }
    "
  
  # Put data into a list format for the model
  dat_meta <- list(logit_outcome = dat$meta$logit_outcome, 
                   logit_var = dat$meta$logit_var,
                   drugB = dat$meta$drugB, 
                   drugA_dose = dat$meta$drugA_dose,
                   drugB_dose = dat$meta$drugB_dose, 
                   trial = dat$meta$trial,
                   cov1 = dat$meta$cov1, cov2 = dat$meta$cov2,
                   n = length(unique(dat$meta$trial)), 
                   nid = length(unique((dat$meta$trial))))
  
  
  ### Full data model
  
  mod_full <-"
    model{
    
    #Likelihood
      for( i in 1:n)
        {
          y[i] ~ dbern(q[i])
        logit(q[i]) <- b[1]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+
          b[5]*cov1[i] + b[6]*cov2[i] + b[7]*cov1[i]*drugA_dose[i] + 
          b[8]*cov1[i]*drugB_dose[i] + u[trial[i]]
        }
    
    for (j in 1:nid)
      {
        u[j]~dnorm(0, prec_u)
      }
    #priors 
    for(j in 1:8) { b[j]~dnorm(0, 1E-6)}
    
    prec_u <- pow(sigma, 2)
    sigma ~ dgamma(0.001, 0.001)
    
    
    }
    "
  
  # Put data into a list format for the model
  dat_full <- list(y = dat$full$y,
                   drugB = dat$full$drugB, 
                   drugA_dose = dat$full$drugA_dose,
                   drugB_dose = dat$full$drugB_dose, 
                   trial = dat$full$trial,
                   cov1 = dat$full$cov1, cov2 = dat$full$cov2,
                   n = length(unique(dat$full$id)), 
                   nid = length(unique((dat$full$trial))))
  
  meta_fit <- fit_mod(model = mod_meta, dat = dat_meta, 
                      burnin = burnin, thinning = thin, iterations = iters)
  full_fit <- fit_mod(model = mod_full, dat = dat_full, 
                      burnin = burnin, thinning = thin, iterations = iters)
  
  
  mod_meta_nocov <-"
    model{
    
    #Likelihood
      for( i in 1:n)
        {
          logit_outcome[i]~dnorm(mu[i], 1/(logit_var[i]))
          mu[i]<-b[1]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+u[trial[i]]
        }
    
    for (j in 1:nid)
      {
        u[j]~dnorm(0, prec_u)
      }
    #priors 
    for(j in 1:4) { b[j]~dnorm(0, 1E-6)}
    
    prec_u <- pow(sigma, 2)
    sigma ~ dgamma(0.001, 0.001)
    
    
    }
    "
  
  # Put data into a list format for the model
  dat_meta_nocov <- list(logit_outcome = dat$meta$logit_outcome, 
                         logit_var = dat$meta$logit_var,
                         drugB = dat$meta$drugB, 
                         drugA_dose = dat$meta$drugA_dose,
                         drugB_dose = dat$meta$drugB_dose, 
                         trial = dat$meta$trial,
                         n = length(unique(dat$meta$trial)), 
                         nid = length(unique((dat$meta$trial))))
  
  
  
  
  meta_fit_nocov <- fit_mod(model = mod_meta_nocov, dat = dat_meta_nocov, 
                            burnin = burnin, thinning = thin, iterations = iters)
  
  
  
  meta_summary_nocov <- summary(meta_fit_nocov[,1:5])[[1]]
  meta_summary_median_nocov <- summary(meta_fit_nocov[,1:5])[[2]]
  
  
  meta_summary <- summary(meta_fit[,1:7])[[1]]
  meta_summary_median <- summary(meta_fit[,1:7])[[2]]
  full_summary <- summary(full_fit[,1:9])[[1]]
  full_summary_median <- summary(full_fit[,1:9])[[2]]
  
  meta_ess <- effectiveSize(meta_fit[,1:7])
  full_ess <- effectiveSize(full_fit[,1:9])
  meta_nocov_ess <- effectiveSize(meta_fit_nocov[,1:5])
  
  gr_meta <- as.data.frame(gelman.diag(meta_fit)$psrf[1:7,])
  gr_meta$type = "meta"
  gr_meta$param <- c("b1", "b2", "b3", "b4", "b5", "b6", "prec")
  gr_full <- as.data.frame(gelman.diag(full_fit)$psrf[1:9,])
  gr_full$type = "full"
  gr_full$param <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "prec")
  gr_meta_nocov <- as.data.frame(gelman.diag(meta_fit_nocov)$psrf[1:5,])
  gr_meta_nocov$type = "meta_nocov"
  gr_meta_nocov$param <- c("b1", "b2", "b3", "b4", "prec")
  
  gr_stats <- rbind(gr_meta, gr_full, gr_meta_nocov)
  rownames(gr_stats) <- NULL
  
  samps_meta <- meta_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:6]
  
  samps_full <- full_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:8]
  
  names(samps_meta) <- c("b1", "b2", "b3", "b4", "b5", "b6")
  
  names(samps_full) <- c("b1", "b2", "b3", "b4", 
                         "b5", "b6", "b7", "b8")
  
  samps_meta_nocov <- meta_fit_nocov %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:4]
  
  names(samps_meta_nocov) <- c("b1", "b2", "b3", "b4")
  
  mad_meta_nocov <- samps_meta_nocov %>%
    mutate(b1_mad = abs(median(b1)-b1),
           b2_mad = abs(median(b2)-b2),
           b3_mad = abs(median(b3)-b3),
           b4_mad = abs(median(b4)-b4)) %>%
    summarise(type = "meta_nocov",
              b1_mad = median(b1_mad)*1.4826,
              b2_mad = median(b2_mad)*1.4826,
              b3_mad = median(b3_mad)*1.4826,
              b4_mad = median(b4_mad)*1.4826,
              b5_mad = NA,
              b6_mad = NA,
              b7_mad = NA,
              b8_mad = NA)
  
  mad_meta <- samps_meta %>%
    mutate(b1_mad = abs(median(b1)-b1),
           b2_mad = abs(median(b2)-b2),
           b3_mad = abs(median(b3)-b3),
           b4_mad = abs(median(b4)-b4),
           b5_mad = abs(median(b5)-b5),
           b6_mad = abs(median(b6)-b6)) %>%
    summarise(type = "meta",
              b1_mad = median(b1_mad)*1.4826,
              b2_mad = median(b2_mad)*1.4826,
              b3_mad = median(b3_mad)*1.4826,
              b4_mad = median(b4_mad)*1.4826,
              b5_mad = median(b5_mad)*1.4826,
              b6_mad = median(b6_mad)*1.4826,
              b7_mad = NA,
              b8_mad = NA)
  
  mad_full <- samps_full %>%
    mutate(b1_mad = abs(median(b1)-b1),
           b2_mad = abs(median(b2)-b2),
           b3_mad = abs(median(b3)-b3),
           b4_mad = abs(median(b4)-b4),
           b5_mad = abs(median(b5)-b5),
           b6_mad = abs(median(b6)-b6),
           b7_mad = abs(median(b7)-b7),
           b8_mad = abs(median(b8)-b8)) %>%
    summarise(type = "full",
              b1_mad = median(b1_mad)*1.4826,
              b2_mad = median(b2_mad)*1.4826,
              b3_mad = median(b3_mad)*1.4826,
              b4_mad = median(b4_mad)*1.4826,
              b5_mad = median(b5_mad)*1.4826,
              b6_mad = median(b6_mad)*1.4826,
              b7_mad = median(b7_mad)*1.4826,
              b8_mad = median(b8_mad)*1.4826)
  
  ratios <- round(mad_meta[,2:9]/mad_full[,2:9], 3)
  names(ratios) <- names(mad_full)[2:9]
  ratios2 <- c(round(mad_meta_nocov[,2:9]/mad_full[,2:9], 3))
  names(ratios2)<- names(mad_full)[2:9]
  rats <- rbind(ratios, ratios2)
  mad <- as.data.frame(rbind(mad_meta[,2:9], mad_meta_nocov[,2:9], mad_full[,2:9], rats))
  mad$type <- c("meta", "meta_nocov", "full", "ratio_meta_full", "ratio_meta_nocov_full")
  
  equiv <- equiv_linear_inter(meta_fit, full_fit, meta_fit_nocov,
                             balanced = balanced_trials)
  
  result <- list("meta_summary" = meta_summary, 
                 "full_summary" = full_summary, 
                 "meta_nocov_summary" = meta_summary_nocov,
                 "meta_summary_median" = meta_summary_median,
                 "full_summary_median" = full_summary_median,
                 "meta_nocov_median" = meta_summary_median_nocov,
                 "equiv" = equiv, 
                 "tau2" = tau2,
                 "mad" = mad)
    

  return(result)
}





fit_one_type_inter <- function(tau2, i, doses = doses, 
                             nA = 75, nB = 75, 
                             nsubjA = 100, nsubjB = 100,
                             minsubj = 50, maxsubj = 200,
                             iters = iters, thin = thin, burnin = burnin, 
                             fit_type = c("lbc", "luc"), seed){
  
  
  if(fit_type == "lbc"){
    set.seed(seed)
    
    run <- full_run_inter(balanced_trials = TRUE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  if(fit_type == "luc"){
    set.seed(seed)
    
    run <- full_run_inter(balanced_trials = FALSE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  
  val <- list("run" = run, "fit_type" = fit_type, "seed" = seed)
  
  return(val)
}
