################################################################################
# Function to fit models and do calculations on posterior samples
## nA, nB : number of studies in drugs A and B
## nsubjA, nsubjB : if studies balanced, number of subjects in each
## minsubj, maxsubj : if studies unbalanced, range of study sizes
## with_covars, balanced_trials : settings for sim condition
## tau2 : level of between-study variability
## doses : pre-specified dose vector, one value for each study
## iters : number of samples per chain
## thin : thinning interval for sampling
## burnin : number of burn-in samples to use per chain
#
# Notes
## 
################################################################################


full_run <- function(with_covars = c(TRUE, FALSE),
                     balanced_trials = c(TRUE, FALSE),
                     tau2 = 0,
                     doses = doses,
                     nA = 75, nB = 75,
                     nsubjA = 100, nsubjB = 100,
                     minsubj = 50, maxsubj = 200,
                     iters = 10000, thin = 1, burnin = 5000){
  
  # Generate dataset
  dat <- sim_linear(with_covars = with_covars, 
                    balanced_trials = balanced_trials, 
                    tau2 = tau2, doses = doses, nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj)
  
  # Model fitting without subject-level covariates
  if(!(with_covars)){
    
    # SL meta specification
    mod_studylevel <-"
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
    # SL meta data list
    dat_studylevel <- list(logit_outcome = dat$studylevel$logit_outcome, 
                     logit_var = dat$studylevel$logit_var,
                     drugB = dat$studylevel$drugB, 
                     drugA_dose = dat$studylevel$drugA_dose,
                     drugB_dose = dat$studylevel$drugB_dose, 
                     trial = dat$studylevel$trial,
                     n = length(unique(dat$studylevel$trial)), 
                     nid = length(unique((dat$studylevel$trial))))
    
    # IPD meta model specification
    mod_ipd <-"
      model{
      
      #Likelihood
        for( i in 1:n)
          {
            y[i] ~ dbern(q[i])
          logit(q[i]) <- b[1]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+u[trial[i]]
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
    
    # IPD meta data list
    dat_ipd <- list(y = dat$ipd$y,
                     drugB = dat$ipd$drugB, 
                     drugA_dose = dat$ipd$drugA_dose,
                     drugB_dose = dat$ipd$drugB_dose, 
                     trial = dat$ipd$trial,
                     n = length(unique(dat$ipd$id)), 
                     nid = length(unique((dat$ipd$trial))))
    
    # Fit models using JAGS
    studylevel_fit <- fit_mod(model = mod_studylevel, dat = dat_studylevel, 
                        burnin = burnin, thinning = thin, iterations = iters)
    ipd_fit <- fit_mod(model = mod_ipd, dat = dat_ipd, 
                        burnin = burnin, thinning = thin, iterations = iters)
    
    # Pull mean and median variable summaries
    studylevel_summary <- summary(studylevel_fit[,1:5])[[1]]
    studylevel_summary_median <- summary(studylevel_fit[,1:5])[[2]]
    ipd_summary <- summary(ipd_fit[,1:5])[[1]]
    ipd_summary_median <- summary(ipd_fit[,1:5])[[2]]
    
    # Calculate effective sample size (ESS) for betas and tau
    studylevel_ess <- c(round(effectiveSize(studylevel_fit[,1:5]), 3))
    ipd_ess <- c(round(effectiveSize(ipd_fit[,1:5]), 3))
    
    # Organize ESS into dataframe
    ess <- rbind(studylevel_ess, ipd_ess)
    ess <- as.data.frame(ess)
    ess <- data.frame(type = c("studylevel_ess", "ipd_ess"), ess)
    colnames(ess) <- c("type", "b1_ess", "b2_ess", 
                       "b3_ess", "b4_ess", "prec_ess")
    rownames(ess) <- NULL
    
    # Calculate Gelman-Rubin statistic on betas and tau, put in dataframe
    gr_studylevel <- as.data.frame(gelman.diag(studylevel_fit)$psrf[1:5,])
    gr_studylevel$type = "studylevel"
    gr_studylevel$param <- c("b1", "b2", "b3", "b4", "prec")
    gr_ipd <- as.data.frame(gelman.diag(ipd_fit)$psrf[1:5,])
    gr_ipd$type = "ipd"
    gr_ipd$param <- c("b1", "b2", "b3", "b4", "prec")
    gr_stats <- rbind(gr_studylevel, gr_ipd)
    rownames(gr_stats) <- NULL
    
    # Combine beta samples from all chains into one dataframe
    samps_studylevel <- studylevel_fit %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      .[,1:4]
    
    samps_ipd <- ipd_fit %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      .[,1:4]
    
    names(samps_studylevel) <- names(samps_ipd) <- c("b1", "b2", "b3", "b4")
    
    # Calculate median absolute deviation (MAD) for betas
    mad_studylevel <- samps_studylevel %>%
      mutate(b1_mad = abs(median(b1)-b1),
             b2_mad = abs(median(b2)-b2),
             b3_mad = abs(median(b3)-b3),
             b4_mad = abs(median(b4)-b4)) %>%
      summarise(type = "studylevel",
                b1_mad = round(median(b1_mad)*1.4826, 3),
                b2_mad = round(median(b2_mad)*1.4826, 3),
                b3_mad = round(median(b3_mad)*1.4826, 3),
                b4_mad = round(median(b4_mad)*1.4826, 3))
    
    mad_ipd <- samps_ipd %>%
      mutate(b1_mad = abs(median(b1)-b1),
             b2_mad = abs(median(b2)-b2),
             b3_mad = abs(median(b3)-b3),
             b4_mad = abs(median(b4)-b4)) %>%
      summarise(type = "ipd", 
                b1_mad = round(median(b1_mad)*1.4826, 3),
                b2_mad = round(median(b2_mad)*1.4826, 3),
                b3_mad = round(median(b3_mad)*1.4826, 3),
                b4_mad = round(median(b4_mad)*1.4826, 3))
    
    # Calculate ratio of MAD between SL and IPD models
    ratios <- c("ratio", round(mad_studylevel[,2:5]/mad_ipd[,2:5], 3))
    names(ratios) <- names(mad_ipd)
    mad <- rbind(mad_studylevel, mad_ipd, ratios)
    
    # Calculate equivalence relationship
    equiv <- equiv_linear_both(studylevel_fit, ipd_fit, 
                               balanced = balanced_trials,
                               covars = with_covars)
  }
  
  # Model fitting with subject-level covariates
  else {
    # SL-C model statement
    mod_studylevel <-"
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
    
    # SL-C model data
    dat_studylevel <- list(logit_outcome = dat$studylevel$logit_outcome, 
                     logit_var = dat$studylevel$logit_var,
                     drugB = dat$studylevel$drugB, 
                     drugA_dose = dat$studylevel$drugA_dose,
                     drugB_dose = dat$studylevel$drugB_dose, 
                     trial = dat$studylevel$trial,
                     cov1 = dat$studylevel$cov1, cov2 = dat$studylevel$cov2,
                     n = length(unique(dat$studylevel$trial)), 
                     nid = length(unique((dat$studylevel$trial))))
    
    # SL-NC model statement
    mod_studylevel_nocov <-"
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
    
    # SL-NC model data
    dat_studylevel_nocov <- list(logit_outcome = dat$studylevel$logit_outcome, 
                           logit_var = dat$studylevel$logit_var,
                           drugB = dat$studylevel$drugB, 
                           drugA_dose = dat$studylevel$drugA_dose,
                           drugB_dose = dat$studylevel$drugB_dose, 
                           trial = dat$studylevel$trial,
                           n = length(unique(dat$studylevel$trial)), 
                           nid = length(unique((dat$studylevel$trial))))
    
    
    # IPD meta model statement
    
    mod_ipd <-"
      model{
      
      #Likelihood
        for( i in 1:n)
          {
            y[i] ~ dbern(q[i])
          logit(q[i]) <- b[1]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+
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
    
    # IPD model data
    dat_ipd <- list(y = dat$ipd$y,
                     drugB = dat$ipd$drugB, 
                     drugA_dose = dat$ipd$drugA_dose,
                     drugB_dose = dat$ipd$drugB_dose, 
                     trial = dat$ipd$trial,
                     cov1 = dat$ipd$cov1, cov2 = dat$ipd$cov2,
                     n = length(unique(dat$ipd$id)), 
                     nid = length(unique((dat$ipd$trial))))
    
    # All model fits
    studylevel_fit <- fit_mod(model = mod_studylevel, dat = dat_studylevel, 
                        burnin = burnin, thinning = thin, iterations = iters)
    studylevel_fit_nocov <- fit_mod(model = mod_studylevel_nocov, dat = dat_studylevel_nocov, 
                              burnin = burnin, thinning = thin, iterations = iters)
    ipd_fit <- fit_mod(model = mod_ipd, dat = dat_ipd, 
                        burnin = burnin, thinning = thin, iterations = iters)
    
    # Pull model summaries
    studylevel_summary <- summary(studylevel_fit[,1:7])[[1]]
    studylevel_summary_median <- summary(studylevel_fit[,1:7])[[2]]
    studylevel_summary_nocov <- summary(studylevel_fit_nocov[,1:5])[[1]]
    studylevel_summary_median_nocov <- summary(studylevel_fit_nocov[,1:5])[[2]]
    ipd_summary <- summary(ipd_fit[,1:7])[[1]]
    ipd_summary_median <- summary(ipd_fit[,1:7])[[2]]
    
    # Calculate ESS for all models
    studylevel_ess <- effectiveSize(studylevel_fit[,1:7])
    studylevel_nocov_ess <- effectiveSize(studylevel_fit[,1:5])
    ipd_ess <- effectiveSize(ipd_fit[,1:7])
    ess <- rbind(studylevel_ess, ipd_ess)
    ess <- as.data.frame(ess)
    ess <- data.frame(type = c("studylevel_ess", "ipd_ess"), ess)
    colnames(ess) <- c("type", "b1_ess", "b2_ess", "b3_ess", 
                       "b4_ess", "b5_ess", "b6_ess", "prec_ess")
    rownames(ess) <- NULL
    
    # Calculate Gelman-Rubin for all models
    gr_studylevel <- as.data.frame(gelman.diag(studylevel_fit)$psrf[1:7,])
    gr_studylevel$type = "studylevel"
    gr_studylevel$param <- c("b1", "b2", "b3", "b4", "b5", "b6", "prec")
    gr_studylevel_nocov <- as.data.frame(gelman.diag(studylevel_fit)$psrf[1:5,])
    gr_studylevel_nocov$type = "studylevel_nocov"
    gr_studylevel_nocov$param <- c("b1", "b2", "b3", "b4", "prec")
    gr_ipd <- as.data.frame(gelman.diag(ipd_fit)$psrf[1:7,])
    gr_ipd$type = "ipd"
    gr_ipd$param <- c("b1", "b2", "b3", "b4", "b5", "b6", "prec")
    gr_stats <- rbind(gr_studylevel, gr_studylevel_nocov, gr_ipd)
    rownames(gr_stats) <- NULL
    
    # Pull and bind posterior samples for all fits
    samps_studylevel <- studylevel_fit %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      .[,1:6]
    
    samps_studylevel_nocov <- studylevel_fit_nocov %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      .[,1:4]
    
    samps_ipd <- ipd_fit %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      .[,1:6]
    
    names(samps_studylevel) <- names(samps_ipd) <- c("b1", "b2", "b3", 
                                                "b4", "b5", "b6")
    names(samps_studylevel_nocov) <- c("b1", "b2", "b3", "b4")
    
    # Calculate median absolute deviation for all parameters
    mad_studylevel <- samps_studylevel %>%
      mutate(b1_mad = abs(median(b1)-b1),
             b2_mad = abs(median(b2)-b2),
             b3_mad = abs(median(b3)-b3),
             b4_mad = abs(median(b4)-b4),
             b5_mad = abs(median(b5)-b5),
             b6_mad = abs(median(b6)-b6)) %>%
      summarise(type = "studylevel",
                b1_mad = median(b1_mad)*1.4826,
                b2_mad = median(b2_mad)*1.4826,
                b3_mad = median(b3_mad)*1.4826,
                b4_mad = median(b4_mad)*1.4826,
                b5_mad = median(b5_mad)*1.4826,
                b6_mad = median(b6_mad)*1.4826)
    
    mad_studylevel_nocov <- samps_studylevel_nocov %>%
      mutate(b1_mad = abs(median(b1)-b1),
             b2_mad = abs(median(b2)-b2),
             b3_mad = abs(median(b3)-b3),
             b4_mad = abs(median(b4)-b4)) %>%
      summarise(type = "studylevel_nocov",
                b1_mad = median(b1_mad)*1.4826,
                b2_mad = median(b2_mad)*1.4826,
                b3_mad = median(b3_mad)*1.4826,
                b4_mad = median(b4_mad)*1.4826,
                b5_mad = NA,
                b6_mad = NA)
    
    mad_ipd <- samps_ipd %>%
      mutate(b1_mad = abs(median(b1)-b1),
             b2_mad = abs(median(b2)-b2),
             b3_mad = abs(median(b3)-b3),
             b4_mad = abs(median(b4)-b4),
             b5_mad = abs(median(b5)-b5),
             b6_mad = abs(median(b6)-b6)) %>%
      summarise(type = "ipd",
                b1_mad = median(b1_mad)*1.4826,
                b2_mad = median(b2_mad)*1.4826,
                b3_mad = median(b3_mad)*1.4826,
                b4_mad = median(b4_mad)*1.4826,
                b5_mad = median(b5_mad)*1.4826,
                b6_mad = median(b6_mad)*1.4826)
    
    # MAD ratios comparing SL-C to IPD and SL-NC to IPD
    ratios <- c("ratio_studylevel_ipd", round(mad_studylevel[,2:7]/mad_ipd[,2:7], 3))
    ratios2 <- c("ratio_studylevel_nocov_ipd", round(mad_studylevel_nocov[,2:5]/mad_ipd[,2:5], 3), NA, NA)
    names(ratios) <- names(ratios2) <- names(mad_ipd)
    mad <- rbind(mad_studylevel, mad_studylevel_nocov, mad_ipd, ratios, ratios2)
    
    # Calculate equivalence relationship
    equiv <- equiv_linear_both(studylevel_fit, ipd_fit, 
                               balanced = balanced_trials,
                               covars = with_covars,
                               studylevel_nocov_samples = studylevel_fit_nocov)
    
  }
  
  
  
  result <- list("studylevel_summary" = studylevel_summary, 
                 "studylevel_nocov_summary" = studylevel_summary_nocov, 
                 "ipd_summary" = ipd_summary, 
                 "studylevel_summary_median" = studylevel_summary_median,
                 "studylevel_summary_median_nocov" = studylevel_summary_median_nocov,
                 "ipd_summary_median" = ipd_summary_median,
                 "ess" = ess,
                 "gr_stats" = gr_stats,
                 "equiv" = equiv, 
                 "tau2" = tau2,
                 "mad" = mad)
  
  return(result)
}


################################################################################
# Generic wrapper to repeatedly fit models in specified settings
## nA, nB : number of studies in drugs A and B
## nsubjA, nsubjB : if studies balanced, number of subjects in each
## minsubj, maxsubj : if studies unbalanced, range of study sizes
## with_covars, balanced_trials : settings for sim condition
## tau2 : level of between-study variability
## doses : pre-specified dose vector, one value for each study
## iters : number of samples per chain
## thin : thinning interval for sampling
## burnin : number of burn-in samples to use per chain
## fit_type : "l" = linear, "b" = balanced, "u" = unbalanced, 
##            "n" = no extra covs, "c" = extra covs
## seed : specify seed for reproducibility
#
# Notes
## 
################################################################################



fit_one_type_tau <- function(tau2, i, doses = doses, 
                             nA = 75, nB = 75, 
                             nsubjA = 100, nsubjB = 100,
                             minsubj = 50, maxsubj = 200,
                             iters = iters, thin = thin, burnin = burnin, 
                             fit_type = c("lbn", "lun", "lbc", "luc"), seed){
  
  
  
  if(fit_type == "lbn"){
    set.seed(seed)
    
    run <- full_run(with_covars = FALSE, balanced_trials = TRUE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  
  if(fit_type == "lun"){
    set.seed(seed)
    
    run <- full_run(with_covars = FALSE, balanced_trials = FALSE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  if(fit_type == "lbc"){
    set.seed(seed)
    
    run <- full_run(with_covars = TRUE, balanced_trials = TRUE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  if(fit_type == "luc"){
    set.seed(seed)
    
    run <- full_run(with_covars = TRUE, balanced_trials = FALSE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  
  val <- list("run" = run, "fit_type" = fit_type, "seed" = seed)
  
  return(val)
}
