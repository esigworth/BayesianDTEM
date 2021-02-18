################################################################################
# Function to fit models and do calculations on posterior samples
## nA, nB : number of studies in drugs A and B
## nsubjA, nsubjB : if studies balanced, number of subjects in each
## minsubj, maxsubj : if studies unbalanced, range of study sizes
## balanced_trials : settings for sim condition
## tau2 : level of between-study variability
## doses : pre-specified dose vector, one value for each study
## iters : number of samples per chain
## thin : thinning interval for sampling
## burnin : number of burn-in samples to use per chain

################################################################################

full_run_inter <- function(balanced_trials = c(TRUE, FALSE),
                           tau2 = 0,
                           doses = doses,
                           nA = 75, nB = 75,
                           nsubjA = 100, nsubjB = 100,
                           minsubj = 50, maxsubj = 200,
                           iters = 10000, thin = 1, burnin = 5000){
  
  dat <- sim_linear_slope(balanced_trials = balanced_trials, 
                        tau2 = tau2, doses = doses, nA = nA, nB = nB,
                        nsubjA = nsubjA, nsubjB = nsubjB,
                        minsubj = minsubj, maxsubj = maxsubj)
  
  
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
  
  # Put data into a list format for the model
  dat_studylevel <- list(logit_outcome = dat$studylevel$logit_outcome, 
                   logit_var = dat$studylevel$logit_var,
                   drugB = dat$studylevel$drugB, 
                   drugA_dose = dat$studylevel$drugA_dose,
                   drugB_dose = dat$studylevel$drugB_dose, 
                   trial = dat$studylevel$trial,
                   cov1 = dat$studylevel$cov1, cov2 = dat$studylevel$cov2,
                   n = length(unique(dat$studylevel$trial)), 
                   nid = length(unique((dat$studylevel$trial))))
  
  
  ### Full data model
  
  mod_ipd <-"
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
  dat_ipd <- list(y = dat$ipd$y,
                   drugB = dat$ipd$drugB, 
                   drugA_dose = dat$ipd$drugA_dose,
                   drugB_dose = dat$ipd$drugB_dose, 
                   trial = dat$ipd$trial,
                   cov1 = dat$ipd$cov1, cov2 = dat$ipd$cov2,
                   n = length(unique(dat$ipd$id)), 
                   nid = length(unique((dat$ipd$trial))))
  
  studylevel_fit <- fit_mod(model = mod_studylevel, dat = dat_studylevel, 
                      burnin = burnin, thinning = thin, iterations = iters)
  ipd_fit <- fit_mod(model = mod_ipd, dat = dat_ipd, 
                      burnin = burnin, thinning = thin, iterations = iters)
  
  
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
  
  # Put data into a list format for the model
  dat_studylevel_nocov <- list(logit_outcome = dat$studylevel$logit_outcome, 
                         logit_var = dat$studylevel$logit_var,
                         drugB = dat$studylevel$drugB, 
                         drugA_dose = dat$studylevel$drugA_dose,
                         drugB_dose = dat$studylevel$drugB_dose, 
                         trial = dat$studylevel$trial,
                         n = length(unique(dat$studylevel$trial)), 
                         nid = length(unique((dat$studylevel$trial))))
  
  
  
  
  studylevel_fit_nocov <- fit_mod(model = mod_studylevel_nocov, dat = dat_studylevel_nocov, 
                            burnin = burnin, thinning = thin, iterations = iters)
  
  
  
  studylevel_summary_nocov <- summary(studylevel_fit_nocov[,1:5])[[1]]
  studylevel_summary_median_nocov <- summary(studylevel_fit_nocov[,1:5])[[2]]
  
  
  studylevel_summary <- summary(studylevel_fit[,1:7])[[1]]
  studylevel_summary_median <- summary(studylevel_fit[,1:7])[[2]]
  ipd_summary <- summary(ipd_fit[,1:9])[[1]]
  ipd_summary_median <- summary(ipd_fit[,1:9])[[2]]

  
  samps_studylevel <- studylevel_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:6]
  
  samps_ipd <- ipd_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:8]
  
  names(samps_studylevel) <- c("b1", "b2", "b3", "b4", "b5", "b6")
  
  names(samps_ipd) <- c("b1", "b2", "b3", "b4", 
                         "b5", "b6", "b7", "b8")
  
  samps_studylevel_nocov <- studylevel_fit_nocov %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:4]
  
  names(samps_studylevel_nocov) <- c("b1", "b2", "b3", "b4")
  
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
              b6_mad = NA,
              b7_mad = NA,
              b8_mad = NA)
  
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
              b6_mad = median(b6_mad)*1.4826,
              b7_mad = NA,
              b8_mad = NA)
  
  mad_ipd <- samps_ipd %>%
    mutate(b1_mad = abs(median(b1)-b1),
           b2_mad = abs(median(b2)-b2),
           b3_mad = abs(median(b3)-b3),
           b4_mad = abs(median(b4)-b4),
           b5_mad = abs(median(b5)-b5),
           b6_mad = abs(median(b6)-b6),
           b7_mad = abs(median(b7)-b7),
           b8_mad = abs(median(b8)-b8)) %>%
    summarise(type = "ipd",
              b1_mad = median(b1_mad)*1.4826,
              b2_mad = median(b2_mad)*1.4826,
              b3_mad = median(b3_mad)*1.4826,
              b4_mad = median(b4_mad)*1.4826,
              b5_mad = median(b5_mad)*1.4826,
              b6_mad = median(b6_mad)*1.4826,
              b7_mad = median(b7_mad)*1.4826,
              b8_mad = median(b8_mad)*1.4826)
  
  ratios <- round(mad_studylevel[,2:9]/mad_ipd[,2:9], 3)
  names(ratios) <- names(mad_ipd)[2:9]
  ratios2 <- c(round(mad_studylevel_nocov[,2:9]/mad_ipd[,2:9], 3))
  names(ratios2)<- names(mad_ipd)[2:9]
  rats <- rbind(ratios, ratios2)
  mad <- as.data.frame(rbind(mad_studylevel[,2:9], mad_studylevel_nocov[,2:9], mad_ipd[,2:9], rats))
  mad$type <- c("studylevel", "studylevel_nocov", "ipd", "ratio_studylevel_ipd", "ratio_studylevel_nocov_ipd")
  
  equiv <- equiv_linear_inter(studylevel_fit, ipd_fit, studylevel_fit_nocov,
                              balanced = balanced_trials)
  
  result <- list("studylevel_summary" = studylevel_summary, 
                 "ipd_summary" = ipd_summary, 
                 "studylevel_nocov_summary" = studylevel_summary_nocov,
                 "studylevel_summary_median" = studylevel_summary_median,
                 "ipd_summary_median" = ipd_summary_median,
                 "studylevel_nocov_median" = studylevel_summary_median_nocov,
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
## balanced_trials : settings for sim condition
## tau2 : level of between-study variability
## doses : pre-specified dose vector, one value for each study
## iters : number of samples per chain
## thin : thinning interval for sampling
## burnin : number of burn-in samples to use per chain
## fit_type : "l" = linear, "b" = balanced, "u" = unbalanced, "c" = extra covs
## seed : specify seed for reproducibility
#
# Notes
## 
################################################################################


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
