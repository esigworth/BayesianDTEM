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

################################################################################

full_run_varied <- function(balanced_trials = c(TRUE, FALSE),
                     tau2 = 0,
                     doses = doses,
                     nA = 75, nB = 75,
                     nsubjA = 100, nsubjB = 100,
                     minsubj = 50, maxsubj = 200,
                     iters = 10000, thin = 1, burnin = 5000){
  
  dat <- sim_varied(balanced_trials = balanced_trials, 
                    tau2 = tau2, doses = doses, nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj)
  
  mod_studylevel <-"
    model{
    
    #Likelihood
      for( i in 1:n)
        {
          logit_outcome[i]~dnorm(mu[i], 1/(logit_var[i]))
          mu[i]<-b[1]*drugA[i]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+u[trial[i]]
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
  dat_studylevel <- list(logit_outcome = dat$studylevel$logit_outcome, 
                   logit_var = dat$studylevel$logit_var,
                   drugA=dat$studylevel$drugA, drugB = dat$studylevel$drugB, 
                   drugA_dose = dat$studylevel$drugA_dose,
                   drugB_dose = dat$studylevel$drugB_dose, 
                   trial = dat$studylevel$trial,
                   n = length(unique(dat$studylevel$trial)), 
                   nid = length(unique((dat$studylevel$trial))))
  
  
  ### ipd data model
  
  mod_ipd <-"
    model{
    
    #Likelihood
      for( i in 1:n)
        {
          y[i] ~ dbern(q[i])
        logit(q[i]) <- b[1]*drugA[i]+b[2]*drugB[i]+b[3]*drugA_dose[i]+b[4]*drugB_dose[i]+u[trial[i]]
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
  dat_ipd <- list(y = dat$ipd$y,
                   drugA=dat$ipd$drugA, drugB = dat$ipd$drugB, 
                   drugA_dose = dat$ipd$drugA_dose,
                   drugB_dose = dat$ipd$drugB_dose, 
                   trial = dat$ipd$trial,
                   n = length(unique(dat$ipd$id)), 
                   nid = length(unique((dat$ipd$trial))))
  
  studylevel_fit <- fit_mod(model = mod_studylevel, dat = dat_studylevel, burnin = burnin, thinning = thin, iterations = iters)
  ipd_fit <- fit_mod(model = mod_ipd, dat = dat_ipd, burnin = burnin, thinning = thin, iterations = iters)
  
  studylevel_summary <- summary(studylevel_fit[,1:5])[[1]]
  studylevel_summary_median <- summary(studylevel_fit[,1:5])[[2]]
  ipd_summary <- summary(ipd_fit[,1:5])[[1]]
  ipd_summary_median <- summary(ipd_fit[,1:5])[[2]]
  
  studylevel_ess <- c(round(effectiveSize(studylevel_fit[,1:5]), 3))
  ipd_ess <- c(round(effectiveSize(ipd_fit[,1:5]), 3))
  
  ess <- rbind(studylevel_ess, ipd_ess)
  ess <- as.data.frame(ess)
  ess <- data.frame(type = c("studylevel_ess", "ipd_ess"), ess)
  colnames(ess) <- c("type", "b1_ess", "b2_ess", "b3_ess", "b4_ess", "prec_ess")
  rownames(ess) <- NULL
  
  gr_studylevel <- as.data.frame(gelman.diag(studylevel_fit)$psrf[1:5,])
  gr_studylevel$type = "studylevel"
  gr_studylevel$param <- c("b1", "b2", "b3", "b4", "prec")
  gr_ipd <- as.data.frame(gelman.diag(ipd_fit)$psrf[1:5,])
  gr_ipd$type = "ipd"
  gr_ipd$param <- c("b1", "b2", "b3", "b4", "prec")
  
  gr_stats <- rbind(gr_studylevel, gr_ipd)
  rownames(gr_stats) <- NULL
  
  samps_studylevel <- studylevel_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:4]
  
  
  
  samps_ipd <- ipd_fit %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    .[,1:4]
  
  names(samps_studylevel) <- names(samps_ipd) <- c("b1", "b2", "b3", "b4")
  
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
  
  ratios <- c("ratio", round(mad_studylevel[,2:5]/mad_ipd[,2:5], 3))
  names(ratios) <- names(mad_ipd)
  mad <- rbind(mad_studylevel, mad_ipd, ratios)
  
  equiv <- equiv_varied(studylevel_fit, ipd_fit, 
                             balanced = balanced_trials)
  
  result <- list("studylevel_summary" = studylevel_summary, 
                 "ipd_summary" = ipd_summary, 
                 "studylevel_summary_median" = studylevel_summary_median,
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



fit_one_varied <- function(tau2, i, doses = doses, 
                             nA = 75, nB = 75, 
                             nsubjA = 100, nsubjB = 100,
                             minsubj = 50, maxsubj = 200,
                             iters = iters, thin = thin, burnin = burnin, 
                             fit_type = c("lbn", "lun"), seed){
  
  
  
  if(fit_type == "lbn"){
    set.seed(seed)
    
    run <- full_run_varied(balanced_trials = TRUE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  
  
  if(fit_type == "lun"){
    set.seed(seed)
    
    run <- full_run_varied(balanced_trials = FALSE, 
                    tau2 = tau2, doses = doses,
                    nA = nA, nB = nB,
                    nsubjA = nsubjA, nsubjB = nsubjB,
                    minsubj = minsubj, maxsubj = maxsubj,
                    iters = iters, thin = thin, burnin = burnin)
  }
  val <- list("run" = run, "fit_type" = fit_type, "seed" = seed)
  
  return(val)
}
