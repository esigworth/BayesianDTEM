library(bridgesampling)
library(R2jags)
library(Hmisc)
library(rms)
library(readr)
library(tidyr)
library(dplyr)
library(magrittr)
library(knitr)
library(BayesianTools)
library(matrixStats)



logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  1/(1+exp(-x))
}

################################################################################
# Data generation
## nA, nB : number of studies in drugs A and B
## nsubjA, nsubjB : if studies balanced, number of subjects in each
## minsubj, maxsubj : if studies unbalanced, range of study sizes
## betas : coefficients for study-level covariates
## with_covars, balanced_trials : settings for sim condition
## ind_covars : coefficients for subject-level covariates
## tau2 : level of between-study variability
## doses : pre-specified dose vector, one value for each study
#
#
#
# Notes
## Code assumes 75 studies in each drug
## doses vector is a random sample of 150 from a Unif(-1, 1) distribution
################################################################################


sim_linear <- function(nA = 75, nB = 75, 
                       nsubjA = 100, nsubjB = 100, 
                       minsubj = 50,
                       maxsubj = 200, 
                       betas = c(-.6, -.2, .5, .9), 
                       with_covars = TRUE,
                       balanced_trials = TRUE,
                       ind_covars = c(.2, .5),
                       tau2 = 0,
                       doses){
  
  dat_valid <- FALSE
  
  # Split doses between studies
  doses_A <- doses[1:75]
  doses_B <- doses[76:150]
  
  # Generate data without additional covariates
  if(with_covars == FALSE){
    drugA <- data.frame()
    drugB <- data.frame()
    
    while(dat_valid == FALSE){
      
      # Generate all data for studies in drug A
      for(i in 1:nA){
        # Random intercept for trial i
        rand_int <- rnorm(1, 0, sqrt(tau2))
        # Set number of subjects in study
        n_subj <- ifelse(balanced_trials, nsubjA, 
                         sample(seq(minsubj, maxsubj, by = 10), 1))
        
        # Pull specific dose for study i
        mnA <- doses_A[i]
        # Generate log-odds per subject
        zA <- betas[1] + betas[3]*mnA + rand_int
        # Transform to probability
        prA <- expit(zA)
        # Sample outcome based on probability
        yA <- rbinom(n_subj, 1, prA)
        # Put all subject-level data in dataframe
        dA <- data.frame(trial = i, 
                         dose = mnA, 
                         z = zA, 
                         pr = prA, 
                         y = yA, 
                         drug = "A",
                         rand_int = rand_int)
        # Add all rows from trial i to overall drugA dataframe
        drugA <- rbind(drugA, dA)
      }
      # Repeat process for drug B
      for(i in 1:nB){
        rand_int <- rnorm(1, 0, sqrt(tau2))
        n_subj <- ifelse(balanced_trials, nsubjB, 
                         sample(seq(minsubj, maxsubj, by = 10), 1))
        
        mnB <- doses_B[i]
        
        zB <- betas[1] + betas[2] + betas[4]*mnB + rand_int
        prB <- expit(zB)
        yB <- rbinom(n_subj, 1, prB)
        dB <- data.frame(trial = i+nA, 
                         dose = mnB, 
                         z = zB, 
                         pr = prB, 
                         y = yB, 
                         drug = "B", rand_int = rand_int)
        drugB <- rbind(drugB, dB)
      }
      
      # Combine data from all trials
      all_trials <- rbind(drugA, drugB)
      all_trials$id <- 1:nrow(all_trials)
      
      # Create indicator variables
      all_trials <- all_trials %>%
        mutate(drugA = drug == "A",
               drugB = drug == "B",
               drugA_dose = (drug == "A")*dose,
               drugB_dose = (drug == "B")*dose)
      
      # Summarize data to the study level and calculate model variables
      sim_dat <- all_trials %>%
        group_by(drug, trial) %>%
        summarise(median_dose = median(dose), 
                  event_rate = mean(y),
                  size = n(),
                  n_obs = sum(y)) %>%
        mutate(logit_outcome = logit(event_rate),
               logit_var = (1/n_obs) + (1/(size - n_obs)))
      
      # Redefine drug and dose variables for median dose
      sim_dat <- sim_dat %>%
        mutate(drugA = drug == "A",
               drugB = drug == "B",
               drugA_dose = (drug == "A")*median_dose,
               drugB_dose = (drug == "B")*median_dose)
      
      # Confirm event rate != 0% or 100% for all studies, if it was then regenerate dataset
      dat_valid <- ifelse(min(sim_dat$event_rate) > 0 & max(sim_dat$event_rate) < 1, TRUE, FALSE)
      
    }
    
  }
  
  # Generate data with additional subject-level covariates
  else {
    drugA <- data.frame()
    drugB <- data.frame()
    
    # Confirm number of desired subject-level covariates
    n_indcovars <- length(ind_covars)
    
    while(dat_valid == FALSE){
      
      # Generate data for studies in drug A
      for(i in 1:nA){
        # Generate random intercept
        rand_int <- rnorm(1, 0, sqrt(tau2))
        # Set number of subjects based on sim setting
        n_subj <- ifelse(balanced_trials, nsubjA, 
                         sample(seq(minsubj, maxsubj, by = 10), 1))
        # Pull dose
        mnA <- doses_A[i]
        # Generate binary covariate at subject-level
        cov1 <- rbinom(n_subj, 1, runif(1, 0.2, 0.5))
        # Generate continuous covariate at subject-level
        cov2 <- rnorm(n_subj, 0, 1)
        # Calculate log-odds for each subject
        zA <- betas[1] + betas[3]*mnA + ind_covars[1]*cov1 + ind_covars[2]*cov2 + rand_int
        # Get probability of outcome
        prA <- expit(zA)
        # Draw outcome for each subject
        yA <- rbinom(n_subj, 1, prA)
        # Put all subject data into a dataframe
        dA <- data.frame(trial = i, 
                         dose = mnA, 
                         cov1 = cov1,
                         cov2 = cov2,
                         z = zA, 
                         pr = prA, 
                         y = yA, 
                         drug = "A",
                         rand_int = rand_int)
        # Add study data to full dataframe
        drugA <- rbind(drugA, dA)
      }
      # Repeat for studies in drug B
      for(i in 1:nB){
        rand_int <- rnorm(1, 0, sqrt(tau2))
        n_subj <- ifelse(balanced_trials, nsubjB, 
                         sample(seq(minsubj, maxsubj, by = 10), 1))
        mnB <- doses_B[i]
        cov1 <- rbinom(n_subj, 1, runif(1, 0.2, 0.5))
        cov2 <- rnorm(n_subj, 0, 1)
        zB <- betas[1] + betas[2] + betas[4]*mnB + ind_covars[1]*cov1 + ind_covars[2]*cov2 + rand_int
        prB <- expit(zB)
        yB <- rbinom(n_subj, 1, prB)
        dB <- data.frame(trial = i+nA, 
                         dose = mnB, 
                         cov1 = cov1,
                         cov2 = cov2,
                         z = zB, 
                         pr = prB, 
                         y = yB, 
                         drug = "B", 
                         rand_int = rand_int)
        drugB <- rbind(drugB, dB)
      }
      
      # Organize and summarize all data as above
      all_trials <- rbind(drugA, drugB)
      all_trials$id <- 1:nrow(all_trials)
      
      all_trials <- all_trials %>%
        mutate(drugA = drug == "A",
               drugB = drug == "B",
               drugA_dose = (drug == "A")*dose,
               drugB_dose = (drug == "B")*dose)
      
      sim_dat <- all_trials %>%
        group_by(drug, trial) %>%
        summarise(median_dose = median(dose), 
                  event_rate = mean(y),
                  size = n(),
                  n_obs = sum(y),
                  cov1 = mean(cov1),
                  cov2 = median(cov2)) %>%
        mutate(logit_outcome = logit(event_rate),
               logit_var = (1/n_obs) + (1/(size - n_obs)))
      
      
      sim_dat <- sim_dat %>%
        mutate(drugA = drug == "A",
               drugB = drug == "B",
               drugA_dose = (drug == "A")*median_dose,
               drugB_dose = (drug == "B")*median_dose)
      
      # Check for valid datasets
      dat_valid <- ifelse(min(sim_dat$event_rate) > 0 & max(sim_dat$event_rate) < 1, TRUE, FALSE)
    }
  }
  
  
  return(list("studylevel" = sim_dat, "ipd" = all_trials))
  
}




################################################################################
# Model fitting 
## model : jags model statement
## dat : jags data list
## chains : number of desired sampling chains
## iterations : number of samples per chain
## thin : thinning interval for sampling
## burnin : number of burn-in samples to use per chain
## varnames : desired variables for model to output estimates of
#
# Notes
## 
################################################################################

fit_mod <- function(model, dat, chains = 4, 
                    burnin = 5000, iterations = 20000, thinning = 2,
                    varnames = c("b", "prec_u", "u")) {
  
  mod <- jags.model(file=textConnection(model), data = dat, n.chains = chains)
  
  update(mod, burnin)
  
  samples <- coda.samples(mod, variable.names = varnames, 
                          n.iter = iterations, thin = thinning)
  
  return(samples)
}


################################################################################
# Equivalence calculation
## studylevel_samples : posterior samples for SL (or SL-C) model
## ipd_samples : posterior samples for IPD model
## dose_range : min and max allowed doses
## betas : pre-specified study-level coefficients
## balanced, covars : specifications for study settings and 
##                    presence of additional covariates
## studylevel_nocov_samples : posterior samples for SL-NC model
#
# Notes
## 
################################################################################



equiv_linear_both <- function(studylevel_samples, ipd_samples, 
                              dose_range = c(-1, 1), 
                              betas = c(-.6, -.2, .5, .9),
                              balanced = c(TRUE, FALSE), 
                              covars = c(TRUE, FALSE),
                              studylevel_nocov_samples = NULL){
  
  if(is.null(studylevel_nocov_samples)){
    all_samps_studylevel <- studylevel_samples %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows()
    
    possible_doses <- seq(dose_range[1], dose_range[2], by = .2)
    n <- length(possible_doses)
    
    vars_studylevel <- all_samps_studylevel[,1:4]
    
    vals_studylevel <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_studylevel)),
                            B2 = rep(vars_studylevel[,2], n),
                            B1 = rep(vars_studylevel[,1], n),
                            B3 = rep(vars_studylevel[,3], n),
                            B4 = rep(vars_studylevel[,4], n))
    
    vals_studylevel <- vals_studylevel %>%
      mutate(d1_studylevel = (B3*d0 - B2)/B4)
    
    
    
    
    vals2_studylevel <- vals_studylevel %>% 
      group_by(d0) %>% 
      summarize(lower_studylevel = quantile(d1_studylevel, probs = .025), 
                upper_studylevel = quantile(d1_studylevel, probs = 0.975),
                median_studylevel = median(d1_studylevel)) %>%
      mutate(vals_true = (betas[3]*d0 - betas[2])/betas[4])
    
    all_samps_ipd <- ipd_samples %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows()
    
    vars_ipd <- all_samps_ipd[,1:4]
    
    vals_ipd <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_ipd)),
                            B2 = rep(vars_ipd[,2], n),
                            B1 = rep(vars_ipd[,1], n),
                            B3 = rep(vars_ipd[,3], n),
                            B4 = rep(vars_ipd[,4], n))
    
    vals_ipd <- vals_ipd %>%
      mutate(d1_ipd = (B3*d0 - B2)/B4)
    
    
    vals2_ipd <- vals_ipd %>% 
      group_by(d0) %>% 
      summarize(lower_ipd = quantile(d1_ipd, probs = .025), 
                upper_ipd = quantile(d1_ipd, probs = 0.975),
                median_ipd = median(d1_ipd)) 
    
    vals2_all <- full_join(vals2_studylevel, vals2_ipd, by = "d0")
  }
  
  else{
    all_samps_studylevel <- studylevel_samples %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows()
    
    possible_doses <- seq(dose_range[1], dose_range[2], by = .2)
    n <- length(possible_doses)
    
    vars_studylevel <- all_samps_studylevel[,1:4]
    
    vals_studylevel <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_studylevel)),
                            B2 = rep(vars_studylevel[,2], n),
                            B1 = rep(vars_studylevel[,1], n),
                            B3 = rep(vars_studylevel[,3], n),
                            B4 = rep(vars_studylevel[,4], n))
    
    vals_studylevel <- vals_studylevel %>%
      mutate(d1_studylevel = (B3*d0 - B2)/B4)
    
    
    vals2_studylevel <- vals_studylevel %>% 
      group_by(d0) %>% 
      summarize(lower_studylevel = quantile(d1_studylevel, probs = .025), 
                upper_studylevel = quantile(d1_studylevel, probs = 0.975),
                median_studylevel = median(d1_studylevel)) %>%
      mutate(vals_true = (betas[3]*d0 - betas[2])/betas[4])
    
    all_samps_studylevel_nocov <- studylevel_nocov_samples %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows()
    
    
    vars_studylevel_nocov <- all_samps_studylevel_nocov[,1:4]
    
    vals_studylevel_nocov <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_studylevel_nocov)),
                                  B2 = rep(vars_studylevel_nocov[,2], n),
                                  B1 = rep(vars_studylevel_nocov[,1], n),
                                  B3 = rep(vars_studylevel_nocov[,3], n),
                                  B4 = rep(vars_studylevel_nocov[,4], n))
    
    vals_studylevel_nocov <- vals_studylevel_nocov %>%
      mutate(d1_studylevel_nocov = (B3*d0 - B2 )/B4)
    
    
    vals2_studylevel_nocov <- vals_studylevel_nocov %>% 
      group_by(d0) %>% 
      summarize(lower_studylevel_nocov = quantile(d1_studylevel_nocov, probs = .025), 
                upper_studylevel_nocov = quantile(d1_studylevel_nocov, probs = 0.975),
                median_studylevel_nocov = median(d1_studylevel_nocov))
    
    
    
    all_samps_ipd <- ipd_samples %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows()
    
    vars_ipd <- all_samps_ipd[,1:4]
    
    vals_ipd <- data.frame(d0 = rep(possible_doses, each = nrow(all_samps_ipd)),
                            B2 = rep(vars_ipd[,2], n),
                            B1 = rep(vars_ipd[,1], n),
                            B3 = rep(vars_ipd[,3], n),
                            B4 = rep(vars_ipd[,4], n))
    
    vals_ipd <- vals_ipd %>%
      mutate(d1_ipd = (B3*d0 - B2)/B4)
    
    
    vals2_ipd <- vals_ipd %>% 
      group_by(d0) %>% 
      summarize(lower_ipd = quantile(d1_ipd, probs = .025), 
                upper_ipd = quantile(d1_ipd, probs = 0.975),
                median_ipd = median(d1_ipd)) 
    
    vals2_all <- full_join(vals2_studylevel, vals2_studylevel_nocov, by = "d0") %>%
      full_join(vals2_ipd, by = "d0")
  }
  
  
  return(vals2_all)
}

