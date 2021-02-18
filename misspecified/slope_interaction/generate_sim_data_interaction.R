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



logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  1/(1+exp(-x))
}


sim_linear_int <- function(nA = 75, nB = 75, 
                       nsubjA = 100, nsubjB = 100, 
                       minsubj = 50,
                       maxsubj = 200, 
                       betas = c(-.6, -.2, .5, .9, .2, .5, -.3, .4), 
                       balanced_trials = TRUE,
                       tau2 = 0,
                       doses){
  
  dat_valid <- FALSE
  
  doses_A <- doses[1:75]
  doses_B <- doses[76:150]
  
  drugA <- data.frame()
  drugB <- data.frame()
  
  
  while(dat_valid == FALSE){
    
    for(i in 1:nA){
      rand_int <- rnorm(1, 0, sqrt(tau2))
      n_subj <- ifelse(balanced_trials, nsubjA, 
                       sample(seq(minsubj, maxsubj, by = 10), 1))
      mnA <- doses_A[i]
      cov1 <- rbinom(n_subj, 1, runif(1, 0.2, 0.5))
      cov2 <- rnorm(n_subj, 0, 1)
      zA <- betas[1] + betas[3]*mnA + betas[5]*cov1 + betas[6]*cov2 + betas[7]*cov1*mnA + rand_int
      prA <- expit(zA)
      yA <- rbinom(n_subj, 1, prA)
      dA <- data.frame(trial = i, 
                       dose = mnA, 
                       cov1 = cov1,
                       cov2 = cov2,
                       z = zA, 
                       pr = prA, 
                       y = yA, 
                       drug = "A",
                       rand_int = rand_int)
      drugA <- rbind(drugA, dA)
    }
    for(i in 1:nB){
      rand_int <- rnorm(1, 0, sqrt(tau2))
      n_subj <- ifelse(balanced_trials, nsubjB, 
                       sample(seq(minsubj, maxsubj, by = 10), 1))
      mnB <- doses_B[i]
      cov1 <- rbinom(n_subj, 1, runif(1, 0.2, 0.5))
      cov2 <- rnorm(n_subj, 0, 1)
      zB <- betas[1] + betas[2] + betas[4]*mnB + betas[5]*cov1 + betas[6]*cov2 + betas[8]*cov1*mnB + rand_int
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
    
    dat_valid <- ifelse(min(sim_dat$event_rate) > 0 & max(sim_dat$event_rate) < 1, TRUE, FALSE)
  }
  
  
  return(list("meta" = sim_dat, "full" = all_trials))
  
}


