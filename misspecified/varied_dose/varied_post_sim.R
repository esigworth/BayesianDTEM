eval_sim_varied <- function(all_sims, tau2, total, betas = c(-.6, -.8, .5, .9), 
                             bal = c("balanced", "unbalanced")){
  
  expr1 <- expression(beta[1])
  expr2 <- expression(beta[2])
  expr3 <- expression(beta[3])
  expr4 <- expression(beta[4])
  
  prec <- 1/tau2
  
  
  studylevel <- ipd <- studylevel_median <- ipd_median <- vector(mode = "list", length = total)
  equiv <- mad <- ess <- gr_stats <- vector(mode = "list", length = total)
  
  
  for(m in 1:total){
    
    studylevel[[m]] <- all_sims[[m]]$run$studylevel_summary
    ipd[[m]] <- all_sims[[m]]$run$ipd_summary
    
    studylevel_median[[m]] <- all_sims[[m]]$run$studylevel_summary_median
    ipd_median[[m]] <- all_sims[[m]]$run$ipd_summary_median
    
    equiv[[m]] <- all_sims[[m]]$run$equiv
    
    mad[[m]] <- all_sims[[m]]$run$mad
    
    ess[[m]] <- all_sims[[m]]$run$ess
    
    gr_stats[[m]] <- all_sims[[m]]$run$gr_stats
  }
  
  betas <- betas[1:4]
  mad_studylevel <- mad_ipd <- mad_ratios <- data.frame(type = rep(NA, total),
                                                   b1_mad = rep(NA, total), 
                                                   b2_mad = rep(NA, total), 
                                                   b3_mad = rep(NA, total), 
                                                   b4_mad = rep(NA, total))
  
  
  for(j in 1:total){
    mad_studylevel[j,] <- mad[[j]][1,]
    mad_ipd[j,] <- mad[[j]][2,]
    mad_ratios[j,] <- mad[[j]][3,]
  }
  
  
  gr_stats2 <- gr_stats %>%
    lapply(., function(x) as.data.frame(x)) %>%
    bind_rows() %>%
    group_by(type, param) %>%
    summarise(med_gr = median(`Point est.`), med_upper_gr = median(`Upper C.I.`)) 
  
  mad <- rbind(mad_studylevel, mad_ipd)
  
  mad2 <- mad %>%
    pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad)) %>%
    group_by(type, name) %>%
    summarise(low = quantile(value, probs = .025),
              high = quantile(value, probs = .975),
              med = median(value)) 
  
  mad2$param <- rep(c("B1", "B2", "B3", "B4"),2)
  mad2$fit <- paste(mad2$param, mad2$type, sep = "_")
  mad2$pos <- c(seq(1.25, 7.25, by = 2), seq(1.75, 7.75, by = 2))
  
  mad_ratios2 <- mad_ratios %>%
    pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad)) %>%
    group_by(name) %>%
    summarise(low = quantile(value, probs = .025),
              high = quantile(value, probs = .975),
              med = median(value))
  
  mad_ratios2$param <- c("B1", "B2", "B3", "B4")
  
  
  
  
  
  
  eff <- studylevel_bias <- ipd_bias <- studylevel_se <- ipd_se <- studylevel_cov_prob <- ipd_cov_prob <- vector(mode = "list", length = total)
  
  betas_n <- c(betas, prec)
  
  for(i in 1:total){
    
    studylevel_cov_prob[[i]] <- c((betas_n > studylevel_median[[i]][,1]) & (betas_n < studylevel_median[[i]][,5]))
    ipd_cov_prob[[i]] <- c((betas_n > ipd_median[[i]][,1]) & (betas_n < ipd_median[[i]][,5]))
    
    eff[[i]] <- ipd[[i]][,2]^2 / studylevel[[i]][,2]^2
    
    studylevel_bias[[i]] <- round((studylevel[[i]][,1] - betas_n)/betas_n, 2)
    
    ipd_bias[[i]] <- round((ipd[[i]][,1] - betas_n)/betas_n, 2)
    
    studylevel_se[[i]] <- round((studylevel[[i]][,1] - betas_n)^2, 2)
    ipd_se[[i]] <- round((ipd[[i]][,1] - betas_n)^2, 2)
  }
  
  studylevel_beta <- ipd_beta <- data.frame(b1=rep(NA, total), 
                                       b2=rep(NA, total), 
                                       b3=rep(NA, total), 
                                       b4=rep(NA, total))
  
  for(i in 1:total){
    studylevel_beta[i,] <- studylevel[[i]][1:4,1]
    ipd_beta[i,] <- ipd[[i]][1:4,1]
  }
  
  
  studylevel_avg_bias <- rep(NA, 5)
  ipd_avg_bias <- rep(NA, 5)
  
  
  studylevel_low_high_bias <- matrix(NA, ncol = 2, nrow = 5)
  ipd_low_high_bias <- matrix(NA, ncol = 2, nrow = 5)
  
  avg_eff <- rep(NA,5)
  low_high_eff <- matrix(NA, ncol = 2, nrow = 5)
  
  studylevel_mse <- ipd_mse <- rep(NA, 5)
  
  for(i in 1:5){
    studylevel_avg_bias[i] <- mean(sapply(studylevel_bias, "[[", i), na.rm = TRUE) 
    ipd_avg_bias[i] <- mean(sapply(ipd_bias, "[[", i), na.rm = TRUE)
    
    studylevel_low_high_bias[i,] <- quantile(sapply(studylevel_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    ipd_low_high_bias[i,] <- quantile(sapply(ipd_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    avg_eff[i] <- mean(sapply(eff, "[[", i))
    
    low_high_eff[i,] <- quantile(sapply(eff, "[[", i), probs = c(.025, .975))
    
    studylevel_mse[i] <- mean(sapply(studylevel_se, "[[", i))
    ipd_mse[i] <- mean(sapply(ipd_se, "[[", i))
    
  }
  
  mse_ratio <- round(c(studylevel_mse/ipd_mse),3)
  mse <- rbind(studylevel_mse, ipd_mse, mse_ratio)
  mse <- cbind(rep("mse", 3), c("studylevel", "ipd", "ratio"), mse)
  
  
  colnames(mse) <- c("stat", "type", "beta1", "beta2", "beta3", "beta4", "prec")
  mse <- as.data.frame(mse)
  rownames(mse) <- NULL
  
  # Percent bias figures
  
  Type <- rep(c("studylevel", "ipd"),5)
  
  
  bias <- data.frame(studylevel = c(studylevel_avg_bias), ipd = c(ipd_avg_bias), 
                     Beta = c(paste0("B", 1:4, sep = ""), "prec")) %>%
    pivot_longer(cols = c(studylevel, ipd),
                 names_to = "Model", values_to = "Percent Bias") 
  
  bias$Type <- Type
  
  err_bars <- data.frame(rbind(studylevel_low_high_bias, ipd_low_high_bias))
  names(err_bars) <- c("Low", "High")
  err_bars$Type <- rep(c(rep("studylevel", 5), rep("ipd", 5)))
  err_bars$Beta <- rep(c(paste0("B", 1:4, sep = ""), "prec"), 2)
  
  bias <- full_join(bias, err_bars, by = c("Type", "Beta"))
  
  bias$Model_Coef <- paste(bias$Beta, bias$Type, sep = "_")
  
  fig_dat <- filter(bias, Beta != "prec")
  
  fig_dat$pos <- c(1.25, 1.75, 3.25, 3.75, 5.25, 5.75, 7.25, 7.75)
  
  min_pb <- min(fig_dat$Low) - .1
  max_pb <- max(fig_dat$High) + .1
  
  
  avg_eff <- data.frame(Efficiency = avg_eff,
                        Beta = c(paste0("B", 1:4, sep=""), "prec")) 
  
  eff_bars <- data.frame(low_high_eff)
  names(eff_bars) <- c("Low", "High")
  eff_bars$Beta <- c(paste0("B", 1:4, sep=""), "prec")
  
  eff_dat <- full_join(avg_eff, eff_bars, by = c("Beta"))
  
  
  avg_equiv <- data.frame(d0 = seq(-1, 1, by = .2),
                          lower_studylevel = NA, 
                          upper_studylevel = NA, 
                          median_studylevel = NA,
                          vals_true = equiv[[1]][,5],
                          lower_ipd = NA, 
                          upper_ipd = NA, 
                          median_ipd = NA
  )
  
  for(i in c(2:4, 6:8)){
    avg_equiv[,i] <- rowMedians(sapply(equiv, "[[", i))
  }
  
 
  # Coverage probability using credible intervals  
  
  studylevel_coverage <- rep(NA, 5)
  ipd_coverage <- rep(NA, 5)
  
  for(i in 1:5){
    studylevel_coverage[i] <- mean(sapply(studylevel_cov_prob, "[[", i)) 
    ipd_coverage[i] <- mean(sapply(ipd_cov_prob, "[[", i))
  }
  
  coverage <- data.frame(type = c("studylevel", "ipd"), b1 = NA, b2 = NA, b3 = NA, b4 = NA, prec = NA)
  coverage[,2:6] <- rbind(studylevel_coverage, ipd_coverage)  
  
  
  all_vals <- list("bias" = bias, 
                   "eff_dat" = eff_dat, 
                   "studylevel" = studylevel, 
                   "ipd" = ipd, 
                   "mad" = mad2, 
                   "mad_ratios" = mad_ratios2,
                   "mse" = mse,
                   "coverage" = coverage,
                   "gr_stats" = gr_stats2,
                   "equiv_dat" = avg_equiv)
  
  return(all_vals)
}

