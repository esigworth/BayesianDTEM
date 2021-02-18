source("slope_component_functions.R")
source("slope_full_run.R")

eval_sim_results_slope <- function(all_sims, tau2, total, betas = c(-.6, -.2, .5, .9, .2, .5, -.3, .4), 
                                   bal = c("balanced", "unbalanced"), covars = "with covars"){
  
  expr1 <- expression(beta[1])
  expr2 <- expression(beta[2])
  expr3 <- expression(beta[3])
  expr4 <- expression(beta[4])
  expr5 <- expression(beta[5])
  expr6 <- expression(beta[6])
  expr7 <- expression(beta[7])
  expr8 <- expression(beta[8])
  
  prec <- 1/tau2
  
  
  studylevel <- studylevel_nocov <- ipd <- studylevel_median <- studylevel_nocov_median <- ipd_median <- vector(mode = "list", length = total)
  equiv <- mad <- vector(mode = "list", length = total)
  
  
  for(m in 1:total){
    
    studylevel[[m]] <- all_sims[[m]]$run$studylevel_summary
    studylevel_nocov[[m]] <- all_sims[[m]]$run$studylevel_nocov_summary
    ipd[[m]] <- all_sims[[m]]$run$ipd_summary
    
    studylevel_median[[m]] <- all_sims[[m]]$run$studylevel_summary_median
    studylevel_nocov_median[[m]] <- all_sims[[m]]$run$studylevel_nocov_median
    ipd_median[[m]] <- all_sims[[m]]$run$ipd_summary_median
    
    equiv[[m]] <- all_sims[[m]]$run$equiv
    
    mad[[m]] <- all_sims[[m]]$run$mad
    
    
  }
  
  
  mad_studylevel <- mad_ipd <- mad_ratios <- mad_studylevel_nocov <- mad_ratios_nocov <- data.frame(b1_mad = rep(NA, total), 
                                                                                         b2_mad = rep(NA, total), 
                                                                                         b3_mad = rep(NA, total), 
                                                                                         b4_mad = rep(NA, total),
                                                                                         b5_mad = rep(NA, total),
                                                                                         b6_mad = rep(NA, total),
                                                                                         b7_mad = rep(NA, total),
                                                                                         b8_mad = rep(NA, total),
                                                                                         type = rep(NA, total))
  
  
  for(j in 1:total){
    mad_studylevel[j,] <- mad[[j]][1,]
    mad_studylevel_nocov[j,] <- mad[[j]][2,]
    mad_ipd[j,] <- mad[[j]][3,]
    mad_ratios[j,] <- mad[[j]][4,]
    mad_ratios_nocov[j,] <- mad[[j]][5,]
  }
  
  
  mad <- rbind(mad_studylevel, mad_ipd, mad_studylevel_nocov)
  mad_ratios <- rbind(mad_ratios, mad_ratios_nocov)
  
  
  mad2 <- mad %>%
    pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad, b5_mad, b6_mad, b7_mad, b8_mad)) %>%
    mutate(value = as.double(unlist(value))) %>%
    group_by(type, name) %>%
    summarise(low = quantile(value, probs = .025, na.rm = TRUE),
              high = quantile(value, probs = .975, na.rm = TRUE),
              med = median(value, na.rm = TRUE)) 
  
  mad2$pos <- c(seq(1.25, 15.25, by = 2), seq(1.75, 15.75, by = 2), seq(1.5, 15.5, by = 2))
  
  
  mad2$param <- rep(c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8"),3)
  mad2$fit <- paste(mad2$param, mad2$type, sep = "_")
  
  mad_ratios2 <- mad_ratios %>%
    pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad, b5_mad, b6_mad, b7_mad, b8_mad)) %>%
    group_by(type, name) %>%
    summarise(low = quantile(value, probs = .025, na.rm = TRUE),
              high = quantile(value, probs = .975, na.rm = TRUE),
              med = median(value, na.rm = TRUE))
  
  
  mad_ratios2$param <- rep(c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8"), 2)
  mad_ratios2$pos <- c(seq(1.25, 15.25, by = 2), seq(1.75, 15.75, by = 2))
  
  
  eff <- eff_nocov <- studylevel_bias <- studylevel_nocov_bias <- ipd_bias <- studylevel_se <- studylevel_nocov_se <- ipd_se <- studylevel_cov_prob <- studylevel_nocov_cov_prob <- ipd_cov_prob <- vector(mode = "list", length = total)
  
  betas_n <- c(betas, prec)
  
  for(i in 1:total){
    
    studylevel_cov_prob[[i]] <- c((betas_n[c(1:6,9)] > studylevel_median[[i]][,1]) & (betas_n[c(1:6,9)] < studylevel_median[[i]][,5]))
    studylevel_nocov_cov_prob[[i]] <- c((betas_n[c(1:4,9)] > studylevel_nocov_median[[i]][,1]) & (betas_n[c(1:4,9)] < studylevel_nocov_median[[i]][,5]))
    ipd_cov_prob[[i]] <- c((betas_n > ipd_median[[i]][,1]) & (betas_n < ipd_median[[i]][,5]))
    
    eff[[i]] <- ipd[[i]][c(1:6,9),2]^2 / studylevel[[i]][,2]^2
    eff_nocov[[i]] <- ipd[[i]][c(1:4,9),2]^2 / studylevel_nocov[[i]][,2]^2
    
    studylevel_bias[[i]] <- round((studylevel[[i]][,1] - betas_n[c(1:6,9)])/betas_n[c(1:6,9)], 2)
    studylevel_nocov_bias[[i]] <- round((studylevel_nocov[[i]][,1] - betas_n[c(1:4,9)])/betas_n[c(1:4,9)], 2)
    ipd_bias[[i]] <- round((ipd[[i]][,1] - betas_n)/betas_n, 2)
    
    studylevel_se[[i]] <- round((studylevel[[i]][,1] - betas_n[c(1:6,9)])^2, 2)
    studylevel_nocov_se[[i]] <- round((studylevel_nocov[[i]][,1] - betas_n[c(1:4,9)])^2, 2)
    ipd_se[[i]] <- round((ipd[[i]][,1] - betas_n)^2, 2)
  }
  
  studylevel_beta <- studylevel_nocov_beta <- ipd_beta <- data.frame(b1=rep(NA, total), 
                                                          b2=rep(NA, total), 
                                                          b3=rep(NA, total), 
                                                          b4=rep(NA, total),
                                                          b5=rep(NA, total),
                                                          b6=rep(NA, total),
                                                          b7=rep(NA, total),
                                                          b8=rep(NA, total))
  
  for(i in 1:total){
    studylevel_beta[i,] <- c(studylevel[[i]][1:6,1], NA, NA)
    studylevel_nocov_beta[i,] <- c(studylevel_nocov[[i]][1:4,1], NA, NA, NA, NA)
    ipd_beta[i,] <- ipd[[i]][1:8,1]
  }
  
  studylevel_avg_bias <- rep(NA, 7)
  studylevel_nocov_avg_bias <- rep(NA, 5)
  ipd_avg_bias <- rep(NA, 9)
  
  studylevel_low_high_bias <- matrix(NA, ncol = 2, nrow = 7)
  studylevel_nocov_low_high_bias <- matrix(NA, ncol = 2, nrow = 5)
  ipd_low_high_bias <- matrix(NA, ncol = 2, nrow = 9)
  
  avg_eff <- rep(NA,7)
  avg_nocov_eff <- rep(NA, 5)
  low_high_eff <- matrix(NA, ncol = 2, nrow = 7)
  low_high_nocov_eff <- matrix(NA, ncol = 2, nrow = 5)
  
  studylevel_mse <- ipd_mse <- rep(NA, 7)
  studylevel_nocov_mse <- rep(NA, 5)
  
  for(i in 1:9){
    ipd_avg_bias[i] <- mean(sapply(ipd_bias, "[[", i), na.rm = TRUE)
    
    ipd_low_high_bias[i,] <- quantile(sapply(ipd_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    ipd_mse[i] <- mean(sapply(ipd_se, "[[", i))
    
  }
  
  for(i in 1:7){
    studylevel_avg_bias[i] <- mean(sapply(studylevel_bias, "[[", i), na.rm = TRUE) 
    
    studylevel_low_high_bias[i,] <- quantile(sapply(studylevel_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    avg_eff[i] <- mean(sapply(eff, "[[", i))
    
    low_high_eff[i,] <- quantile(sapply(eff, "[[", i), probs = c(.025, .975))
    
    studylevel_mse[i] <- mean(sapply(studylevel_se, "[[", i))
    
  }
  
  
  for(i in 1:5){
    studylevel_nocov_avg_bias[i] <- mean(sapply(studylevel_nocov_bias, "[[", i), na.rm = TRUE) 
    
    studylevel_nocov_low_high_bias[i,] <- quantile(sapply(studylevel_nocov_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    avg_nocov_eff[i] <- mean(sapply(eff_nocov, "[[", i))
    
    low_high_nocov_eff[i,] <- quantile(sapply(eff_nocov, "[[", i), probs = c(.025, .975))
    
    studylevel_nocov_mse[i] <- mean(sapply(studylevel_nocov_se, "[[", i))
    
  }
  
  mse_ratio <- round(c(studylevel_mse[1:6]/ipd_mse[1:6]),3)
  mse_ratio <-c(mse_ratio, NA, NA)
  mse_nocov_ratio <- round(c(studylevel_nocov_mse[1:4]/ipd_mse[1:4]),3)
  mse_nocov_ratio <- c(mse_nocov_ratio[1:4], NA, NA, NA, NA)
  studylevel_nocov_mse <- c(studylevel_nocov_mse[1:4], NA, NA, NA, NA)
  studylevel_mse <- c(studylevel_mse[1:6], NA, NA)
  mse <- data.frame(rbind(studylevel_mse, ipd_mse[1:8], studylevel_nocov_mse, mse_ratio, mse_nocov_ratio))
  mse <- data.frame(stat = "mse", type = c("studylevel", "ipd", "studylevel_nocov", "ratio_slcipd", "ratio_slncipd"), mse)
  
  colnames(mse) <- c("stat", "type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8")
  mse <- as.data.frame(mse)
  rownames(mse) <- NULL
  
  mse <- do.call(data.frame,lapply(mse, function(x) replace(x, x == "Inf" | x == "NaN",NA)))
  
  # Percent bias figures
  
  Type <- rep(c( "sl_nc", "sl_c", "ipd"),8)
  
  
  bias <- data.frame(studylevel = c(studylevel_avg_bias[1:6], NA, NA), 
                     ipd = c(ipd_avg_bias[1:8]), 
                     studylevel_nocov = c(studylevel_nocov_avg_bias[1:4], NA, NA, NA, NA),
                     Beta = c(paste0("B", 1:8, sep = ""))) %>%
    pivot_longer(cols = c(studylevel, studylevel_nocov, ipd),
                 names_to = "Model", values_to = "Percent Bias") 
  
  bias$Type <- Type
  
  err_bars <- data.frame(rbind(studylevel_low_high_bias[1:6,], studylevel_nocov_low_high_bias[1:4,], ipd_low_high_bias[1:8,]))
  names(err_bars) <- c("Low", "High")
  err_bars$Type <- rep(c(rep("studylevel", 6), rep("studylevel NC", 4), rep("ipd", 8)))
  err_bars$Beta <- c(paste0("B", 1:6, sep = ""), paste0("B", 1:4, sep = ""), paste0("B", 1:8, sep = ""))
  
  bias <- full_join(bias, err_bars, by = c("Type", "Beta")) 
  
  
  bias$Model_Coef <- paste(bias$Beta, bias$Type, sep = "_")

  avg_nocov_eff <- c(avg_nocov_eff[1:4], NA, NA)
  avg_eff <- data.frame(type = c("slcipd", "slncipd"), rbind(avg_eff[1:6], avg_nocov_eff))
  colnames(avg_eff) <- c("type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6")
  rownames(avg_eff) <- NULL
  
  low_high_eff <- data.frame(rbind(low_high_eff[1:6,], low_high_nocov_eff[1:4,]))
  low_high_eff$Type <- c(rep("slcipd", 6), rep("slncipd", 4))
  colnames(low_high_eff)[1:2] <- c("Low", "High")
  low_high_eff$Beta <- c(c(paste0("B", 1:6, sep="")), c(paste0("B", 1:4, sep="")))
  
  avg_eff <- pivot_longer(avg_eff, cols = 2:7, names_to = "Beta", values_to = "Efficiency")
  avg_eff$Beta <- c(c(paste0("B", 1:6, sep="")), c(paste0("B", 1:6, sep="")))
  
  
  eff_dat <- full_join(avg_eff, low_high_eff, by = c("Beta", "type" = "Type")) %>% arrange(Beta)

  
  d0 <- seq(-1, 1, by = .2)

  vals_true_female <- ((betas[3]+betas[7])*d0 - betas[2])/(betas[4]+betas[8])
  vals_true_male <- ((betas[3])*d0 - betas[2])/(betas[4])
  
  
  avg_equiv <- data.frame(d0 = seq(-1, 1, by = .2),
                          lower_studylevel = NA,
                          upper_studylevel = NA,
                          median_studylevel = NA,
                          vals_true_female = vals_true_female,
                          lower_ipd_female = NA, 
                          upper_ipd_female = NA, 
                          median_ipd_female = NA,
                          lower_ipd_male = NA, 
                          upper_ipd_male = NA, 
                          median_ipd_male = NA,
                          lower_studylevel_nocov = NA,
                          upper_studylevel_nocov = NA,
                          median_studylevel_nocov = NA,
                          vals_true_male = vals_true_male
  )
  
  
  for(i in c(2:4, 6:14)){
    avg_equiv[,i] <- rowMedians(sapply(equiv, "[[", i))
  }

  # Coverage probability using credible intervals  
  
  studylevel_coverage <- rep(NA, 6)
  ipd_coverage <- rep(NA, 8)
  studylevel_nocov_coverage <- rep(NA, 4)
  
  for(i in 1:8){
    ipd_coverage[i] <- mean(sapply(ipd_cov_prob, "[[", i))
  }
  
  for(i in 1:6){
    studylevel_coverage[i] <- mean(sapply(studylevel_cov_prob, "[[", i)) 
  }
  
  for(i in 1:4){
    studylevel_nocov_coverage[i] <- mean(sapply(studylevel_nocov_cov_prob, "[[", i)) 
  }
  
  studylevel_nocov_coverage <- c(studylevel_nocov_coverage[1:4], NA, NA, NA, NA)
  studylevel_coverage <- c(studylevel_coverage, NA, NA)
  
  
  
  coverage <- data.frame(type = c("studylevel", "studylevel_nocov", "ipd"), b1 = NA, b2 = NA, b3 = NA, b4 = NA, b5 = NA, b6 = NA, b7 = NA, b8 = NA)
  coverage[,2:9] <- rbind(studylevel_coverage, studylevel_nocov_coverage, ipd_coverage)  
  
  all_vals <- list("bias" = bias, 
                   "eff_dat" = eff_dat, 
                   "studylevel" = studylevel, 
                   "studylevel_nocov" = studylevel_nocov,
                   "ipd" = ipd, 
                   "mad" = mad2, 
                   "mad_ratios" = mad_ratios2,
                   "mse" = mse,
                   "coverage" = coverage,
                   "equiv_dat" = avg_equiv)
  
  
  
  
  return(all_vals)
}

