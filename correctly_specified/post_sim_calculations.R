source("component_functions.R")
source("full_run.R")

eval_sim_results_ref <- function(all_sims, tau2, total, betas = c(-.6, -.2, .5, .9, .2, .5), 
                             bal = c("balanced", "unbalanced"), covars = c("with covars", "no covars")){
  
  expr1 <- expression(beta[1])
  expr2 <- expression(beta[2])
  expr3 <- expression(beta[3])
  expr4 <- expression(beta[4])
  expr5 <- expression(beta[5])
  expr6 <- expression(beta[6])
  
  prec <- 1/tau2
  
  
  studylevel <- studylevel_nocov <- ipd <- studylevel_median <- studylevel_nocov_median <- ipd_median <- vector(mode = "list", length = total)
  equiv <- mad <- ess <- gr_stats <- vector(mode = "list", length = total)
  
  
  for(m in 1:total){
    
    studylevel[[m]] <- all_sims[[m]]$run$studylevel_summary
    studylevel_nocov[[m]] <- all_sims[[m]]$run$studylevel_nocov_summary
    ipd[[m]] <- all_sims[[m]]$run$ipd_summary
    
    studylevel_median[[m]] <- all_sims[[m]]$run$studylevel_summary_median
    studylevel_nocov_median[[m]] <- all_sims[[m]]$run$studylevel_summary_median_nocov
    ipd_median[[m]] <- all_sims[[m]]$run$ipd_summary_median
    
    equiv[[m]] <- all_sims[[m]]$run$equiv
    
    mad[[m]] <- all_sims[[m]]$run$mad
    
    ess[[m]] <- all_sims[[m]]$run$ess
    
    gr_stats[[m]] <- all_sims[[m]]$run$gr_stats
  }
  
  if(covars == "no covars"){
    betas <- betas[1:4]
    mad_studylevel <- mad_ipd <- mad_ratios <- data.frame(type = rep(NA, total),
                                                     b1_mad = rep(NA, total), 
                                                     b2_mad = rep(NA, total), 
                                                     b3_mad = rep(NA, total), 
                                                     b4_mad = rep(NA, total))
    
    ess_studylevel <- ess_ipd <- data.frame(type = rep("meh", total),
                                       b1_ess = rep(NA, total), 
                                       b2_ess = rep(NA, total), 
                                       b3_ess = rep(NA, total), 
                                       b4_ess = rep(NA, total))
    
    
    for(j in 1:total){
      mad_studylevel[j,] <- mad[[j]][1,]
      mad_ipd[j,] <- mad[[j]][2,]
      mad_ratios[j,] <- mad[[j]][3,]
      
      ess_studylevel[j,] <- ess[[j]][1,1:5]
      ess_ipd[j,] <- ess[[j]][2,1:5]
    }
    
    ess_studylevel$type <- "studylevel"
    ess_ipd$type <- "ipd"
    
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
    
    mad_ratios2 <- mad_ratios %>%
      pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad)) %>%
      group_by(name) %>%
      summarise(low = quantile(value, probs = .025),
                high = quantile(value, probs = .975),
                med = median(value))
    
    mad_ratios2$param <- c("B1", "B2", "B3", "B4")
    
    ess <- rbind(ess_studylevel, ess_ipd)
    
    ess2 <- ess %>%
      pivot_longer(cols = c(b1_ess, b2_ess, b3_ess, b4_ess)) %>%
      group_by(type, name) %>%
      summarise(low = quantile(value, probs = .025),
                high = quantile(value, probs = .975),
                med = median(value)) 
    
    ess2$param <- rep(c("B1", "B2", "B3", "B4"),2)
    ess2$fit <- paste(ess2$param, ess2$type, sep = "_")
    
    
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
                     "ess" = ess2,
                     "gr_stats" = gr_stats2,
                     "equiv_dat" = avg_equiv)
    
  }
  
  else {
    mad_studylevel <- mad_ipd <- mad_ratios <- mad_studylevel_nocov <- mad_ratios_nocov <- data.frame(type = rep(NA, total),
                                                                                           b1_mad = rep(NA, total), 
                                                                                           b2_mad = rep(NA, total), 
                                                                                           b3_mad = rep(NA, total), 
                                                                                           b4_mad = rep(NA, total),
                                                                                           b5_mad = rep(NA, total),
                                                                                           b6_mad = rep(NA, total))
    
    ess_studylevel <- ess_ipd <- data.frame(type = rep("mod", total),
                                       b1_ess = rep(NA, total), 
                                       b2_ess = rep(NA, total), 
                                       b3_ess = rep(NA, total), 
                                       b4_ess = rep(NA, total),
                                       b5_ess = rep(NA, total),
                                       b6_ess = rep(NA, total))
    
    
    for(j in 1:total){
      mad_studylevel[j,] <- mad[[j]][1,]
      mad_studylevel_nocov[j,] <- mad[[j]][2,]
      mad_ipd[j,] <- mad[[j]][3,]
      mad_ratios[j,] <- mad[[j]][4,]
      mad_ratios_nocov[j,] <- mad[[j]][5,]
      
      ess_studylevel[j,] <- ess[[j]][1,1:7]
      ess_ipd[j,] <- ess[[j]][2,1:7]
    }
    
    ess_studylevel$type <- "studylevel"
    ess_ipd$type <- "ipd"
    
    gr_stats2 <- gr_stats %>%
      lapply(., function(x) as.data.frame(x)) %>%
      bind_rows() %>%
      group_by(type, param) %>%
      summarise(med_gr = median(`Point est.`), med_upper_gr = median(`Upper C.I.`)) 
    
    mad <- rbind(mad_studylevel, mad_ipd, mad_studylevel_nocov)
    mad_ratios <- rbind(mad_ratios, mad_ratios_nocov)
    
    mad2 <- mad %>%
      pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad, b5_mad, b6_mad)) %>%
      group_by(type, name) %>%
      summarise(low = quantile(value, probs = .025, na.rm = TRUE),
                high = quantile(value, probs = .975, na.rm = TRUE),
                med = median(value, na.rm = TRUE)) 
    
    
    mad2$param <- rep(c("B1", "B2", "B3", "B4", "B5", "B6"),3)
    mad2$fit <- paste(mad2$param, mad2$type, sep = "_")
    
    mad_ratios2 <- mad_ratios %>%
      pivot_longer(cols = c(b1_mad, b2_mad, b3_mad, b4_mad, b5_mad, b6_mad)) %>%
      group_by(type, name) %>%
      summarise(low = quantile(value, probs = .025, na.rm = TRUE),
                high = quantile(value, probs = .975, na.rm = TRUE),
                med = median(value, na.rm = TRUE))
    
    
    mad_ratios2$param <- rep(c("B1", "B2", "B3", "B4", "B5", "B6"), 2)


    
    eff <- eff_nocov <- studylevel_bias <- studylevel_nocov_bias <- ipd_bias <- studylevel_se <- studylevel_nocov_se <- ipd_se <- studylevel_cov_prob <- studylevel_nocov_cov_prob <- ipd_cov_prob <- vector(mode = "list", length = total)
    
    betas_n <- c(betas, prec)
    
    for(i in 1:total){
      
      studylevel_cov_prob[[i]] <- c((betas_n > studylevel_median[[i]][,1]) & (betas_n < studylevel_median[[i]][,5]))
      studylevel_nocov_cov_prob[[i]] <- c((betas_n[c(1:4,7)] > studylevel_nocov_median[[i]][,1]) & (betas_n[c(1:4,7)] < studylevel_nocov_median[[i]][,5]))
      ipd_cov_prob[[i]] <- c((betas_n > ipd_median[[i]][,1]) & (betas_n < ipd_median[[i]][,5]))
      
      eff[[i]] <- ipd[[i]][,2]^2 / studylevel[[i]][,2]^2
      eff_nocov[[i]] <- ipd[[i]][c(1:4,7),2]^2 / studylevel_nocov[[i]][,2]^2
      
      studylevel_bias[[i]] <- round((studylevel[[i]][,1] - betas_n)/betas_n, 2)
      studylevel_nocov_bias[[i]] <- round((studylevel_nocov[[i]][,1] - betas_n[c(1:4,7)])/betas_n[c(1:4,7)], 2)
      ipd_bias[[i]] <- round((ipd[[i]][,1] - betas_n)/betas_n, 2)
      
      studylevel_se[[i]] <- round((studylevel[[i]][,1] - betas_n)^2, 2)
      studylevel_nocov_se[[i]] <- round((studylevel_nocov[[i]][,1] - betas_n[c(1:4,7)])^2, 2)
      ipd_se[[i]] <- round((ipd[[i]][,1] - betas_n)^2, 2)
    }
    
    studylevel_beta <- studylevel_nocov_beta <- ipd_beta <- data.frame(b1=rep(NA, total), 
                                                            b2=rep(NA, total), 
                                                            b3=rep(NA, total), 
                                                            b4=rep(NA, total),
                                                            b5=rep(NA, total),
                                                            b6=rep(NA, total))
    
    for(i in 1:total){
      studylevel_beta[i,] <- studylevel[[i]][1:6,1]
      studylevel_nocov_beta[i,] <- c(studylevel_nocov[[i]][1:4,1], NA, NA)
      ipd_beta[i,] <- ipd[[i]][1:6,1]
    }
    
    
    studylevel_avg_bias <- rep(NA, 7)
    studylevel_nocov_avg_bias <- rep(NA, 5)
    ipd_avg_bias <- rep(NA, 7)
    
    studylevel_low_high_bias <- matrix(NA, ncol = 2, nrow = 7)
    studylevel_nocov_low_high_bias <- matrix(NA, ncol = 2, nrow = 5)
    ipd_low_high_bias <- matrix(NA, ncol = 2, nrow = 7)
    
    avg_eff <- rep(NA,7)
    avg_nocov_eff <- rep(NA, 5)
    low_high_eff <- matrix(NA, ncol = 2, nrow = 7)
    low_high_nocov_eff <- matrix(NA, ncol = 2, nrow = 5)
    
    studylevel_mse <- ipd_mse <- rep(NA, 7)
    studylevel_nocov_mse <- rep(NA, 5)
    
    for(i in 1:7){
      studylevel_avg_bias[i] <- mean(sapply(studylevel_bias, "[[", i), na.rm = TRUE) 
      ipd_avg_bias[i] <- mean(sapply(ipd_bias, "[[", i), na.rm = TRUE)
      
      studylevel_low_high_bias[i,] <- quantile(sapply(studylevel_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
      ipd_low_high_bias[i,] <- quantile(sapply(ipd_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
      
      avg_eff[i] <- mean(sapply(eff, "[[", i))
      
      low_high_eff[i,] <- quantile(sapply(eff, "[[", i), probs = c(.025, .975))
      
      studylevel_mse[i] <- mean(sapply(studylevel_se, "[[", i))
      ipd_mse[i] <- mean(sapply(ipd_se, "[[", i))
      
    }
    
    studylevel_nocov_bias <- studylevel_nocov_bias[sapply(studylevel_nocov_bias, length) > 0]
    eff_nocov <- eff_nocov[sapply(eff_nocov, length) > 0]
    studylevel_nocov_se <- studylevel_nocov_se[sapply(studylevel_nocov_se, length) > 0]
    
    for(i in 1:5){
      studylevel_nocov_avg_bias[i] <- mean(sapply(studylevel_nocov_bias, "[[", i), na.rm = TRUE) 
      
      studylevel_nocov_low_high_bias[i,] <- quantile(sapply(studylevel_nocov_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
      
      avg_nocov_eff[i] <- mean(sapply(eff_nocov, "[[", i))
      
      low_high_nocov_eff[i,] <- quantile(sapply(eff_nocov, "[[", i), probs = c(.025, .975))
      
      studylevel_nocov_mse[i] <- mean(sapply(studylevel_nocov_se, "[[", i))
      
    }
    
    mse_ratio <- round(c(studylevel_mse/ipd_mse),3)
    mse_nocov_ratio <- round(c(studylevel_nocov_mse[1:4]/ipd_mse[1:4]),3)
    mse_nocov_ratio <- c(mse_nocov_ratio[1:4], NA, NA, NA)
    studylevel_nocov_mse <- c(studylevel_nocov_mse[1:4], NA, NA, NA)
    mse <- rbind(studylevel_mse, ipd_mse, studylevel_nocov_mse, mse_ratio, mse_nocov_ratio)
    mse <- cbind(rep("mse", 5), c("studylevel", "ipd", "studylevel_nocov", "ratio_slc_ipd", "ratio_slnc_ipd"), mse)
    
    
    colnames(mse) <- c("stat", "type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "prec")
    mse <- as.data.frame(mse)
    rownames(mse) <- NULL
    
    mse <- do.call(data.frame,lapply(mse, function(x) replace(x, x == "Inf" | x == "NaN",NA)))
    
    # Percent bias figures
    
    Type <- rep(c("studylevel C", "studylevel NC", "ipd"),7)
    
    
    bias <- data.frame(studylevel = c(studylevel_avg_bias), ipd = c(ipd_avg_bias), studylevel_nocov = c(studylevel_nocov_avg_bias[1:4], NA, NA, studylevel_nocov_avg_bias[5]),
                       Beta = c(paste0("B", 1:6, sep = ""), "prec")) %>%
      pivot_longer(cols = c(studylevel, studylevel_nocov, ipd),
                   names_to = "Model", values_to = "Percent Bias") 
    
    bias$Type <- Type
    
    err_bars <- data.frame(rbind(studylevel_low_high_bias, studylevel_nocov_low_high_bias, ipd_low_high_bias))
    names(err_bars) <- c("Low", "High")
    err_bars$Type <- rep(c(rep("studylevel C", 7), rep("studylevel NC", 5), rep("ipd", 7)))
    err_bars$Beta <- c(paste0("B", 1:6, sep = ""), "prec", paste0("B", 1:4, sep = ""), "prec", paste0("B", 1:6, sep = ""), "prec")
    
    bias <- full_join(bias, err_bars, by = c("Type", "Beta"))
    
    bias$Model_Coef <- paste(bias$Beta, bias$Type, sep = "_")

    
    avg_nocov_eff <- c(avg_nocov_eff[1:4], NA, NA, avg_nocov_eff[5])
    avg_eff <- data.frame(type = c("slcipd", "slncipd"), rbind(avg_eff, avg_nocov_eff))
    colnames(avg_eff) <- c("type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "prec")
    rownames(avg_eff) <- NULL
    
    low_high_eff <- data.frame(rbind(low_high_eff, low_high_nocov_eff))
    low_high_eff$Type <- c(rep("slcipd", 7), rep("slncipd", 5))
    colnames(low_high_eff)[1:2] <- c("Low", "High")
    low_high_eff$Beta <- c(c(paste0("B", 1:6, sep=""), "prec"), c(paste0("B", 1:4, sep=""), "prec"))
    
    avg_eff <- pivot_longer(avg_eff, cols = 2:8, names_to = "Beta", values_to = "Efficiency")
    avg_eff$Beta <- c(c(paste0("B", 1:6, sep=""), "prec"), c(paste0("B", 1:6, sep=""), "prec"))
    
    
    eff_dat <- full_join(avg_eff, low_high_eff, by = c("Beta", "type" = "Type")) %>% arrange(Beta)
    eff_dat <- filter(eff_dat, Beta != "prec")

    avg_equiv <- data.frame(d0 = seq(-1, 1, by = .2),
                            lower_studylevel = NA, 
                            upper_studylevel = NA, 
                            median_studylevel = NA,
                            vals_true = equiv[[1]][,5],
                            lower_studylevel_nocov = NA,
                            upper_studylevel_nocov = NA,
                            median_studylevel_nocov = NA,
                            lower_ipd = NA, 
                            upper_ipd = NA, 
                            median_ipd = NA
    )
    
    for(i in c(2:4, 6:8, 9:11)){
      avg_equiv[,i] <- rowMedians(sapply(equiv, "[[", i))
    }
    
    # Coverage probability using credible intervals  
    
    studylevel_coverage <- rep(NA, 7)
    ipd_coverage <- rep(NA, 7)
    studylevel_nocov_coverage <- rep(NA, 5)
    
    for(i in 1:7){
      studylevel_coverage[i] <- mean(sapply(studylevel_cov_prob, "[[", i)) 
      ipd_coverage[i] <- mean(sapply(ipd_cov_prob, "[[", i))
    }
    
    studylevel_nocov_cov_prob <- studylevel_nocov_cov_prob[sapply(studylevel_nocov_cov_prob, length) > 0]
    
    for(i in 1:5){
      studylevel_nocov_coverage[i] <- mean(sapply(studylevel_nocov_cov_prob, "[[", i)) 
    }
    
    studylevel_nocov_coverage <- c(studylevel_nocov_coverage[1:4], NA, NA, studylevel_nocov_coverage[5])
    
    
    
    coverage <- data.frame(type = c("studylevel", "studylevel_nocov", "ipd"), b1 = NA, b2 = NA, b3 = NA, b4 = NA, b5 = NA, b6 = NA, prec = NA)
    coverage[,2:8] <- rbind(studylevel_coverage, studylevel_nocov_coverage, ipd_coverage)  
    
    all_vals <- list("bias" = bias, 
                     "eff_dat" = eff_dat, 
                     "studylevel" = studylevel,
                     "studylevel_nocov" = studylevel_nocov,
                     "ipd" = ipd, 
                     "mad" = mad2, 
                     "mad_ratios" = mad_ratios2,
                     "mse" = mse,
                     "coverage" = coverage,
                     "gr_stats" = gr_stats2,
                     "equiv_dat" = avg_equiv)
    
  }
  
  
  
  return(all_vals)
}

