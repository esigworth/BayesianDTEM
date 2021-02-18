setwd("~/Box/Time-to-Event -- Elizabeth Sigworth PhD Dissertation/Meta-analysis/Simulation/accre_bundles/interaction_term")

source("generate_sim_data_interaction.R")
source("fit_sim_models_accre.R")
source("equiv_functions_interaction.R")

library(kableExtra)
detach(package:plyr)

colors <- c("Median Full Model" = "#00BFC4", 
            "Credible interval Full Model" = "#00BFC4", 
            "Median Meta Model" = "#F8766D",
            "Credible interval Meta Model" = "#F8766D",
            "True relationship" = "#7CAE00")


make_equiv_int <- function(balanced, vals2_all, tau2){
  
  
    bal <- ifelse(balanced, "balanced", "unbalanced")
    main_title <- paste0("Average toxo-equivalence, drug A and B, ", 
                         bal, ", tau2 = ", tau2, collapse = "")
    
    
    colors <- c("Median Full Model - male" = "#0fd960", 
                "Credible interval Full Model - male" = "#0fd960", 
                "Median Full Model - female" = "#eddb11", 
                "Credible interval Full Model - female" = "#eddb11", 
                "Median Meta Model" = "#F8766D",
                "Credible interval Meta Model" = "#F8766D",
                "Median Meta Nocov Model" = "#C77CFF",
                "Credible interval Meta Nocov Model" = "#C77CFF",
                "True relationship - male" = "#080ba3",
                "True relationship - female" = "#a30808")
    
    # Plot equivalence relationship with 95% credible intervals
    fg <- ggplot(aes(x = d0), data = vals2_all) +
      ylim(-3.2, 3.4) +
      geom_line(aes(x = d0, y = vals_true_male, color = "True relationship - male"), 
                data = vals2_all, size = 2, linetype = 1) +
      geom_line(aes(x = d0, y = vals_true_female, color = "True relationship - female"), 
                data = vals2_all, size = 2, linetype = 1) +
      # geom_point(aes(x = d0, y = median_full_male, color = "Median Full Model - male"), 
      #            data = vals2_all, size =2) +
      geom_line(aes(x = d0, y = median_full_male, color = "Median Full Model - male"), 
                data = vals2_all, size = 2, linetype = 2) +
      geom_ribbon(aes(ymin=lower_full_male, ymax = upper_full_male, fill = "Credible interval Full Model - male"), 
                  alpha = .2, data = vals2_all) +
      # geom_point(aes(x = d0, y = median_full_female, color = "Median Full Model - female"), 
      #            data = vals2_all, size =2) +
      geom_line(aes(x = d0, y = median_full_female, color = "Median Full Model - female"), 
                data = vals2_all, size = 2, linetype = 2) +
      geom_ribbon(aes(ymin=lower_full_female, ymax = upper_full_female, fill = "Credible interval Full Model - female"), 
                  alpha = .2, data = vals2_all) +
      # geom_point(aes(x = d0, y = median_meta_nocov, color = "Median Meta Nocov Model"), 
      #            data = vals2_all, size =2) +
      geom_line(aes(x = d0, y = median_meta_nocov, color = "Median Meta Nocov Model"), 
                data = vals2_all, size = 2, linetype = 1) +
      geom_ribbon(aes(ymin=lower_meta_nocov, ymax = upper_meta_nocov, fill = "Credible interval Meta Nocov Model"), 
                  alpha = .2, data = vals2_all) +
      # geom_point(aes(x = d0, y = median_meta, color = "Median Meta Model"), 
      #            data = vals2_all, size =2) +
      geom_line(aes(x = d0, y = median_meta, color = "Median Meta Model"), 
                data = vals2_all, size = 2, linetype = 2) +
      geom_ribbon(aes(ymin=lower_meta, ymax = upper_meta, fill = "Credible interval Meta Model"), 
                  alpha = .1, data = vals2_all) +
      xlab("Delivered drug A") +
      ylab("Delivered drug B") +
      ggtitle(main_title) +
      theme_light() +
      theme(text = element_text(size=12),
            axis.text.x = element_text(size = 10),
            legend.position = c(.5, .95),
            legend.justification = c("left", "top"),
            legend.box.just = "left",
            legend.margin = margin(0, 0, 0, 0),
            legend.spacing = unit(0, 'cm'),
            legend.key = element_rect(colour = NA)
      ) +
      guides(fill=guide_legend(ncol=2),
             color = guide_legend(ncol=2)) +
      labs(fill = NULL,
           color = NULL) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values=colors)
  
   
  
  return(fg)
}


eval_sim_results_inter <- function(all_sims, tau2, total, betas = c(-.6, -.2, .5, .9, .2, .5, -.3, .4), 
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
  
  
  meta <- meta_nocov <- full <- meta_median <- meta_nocov_median <- full_median <- vector(mode = "list", length = total)
  equiv <- mad <- vector(mode = "list", length = total)
  
  
  for(m in 1:total){
    
    meta[[m]] <- all_sims[[m]]$run$meta_summary
    meta_nocov[[m]] <- all_sims[[m]]$run$meta_nocov_summary
    full[[m]] <- all_sims[[m]]$run$full_summary
    
    meta_median[[m]] <- all_sims[[m]]$run$meta_summary_median
    meta_nocov_median[[m]] <- all_sims[[m]]$run$meta_nocov_median
    full_median[[m]] <- all_sims[[m]]$run$full_summary_median
    
    equiv[[m]] <- all_sims[[m]]$run$equiv
    
    mad[[m]] <- all_sims[[m]]$run$mad
    
  
  }


  mad_meta <- mad_full <- mad_ratios <- mad_meta_nocov <- mad_ratios_nocov <- data.frame(b1_mad = rep(NA, total), 
                                                                                         b2_mad = rep(NA, total), 
                                                                                         b3_mad = rep(NA, total), 
                                                                                         b4_mad = rep(NA, total),
                                                                                         b5_mad = rep(NA, total),
                                                                                         b6_mad = rep(NA, total),
                                                                                         b7_mad = rep(NA, total),
                                                                                         b8_mad = rep(NA, total),
                                                                                         type = rep(NA, total))

  
  for(j in 1:total){
    mad_meta[j,] <- mad[[j]][1,]
    mad_meta_nocov[j,] <- mad[[j]][2,]
    mad_full[j,] <- mad[[j]][3,]
    mad_ratios[j,] <- mad[[j]][4,]
    mad_ratios_nocov[j,] <- mad[[j]][5,]
  }
  
  
  mad <- rbind(mad_meta, mad_full, mad_meta_nocov)
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
  
  
  eff <- eff_nocov <- meta_bias <- meta_nocov_bias <- full_bias <- meta_se <- meta_nocov_se <- full_se <- meta_cov_prob <- meta_nocov_cov_prob <- full_cov_prob <- vector(mode = "list", length = total)
  
  betas_n <- c(betas, prec)
  
  for(i in 1:total){
    
    meta_cov_prob[[i]] <- c((betas_n[c(1:6,9)] > meta_median[[i]][,1]) & (betas_n[c(1:6,9)] < meta_median[[i]][,5]))
    meta_nocov_cov_prob[[i]] <- c((betas_n[c(1:4,9)] > meta_nocov_median[[i]][,1]) & (betas_n[c(1:4,9)] < meta_nocov_median[[i]][,5]))
    full_cov_prob[[i]] <- c((betas_n > full_median[[i]][,1]) & (betas_n < full_median[[i]][,5]))
    
    eff[[i]] <- full[[i]][c(1:6,9),2]^2 / meta[[i]][,2]^2
    eff_nocov[[i]] <- full[[i]][c(1:4,9),2]^2 / meta_nocov[[i]][,2]^2
    
    meta_bias[[i]] <- round((meta[[i]][,1] - betas_n[c(1:6,9)])/betas_n[c(1:6,9)], 2)
    meta_nocov_bias[[i]] <- round((meta_nocov[[i]][,1] - betas_n[c(1:4,9)])/betas_n[c(1:4,9)], 2)
    full_bias[[i]] <- round((full[[i]][,1] - betas_n)/betas_n, 2)
    
    meta_se[[i]] <- round((meta[[i]][,1] - betas_n[c(1:6,9)])^2, 2)
    meta_nocov_se[[i]] <- round((meta_nocov[[i]][,1] - betas_n[c(1:4,9)])^2, 2)
    full_se[[i]] <- round((full[[i]][,1] - betas_n)^2, 2)
  }
  
  meta_beta <- meta_nocov_beta <- full_beta <- data.frame(b1=rep(NA, total), 
                                                          b2=rep(NA, total), 
                                                          b3=rep(NA, total), 
                                                          b4=rep(NA, total),
                                                          b5=rep(NA, total),
                                                          b6=rep(NA, total),
                                                          b7=rep(NA, total),
                                                          b8=rep(NA, total))
  
  for(i in 1:total){
    meta_beta[i,] <- c(meta[[i]][1:6,1], NA, NA)
    meta_nocov_beta[i,] <- c(meta_nocov[[i]][1:4,1], NA, NA, NA, NA)
    full_beta[i,] <- full[[i]][1:8,1]
  }
  
  median_meta <- apply(meta_beta, 2, median)
  median_meta_nocov <- apply(meta_nocov_beta, 2, median)
  median_full <- apply(full_beta, 2, median)
  
  d0 <- seq(-1, 1, by = .2)
  
 
  vals_meta <- (median_meta[3]*d0 - median_meta[2])/median_meta[4]
  vals_meta_nocov <- (median_meta_nocov[3]*d0 - median_meta_nocov[2])/median_meta[4]
  vals_full_female <- ((median_full[3]+median_full[7])*d0 - median_full[2])/(median_full[4]+median_full[8])
  vals_full_male <- ((median_full[3])*d0 - median_full[2])/(median_full[4])
  vals_true_female <- ((betas[3]+betas[7])*d0 - betas[2])/(betas[4]+betas[8])
  vals_true_male <- ((betas[3])*d0 - betas[2])/(betas[4])
  
  median_param_equiv <- cbind(d0, vals_meta, vals_meta_nocov, vals_full_female, vals_full_male, vals_true_female, vals_true_male)
  
  
  
  
  meta_avg_bias <- rep(NA, 7)
  meta_nocov_avg_bias <- rep(NA, 5)
  full_avg_bias <- rep(NA, 9)
  
  meta_low_high_bias <- matrix(NA, ncol = 2, nrow = 7)
  meta_nocov_low_high_bias <- matrix(NA, ncol = 2, nrow = 5)
  full_low_high_bias <- matrix(NA, ncol = 2, nrow = 9)
  
  avg_eff <- rep(NA,7)
  avg_nocov_eff <- rep(NA, 5)
  low_high_eff <- matrix(NA, ncol = 2, nrow = 7)
  low_high_nocov_eff <- matrix(NA, ncol = 2, nrow = 5)
  
  meta_mse <- full_mse <- rep(NA, 7)
  meta_nocov_mse <- rep(NA, 5)
  
  for(i in 1:9){
    full_avg_bias[i] <- mean(sapply(full_bias, "[[", i), na.rm = TRUE)
    
    full_low_high_bias[i,] <- quantile(sapply(full_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    full_mse[i] <- mean(sapply(full_se, "[[", i))
    
  }
  
  for(i in 1:7){
    meta_avg_bias[i] <- mean(sapply(meta_bias, "[[", i), na.rm = TRUE) 
    
    meta_low_high_bias[i,] <- quantile(sapply(meta_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    avg_eff[i] <- mean(sapply(eff, "[[", i))
    
    low_high_eff[i,] <- quantile(sapply(eff, "[[", i), probs = c(.025, .975))
    
    meta_mse[i] <- mean(sapply(meta_se, "[[", i))
    
  }

  
  for(i in 1:5){
    meta_nocov_avg_bias[i] <- mean(sapply(meta_nocov_bias, "[[", i), na.rm = TRUE) 
    
    meta_nocov_low_high_bias[i,] <- quantile(sapply(meta_nocov_bias, "[[", i), probs = c(.025, .975), na.rm = TRUE)
    
    avg_nocov_eff[i] <- mean(sapply(eff_nocov, "[[", i))
    
    low_high_nocov_eff[i,] <- quantile(sapply(eff_nocov, "[[", i), probs = c(.025, .975))
    
    meta_nocov_mse[i] <- mean(sapply(meta_nocov_se, "[[", i))
    
  }
  
  mse_ratio <- round(c(meta_mse[1:6]/full_mse[1:6]),3)
  mse_ratio <-c(mse_ratio, NA, NA)
  mse_nocov_ratio <- round(c(meta_nocov_mse[1:4]/full_mse[1:4]),3)
  mse_nocov_ratio <- c(mse_nocov_ratio[1:4], NA, NA, NA, NA)
  meta_nocov_mse <- c(meta_nocov_mse[1:4], NA, NA, NA, NA)
  meta_mse <- c(meta_mse[1:6], NA, NA)
  mse <- data.frame(rbind(meta_mse, full_mse[1:8], meta_nocov_mse, mse_ratio, mse_nocov_ratio))
  mse <- data.frame(stat = "mse", type = c("meta", "full", "meta_nocov", "ratio_mf", "ratio_mncf"), mse)
  
  colnames(mse) <- c("stat", "type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8")
  mse <- as.data.frame(mse)
  rownames(mse) <- NULL
  
  mse <- do.call(data.frame,lapply(mse, function(x) replace(x, x == "Inf" | x == "NaN",NA)))
  
  # Percent bias figures
  
  Type <- rep(c( "Meta NC", "Meta", "Full"),8)
  
  
  bias <- data.frame(meta = c(meta_avg_bias[1:6], NA, NA), 
                     full = c(full_avg_bias[1:8]), 
                     meta_nocov = c(meta_nocov_avg_bias[1:4], NA, NA, NA, NA),
                     Beta = c(paste0("B", 1:8, sep = ""))) %>%
    pivot_longer(cols = c(meta, meta_nocov, full),
                 names_to = "Model", values_to = "Percent Bias") 
  
  bias$Type <- Type
  
  err_bars <- data.frame(rbind(meta_low_high_bias[1:6,], meta_nocov_low_high_bias[1:4,], full_low_high_bias[1:8,]))
  names(err_bars) <- c("Low", "High")
  err_bars$Type <- rep(c(rep("Meta", 6), rep("Meta NC", 4), rep("Full", 8)))
  err_bars$Beta <- c(paste0("B", 1:6, sep = ""), paste0("B", 1:4, sep = ""), paste0("B", 1:8, sep = ""))
  
  bias <- full_join(bias, err_bars, by = c("Type", "Beta")) 
  
  
  bias$Model_Coef <- paste(bias$Beta, bias$Type, sep = "_")
  
  fig_dat <- filter(bias, Beta != "B1" & Beta != "B6" & Beta != "B5")
  
  fig_dat$pos <- c(3.25, 3.5, 3.75, 5.25, 5.5, 5.75, 7.25, 7.5, 7.75, 9.25, 9.5, 9.75, 11.25, 11.5, 11.75)
  
  min_pb <- min(fig_dat$Low) - .1
  max_pb <- max(fig_dat$High) + .1
  
  avg_nocov_eff <- c(avg_nocov_eff[1:4], NA, NA)
  avg_eff <- data.frame(type = c("MF", "MNCF"), rbind(avg_eff[1:6], avg_nocov_eff))
  colnames(avg_eff) <- c("type", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6")
  rownames(avg_eff) <- NULL
  
  low_high_eff <- data.frame(rbind(low_high_eff[1:6,], low_high_nocov_eff[1:4,]))
  low_high_eff$Type <- c(rep("MF", 6), rep("MNCF", 4))
  colnames(low_high_eff)[1:2] <- c("Low", "High")
  low_high_eff$Beta <- c(c(paste0("B", 1:6, sep="")), c(paste0("B", 1:4, sep="")))
  
  avg_eff <- pivot_longer(avg_eff, cols = 2:7, names_to = "Beta", values_to = "Efficiency")
  avg_eff$Beta <- c(c(paste0("B", 1:6, sep="")), c(paste0("B", 1:6, sep="")))
  
  
  eff_dat <- full_join(avg_eff, low_high_eff, by = c("Beta", "type" = "Type")) %>% arrange(Beta)
  eff_dat$pos <- c(1.25, 1.75, 3.25, 3.75, 5.25, 5.75, 7.25, 7.75, 13.25, 13.75, 15.25, 15.75)
  
  
  
  
  eff_dat_trial <- filter(eff_dat, Beta != "B5" & Beta != "B6")
  eff_dat_extra <- filter(eff_dat, (Beta == "B5" | Beta == "B6") & type == "MF")
  
  
  avg_equiv <- data.frame(d0 = seq(-1, 1, by = .2),
                          lower_meta = NA,
                          upper_meta = NA,
                          median_meta = NA,
                          vals_true_female = vals_true_female,
                          lower_full_female = NA, 
                          upper_full_female = NA, 
                          median_full_female = NA,
                          lower_full_male = NA, 
                          upper_full_male = NA, 
                          median_full_male = NA,
                          lower_meta_nocov = NA,
                          upper_meta_nocov = NA,
                          median_meta_nocov = NA,
                          vals_true_male = vals_true_male
  )
  
  
  for(i in c(2:4, 6:14)){
    avg_equiv[,i] <- rowMedians(sapply(equiv, "[[", i))
  }
  
  equiv_fig <- make_equiv_int((bal == "balanced"), avg_equiv, tau2)
  
  
  # Coverage probability using credible intervals  
  
  meta_coverage <- rep(NA, 6)
  full_coverage <- rep(NA, 8)
  meta_nocov_coverage <- rep(NA, 4)
  
  for(i in 1:8){
    full_coverage[i] <- mean(sapply(full_cov_prob, "[[", i))
  }
  
  for(i in 1:6){
    meta_coverage[i] <- mean(sapply(meta_cov_prob, "[[", i)) 
  }
  
  for(i in 1:4){
    meta_nocov_coverage[i] <- mean(sapply(meta_nocov_cov_prob, "[[", i)) 
  }
  
  meta_nocov_coverage <- c(meta_nocov_coverage[1:4], NA, NA, NA, NA)
  meta_coverage <- c(meta_coverage, NA, NA)
  
  
  
  coverage <- data.frame(type = c("meta", "meta_nocov", "full"), b1 = NA, b2 = NA, b3 = NA, b4 = NA, b5 = NA, b6 = NA, b7 = NA, b8 = NA)
  coverage[,2:9] <- rbind(meta_coverage, meta_nocov_coverage, full_coverage)  
  
  all_vals <- list("bias" = bias, 
                   "eff_dat" = eff_dat, 
                   "meta" = meta, 
                   "full" = full, 
                   "mad" = mad2, 
                   "mad_ratios" = mad_ratios2,
                   "mse" = mse,
                   "coverage" = coverage,
                   "equiv_dat" = avg_equiv,
                   "median_equiv" = median_param_equiv)

  
  
  
  return(all_vals)
}


# un_0 <- eval_sim_results(un_0_all, tau2 = 0, total = 1000, bal = "unbalanced", covars = "no covars")
# un_half <- eval_sim_results(un_half_all, tau2 = .5, total = 1000, bal = "unbalanced", covars = "no covars")
# un_1 <- eval_sim_results(un_1_all, tau2 = 1, total = 1000, bal = "unbalanced", covars = "no covars")
# 
# bn_0 <- eval_sim_results(bn_0_all, tau2 = 0, total = 1000, bal = "balanced", covars = "no covars")
# bn_half <- eval_sim_results(bn_half_all, tau2 = 0.5, total = 999, bal = "balanced", covars = "no covars")
# bn_1 <- eval_sim_results(bn_1_all, tau2 = 1, total = 1000, bal = "balanced", covars = "no covars")
# 
# bc_0 <- eval_sim_results(bc_0_all, tau2 = 0, total = 250, bal = "balanced", covars = "with covars")
# bc_half <- eval_sim_results(bc_half_all, tau2 = .5, total = 250, bal = "balanced", covars = "with covars")
# bc_1 <- eval_sim_results(bc_1_all, tau2 = 1, total = 248, bal = "balanced", covars = "with covars")
# 
# uc_0 <- eval_sim_results(uc_0_all, tau2 = 0, total = 250, bal = "unbalanced", covars = "with covars")
# uc_half <- eval_sim_results(uc_half_all, tau2 = 0.5, total = 250, bal = "unbalanced", covars = "with covars")
# uc_1 <- eval_sim_results(uc_1_all, tau2 = 1, total = 249, bal = "unbalanced", covars = "with covars")
# 
# 
# bc_chi <- eval_sim_results(bc_0_chi, tau2 = 0, total = 147, bal = "balanced", covars = "with covars")


lbc_0_interaction <- create_full("lbc_0_interaction")
lbc_half_interaction <- create_full("lbc_half_interaction")
lbc_1_interaction <- create_full("lbc_1_interaction")

luc_0_interaction <- create_full("luc_0_interaction")
luc_half_interaction <- create_full("luc_half_interaction")
luc_1_interaction <- create_full("luc_1_interaction")

save(lbc_0_interaction, file = "lbc_0_interaction.Rda")
save(lbc_half_interaction, file = "lbc_half_interaction.Rda")
save(lbc_1_interaction, file = "lbc_1_interaction.Rda")

save(luc_0_interaction, file = "luc_0_interaction.Rda")
save(luc_half_interaction, file = "luc_half_interaction.Rda")
save(luc_1_interaction, file = "luc_1_interaction.Rda")


load("lbc_0_interaction.Rda")
load("lbc_half_interaction.Rda")
load("lbc_1_interaction.Rda")
load("luc_0_interaction.Rda")
load("luc_half_interaction.Rda")
load("luc_1_interaction.Rda")



# Create summary figures in R

bc_0 <- eval_sim_results_inter(lbc_0_interaction, tau2 = 0, total = 100, bal = "balanced", covars = "with covars")
bc_half <- eval_sim_results_inter(lbc_half_interaction, tau2 = 0.5, total = 99, bal = "balanced", covars = "with covars")
bc_1 <- eval_sim_results_inter(lbc_1_interaction, tau2 = 1, total = 98, bal = "balanced", covars = "with covars")

uc_0 <- eval_sim_results_inter(luc_0_interaction, tau2 = 0, total = 100, bal = "unbalanced", covars = "with covars")
uc_half <- eval_sim_results_inter(luc_half_interaction, tau2 = 0.5, total = 96, bal = "unbalanced", covars = "with covars")
uc_1 <- eval_sim_results_inter(luc_1_interaction, tau2 = 1, total = 96, bal = "unbalanced", covars = "with covars")



  expr2 <- expression(beta[2])
  expr3 <- expression(beta[3])
  expr4 <- expression(beta[4])
  expr5 <- expression(beta[5])
  expr6 <- expression(beta[6])
  
  
  
  bc_bias <- rbind(bc_0$bias, bc_half$bias, bc_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 24),
           fit = "Balanced") %>%
    mutate(Type = case_when(Type == "Meta NC" ~ "SL-NC",
                            Type == "Meta" ~ "SL-C",
                            Type == "Full" ~ "IPD",
                            TRUE ~ "IPD")) 
  uc_bias <- rbind(uc_0$bias, uc_half$bias, uc_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 24),
           fit = "Unbalanced") %>%
    mutate(Type = case_when(Type == "Meta NC" ~ "SL-NC",
                            Type == "Meta" ~ "SL-C",
                            Type == "Full" ~ "IPD",
                            TRUE ~ "IPD")) 
  
  all_bias <- rbind(bc_bias, uc_bias) %>% 
    filter(Beta %in% c("B2", "B3", "B4"))
  
  fig_dat <- all_bias

  fig_dat_withcov <- fig_dat

  colrs <- c("IPD" = "#fa0fe2", 
             "SL-C" = "#88C6ED",
             "SL-NC" = "#FAA31B")


  
  fig_dat_withcov$pos <- rep(c(1.8, 1.5, 1.2, 3.8, 3.5, 3.2, 5.8, 5.5, 5.2), 6)
  
  bias_min_withcov <- min(-1, floor(min(fig_dat_withcov$Low, na.rm=T)*10)/10)
  bias_max_withcov <- max(1,ceiling(max(fig_dat_withcov$High,na.rm=T)*10)/10)
  
  
  ggplot(data = fig_dat_withcov, mapping = aes(x = pos, 
                                                                y = `Percent Bias`, 
                                                                color = Type)) +
    coord_cartesian(ylim=c(bias_min_withcov, bias_max_withcov), expand = FALSE, clip = "off", xlim = c(0.5, 6.5)) +
    scale_x_continuous(breaks = seq(1.5, 5.5, by = 2)) +
    geom_hline(yintercept = 0, color = "grey", size = .7) +
    geom_point(size = 2.5, stroke = 1) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = .3, size = 1.1) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size=12),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.margin = unit(c(1,1,4,1), "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = c(0.16, 0.03),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black")) + 
    scale_color_manual(values = colrs) +
    labs(color = "Model Fit") + 
    geom_text(x = 1.5, y = bias_min_withcov - .425, label = expr2, size = 5, color = "black") +
    geom_text(x = 3.5, y = bias_min_withcov - .425, label = expr3, size = 5, color = "black") +
    geom_text(x = 5.5, y = bias_min_withcov - .425, label = expr4, size = 5, color = "black") +
    facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2))) +
    ylab("Percent Bias") 
  
  
  bc_eff <- rbind(bc_0$eff_dat, bc_half$eff_dat, bc_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 12),
           fit = "Balanced",
           type = ifelse(type == "MF", "SL-C:IPD", "SL-NC:IPD")) %>%
    dplyr::select(-pos) 
    relocate(Efficiency, Beta, Low, High, tau2, fit, type)
  uc_eff <- rbind(uc_0$eff_dat, uc_half$eff_dat, uc_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 12),
           fit = "Unbalanced",
           type = ifelse(type == "MF", "SL-C:IPD", "SL-NC:IPD")) %>%
    dplyr::select(-pos) %>%
    relocate(Efficiency, Beta, Low, High, tau2, fit, type)
  
  all_eff <- rbind(bc_eff, uc_eff) %>%
    arrange(tau2, fit)
  
  eff_dat_nocov <- all_eff %>% filter(Beta %in% c("B2", "B3", "B4"))
  eff_dat_cov <- all_eff %>% filter((Beta == "B5" | Beta == "B6") & !is.na(Efficiency))
  
  eff_dat2 <- eff_dat_nocov
  
  eff_min2 <- floor(min(eff_dat2$Low,na.rm=T)*10)/10
  eff_max2 <- ceiling(max(eff_dat2$High, na.rm=T)*10)/10
  
  eff_dat2$pos <- rep(c(.8, 1.2, 1.8, 2.2, 2.8, 3.2), 6)
  
  colrs <- c("SL-C:IPD" = "#88C6ED",
             "SL-NC:IPD" = "#FAA31B")
  
  ggplot(data = eff_dat2, mapping = aes(x = pos, y = Efficiency, color = type)) +
    coord_cartesian(ylim=c(eff_min2, eff_max2), expand = FALSE, clip = "off", xlim = c(0, 4)) +
    geom_hline(yintercept = 1, color = "grey", size = .7) +
    geom_point(size = 2.5, stroke = 1) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = .1, size = 1.1) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          axis.text.y = element_text(size = 14),
          plot.margin = unit(c(1,1,4,1), "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = c(0.16, 0.03),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black")) + 
    scale_color_manual(values = colrs) +
    ylab("Relative efficiency") +
    labs(color = "Ratio") + 
    geom_text(x = 1, y = eff_min2 - .4, label = expr2, size = 5, colour = "black") +
    geom_text(x = 2, y = eff_min2 - .4, label = expr3, size = 5, colour = "black") +
    geom_text(x = 3, y = eff_min2 - .4, label = expr4, size = 5, colour = "black") + 
    facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2)))
  
  eff_min_cov <- floor(min(eff_dat_cov$Low)*100)/100
  eff_max_cov <- ceiling(max(eff_dat_cov$High)*100)/100
  
  ggplot(data = eff_dat_cov, mapping = aes(x = Beta, y = Efficiency)) +
    coord_cartesian(ylim=c(eff_min_cov-.005, eff_max_cov), 
                    expand = FALSE, clip = "off", xlim = c(0.5, 2.5)) +
    geom_hline(yintercept = 1, color = "grey", size = .7) +
    geom_point(size = 2.5, stroke = 1, color = "#88C6ED") +
    geom_errorbar(aes(ymin = Low, ymax = High), width = .1, size = 1.1, color = "#88C6ED") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          axis.text.y = element_text(size = 14),
          plot.margin = unit(c(1,1,4,1), "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none",
          panel.spacing = unit(1, "lines"),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black")) + 
    ylab("Relative efficiency") +
    geom_text(x = 1, y = eff_min_cov - .015, label = expr5, size = 5) +
    geom_text(x = 2, y = eff_min_cov - .015, label = expr6, size = 5) +
    facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2)))

  
  colors <- c("Median IPD Model - male" = "#FFF000", 
              "Credible interval IPD Model - male" = "#FFF000", 
              "Median IPD Model - female" = "#009F75", 
              "Credible interval IPD Model - female" = "#009F75", 
              "Median SL-C Model" = "#88C6ED",
              "Credible interval SL-C Model" = "#88C6ED",
              "Median SL-NC Model" = "#FAA31B",
              "Credible interval SL-NC Model" = "#FAA31B",
              "True relationship - male" = "#394BA0",
              "True relationship - female" = "#EF4444")
  
  bc_equiv <- rbind(bc_0$equiv_dat, bc_half$equiv_dat, bc_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Balanced")
  uc_equiv <- rbind(uc_0$equiv_dat, uc_half$equiv_dat, uc_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Unbalanced")
  
  equiv_withcov <- rbind(bc_equiv, uc_equiv)
  
  
  # Plot equivalence relationship with 95% credible intervals
  ggplot(aes(x = d0), data = equiv_withcov) +
    ylim(-.7, 1) +
    geom_line(aes(x = d0, y = vals_true_male, color = "True relationship - male"), 
              data = equiv_withcov, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = vals_true_female, color = "True relationship - female"), 
              data = equiv_withcov, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = median_full_male, color = "Median IPD Model - male"), 
              data = equiv_withcov, size = 2, linetype = 2) +
    geom_line(aes(x = d0, y = median_full_female, color = "Median IPD Model - female"), 
              data = equiv_withcov, size = 2, linetype = 2) +
    geom_line(aes(x = d0, y = median_meta_nocov, color = "Median SL-NC Model"), 
              data = equiv_withcov, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = median_meta, color = "Median SL-C Model"), 
              data = equiv_withcov, size = 2, linetype = 2) +
    xlab("Delivered drug A") +
    ylab("Delivered drug B") +
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "bottom",
          legend.justification = c("center"),
          legend.box="vertical",
          legend.box.just = "center",
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing = unit(0, 'cm'),
          legend.key = element_rect(colour = NA),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black")
    ) +
    labs(fill = NULL,
         color = NULL) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values=colors) + 
    facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2)))
  
  
  colors <- c("Median Mega Model - male" = "#FFF000", 
              "Credible interval Mega Model - male" = "#FFF000", 
              "Median Mega Model - female" = "#009F75", 
              "Credible interval Mega Model - female" = "#009F75", 
              "Median Meta-C Model" = "#88C6ED",
              "Credible interval Meta-C Model" = "#88C6ED",
              "Median Meta-NC Model" = "#FAA31B",
              "Credible interval Meta-NC Model" = "#FAA31B",
              "True relationship - male" = "#394BA0",
              "True relationship - female" = "#EF4444")
  
  bc_equiv_med <- as_tibble(rbind(bc_0$median_equiv, bc_half$median_equiv, bc_1$median_equiv)) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Balanced")
  uc_equiv_med <- as_tibble(rbind(uc_0$median_equiv, uc_half$median_equiv, uc_1$median_equiv)) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Unbalanced")
  
  equiv_med <- rbind(bc_equiv_med, uc_equiv_med)
  
  
  # Plot equivalence relationship with 95% credible intervals
  ggplot(aes(x = d0), data = equiv_med) +
    ylim(-.5, 1) +
    geom_line(aes(x = d0, y = vals_true_male, color = "True relationship - male"), 
              data = equiv_med, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = vals_true_female, color = "True relationship - female"), 
              data = equiv_med, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = vals_full_male, color = "Median Mega Model - male"), 
              data = equiv_med, size = 2, linetype = 2) +
    geom_line(aes(x = d0, y = vals_full_female, color = "Median Mega Model - female"), 
              data = equiv_med, size = 2, linetype = 2) +
    geom_line(aes(x = d0, y = vals_meta_nocov, color = "Median Meta-NC Model"), 
              data = equiv_med, size = 2, linetype = 1) +
    geom_line(aes(x = d0, y = vals_meta, color = "Median Meta-C Model"), 
              data = equiv_med, size = 2, linetype = 2) +
    xlab("Delivered drug A") +
    ylab("Delivered drug B") +
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "bottom",
          legend.justification = c("center"),
          legend.box="vertical",
          legend.box.just = "center",
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing = unit(0, 'cm'),
          legend.key = element_rect(colour = NA),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size = 14, colour = "black")
    ) +
    labs(fill = NULL,
         color = NULL) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values=colors) + 
    facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2)))
  
  
  
 
  return(list("bias_plot" = bias_plot_withcov,
              "eff_plot" = eff_plot2,
              "eff_plot_extra" = eff_cov_fig,
              "equiv_plot_male" = male,
              "equiv_plot_female" = female))
  


# Create summary tables


latex_tabs <- function(bc_0, bc_half, bc_1, 
                       uc_0, uc_half, uc_1){
  
  
  # Put together all coverage values

  bc_cov <- rbind(bc_0$coverage, bc_half$coverage, bc_1$coverage) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 3),
           fit = "Balanced, With Covars") %>%
    dplyr::select(-c(b1, b5, b6, b7, b8)) %>%
    pivot_longer(cols = c(b2, b3, b4),
                 names_to = "parameter",
                 values_to = "coverage")
  uc_cov <- rbind(uc_0$coverage, uc_half$coverage, uc_1$coverage) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 3),
           fit = "Unbalanced, With Covars") %>%
    dplyr::select(-c(b1, b5, b6, b7, b8)) %>%
    pivot_longer(cols = c(b2, b3, b4),
                 names_to = "parameter",
                 values_to = "coverage")
  
  # Put together all MAD ratios

  bc_mad <- rbind(bc_0$mad_ratios, bc_half$mad_ratios, bc_1$mad_ratios) %>% 
    filter(param %in% c("B2", "B3", "B4")) %>%
    mutate(tau2 = rep(c(0, .5, 1), each = 3),
           fit = "Balanced, With Covars",
           width = high - low) %>%
    dplyr::select(-c(low, high))
  uc_mad <- rbind(uc_0$mad_ratios, uc_half$mad_ratios, uc_1$mad_ratios) %>% 
    filter(param %in% c("B2", "B3", "B4")) %>%
    mutate(tau2 = rep(c(0, .5, 1), each = 3),
           fit = "Unbalanced, With Covars",
           width = high - low) %>%
    dplyr::select(-c(low, high))
  
  all_cov <- rbind(bc_cov, uc_cov) %>%
    pivot_wider(id_cols = c(tau2, fit, parameter), names_from = type, values_from = coverage) %>%
    pivot_wider(id_cols = c(fit, parameter), names_from = tau2, values_from = c(meta, meta_nocov, full)) %>%
    dplyr::select(fit, parameter, meta_0, meta_nocov_0, full_0, meta_0.5, meta_nocov_0.5, full_0.5, meta_1, meta_nocov_1, full_1) %>%
    mutate_at(3:11, round, 2)
  
  all_mad <- rbind(bc_mad, uc_mad) %>%
    pivot_wider(id_cols = c(tau2, fit, param), names_from = type, values_from = c(med, width)) %>%
    pivot_wider(id_cols = c(fit, param), names_from = tau2, values_from = c(med_ratio_meta_full, med_ratio_meta_nocov_full, width_ratio_meta_full, width_ratio_meta_nocov_full)) %>%
    mutate(parameter = tolower(param)) %>%
    dplyr::select(-param) %>%
    mutate_at(2:13, round, 2) %>%
    rename(med_mf_0 = med_ratio_meta_full_0,
           med_mf_half = med_ratio_meta_full_0.5,
           med_mf_1 = med_ratio_meta_full_1,
           width_mf_0 = width_ratio_meta_full_0,
           width_mf_half = width_ratio_meta_full_0.5,
           width_mf_1 = width_ratio_meta_full_1,
           med_mncf_0 = med_ratio_meta_nocov_full_0,
           med_mncf_half = med_ratio_meta_nocov_full_0.5,
           med_mncf_1 = med_ratio_meta_nocov_full_1,
           width_mncf_0 = width_ratio_meta_nocov_full_0,
           width_mncf_half = width_ratio_meta_nocov_full_0.5,
           width_mncf_1 = width_ratio_meta_nocov_full_1)
  
  all_vals <- inner_join(all_cov, all_mad, by = c("fit" = "fit", "parameter" = "parameter")) %>%
    dplyr::select(fit, parameter, meta_0, meta_nocov_0, full_0, med_mf_0, med_mncf_0, width_mf_0, width_mncf_0, 
                  meta_0.5, meta_nocov_0.5, full_0.5, med_mf_half, med_mncf_half, width_mf_half, width_mncf_half,
                  meta_1, meta_nocov_1, full_1, med_mf_1, med_mncf_1, width_mf_1, width_mncf_1)
  
  kbl(all_vals[,2:23], booktabs = TRUE, col.names = c("Fit", rep(c("M", "MNC", "F", "MF", "MNCF", "MF", "MNCF"), 3)), format = "html", align = "c", digits = 2, caption = "Coverage, Median Absolute Deviation Ratios (MADR), and 95% credible interval width (CIW) of MADR", padding=2L) %>%
    kable_styling("striped") %>%
    add_header_above(c("Statistic" = 1,
                       "Coverage" = 3,
                       "MAD ratio" = 2,
                       "MAD ratio CIW" = 2,
                       "Coverage" = 3,
                       "MAD ratio" = 2,
                       "MAD ratio CIW" = 2,
                       "Coverage" = 3,
                       "MAD ratio" = 2,
                       "MAD ratio CIW" = 2)) %>%
    add_header_above(c(`$T^2$` = 1, 
                       "0" = 7, 
                       "0.5" = 7,
                       "1" = 7)) %>%
    pack_rows("Balanced, With Covars", 1, 3) %>%
    pack_rows("Unbalanced, With Covars", 4, 6) 
}



kbl(MCMCsummary(samples[,1:7]), booktabs = TRUE, format = "latex", digits = 2, caption = "Mean, standard deviation, and quantile estimates for the final model fit to the taxane data, along with the Gelman-Rubin (GR) statistics and effective sample size (ESS) from the sampling chain of each parameter.") %>%
  kable_styling(latex_options = "scale_down") 












summarize_setting <- function(bn_0, bn_half, bn_1, 
                              un_0, un_half, un_1, 
                              bc_0, bc_half, bc_1, 
                              uc_0, uc_half, uc_1,
                              type, bal, covars){
  
  
  expr1 <- expression(beta[1])
  expr2 <- expression(beta[2])
  expr3 <- expression(beta[3])
  expr4 <- expression(beta[4])
  expr5 <- expression(beta[5])
  expr6 <- expression(beta[6])
  
  tau_labs <- c(expression(T^2==0), expression(T^2==0.5), expression(T^2==1))
  names(tau_labs) <- c("0", "0.5", "1")
  
  bn_bias <- rbind(bn_0$bias, bn_half$bias, bn_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 10),
           fit = "Balanced, No Covars")
  un_bias <- rbind(un_0$bias, un_half$bias, un_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 10),
           fit = "Unbalanced, No Covars")
  bc_bias <- rbind(bc_0$bias, bc_half$bias, bc_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 14),
           fit = "Balanced, With Covars")
  uc_bias <- rbind(uc_0$bias, uc_half$bias, uc_1$bias) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 14),
           fit = "Unbalanced, With Covars")
  
  all_bias <- rbind(bn_bias, un_bias, bc_bias, uc_bias) %>% 
    filter(Beta != "prec" & Beta != "B5" & Beta != "B6")
  
  fig_dat <- all_bias
  
  fig_dat$pos <- rep(c(1.75, 1.25, 3.75, 3.25, 5.75, 5.25, 7.75, 7.25), 12)
  
  bias_min <- min(-1, floor(min(fig_dat$Low)*10)/10)
  bias_max <- max(1,ceiling(max(fig_dat$High)*10)/10)
  
  bias_plot <- ggplot(data = fig_dat, mapping = aes(x = pos, 
                                                    y = `Percent Bias`, 
                                                    color = Type)) +
    coord_cartesian(ylim=c(bias_min, bias_max), expand = FALSE, clip = "off", xlim = c(0.5, 8.5)) +
    scale_x_continuous(breaks = seq(1.5, 7.5, by = 2)) +
    geom_hline(yintercept = 0, color = "grey", size = .7) +
    geom_point(size = 2.5, stroke = 1) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = .3, size = 1.1) +
    theme_light() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(1,1,4,1), "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = c(0.12, 0.045),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.key = element_blank()) + 
    labs(color = "Model Fit") + 
    geom_text(x = 1.5, y = bias_min - .3, label = expr1, size = 5, color = "black") +
    geom_text(x = 3.5, y = bias_min - .3, label = expr2, size = 5, color = "black") +
    geom_text(x = 5.5, y = bias_min - .3, label = expr3, size = 5, color = "black") +
    geom_text(x = 7.5, y = bias_min - .3, label = expr4, size = 5, color = "black") +
    facet_grid(fit ~ tau2, labeller = labeller(tau2 = tau_labs)) +
    ylab("Percent Bias") 
  
  
  bn_eff <- rbind(bn_0$eff_dat, bn_half$eff_dat, bn_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 4),
           fit = "Balanced, No Covars")
  un_eff <- rbind(un_0$eff_dat, un_half$eff_dat, un_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 4),
           fit = "Unbalanced, No Covars")
  bc_eff <- rbind(bc_0$eff_dat, bc_half$eff_dat, bc_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 6),
           fit = "Balanced, With Covars")
  uc_eff <- rbind(uc_0$eff_dat, uc_half$eff_dat, uc_1$eff_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 6),
           fit = "Unbalanced, With Covars")
  
  all_eff <- rbind(bn_eff, un_eff, bc_eff, uc_eff) 
  
  eff_dat_trial <- filter(all_eff, Beta != "B5" & Beta != "B6")
  
  eff_min <- floor(min(eff_dat_trial$Low)*10)/10
  eff_max <- ceiling(max(eff_dat_trial$High)*10)/10
  
  eff_plot <- ggplot(data = eff_dat_trial, mapping = aes(x = Beta, y = Efficiency)) +
    coord_cartesian(ylim=c(eff_min, eff_max), expand = FALSE, clip = "off", xlim = c(0.5, 4.5)) +
    geom_hline(yintercept = 1, color = "grey", size = .7) +
    geom_point(size = 2.5, stroke = 1) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = .2) +
    theme_light() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(1,1,4,1), "lines"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none") + 
    ylab("Relative efficiency") +
    geom_text(x = 1, y = eff_min - .04, label = expr1, size = 5) +
    geom_text(x = 2, y = eff_min - .04, label = expr2, size = 5) +
    geom_text(x = 3, y = eff_min - .04, label = expr3, size = 5) +
    geom_text(x = 4, y = eff_min - .04, label = expr4, size = 5) + 
    facet_grid(fit ~ tau2)
  
  
  colors <- c("Median Full Model" = "#00BFC4", 
              "Credible interval Full Model" = "#00BFC4", 
              "Median Meta Model" = "#C41467",
              "Credible interval Meta Model" = "#C41467",
              "True relationship" = "#67C414")
  
  bn_equiv <- rbind(bn_0$equiv_dat, bn_half$equiv_dat, bn_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Balanced, No Covars")
  un_equiv <- rbind(un_0$equiv_dat, un_half$equiv_dat, un_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Unbalanced, No Covars")
  bc_equiv <- rbind(bc_0$equiv_dat, bc_half$equiv_dat, bc_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Balanced, With Covars")
  uc_equiv <- rbind(uc_0$equiv_dat, uc_half$equiv_dat, uc_1$equiv_dat) %>% 
    mutate(tau2 = rep(c(0, .5, 1), each = 11),
           fit = "Unbalanced, With Covars")
  
  all_equiv <- rbind(bn_equiv, un_equiv, bc_equiv, uc_equiv) 
  
  vals2_all <- all_equiv
  
  # Plot equivalence relationship with 95% credible intervals
  fg <- ggplot(aes(x = d0), data = vals2_all) +
    ylim(-1.5, 2.2) +
    geom_point(aes(x = d0, y = median_full, color = "Median Full Model"), 
               data = vals2_all, size =2) +
    geom_line(aes(x = d0, y = median_full, color = "Median Full Model"), 
              data = vals2_all, size = 1, linetype = 4) +
    geom_ribbon(aes(ymin=lower_full, ymax = upper_full, fill = "Credible interval Full Model"), 
                alpha = .2, data = vals2_all) +
    geom_ribbon(aes(ymin=lower_full, ymax = upper_full), 
                data = vals2_all, 
                color = colors["Credible interval Full Model"], 
                show.legend = FALSE, fill = NA, size = 1, linetype = 4) +
    geom_point(aes(x = d0, y = median_meta, color = "Median Meta Model"), 
               data = vals2_all, size =2) +
    geom_line(aes(x = d0, y = median_meta, color = "Median Meta Model"), 
              data = vals2_all, size = 1, linetype = 3) +
    geom_ribbon(aes(ymin=lower_meta, ymax = upper_meta, fill = "Credible interval Meta Model"), 
                alpha = .1, data = vals2_all) +
    geom_ribbon(aes(ymin=lower_meta, ymax = upper_meta), 
                data = vals2_all, 
                color = colors["Credible interval Meta Model"], 
                show.legend=FALSE, fill = NA, size = 1, linetype = 3) +
    geom_line(aes(x = d0, y = vals_true, color = "True relationship"), 
              data = vals2_all, size = 1, linetype = 2) +
    xlab("Delivered drug A") +
    ylab("Delivered drug B") +
    theme_light() +
    theme(text = element_text(size=12),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom",
          legend.justification = c("center"),
          legend.box.just = "left",
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing = unit(0, 'cm'),
          legend.key = element_rect(colour = NA)
    ) +
    labs(fill = NULL,
         color = NULL) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values=colors) + 
    facet_grid(fit ~ tau2)
  
  # Put together into one big grid
  
  
}



# Adjust gap between different params to be larger, between same smaller
# For simulations, only need to show performance for parameters
# For real data, demonstrate fit on real data b/c only have one context
# Should discuss performance of beta_5 and beta_6 separate from the rest since they are different parameters than the ones from the full model
# Double check coverage of lbc_0
# Point to make: when there are covariates, the way we generate them they are not confounders (we adjust for that to reduce trial-to-trial heterogenicity, i.e. don't impact treatment decision); that makes sense for our motivation, pulled clinical trial data where they assign treatment at random