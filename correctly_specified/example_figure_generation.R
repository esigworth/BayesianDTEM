source("post_sim_calculations.R")

load("uc_half.Rda")
load("uc_1.Rda")
load("uc_0.Rda")
load("bc_half.Rda")
load("bc_1.Rda")
load("bc_0.Rda")

bc_0 <- eval_sim_results_ref(ref_bc_0, tau2 = 0, total = length(ref_bc_0), bal = "balanced", covars = "with covars")
bc_half <- eval_sim_results_ref(ref_bc_half, tau2 = 0.5, total = length(ref_bc_half), bal = "balanced", covars = "with covars")
bc_1 <- eval_sim_results_ref(ref_bc_1, tau2 = 1, total = length(ref_bc_1), bal = "balanced", covars = "with covars")

uc_0 <- eval_sim_results_ref(rec_uc_0, tau2 = 0, total = length(ref_uc_0), bal = "unbalanced", covars = "with covars")
uc_half <- eval_sim_results_ref(ref_uc_half, tau2 = 0.5, total = length(ref_uc_half), bal = "unbalanced", covars = "with covars")
uc_1 <- eval_sim_results_ref(ref_uc_1, tau2 = 1, total = length(ref_uc_1), bal = "unbalanced", covars = "with covars")

expr1 <- expression(beta[1])
expr2 <- expression(beta[2])
expr3 <- expression(beta[3])
expr4 <- expression(beta[4])
expr5 <- expression(beta[5])
expr6 <- expression(beta[6])



bc_bias <- rbind(bc_0$bias, bc_half$bias, bc_1$bias) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 21),
         fit = "Balanced")
uc_bias <- rbind(uc_0$bias, uc_half$bias, uc_1$bias) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 21),
         fit = "Unbalanced")

all_bias <- rbind(bc_bias, uc_bias) %>% 
  filter(Beta != "prec" & Beta != "B5" & Beta != "B6" & Beta != "B1") 

fig_dat_withcov <- all_bias %>%
  mutate(Model = ifelse(Model == "studylevel", "SL-C", ifelse(Model == "ipd", "IPD", "SL-NC")),
         Type = Model)

fig_dat_withcov$pos <- rep(c(1.5, 1.8, 1.2, 3.5, 3.8, 3.2, 5.5, 5.8, 5.2), 6)

bias_min_withcov <- min(-1, floor(min(fig_dat_withcov$Low)*10)/10)
bias_max_withcov <- max(1,ceiling(max(fig_dat_withcov$High)*10)/10)

colrs <- c("IPD" = "#fa0fe2", 
           "SL-C" = "#88C6ED",
           "SL-NC" = "#FAA31B")

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
        legend.position = c(0.16, 0.04),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black")) + 
  scale_color_manual(values = colrs) +
  labs(color = "Model Fit") + 
  geom_text(x = 1.5, y = bias_min_withcov - .44, label = expr2, size = 5, color = "black") +
  geom_text(x = 3.5, y = bias_min_withcov - .44, label = expr3, size = 5, color = "black") +
  geom_text(x = 5.5, y = bias_min_withcov - .44, label = expr4, size = 5, color = "black") +
  facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2))) +
  ylab("Percent Bias") 


bc_eff <- rbind(bc_0$eff_dat, bc_half$eff_dat, bc_1$eff_dat) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 12),
         fit = "Balanced",
         type = rep(c("SL-C:IPD", "SL-NC:IPD"), 18)) %>%
  relocate(Efficiency, Beta, Low, High, tau2, fit, type)
uc_eff <- rbind(uc_0$eff_dat, uc_half$eff_dat, uc_1$eff_dat) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 12),
         fit = "Unbalanced",
         type = rep(c("SL-C:IPD", "SL-NC:IPD"), 18)) %>%
  relocate(Efficiency, Beta, Low, High, tau2, fit, type)

all_eff <- rbind(bc_eff, uc_eff) 

eff_dat2 <- all_eff %>% filter(Beta != "B5" & Beta != "B6" & Beta != "B1")

colrs <- c("SL:IPD" = "#88C6ED",
           "SL-C:IPD" = "#88C6ED",
           "SL-NC:IPD" = "#FAA31B")

eff_min2 <- floor(min(eff_dat2$Low)*10)/10
eff_max2 <- ceiling(max(eff_dat2$High)*10)/10

eff_dat2$pos <- rep(c(.8, 1.2, 1.8, 2.2, 2.8, 3.2), 6)


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
        legend.position = c(0.16, 0.033),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 14, colour = "black")) + 
  scale_color_manual(values = colrs) +
  ylab("Relative efficiency") +
  labs(color = "Ratio") + 
  geom_text(x = 1, y = eff_min2 - .05, label = expr2, size = 5, colour = "black") +
  geom_text(x = 2, y = eff_min2 - .05, label = expr3, size = 5, colour = "black") +
  geom_text(x = 3, y = eff_min2 - .05, label = expr4, size = 5, colour = "black") + 
  facet_grid(fit ~ tau2, labeller = label_bquote(col = tau^2 : .(tau2)))


colors <- c("Median IPD Model" = "#fa0fe2", 
            "Credible interval IPD Model" = "#fa0fe2", 
            "Median SL Model" = "#88C6ED",
            "Credible interval SL Model" = "#88C6ED",
            "Median SL-C Model" = "#88C6ED",
            "Credible interval SL-C Model" = "#88C6ED",
            "Median SL-NC Model" = "#FAA31B",
            "Credible interval SL-NC Model" = "#FAA31B",
            "True relationship" = "#7CAE00")

bc_equiv <- rbind(bc_0$equiv_dat, bc_half$equiv_dat, bc_1$equiv_dat) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 11),
         fit = "Balanced")
uc_equiv <- rbind(uc_0$equiv_dat, uc_half$equiv_dat, uc_1$equiv_dat) %>% 
  mutate(tau2 = rep(c(0, .5, 1), each = 11),
         fit = "Unbalanced")

equiv_withcov <- rbind(bc_equiv, uc_equiv)

ggplot(aes(x = d0), data = equiv_withcov) +
  ylim(-1.5, 2.3) +
  geom_ribbon(aes(ymin=lower_studylevel, ymax = upper_studylevel, fill = "Credible interval SL-C Model"), 
              alpha = .1, data = equiv_withcov) +
  geom_ribbon(aes(ymin=lower_studylevel_nocov, ymax = upper_studylevel_nocov, fill = "Credible interval SL-NC Model"), 
              alpha = .2, data = equiv_withcov) +
  geom_ribbon(aes(ymin=lower_ipd, ymax = upper_ipd, fill = "Credible interval IPD Model"), 
              alpha = .2, data = equiv_withcov) +
  geom_line(aes(x = d0, y = vals_true, color = "True relationship"), 
            data = equiv_withcov, size = 2, linetype = 1) +
  geom_line(aes(x = d0, y = median_ipd, color = "Median IPD Model"), 
            data = equiv_withcov, size = 2, linetype = 5) +
  geom_ribbon(aes(ymin=lower_ipd, ymax = upper_ipd), 
              data = equiv_withcov, 
              color = colors["Credible interval IPD Model"], 
              show.legend = FALSE, fill = NA, size = 1, linetype = 5) +
  geom_line(aes(x = d0, y = median_studylevel, color = "Median SL-C Model"), 
            data = equiv_withcov, size = 2, linetype = 2) +
  geom_ribbon(aes(ymin=lower_studylevel, ymax = upper_studylevel), 
              data = equiv_withcov, 
              color = colors["Credible interval SL-C Model"], 
              show.legend=FALSE, fill = NA, size = 1, linetype = 2) +
  geom_line(aes(x = d0, y = median_studylevel_nocov, color = "Median SL-NC Model"), 
            data = equiv_withcov, size = 2, linetype = 3) +
  geom_ribbon(aes(ymin=lower_studylevel_nocov, ymax = upper_studylevel_nocov), 
              data = equiv_withcov, 
              color = colors["Credible interval SL-NC Model"], 
              show.legend = FALSE, fill = NA, size = 1, linetype = 3) +
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
