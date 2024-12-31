############################## FRONT MATTER ##############################
#   SCRIPT:     simulation_figs.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      Script to construct simulation figures

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, fredr, tidycensus, haven, xtable, geomtextpath)

# Set project root:
i_am("src/R/s_simulation_figs.R")

# STEP 1: Simulations 1-3 ----

# Import simulated data summary stats:
sim_data <- fread(here("data/simulations/sim_data_sim1to6.csv")) 
sim_data_labels = data.table(
  sim = c("sim1", "sim2", "sim3", "sim4", "sim5", "sim6"),
  label = c(
    "Simulation 1: Random Idiosynratic Pref.",
    "Simulation 2: Negative Correlation",
    "Simulation 3: Positive Correlation",
    "Simulation 4: Additive Msmt. Error",
    "Simulation 5: Nonlinear Msmt. Error",
    "Simulation 6: Negative\nCorrelation + Additive Msmt. Error"
  )
)
sim_data <- merge(sim_data, sim_data_labels)

# Plot example DGP:
fig_dgp1 <- ggplot(sim_data[sim %in% c("sim1", "sim2", "sim3")], aes(xi, cost, color = sim)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", formula = 'y ~ x', color = "black", linewidth = 1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  facet_grid(~label) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "none", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("Idiosyncratic Preference ("*xi*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Travel Cost")
ggsave(here("output/simulations/fig_dgp1.jpg"), plot = fig_dgp1, width = 8, height = 4, units = "in", dpi = 400)
rm(fig_dgp1)

# Import estimated WTP
wtp <- lapply(1:6, function(x){
  fread(here(paste0("data/simulations/est_wtp_sim", x, ".csv")))
}) %>%
  rbindlist(.)
wtp <- merge(wtp, sim_data_labels)
rm(sim_data_labels)

# Plot standard estimates:
fig_wtp1 <- ggplot(wtp[est == "base" & sim %in% c("sim1", "sim2", "sim3")], aes(x = wtp,  color = sim, fill = sim)) +
  geom_vline(xintercept = -0.5) +
  geom_density(aes(y = after_stat(ndensity)), alpha = 0.5) +
  facet_grid(~label) +
  #scale_x_continuous(limits = c(-0.8, -0.2)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "none", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("WTP Estimate ("*beta/alpha*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Scaled Density")
ggsave(here("output/simulations/fig_wtp1.jpg"), plot = fig_wtp1, width = 8, height = 4, units = "in", dpi = 400)
rm(fig_wtp1)

# STEP 2: Simulations 4-6 ----

# Plot example DGP:
fig_dgp2 <- ggplot(sim_data[sim %in% c("sim4", "sim5", "sim6")], aes(cost, cost_obs, color = sim)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", formula = 'y ~ x', color = "black", linewidth = 1) +
  scale_color_manual(values = c("#999999", "#0072B2", "#D55E00")) +
  facet_grid(~label) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "none", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = "True Travel Cost") +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Observed Travel Cost")
ggsave(here("output/simulations/fig_dgp2.jpg"), plot = fig_dgp2, width = 8, height = 4, units = "in", dpi = 400)
rm(fig_dgp2)

# Plot standard estimates:
fig_wtp2 <- ggplot(wtp[est == "base" & sim %in% c("sim4", "sim5", "sim6")], aes(x = wtp,  color = sim, fill = sim)) +
  geom_vline(xintercept = -0.5) +
  geom_density(aes(y = after_stat(ndensity)), alpha = 0.5) +
  facet_grid(~label) +
  #scale_x_continuous(limits = c(-0.8, -0.2)) +
  scale_color_manual(values = c("#999999", "#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#999999", "#0072B2", "#D55E00")) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "none", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("WTP Estimate ("*beta/alpha*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Scaled Density")
ggsave(here("output/simulations/fig_wtp2.jpg"), plot = fig_wtp2, width = 8, height = 4, units = "in", dpi = 400)
rm(fig_wtp2)

# STEP 3: Figure comparing estimators for simulations 1-6 ----

# Add estimation labels:
wtp[ , est_label := ifelse(est == "base", "No Correction", ifelse(est == "groupfe", "Group FE", "Control Function"))]

# Figure comparing estimators:
fig_wtp3 <- ggplot(wtp[est != "groupfe"], aes(x = wtp,  color = sim, fill = sim, linetype = est_label, alpha = est_label)) +
  geom_vline(xintercept = -0.5) +
  geom_density(aes(y = after_stat(ndensity))) +
  facet_wrap(~label,ncol = 3) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#D55E00"), guide = "none") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#D55E00"), guide = "none") +
  scale_alpha_manual(values = c(0.5, 0.0, 0.0)) +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "bottom", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("WTP Estimate ("*beta/alpha*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Scaled Density")
ggsave(here("output/simulations/fig_wtp3.jpg"), plot = fig_wtp3, width = 8, height = 6, units = "in", dpi = 400)
rm(fig_wtp3)

# STEP 4: Table comparing estimators for simulations 1-6 ----

# Calculate summary statistics:
wtp_sumstat <-
  wtp[est != "groupfe", .(
    wtp_mean = mean(wtp),
    wtp_bias = mean(wtp + 0.5),
    wtp_mse = mean((wtp + 0.5) ^ 2),
    wtp_tstat = (mean(wtp) + 0.5) / (sd(wtp))
  ), by = c("sim", "est")]

# Reshape and save as tex:
wtp_sumstat_wide <-
  dcast(
    wtp_sumstat,
    sim ~ est,
    value.var = c("wtp_mean", "wtp_bias", "wtp_mse", "wtp_tstat")
  )
print.xtable(
  xtable(wtp_sumstat_wide[, c(2, 4, 6, 8,3, 5, 7, 9)], digits = 3),
  file = here("output/simulations/tab_sim_est.tex"),
  include.rownames = F,
  floating = F
)
rm(wtp_sumstat, wtp_sumstat_wide)

# STEP 5: Table comparing comparing group fe estimators for simulations 1-6 ----

# Calculate summary statistics:
wtp_sumstat_full <-
  wtp[, .(
    wtp_mean = mean(wtp),
    wtp_bias = mean(wtp + 0.5),
    wtp_mse = mean((wtp + 0.5) ^ 2),
    wtp_tstat = (mean(wtp) + 0.5) / (sd(wtp))
  ), by = c("sim", "est")]

# Reshape and save as tex:
wtp_sumstat_full_wide <-
  dcast(
    wtp_sumstat_full,
    sim ~ est,
    value.var = c("wtp_mean", "wtp_bias", "wtp_mse", "wtp_tstat")
  )
print.xtable(
  xtable(wtp_sumstat_full_wide[, c(2, 5, 8, 11, 4, 7, 10, 13, 3, 6, 9, 12)], digits = 3),
  file = here("output/simulations/tab_sim_est_groupfe.tex"),
  include.rownames = F,
  floating = F
)
rm(wtp_sumstat_full_wide, wtp_sumstat_full)

# STEP 6: Figure showing performance of estimators to different degrees of w/in group correlation in xi ----

# Import data from RC with different w/in group correlations:
rc1_bias <- fread(here("data/simulations/est_wtp_rc1.csv"))

# Plot estimated T-statistics for different degrees of w/in group correlation:
fig_rc1 <- ggplot(rc1_bias, aes(x = corel, y = bias / 0.5, color = est, linetype = est, shape = est)) +
  geom_hline(yintercept = 0.0, color = "#999999", linetype = "solid", linewidth= 0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                     labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_shape_manual(values = c(15, 16, 17),
                     labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                     labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_y_continuous(label = scales::percent) + 
  theme_classic() +
  theme(plot.subtitle = element_text(hjust = -0.1), legend.title = element_blank(), 
        legend.position.inside = c(0.8, 0.8), legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("Degree of Within-group Idiosyncratic Preference Correlation")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Mean Bias (%)")
ggsave(here("output/simulations/fig_rc1.jpg"), plot = fig_rc1, width = 5, height = 4, units = "in", dpi = 400)
rm(rc1_bias, fig_rc1)

# STEP 7: Figure showing performance of estimators to different group sizes ----

# Import data from RC with different group sizes:
rc2_bias <- fread(here("data/simulations/est_wtp_rc2.csv"))

# Plot estimated T-statistics for different degrees of w/in group correlation:
fig_rc2 <- ggplot(rc2_bias, aes(x = n_groups, y = bias / 0.5, color = est, linetype = est, shape = est)) +
  geom_hline(yintercept = 0.0, color = "#999999", linetype = "solid", linewidth= 0.1) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                     labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_shape_manual(values = c(0, 1, 2),
                     labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        labels = c("No Correction", "Control Function", "Origin Group FE")) +
  scale_y_continuous(label = scales::percent) + 
  theme_classic() +
  theme(plot.subtitle = element_text(hjust = -0.1), legend.title = element_blank(), 
        legend.position.inside = c(0.8, 0.7), legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("Number of Origin Groups")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Mean Bias (%)")
ggsave(here("output/simulations/fig_rc2.jpg"), plot = fig_rc2, width = 5, height = 4, units = "in", dpi = 400)
rm(rc2_bias, fig_rc2)

# STEP 8: Figure showing robustness of estimators to nonlinearities in DGP ----

# Import data from RC with different nonlinearities in unobservable:
rc3_wtp <- fread(here("data/simulations/est_wtp_rc3.csv"))

# Add labels:
rc3_wtp[ , est_label := ifelse(est == 1, "No Correction", ifelse(est == 2, "Group FE", "Control Function"))]
rc3_wtp <- merge(rc3_wtp, data.table(
  nonlinear = 1:6,
  nl_type = factor(
    c(rep("Quadratic", 3), rep("Exponential", 3)),
    levels = c("Quadratic", "Exponential")
    ),
  nl_fxn =  factor(
    rep(c("Cost", "Utility", "Cost + Utility"), 2),
    levels = c("Cost", "Utility", "Cost + Utility")
  )
))

# Plot distribution of WTP statistics for different estimators, nonlinearities in unobservable DGP:
fig_rc3 <- ggplot(rc3_wtp[est != 2], aes(x = wtp,  color = as.character(nonlinear), fill = as.character(nonlinear), linetype = est_label, alpha = est_label)) +
  geom_vline(xintercept = -0.5) +
  geom_density(aes(y = after_stat(ndensity))) +
  facet_grid(nl_type~ nl_fxn) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#D55E00"), guide = "none") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#D55E00"), guide = "none") +
  scale_alpha_manual(values = c(0.5, 0.0, 0.0)) +
  theme_bw() +
  scale_x_continuous(breaks = c(-0.55, -0.5, -0.45, -0.4)) +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "bottom", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("WTP Estimate ("*beta/alpha*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Scaled Density")
ggsave(here("output/simulations/fig_rc3.jpg"), plot = fig_rc3, width = 8, height = 6, units = "in", dpi = 400)
rm(fig_rc3, rc3_wtp)

# Import data from RC with different unobservables in true instrument DGP:
rc4_wtp <- fread(here("data/simulations/est_wtp_rc4.csv"))

# Add labels:
rc4_wtp[ , est_label := ifelse(est == 1, "No Correction", ifelse(est == 2, "Group FE", "Control Function"))]
rc4_wtp[ , nl_type := ifelse(nonlinear == 8, "Costs Exponential in Z", "Costs Quadratic in Z")]
rc4_wtp[ , nl_type := factor(nl_type, levels = c("Costs Quadratic in Z", "Costs Exponential in Z"))]

# Plot distribution of WTP statistics for different estimators, nonlinearities in Z DGP:
fig_rc4 <- ggplot(rc4_wtp[est != 2], aes(x = wtp,  color = nl_type, fill = nl_type, linetype = est_label, alpha = est_label)) +
  geom_vline(xintercept = -0.5) +
  geom_density(aes(y = after_stat(ndensity))) +
  facet_wrap(~nl_type, ncol = 2) + 
  scale_color_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
  scale_alpha_manual(values = c(0.5, 0.0)) +
  theme_bw() +
  scale_x_continuous(breaks = c(-0.55, -0.5, -0.45, -0.4)) +
  theme(plot.subtitle = element_text(hjust = -0.01), legend.title = element_blank(), 
        legend.position = "bottom", legend.spacing.y = unit(-0.1, "cm")) +
  xlab(label = expression("WTP Estimate ("*beta/alpha*")")) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Scaled Density")
ggsave(here("output/simulations/fig_rc4.jpg"), plot = fig_rc4, width = 7, height = 4, units = "in", dpi = 400)
  rc4_wtp[, .(
    wtp_mean = mean(wtp),
    wtp_bias = mean(wtp + 0.5),
    wtp_mse = mean((wtp + 0.5) ^ 2),
    wtp_tstat = (mean(wtp) + 0.5) / (sd(wtp))
  ), by = c("nonlinear", "est")]
rm(rc4_wtp, fig_rc4)

# STEP 9: Figure showing robustness of CF to weak instruments ----

# Import calculated F-statistics from weak instrument simulations:
rc5_fstat <- fread(here("data/simulations/est_wtp_rc_weak_fstat.csv"))
rc5_fstat <- rc5_fstat[, .(mean_f_stat = mean(f_stat)), by = z_correl]

# Import calcualted bias from weak instrument simulations:
rc5_bias <- fread(here("data/simulations/est_wtp_rc_weak.csv"))
rc5_bias[ , est_label := ifelse(est == "base", "No Correction", ifelse(est == "groupfe", "Group FE", "Control Function"))]

fig_rc5 <- ggplot(data = rc5_bias[est != "groupfe"]) +
  geom_point(aes(x = z_correl, y=abs(bias)/0.5, color = est, shape = est), show.legend = F) +
  geom_textline( aes(x = z_correl, y=abs(bias)/0.5, color = est, linetype = est, label = est_label), show.legend = F, hjust = 0.85, straight = T) +
  geom_point(data = rc5_fstat, aes(x = z_correl, y =( mean_f_stat/1000)), shape = 4, size = 4, col = "#E69F00") +
  scale_color_manual(values = c("#009E73", "#56B4E9"),) +
  scale_shape_manual(values = c(0, 1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous("Mean Absolute Bias (%)", label = scales::percent, sec.axis = sec_axis(~ .*1000, name = "Mean First-stage F-statistic")) + 
  theme_classic() +
  theme(plot.subtitle = element_text(hjust = -0.1), legend.title = element_blank(), 
        legend.position.inside = c(0.8, 0.7), legend.spacing.y = unit(-0.1, "cm"),
        axis.title.y.right=element_text(color="#E69F00"),
        axis.text.y.right=element_text(color="#E69F00")) +
  xlab(label = "Instrument Strength") +
  ylab(element_blank()) + 
  labs(title = element_blank())
ggsave(here("output/simulations/fig_rc5.jpg"), plot = fig_rc5, width = 6, height = 4, units = "in", dpi = 400)
rm(rc5_fstat, rc5_bias, fig_rc5)
rm(sim_data, wtp)
