library(data.table)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(scales)
library(rstatix)
library(broom)
library(cowplot)
library(multcomp)
library(ggridges)
library(ggpattern)
library(viridisLite)

rm(list = ls())
setwd("/Users/kexinxie/Downloads/GitHub/SpatialSpilloverEffect_Measles/") 
total_cost_inci_all_expt<-fread("data/total_costs_inci_all_expt.csv")


# ------------------------------------------------------------------------------
# Figure 2(A)
# ------------------------------------------------------------------------------
p1 <- ggplot(total_cost_inci_all_expt, aes(x = factor(tau), y = log(incidence),
                                           fill = factor(vhi), color = factor(vhi))) +
  geom_boxplot(position = "dodge2", alpha = 0.8) +
  facet_wrap(~alpha, labeller = label_bquote(alpha == .(alpha))) +
  labs(
    y = "Log of measles incidence",
    x = expression(paste("Transmissibility (", tau, ")")),
    fill = expression(paste("Home isolation/quarantine compliance % (", gamma, ")")),
    color = expression(paste("Home isolation/quarantine compliance % (", gamma, ")"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 15),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    text = element_text(size = 15)
  )

# ------------------------------------------------------------------------------
# Figure 2(B): Total Cost vs Settings
# ------------------------------------------------------------------------------
p2 <- ggplot(total_cost_inci_all_expt, aes(x = factor(tau), y = log(totalCost / 1e6),
                                           fill = factor(vhi), color = factor(vhi))) +
  geom_boxplot(position = "dodge2", alpha = 0.7) +
  facet_wrap(~alpha, labeller = label_bquote(alpha == .(alpha)), scales = "free_y") +
  labs(
    y = "Log of total cost (million US$)",
    x = expression(paste("Transmissibility (", tau, ")")),
    fill = expression(paste("Home isolation/quarantine compliance % (", gamma, ")")),
    color = expression(paste("Home isolation/quarantine compliance % (", gamma, ")"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 15),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    text = element_text(size = 15)
  )

# ------------------------------------------------------------------------------
# Figure 2(C): Total Cost with Outliers Removed
# ------------------------------------------------------------------------------
p3 <- ggplot(total_cost_inci_all_expt, aes(x = factor(tau), y = log(totalCost / 1e6),
                                           fill = factor(vhi), color = factor(vhi))) +
  geom_boxplot(position = "dodge2", alpha = 0.7, outliers=FALSE) +
  facet_wrap(~alpha, labeller = label_bquote(alpha == .(alpha)), scales = "free_y") +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  labs(
    y = "Log of total cost (million US$)",
    x = expression(paste("Transmissibility (", tau, ")")),
    fill = expression(paste("Home isolation/quarantine compliance % (", gamma, ")")),
    color = expression(paste("Home isolation/quarantine compliance % (", gamma, ")"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 15),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    text = element_text(size = 15)
  )

# ------------------------------------------------------------------------------
# Figure 2: Combine all three panels
# ------------------------------------------------------------------------------
pdf("results/plots/incidence and cost vs settings outliers.pdf", width = 12, height = 5.5)
ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom",
          widths = c(4, 5, 5)) +
  draw_plot_label(label = c("A", "B", "C"), size = 15, x = c(0, 0.3, 0.65), y = c(1, 1, 1))
dev.off()

# ------------------------------------------------------------------------------
# Individual Subplots (e.g., for inset or publication)
# ------------------------------------------------------------------------------
p2_1 <- total_cost_inci_all_expt %>%
  filter(alpha == "15" & tau == 0.4) %>%
  ggplot(aes(x = factor(tau), y = log(totalCost / 1e6),
             fill = factor(vhi), color = factor(vhi))) +
  geom_boxplot(position = "dodge2", alpha = 0.7, outliers=FALSE) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  theme_bw() +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

pdf("results/plots/cost tau0.4 alpha15.pdf", width = 1, height = 1.2)
p2_1
dev.off()

# ------------------------------------------------------------------------------
# Figure 3(A): Outbreak Probabilities vs Alpha
# ------------------------------------------------------------------------------
p3 <- total_cost_inci_all_expt %>%
  group_by(alpha) %>%
  summarise(
    `200` = mean(incidence > 200),
    `400` = mean(incidence > 400),
    `800` = mean(incidence > 800),
    `1600` = mean(incidence > 1600),
    `2400` = mean(incidence > 2400)
  ) %>%
  pivot_longer(-alpha, names_to = "threshold", values_to = "prob") %>%
  ggplot(aes(x = factor(alpha), y = prob, group = threshold,
             color = threshold, linetype = threshold, shape = threshold)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  labs(
    y = expression(paste("Probability of outbreak size being > ", kappa)),
    x = expression(paste("Reduction in MMR vaccination level % (", alpha, ")")),
    color = expression(paste("Threshold of outbreak size (", kappa, ")")),
    linetype = expression(paste("Threshold of outbreak size (", kappa, ")")),
    shape = expression(paste("Threshold of outbreak size (", kappa, ")"))
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20), panel.border = element_rect(color = "black"))

# ------------------------------------------------------------------------------
# Figure 3(B): Economic Burden Probabilities vs Alpha
# ------------------------------------------------------------------------------
p4 <- total_cost_inci_all_expt %>%
  group_by(alpha) %>%
  summarise(
    `2000` = mean(totalCost / 1e6 > 2000),
    `2500` = mean(totalCost / 1e6 > 2500),
    `3000` = mean(totalCost / 1e6 > 3000),
    `3500` = mean(totalCost / 1e6 > 3500),
    `4000` = mean(totalCost / 1e6 > 4000)
  ) %>%
  pivot_longer(-alpha, names_to = "threshold", values_to = "prob") %>%
  ggplot(aes(x = factor(alpha), y = prob, group = threshold,
             color = threshold, linetype = threshold, shape = threshold)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  labs(
    y = expression(atop("Probability of economic burden", paste("(million US$) being > ", kappa))),
    x = expression(paste("Reduction in MMR vaccination level % (", alpha, ")")),
    color = expression(paste("Threshold of economic burden (", kappa, ")")),
    linetype = expression(paste("Threshold of economic burden (", kappa, ")")),
    shape = expression(paste("Threshold of economic burden (", kappa, ")"))
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20), panel.border = element_rect(color = "black"))

pdf("results/plots/Probability of outbreak and cost.pdf", width = 16.8, height = 6.8)
ggarrange(p3, p4, ncol = 2, nrow = 1) +
  draw_plot_label(label = c("A", "B"), size = 20, x = c(0, 0.5), y = c(1, 1))
dev.off()

# ------------------------------------------------------------------------------
# Figure 4: Cost Heatmap (Table 1 in appendix)
# ------------------------------------------------------------------------------
mean_total_cost_inci_all_expt <- total_cost_inci_all_expt %>%
  group_by(alpha, tau, vhi) %>%
  summarise(
    totalCost_mean = mean(totalCost) / 1e6,
    totalCost_sd = sd(totalCost) / 1e6,
    directMedicalCost_mean = mean(directMedicalCost) / 1e6,
    directMedicalCost_sd = sd(directMedicalCost) / 1e6,
    prodlosscost_mean = mean(prod_loss_cost) / 1e6,
    prodlosscost_sd = sd(prod_loss_cost) / 1e6,
    mmrcost_mean = mean(mmr_cost) / 1e6
  )

custom_labeller <- labeller(
  CostType = c(
    directMedicalCost_mean = "Treatment~Cost",
    prodlosscost_mean = "Productivity~Loss~Cost",
    totalCost_mean = "Total~Cost"
  ),
  alpha = c(
    `0` = "(MMR~cost==1603.81)",
    `5` = "(MMR~cost==1584.22)",
    `15` = "(MMR~cost==1545.05)",
    `25` = "(MMR~cost==1505.88)"
  ),
  .default = label_parsed
)

p1 <- mean_total_cost_inci_all_expt %>%
  pivot_longer(cols = c(totalCost_mean, directMedicalCost_mean, prodlosscost_mean),
               names_to = "CostType", values_to = "CostValue") %>%
  mutate(cost_group = cut(CostValue,
                          breaks = c(0, 2.5, 5, 10, 50, 100, 500, 1000, 1500, 1600, 1700, 1800, 1900, 2000, 2250, 2500, 3000))) %>%
  ggplot(aes(x = tau, y = vhi, fill = cost_group)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(CostValue, 2)), color = "black", size = 4.5) +
  facet_grid(rows = vars(CostType), cols = vars(alpha), labeller = custom_labeller) +
  scale_fill_manual(values = plasma(20)[7:20], drop = TRUE) +
  labs(
    x = expression(paste("Transmissibility (", tau, ")")),
    y = expression(paste("Home isolation compliance % (", gamma, ")")),
    fill = "Total Cost (million US$)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(color = "black"),
    text = element_text(size = 13)
  )

pdf("results/plots/cost_heatmap.pdf", width = 11, height = 9)
p1
dev.off()
