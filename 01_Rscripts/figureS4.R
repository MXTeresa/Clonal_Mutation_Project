# 1. Source functions ====
source("01_Rscripts/ClonalMuFunctions.R")

# 2. Load packages ====
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(showtext)
library(forcats)
library(parallel)
library(metR)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 14)))
showtext_auto()

# 3. Input data sets ====
sars <- data.frame(ClonalMu = 0:1, freq = c(35, 4)) # R0 = 2.6, 7.4, 14.9

# 4. Input parameters ====
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 3.51, by = 0.05)
lambda_values <- seq(0.01, 4.01, by = 0.1)

# 5. Generate data frame ====
listClonal_sars2.6 <- list_clonal(n_values, 2.6, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sars7.4 <- list_clonal(n_values, 7.4, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sars14.9 <- list_clonal(n_values, 14.9, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_sars2.6 <- LL_meanNb.df(sars, listClonal_sars2.6, seq(0.01, 5.01, by = 0.1), 2.6, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_sars7.4 <- LL_meanNb.df(sars, listClonal_sars7.4, lambda_values, 7.4, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_sars14.9 <- LL_meanNb.df(sars, listClonal_sars14.9, lambda_values, 14.9, mu_values, maxMuGen, maxFS, maxIni)

lambda_sars_hat2.6 <-  df_heatmap_sars2.6$lambda[df_heatmap_sars2.6$prob == max(df_heatmap_sars2.6$prob)]
mu_sars_hat2.6 <-  df_heatmap_sars2.6$mu[df_heatmap_sars2.6$prob == max(df_heatmap_sars2.6$prob)]
min_sars2.6 <- max(df_heatmap_sars2.6$prob)-2.995
max_sars2.6 <- max(df_heatmap_sars2.6$prob)

lambda_sars_hat7.4 <-  df_heatmap_sars7.4$lambda[df_heatmap_sars7.4$prob == max(df_heatmap_sars7.4$prob)]
mu_sars_hat7.4 <-  df_heatmap_sars7.4$mu[df_heatmap_sars7.4$prob == max(df_heatmap_sars7.4$prob)]
min_sars7.4 <- max(df_heatmap_sars7.4$prob)-2.995
max_sars7.4 <- max(df_heatmap_sars7.4$prob)

lambda_sars_hat14.9 <-  df_heatmap_sars14.9$lambda[df_heatmap_sars14.9$prob == max(df_heatmap_sars14.9$prob)]
mu_sars_hat14.9 <-  df_heatmap_sars14.9$mu[df_heatmap_sars14.9$prob == max(df_heatmap_sars14.9$prob)]
min_sars14.9 <- max(df_heatmap_sars14.9$prob)-2.995
max_sars14.9 <- max(df_heatmap_sars14.9$prob)

melt_distri_sars2.6 <- calc_distri_lambda0(sars, 2.6, mu_sars_hat, maxMuGen, maxFS, clonal)
melt_distri_sars7.4 <- calc_distri_lambda0(sars, 7.4, mu_sars_hat, maxMuGen, maxFS, clonal)
melt_distri_sars14.9 <- calc_distri_lambda0(sars, 14.9, mu_sars_hat, maxMuGen, maxFS, clonal)

# 6. Graphs ====
## 6.1 R0=2.6 ====
clonal_distri_sars2.6 <- ggplot(melt_distri_sars2.6, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars2.6 <- ggplot(df_heatmap_sars2.6) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 300) +
  geom_hline(aes(yintercept = mu_sars_hat2.6), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars2.6, max_sars2.6)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3, y = 0.45, label = "mu==0.16", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 2.2, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars2.6 <- ggplot(df_heatmap_sars2.6) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 280) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars2.6, max_sars2.6)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars2.6 <- ggplot(df_heatmap_sars2.6) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars2.6, max_sars2.6)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.2 R0=7.4 ====
clonal_distri_sars7.4 <- ggplot(melt_distri_sars7.4, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars7.4 <- ggplot(df_heatmap_sars7.4) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 300) +
  geom_hline(aes(yintercept = mu_sars_hat7.4), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars7.4, max_sars7.4)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.5, y = 0.7, label = "mu==0.51", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 2.2, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars7.4 <- ggplot(df_heatmap_sars7.4) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 280) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars7.4, max_sars7.4)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars7.4 <- ggplot(df_heatmap_sars7.4) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars7.4, max_sars7.4)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.3 R0=14.9 ====
clonal_distri_sars14.9 <- ggplot(melt_distri_sars14.9, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars14.9 <- ggplot(df_heatmap_sars14.9) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 300) +
  geom_hline(aes(yintercept = mu_sars_hat14.9), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars14.9, max_sars14.9)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.5, y = 0.95, label = "mu==0.76", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 2.7, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars14.9 <- ggplot(df_heatmap_sars14.9) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 280) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars14.9, max_sars14.9)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars14.9 <- ggplot(df_heatmap_sars14.9) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars14.9, max_sars14.9)) +
  xlim(0, 5) +
  ylim(0, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
S4_r1 <- clonal_distri_sars2.6 + lambda_heat_sars2.6 + meanN_heat_sars2.6 + meanNb_heat_sars2.6 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 2.6")), tag_levels = "A")

S4_r2 <- clonal_distri_sars7.4 + lambda_heat_sars7.4 + meanN_heat_sars7.4 + meanNb_heat_sars7.4 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 7.4")), tag_levels = list(c("E", "F", "G", "H")))

S4_r3 <- clonal_distri_sars14.9 + lambda_heat_sars14.9 + meanN_heat_sars14.9 + meanNb_heat_sars14.9 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 14.9")), tag_levels = list(c("I", "J", "K", "L")))

figureS4 <- wrap_elements(S4_r1) / wrap_elements(S4_r2) / wrap_elements(S4_r3)

ggsave(figureS4, file = "figureS4.pdf", path = "02_plots", width = 10, height = 9)
