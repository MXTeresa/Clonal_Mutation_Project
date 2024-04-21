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
sars0.5 <- data.frame(ClonalMu = 0:1, freq = c(37, 2)) # R0 = 7.4
sars3 <- data.frame(ClonalMu = 0:1, freq = c(35, 4))
sars7 <- data.frame(ClonalMu = 0:1, freq = c(30, 9))

# 4. Input parameters ====
R0 <- 7.4
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 3.51, by = 0.05)
lambda_values <- seq(0.01, 4.01, by = 0.1)

# 5. Generate data frame ====
listClonal_sars <- list_clonal(n_values, R0, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_sars_thre0.5 <- LL_meanNb.df(sars0.5, listClonal_sars, seq(0.01, 4.81, by = 0.1), R0, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_sars_thre3 <- LL_meanNb.df(sars3, listClonal_sars, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_sars_thre7 <- LL_meanNb.df(sars7, listClonal_sars, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)

lambda_sars_hat_thre0.5 <-  df_heatmap_sars_thre0.5$lambda[df_heatmap_sars_thre0.5$prob == max(df_heatmap_sars_thre0.5$prob)]
mu_sars_hat_thre0.5 <-  df_heatmap_sars_thre0.5$mu[df_heatmap_sars_thre0.5$prob == max(df_heatmap_sars_thre0.5$prob)]
min_sars_thre0.5 <- max(df_heatmap_sars_thre0.5$prob)-2.995
max_sars_thre0.5 <- max(df_heatmap_sars_thre0.5$prob)

lambda_sars_hat_thre3 <-  df_heatmap_sars_thre3$lambda[df_heatmap_sars_thre3$prob == max(df_heatmap_sars_thre3$prob)]
mu_sars_hat_thre3 <-  df_heatmap_sars_thre3$mu[df_heatmap_sars_thre3$prob == max(df_heatmap_sars_thre3$prob)]
min_sars_thre3 <- max(df_heatmap_sars_thre3$prob)-2.995
max_sars_thre3 <- max(df_heatmap_sars_thre3$prob)

lambda_sars_hat_thre7 <-  df_heatmap_sars_thre7$lambda[df_heatmap_sars_thre7$prob == max(df_heatmap_sars_thre7$prob)]
mu_sars_hat_thre7 <-  df_heatmap_sars_thre7$mu[df_heatmap_sars_thre7$prob == max(df_heatmap_sars_thre7$prob)]
min_sars_thre7 <- max(df_heatmap_sars_thre7$prob)-2.995
max_sars_thre7 <- max(df_heatmap_sars_thre7$prob)

melt_distri_sars_thre0.5 <- calc_distri_lambda0(sars0.5, R0, mu_sars_hat_thre0.5, maxMuGen, maxFS, clonal)
melt_distri_sars_thre3 <- calc_distri_lambda0(sars3, R0, mu_sars_hat_thre3, maxMuGen, maxFS, clonal)
melt_distri_sars_thre7 <- calc_distri_lambda0(sars7, R0, mu_sars_hat_thre7, maxMuGen, maxFS, clonal)

# 6. Graphs ====
## 6.1 Threshold=0.5% ====
clonal_distri_sars_thre0.5 <- ggplot(melt_distri_sars_thre0.5, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars_thre0.5 <- ggplot(df_heatmap_sars_thre0.5) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  geom_hline(aes(yintercept = mu_sars_hat_thre0.5), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre0.5, max_sars_thre0.5)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2.2, y = 0.6, label = "mu==0.31", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 3.2, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars_thre0.5 <- ggplot(df_heatmap_sars_thre0.5) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre0.5, max_sars_thre0.5)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars_thre0.5 <- ggplot(df_heatmap_sars_thre0.5) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre0.5, max_sars_thre0.5)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.2 Threshold=3% ====
clonal_distri_sars_thre3 <- ggplot(melt_distri_sars_thre3, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars_thre3 <- ggplot(df_heatmap_sars_thre3) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 330) +
  geom_hline(aes(yintercept = mu_sars_hat_thre3), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre3, max_sars_thre3)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.5, y = 0.8, label = "mu==0.51", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 2.5, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars_thre3 <- ggplot(df_heatmap_sars_thre3) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre3, max_sars_thre3)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars_thre3 <- ggplot(df_heatmap_sars_thre3) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre3, max_sars_thre3)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.3 Threshold = 7% ====
clonal_distri_sars_thre7 <- ggplot(melt_distri_sars_thre7, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_sars_thre7 <- ggplot(df_heatmap_sars_thre7) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 300) +
  geom_hline(aes(yintercept = mu_sars_hat_thre7), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre7, max_sars_thre7)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.5, y = 1, label = "mu==0.71", parse = TRUE, color = "red") +
  annotate("text", x = 1.2, y = 3.2, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars_thre7 <- ggplot(df_heatmap_sars_thre7) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre7, max_sars_thre7)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars_thre7 <- ggplot(df_heatmap_sars_thre7) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 250) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_thre7, max_sars_thre7)) +
  xlim(0, 5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
S6_r1 <- clonal_distri_sars_thre0.5 + lambda_heat_sars_thre0.5 + meanN_heat_sars_thre0.5 + meanNb_heat_sars_thre0.5 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = "threshold = 0.5%", tag_levels = "A")

S6_r2 <- clonal_distri_sars_thre3 + lambda_heat_sars_thre3 + meanN_heat_sars_thre3 + meanNb_heat_sars_thre3 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = "threshold = 3%", tag_levels = list(c("E", "F", "G", "H")))

S6_r3 <- clonal_distri_sars_thre7 + lambda_heat_sars_thre7 + meanN_heat_sars_thre7 + meanNb_heat_sars_thre7 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = "threshold = 7%", tag_levels = list(c("I", "J", "K", "L")))

figureS6 <- wrap_elements(S6_r1) / wrap_elements(S6_r2) / wrap_elements(S6_r3)

ggsave(figureS6, file = "figureS6.pdf", path = "02_plots", width = 10, height = 9)
