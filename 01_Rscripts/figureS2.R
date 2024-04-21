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
                  axis.title = element_text(size = 14),
                  legend.text = element_text(size = 10)))
showtext_auto()

# 3. Input data sets ====
flu <- data.frame(ClonalMu = 0:3, freq = c(42, 5, 2, 3)) # R0 = 4.4, 11.1, 37.7

# 4. Input parameters ====
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 3.51, by = 0.05)
lambda_values <- seq(0.01, 3.01, by = 0.1)

# 5. Generate data frame ====
listClonal_flu4.4 <- list_clonal(n_values, 4.4, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_flu11.1 <- list_clonal(n_values, 11.1, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_flu37.7 <- list_clonal(n_values, 37.7, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_flu4.4 <- LL_meanNb.df(flu, listClonal_flu4.4, lambda_values, 4.4, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_flu11.1 <- LL_meanNb.df(flu, listClonal_flu11.1, lambda_values, 11.1, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_flu37.7 <- LL_meanNb.df(flu, listClonal_flu37.7, lambda_values, 37.7, mu_values, maxMuGen, maxFS, maxIni)

lambda_flu_hat4.4 <-  df_heatmap_flu4.4$lambda[df_heatmap_flu4.4$prob == max(df_heatmap_flu4.4$prob)]
mu_flu_hat4.4 <-  df_heatmap_flu4.4$mu[df_heatmap_flu4.4$prob == max(df_heatmap_flu4.4$prob)]
min_flu4.4 <- max(df_heatmap_flu4.4$prob)-2.995
max_flu4.4 <- max(df_heatmap_flu4.4$prob)

lambda_flu_hat11.1 <-  df_heatmap_flu11.1$lambda[df_heatmap_flu11.1$prob == max(df_heatmap_flu11.1$prob)]
mu_flu_hat11.1 <-  df_heatmap_flu11.1$mu[df_heatmap_flu11.1$prob == max(df_heatmap_flu11.1$prob)]
min_flu11.1 <- max(df_heatmap_flu11.1$prob)-2.995
max_flu11.1 <- max(df_heatmap_flu11.1$prob)

lambda_flu_hat37.7 <-  df_heatmap_flu37.7$lambda[df_heatmap_flu37.7$prob == max(df_heatmap_flu37.7$prob)]
mu_flu_hat37.7 <-  df_heatmap_flu37.7$mu[df_heatmap_flu37.7$prob == max(df_heatmap_flu37.7$prob)]
min_flu37.7 <- max(df_heatmap_flu37.7$prob)-2.995
max_flu37.7 <- max(df_heatmap_flu37.7$prob)

melt_distri_flu4.4 <- calc_distri_lambda0(flu, 4.4, mu_flu_hat, maxMuGen, maxFS, clonal)
melt_distri_flu11.1 <- calc_distri_lambda0(flu, 11.1, mu_flu_hat, maxMuGen, maxFS, clonal)
melt_distri_flu37.7 <- calc_distri_lambda0(flu, 37.7, mu_flu_hat, maxMuGen, maxFS, clonal)

# 6. Graphs ====
## 6.1 R0=4.4 ====
clonal_distri_flu4.4 <- ggplot(melt_distri_flu4.4, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") + 
  ylim(0, 1)

lambda_heat_flu4.4 <- ggplot(df_heatmap_flu4.4) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 1500) +
  geom_hline(aes(yintercept = mu_flu_hat4.4), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu4.4, max_flu4.4)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2.4, y = 1.4, label = "mu==1.06", parse = TRUE, color = "red") +
  annotate("text", x = 0.8, y = 2.8, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu4.4 <- ggplot(df_heatmap_flu4.4) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 1300) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu4.4, max_flu4.4)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu4.4 <- ggplot(df_heatmap_flu4.4) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 900) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu4.4, max_flu4.4)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.2 R0=11.1 ====
clonal_distri_flu11.1 <- ggplot(melt_distri_flu11.1, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_flu11.1 <- ggplot(df_heatmap_flu11.1) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 1400) +
  geom_hline(aes(yintercept = mu_flu_hat11.1), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu11.1, max_flu11.1)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2, y = 1.9, label = "mu==1.56", parse = TRUE, color = "red") +
  annotate("text", x = 0.8, y = 0.4, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu11.1 <- ggplot(df_heatmap_flu11.1) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 1050) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu11.1, max_flu11.1)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu11.1 <- ggplot(df_heatmap_flu11.1) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 850) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu11.1, max_flu11.1)) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.3 R0=37.7 ====
clonal_distri_flu37.7 <- ggplot(melt_distri_flu37.7, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  ylim(0, 1)

lambda_heat_flu37.7 <- ggplot(df_heatmap_flu37.7) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 1000) +
  geom_hline(aes(yintercept = mu_flu_hat37.7), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu37.7, max_flu37.7),
                                   breaks = c(-52, -53, -54),
                                   labels = c("-52", "-53", "-54")) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2, y = 2.1, label = "mu==1.71", parse = TRUE, color = "red") +
  annotate("text", x = 0.8, y = 0.4, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu37.7 <- ggplot(df_heatmap_flu37.7) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 1000) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu37.7, max_flu37.7),
                                   breaks = c(-52, -53, -54),
                                   labels = c("-52", "-53", "-54")) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu37.7 <- ggplot(df_heatmap_flu37.7) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 1000) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu37.7, max_flu37.7),
                                   breaks = c(-52, -53, -54),
                                   labels = c("-52", "-53", "-54")) +
  xlim(0, 3) +
  ylim(0, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
S2_r1 <- clonal_distri_flu4.4 + lambda_heat_flu4.4 + meanN_heat_flu4.4 + meanNb_heat_flu4.4 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 4.4")), tag_levels = "A")

S2_r2 <- clonal_distri_flu11.1 + lambda_heat_flu11.1 + meanN_heat_flu11.1 + meanNb_heat_flu11.1 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 11.1")), tag_levels = list(c("E", "F", "G", "H")))

S2_r3 <- clonal_distri_flu37.7 + lambda_heat_flu37.7 + meanN_heat_flu37.7 + meanNb_heat_flu37.7 + 
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(title = expression(paste(R[0], " = 37.7")), tag_levels = list(c("I", "J", "K", "L")))

figureS2 <- wrap_elements(S2_r1) / wrap_elements(S2_r2) / wrap_elements(S2_r3)

ggsave(figureS2, file = "figureS2.pdf", path = "02_plots", width = 10, height = 9)
