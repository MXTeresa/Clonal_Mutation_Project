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
flu <- data.frame(ClonalMu = 0:3, freq = c(42, 5, 2, 3)) # mu = 1.75
sars <- data.frame(ClonalMu = 0:3, freq = c(35, 4, 0, 0)) # mu = 0.03

# 4. Input parameters ====
mu_flu <- 1.75
mu_sars <- 0.03
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
R0_values_flu <- seq(1.01, 12.01, by = 0.2)
R0_values_sars <- seq(1.01, 16.01, by = 0.2)
lambda_values <- seq(0.01, 10.01, by = 0.2)

# 5. Generate data frame ====
listClonal_flu_mu <- list_clonal_mu(n_values, R0_values_flu, mu_flu, maxMuGen, maxFS, maxClonal)
listClonal_sars_mu <- list_clonal_mu(n_values, R0_values_sars, mu_sars, maxMuGen, maxFS, maxClonal)

df_heatmap_flu_mu <- LL_meanNb_mu.df(flu, listClonal_flu_mu, lambda_values, R0_values_flu, mu_flu, maxMuGen, maxFS, maxIni)
df_heatmap_sars_mu <- LL_meanNb_mu.df(sars, listClonal_sars_mu, lambda_values, R0_values_sars, mu_sars, maxMuGen, maxFS, maxIni)

lambda_flu_hat_mu <-  df_heatmap_flu_mu$lambda[df_heatmap_flu_mu$prob == max(df_heatmap_flu_mu$prob)]
lambda_sars_hat_mu <-  df_heatmap_sars_mu$lambda[df_heatmap_sars_mu$prob == max(df_heatmap_sars_mu$prob)]
R0_flu_hat <-  df_heatmap_flu_mu$R0[df_heatmap_flu_mu$prob == max(df_heatmap_flu_mu$prob)]
R0_sars_hat <-  df_heatmap_sars_mu$R0[df_heatmap_sars_mu$prob == max(df_heatmap_sars_mu$prob)]
min_flu_mu <- max(df_heatmap_flu_mu$prob)-2.995
max_flu_mu <- max(df_heatmap_flu_mu$prob)
min_sars_mu <- max(df_heatmap_sars_mu$prob)-2.995
max_sars_mu <- max(df_heatmap_sars_mu$prob)

melt_distri_flu_mu <- calc_distri(flu, lambda_flu_hat, R0_flu_hat, mu_flu, maxMuGen, maxFS, clonal, maxIni)
melt_distri_sars_mu <- calc_distri(sars, lambda_sars_hat, R0_sars_hat, mu_sars, maxMuGen, maxFS, clonal, maxIni)

# 6. Graphs ====
## 6.1 Flu ====
clonal_distri_flu_mu <- ggplot(melt_distri_flu_mu, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  ylim(0, 1) +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_flu_mu <- ggplot(df_heatmap_flu_mu) +
  geom_contour_fill(aes(x = lambda, y = R0, z = prob), na.fill = TRUE, bins = 1600) +
  geom_hline(aes(yintercept = R0_flu_hat), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_mu, max_flu_mu)) +
  xlim(0, 5) +
  scale_y_continuous(breaks = seq(1, 50, 2), limits = c(1, 12)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.3, y = 6, label = "R[0]==5.01", parse = TRUE, color = "red") +
  annotate("text", x = 1.3, y = 1.6, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu_mu <- ggplot(df_heatmap_flu_mu) +
  geom_contour_fill(aes(x = meanN, y = R0, z = prob), na.fill = TRUE, bins = 1000) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_mu, max_flu_mu)) +
  xlim(0, 5) +
  scale_y_continuous(breaks = seq(1, 50, 2), limits = c(1, 12)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu_mu <- ggplot(df_heatmap_flu_mu) +
  geom_contour_fill(aes(x = meanNb, y = R0, z = prob), na.fill = TRUE, bins = 1000) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_mu, max_flu_mu)) +
  xlim(0, 5) +
  scale_y_continuous(breaks = seq(1, 50, 2), limits = c(1, 12)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.2 Sars ====
clonal_distri_sars_mu <- ggplot(melt_distri_sars_mu, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                               labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  ylim(0, 1) +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_sars_mu <- ggplot(df_heatmap_sars_mu) +
  geom_contour_fill(aes(x = lambda, y = R0, z = prob), na.fill = TRUE, bins = 300) +
  geom_hline(aes(yintercept = R0_sars_hat), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_mu, max_sars_mu),
                                   breaks = c(-14, -15, -16),
                                   labels = c("-14", "-15", "-16")) +
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0, 5.01)) +
  scale_y_continuous(breaks = seq(1, 10, 0.5), limits = c(1, 3)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2.3, y = 1.4, label = "R[0]==1.21", parse = TRUE, color = "red") +
  annotate("text", x = 1.3, y = 2.7, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars_mu <- ggplot(df_heatmap_sars_mu) +
  geom_contour_fill(aes(x = meanN, y = R0, z = prob), na.fill = TRUE, bins = 300) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_mu, max_sars_mu),
                                   breaks = c(-14, -15, -16),
                                   labels = c("-14", "-15", "-16")) +
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0, 5)) +
  scale_y_continuous(breaks = seq(1, 10, 0.5), limits = c(1, 3)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars_mu <- ggplot(df_heatmap_sars_mu) +
  geom_contour_fill(aes(x = meanNb, y = R0, z = prob), na.fill = TRUE, bins = 300) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_mu, max_sars_mu),
                                   breaks = c(-14, -15, -16),
                                   labels = c("-14", "-15", "-16")) +
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0, 5)) +
  scale_y_continuous(breaks = seq(1, 10, 0.5), limits = c(1, 3)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
figureS3 <- clonal_distri_flu_mu + lambda_heat_flu_mu + meanN_heat_flu_mu + meanNb_heat_flu_mu + 
  clonal_distri_sars_mu + lambda_heat_sars_mu + meanN_heat_sars_mu + meanNb_heat_sars_mu + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A")

ggsave(figureS3, file = "figureS3.pdf", path = "02_plots", width = 10, height = 5)

