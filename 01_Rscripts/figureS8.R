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
sars_new <- data.frame(ClonalMu = 0:3, freq = c(36, 3, 0, 0)) # R0 = 7.4

# 4. Input parameters ====
R0_sars <- 7.4
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 3.51, by = 0.02)
lambda_values <- seq(0.01, 4.01, by = 0.1)

# 5. Generate data frame ====
listClonal_sars_new <- list_clonal(n_values, R0_sars, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_sars_new <- LL_meanNb.df(sars_new, listClonal_sars_new, lambda_values, R0_sars, mu_values, maxMuGen, maxFS, maxIni)

lambda_sars_hat_new <-  df_heatmap_sars_new$lambda[df_heatmap_sars_new$prob == max(df_heatmap_sars_new$prob)]
mu_sars_hat_new <-  df_heatmap_sars_new$mu[df_heatmap_sars_new$prob == max(df_heatmap_sars_new$prob)]
min_sars_new <- max(df_heatmap_sars_new$prob)-2.995
max_sars_new <- max(df_heatmap_sars_new$prob)

melt_distri_sars_new <- calc_distri_lambda0(sars_new, R0_sars, mu_sars_hat_new, maxMuGen, maxFS, clonal)

# 6. Graphs ====
clonal_distri_sars_new <- ggplot(melt_distri_sars_new, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = "none") +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_sars_new <- ggplot(df_heatmap_sars_new) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 150) +
  geom_hline(aes(yintercept = mu_sars_hat_new), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_new, max_sars_new),
                                   breaks = c(-12, -13, -14),
                                   labels = c("-12", "-13", "-14")) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 2, y = 0.7, label = "mu==0.43", parse = TRUE, color = "red") +
  annotate("text", x = 1.1, y = 2.4, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars_new <- ggplot(df_heatmap_sars_new) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 150) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_new, max_sars_new),
                                   breaks = c(-12, -13, -14),
                                   labels = c("-12", "-13", "-14")) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars_new <- ggplot(df_heatmap_sars_new) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 150) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars_new, max_sars_new),
                                   breaks = c(-12, -13, -14),
                                   labels = c("-12", "-13", "-14")) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
figureS8 <- clonal_distri_sars + lambda_heat_sars + meanN_heat_sars + meanNb_heat_sars + 
  clonal_distri_sars_new + lambda_heat_sars_new + meanN_heat_sars_new + meanNb_heat_sars_new + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A")

ggsave(figureS8, file = "figureS8.pdf", path = "02_plots", width = 10, height = 5)
