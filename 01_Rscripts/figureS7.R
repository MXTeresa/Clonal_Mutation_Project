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
flu_new <- data.frame(ClonalMu = 0:3, freq = c(42, 5, 3, 2)) # R0 = 11.1

# 4. Input parameters ====
R0_flu <- 11.1
clonal <- 3
maxFS <- 50
maxMuGen <- 50
maxClonal <- 5
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 3.51, by = 0.02)
lambda_values <- seq(0.01, 4.01, by = 0.1)

# 5. Generate data frame ====
listClonal_flu_new <- list_clonal(n_values, R0_flu, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_flu_new <- LL_meanNb.df(flu_new, listClonal_flu_new, lambda_values, R0_flu, mu_values, maxMuGen, maxFS, maxIni)

lambda_flu_hat_new <-  df_heatmap_flu_new$lambda[df_heatmap_flu_new$prob == max(df_heatmap_flu_new$prob)]
mu_flu_hat_new <-  df_heatmap_flu_new$mu[df_heatmap_flu_new$prob == max(df_heatmap_flu_new$prob)]
min_flu_new <- max(df_heatmap_flu_new$prob)-2.995
max_flu_new <- max(df_heatmap_flu_new$prob)

melt_distri_flu_new <- calc_distri_lambda0(flu_new, R0_flu, mu_flu_hat_new, maxMuGen, maxFS, clonal)

# 6. Graphs ====
clonal_distri_flu_new <- ggplot(melt_distri_flu_new, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_flu_new <- ggplot(df_heatmap_flu_new) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  geom_hline(aes(yintercept = mu_flu_hat_new), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_new, max_flu_new)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3, y = 1.7, label = "mu==1.45", parse = TRUE, color = "red") +
  annotate("text", x = 1.1, y = 0.5, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu_new <- ggplot(df_heatmap_flu_new) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_new, max_flu_new)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu_new <- ggplot(df_heatmap_flu_new) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu_new, max_flu_new)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
figureS7 <- clonal_distri_flu + lambda_heat_flu + meanN_heat_flu + meanNb_heat_flu + 
  clonal_distri_flu_new + lambda_heat_flu_new + meanN_heat_flu_new + meanNb_heat_flu_new + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A")

ggsave(figureS7, file = "figureS7.pdf", path = "02_plots", width = 10, height = 5)
