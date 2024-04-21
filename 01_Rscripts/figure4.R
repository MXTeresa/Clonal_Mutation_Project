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
flu <- data.frame(ClonalMu = 0:3, freq = c(42, 5, 2, 3)) # R0 = 11.1
sars <- data.frame(ClonalMu = 0:3, freq = c(35, 4, 0, 0)) # R0 = 7.4

# 4. Input parameters ====
R0_flu <- 11.1
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
listClonal_flu <- list_clonal(n_values, R0_flu, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sars <- list_clonal(n_values, R0_sars, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_flu <- LL_meanNb.df(flu, listClonal_flu, lambda_values, R0_flu, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_sars <- LL_meanNb.df(sars, listClonal_sars, lambda_values, R0_sars, mu_values, maxMuGen, maxFS, maxIni)

lambda_flu_hat <-  df_heatmap_flu$lambda[df_heatmap_flu$prob == max(df_heatmap_flu$prob)]
lambda_sars_hat <-  df_heatmap_sars$lambda[df_heatmap_sars$prob == max(df_heatmap_sars$prob)]
mu_flu_hat <-  df_heatmap_flu$mu[df_heatmap_flu$prob == max(df_heatmap_flu$prob)]
mu_sars_hat <-  df_heatmap_sars$mu[df_heatmap_sars$prob == max(df_heatmap_sars$prob)]
min_flu <- max(df_heatmap_flu$prob)-2.995
max_flu <- max(df_heatmap_flu$prob)
min_sars <- max(df_heatmap_sars$prob)-2.995
max_sars <- max(df_heatmap_sars$prob)

melt_distri_flu <- calc_distri_lambda0(flu, R0_flu, mu_flu_hat, maxMuGen, maxFS, clonal)
melt_distri_sars <- calc_distri_lambda0(sars, R0_sars, mu_sars_hat, maxMuGen, maxFS, clonal)

## 5.1 lambda approaches 0 ====

df_flu_lambda0 <- df_heatmap_flu_prop[df_heatmap_flu_prop$prop == 1, ]
names(df_flu_lambda0) <- c("lambda", names(df_flu_lambda0)[-1])
df_flu_lambda0$lambda <- 0
df_flu_lambda0$meanN <- 1
df_flu_lambda0$meanNb <- 1
df_heatmap_flu_fin <- rbind(df_heatmap_flu, df_flu_lambda0)

df_sars_lambda0 <- df_heatmap_sars_prop[df_heatmap_sars_prop$prop == 1, ]
names(df_sars_lambda0) <- c("lambda", names(df_sars_lambda0)[-1])
df_sars_lambda0$lambda <- 0
df_sars_lambda0$meanN <- 1
df_sars_lambda0$meanNb <- 1
df_heatmap_sars_fin <- rbind(df_heatmap_sars, df_sars_lambda0)

# 6. Graphs ====
## 6.1 Flu ====
clonal_distri_flu <- ggplot(melt_distri_flu, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_flu <- ggplot(df_heatmap_flu_fin) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  geom_hline(aes(yintercept = mu_flu_hat), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu, max_flu)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3, y = 1.8, label = "mu==1.55", parse = TRUE, color = "red") +
  annotate("text", x = 1.1, y = 0.5, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_flu <- ggplot(df_heatmap_flu_fin) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 400) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu, max_flu)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_flu <- ggplot(df_heatmap_flu_fin) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 400) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_flu, max_flu)) +
  xlim(0, 4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.2 Sars ====
clonal_distri_sars <- ggplot(melt_distri_sars, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = c(0.75, 0.93),
        legend.key.size = unit(0.8,"line"))

lambda_heat_sars <- ggplot(df_heatmap_sars_fin) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 350) +
  geom_hline(aes(yintercept = mu_sars_hat), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.03), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars, max_sars)) +
  xlim(0,4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 3.2, y = 0.8, label = "mu==0.51", parse = TRUE, color = "red") +
  annotate("text", x = 1.1, y = 2.4, label = "lambda %->% 0", parse = TRUE, color = "red")

meanN_heat_sars <- ggplot(df_heatmap_sars_fin) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 200) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars, max_sars)) +
  xlim(0,4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_sars <- ggplot(df_heatmap_sars_fin) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 200) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_sars, max_sars)) +
  xlim(0,4) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
figure4 <- clonal_distri_flu + lambda_heat_flu + meanN_heat_flu + meanNb_heat_flu + 
  clonal_distri_sars + lambda_heat_sars + meanN_heat_sars + meanNb_heat_sars + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A")

ggsave(figure4, file = "figure4.pdf", path = "02_plots", width = 10, height = 5)

