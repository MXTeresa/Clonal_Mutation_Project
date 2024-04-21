# 1. Source functions ====
source("01_Rscripts/ClonalMuFunctions.R")

# 2. Load packages ====
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(showtext)
library(metR)
theme_set(theme_minimal() +
          theme(axis.ticks = element_line(colour = "grey50"),
                axis.line = element_line(colour = "grey50"),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 15),
                legend.text = element_text(size = 12)))
showtext_auto()

# 3. Input data sets ====
sim <- data.frame(ClonalMu = 0:8, freq = c(78, 12, 7, 1, 0, 1, 0, 1, 0)) # lambda = 2.1, R0 = 1.6, mu = 0.4

# 4. Input parameters ====
lambda <- 2.1
R0 <- 1.6
mu <- 0.4
maxFS <- 50
maxMuGen <- 50
maxClonal <- 8
maxIni <- 15
n_values <- 1:maxIni
mu_values <- seq(0.01, 1.2, by = 0.01)
R0_values <- seq(1.2, 2.8, by = 0.02)
lambda_values <- seq(0.01, 7.01, by = 0.1)

# 5. Generate and modify data frame ====
## 5.1 List needed ====
listClonal_sim <- list_clonal(n_values, R0, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sim_mu <- list_clonal_mu(n_values, R0_values, mu, maxMuGen, maxFS, maxClonal)
## 5.2 Clonal mutation heat map ====
clonal_0 <- pmf.clonal_params(0, listClonal_sim, n_values, R0, mu_values, maxMuGen, maxFS)
clonal_1 <- pmf.clonal_params(1, listClonal_sim, n_values, R0, mu_values, maxMuGen, maxFS)
clonal_2 <- pmf.clonal_params(2, listClonal_sim, n_values, R0, mu_values, maxMuGen, maxFS)
## 5.3 Adjusted Poisson ====
df_pois <- pmf.ini_size(lambda, R0, 8)
## 5.4 Log likelihood heat map ====
df_heatmap <- LL_meanNb.df(sim, listClonal_sim, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)
df_heatmap_mu <- LL_meanNb_mu.df(sim, listClonal_sim_mu,  seq(0.01, 14.01, by = 0.4), R0_values, mu, maxMuGen, maxFS, maxIni)
### 5.4.1 MLE ====
lambda_hat <- df_heatmap$lambda[df_heatmap$prob == max(df_heatmap$prob)]
mu_hat <- df_heatmap$mu[df_heatmap$prob == max(df_heatmap$prob)]
lambda_hat_mu <- df_heatmap_mu$lambda[df_heatmap_mu$prob == max(df_heatmap_mu$prob)]
R0_hat_mu <- df_heatmap_mu$R0[df_heatmap_mu$prob == max(df_heatmap_mu$prob)]
### 5.4.2 95% confidence interval ====
max <- max(df_heatmap$prob)
min <- max(df_heatmap$prob)-2.995 #chi square degrees of freedom 1
max_mu <- max(df_heatmap_mu$prob)
min_mu <- max(df_heatmap_mu$prob)-2.995
## 5.5 Simulated and calculated distribution ====
melt_distri <- calc_distri(sim, lambda_hat, R0, mu_hat, maxMuGen, maxFS, maxClonal, maxIni)

# 6. Graphs ====
pois_adjust <- ggplot(df_pois, aes(N, prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("Poisson", "Poisson | infection")) +
  labs(x = expression(paste("initial viral population size ", Nu)), y = "probability", fill = "") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  theme(legend.position = c(0.82, 0.93)) +
  theme(legend.text = element_text(size = 14))

clonal_distri <- ggplot(melt_distri, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"),
                    labels = list("observed", "estimated")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  theme(legend.position = c(0.82, 0.93)) +
  theme(legend.text = element_text(size = 14))

## 6.1 Clonal heatmap ====
clonal_0_heat <- ggplot(clonal_0, aes(x = n, y = mu)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE) +
  labs(x = expression(Nu), y = expression(mu), fill = "") +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  theme(legend.position = "right", legend.direction="vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm"))

clonal_1_heat <- ggplot(clonal_1, aes(x = n, y = mu)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE) +
  labs(x = expression(Nu), y = expression(mu), fill = "") +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  theme(legend.position = "right", legend.direction="vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm"))

clonal_2_heat <- ggplot(clonal_2, aes(x = n, y = mu)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE) +
  labs(x = expression(Nu), y = expression(mu), fill = "") +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm"))

## 6.2 Estimate mu ====
lambda_heat <- ggplot(df_heatmap) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 850) +
  geom_hline(aes(yintercept = 0.4), color = "black") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  geom_hline(aes(yintercept = mu_hat), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(aes(xintercept = lambda_hat), linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min, max)) +
  xlim(0,7) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 6, y = 0.3, label = "mu==0.39", parse = TRUE, color = "red") +
  annotate("text", x = 3.4, y = 1.1, label = "lambda==2.3", parse = TRUE, color = "red")
  

meanN_heat <- ggplot(df_heatmap) +
  geom_contour_fill(aes(x = meanN, y = mu, z = prob), na.fill = TRUE, bins = 850) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min, max)) +
  xlim(0,7) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat <- ggplot(df_heatmap) +
  geom_contour_fill(aes(x = meanNb, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min, max)) +
  xlim(0,7) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

## 6.3 Estimate R0 ====
lambda_heat_mu <- ggplot(df_heatmap_mu) +
  geom_contour_fill(aes(x = lambda, y = R0, z = prob), na.fill = TRUE, bins = 850) +
  geom_hline(aes(yintercept = 1.6), color = "black") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  geom_hline(aes(yintercept = R0_hat_mu), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(aes(xintercept = lambda_hat_mu), linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_mu, max_mu)) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
  ylim(1, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 7.5, y = 1.7, label = "R[0]==1.6", parse = TRUE, color = "red") +
  annotate("text", x = 4, y = 2.5, label = "lambda==2.4", parse = TRUE, color = "red")


meanN_heat_mu <- ggplot(df_heatmap_mu) +
  geom_contour_fill(aes(x = meanN, y = R0, z = prob), na.fill = TRUE, bins = 500) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_mu, max_mu)) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
  ylim(1, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

meanNb_heat_mu <- ggplot(df_heatmap_mu) +
  geom_contour_fill(aes(x = meanNb, y = R0, z = prob), na.fill = TRUE, bins = 500) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(bar(Nu)[italic(b)]), y = expression(R[0]), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(min_mu, max_mu)) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
  ylim(1, 3) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) 

# 7. Put graphs together ====
figure3 <- (pois_adjust + clonal_distri)/
  (clonal_0_heat + clonal_1_heat + clonal_2_heat)/
  (lambda_heat + meanN_heat + meanNb_heat) /
  (lambda_heat_mu + meanN_heat_mu + meanNb_heat_mu) +
  plot_annotation(tag_levels = "A")

ggsave(figure3, file = "figure3.pdf", path = "02_plots", width = 10, height = 12)

