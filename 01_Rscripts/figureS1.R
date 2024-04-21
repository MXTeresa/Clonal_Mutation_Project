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
                  axis.title = element_text(size = 15)))
showtext_auto()

# 3. Input data sets ====
sim <- data.frame(ClonalMu = 0:8, freq = c(78, 12, 7, 1, 0, 1, 0, 1, 0)) # lambda = 2.1, R0 = 1.6, mu = 0.4

# 4. Input parameters ====
lambda <- 2.1
mu <- 0.4
maxFS <- 50
maxMuGen <- 50
maxClonal <- 8
maxIni <- 10
n_values <- 1:maxIni
mu_values <- seq(0.01, 1.2, by = 0.01)
lambda_values <- seq(0.01, 7.01, by = 0.1)

# 5. Generate and modify data frame ====
## 5.1 Profile likelihood ====
df_heatmap_profile <- read.csv("03_data/df_heatmap_profile.csv")[,-1]
lambda_profile_hat <- df_heatmap_profile$lambda[df_heatmap_profile$prob == max(df_heatmap_profile$prob)]
R0_profile_hat <- df_heatmap_profile$R0[df_heatmap_profile$prob == max(df_heatmap_profile$prob)]
mu_profile_hat <- df_heatmap_profile$mu[df_heatmap_profile$prob == max(df_heatmap_profile$prob)]

df_lambda_profile <- df_heatmap_profile %>%
  group_by(lambda) %>%
  filter(prob == max(prob))

df_R0_profile <- df_heatmap_profile %>%
  group_by(R0) %>%
  filter(prob == max(prob))

df_mu_profile <- df_heatmap_profile %>%
  group_by(mu) %>%
  filter(prob == max(prob))

## 5.2 Heat map ====
listClonal_sim1.3 <- list_clonal(n_values, 1.3, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sim1.6 <- list_clonal(n_values, 1.6, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sim3 <- list_clonal(n_values, 3, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap1.3 <- LL_meanNb.df(sim, listClonal_sim1.3, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)
lambda_hat1.3 <- df_heatmap1.3$lambda[df_heatmap1.3$prob == max(df_heatmap1.3$prob)]
mu_hat1.3 <- df_heatmap1.3$mu[df_heatmap1.3$prob == max(df_heatmap1.3$prob)]

df_heatmap1.6 <- LL_meanNb.df(sim, listClonal_sim1.6, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)
lambda_hat1.6 <- df_heatmap1.6$lambda[df_heatmap1.6$prob == max(df_heatmap1.6$prob)]
mu_hat1.6 <- df_heatmap1.6$mu[df_heatmap1.6$prob == max(df_heatmap1.6$prob)]

df_heatmap3 <- LL_meanNb.df(sim, listClonal_sim3, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni)
lambda_hat3 <- df_heatmap3$lambda[df_heatmap3$prob == max(df_heatmap3$prob)]
mu_hat3 <- df_heatmap3$mu[df_heatmap3$prob == max(df_heatmap3$prob)]

# 6. Graphs ====
## 6.1 Profile likelihood ====
lambda_profile <- ggplot(df_lambda_profile, aes(x = lambda, y = prob)) +
  geom_line(color = "blue") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  labs(x = expression(lambda), y = "log-likelihood") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 10, 1)) +
  scale_y_continuous(limits = c(-84.01628, -80.8), breaks = seq(-100, 0, 1))

R0_profile <- ggplot(df_R0_profile, aes(x = R0, y = prob)) +
  geom_line(color = "blue") +
  geom_vline(aes(xintercept = 1.6), color = "black") +
  labs(x = expression(R[0]), y = "log-likelihood") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(limits = c(1.3, 3), breaks = seq(1, 5, 0.3)) +
  scale_y_continuous(limits = c(-84.01628, -80.8), breaks = seq(-100, 0, 1))

mu_profile <- ggplot(df_mu_profile, aes(x = mu, y = prob)) +
  geom_line(color = "blue") +
  geom_vline(aes(xintercept = 0.4), color = "black") +
  labs(x = expression(mu), y = "log-likelihood") +
  coord_cartesian(expand = FALSE) +
  scale_x_continuous(limits = c(0.2, 0.9), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(-84.01628, -80.8), breaks = seq(-100, 0, 1))

## 6.2 Heat map ====
heatmap_sim1.3 <- ggplot(df_heatmap1.3) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  geom_hline(aes(yintercept = 0.4), color = "black") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  geom_hline(aes(yintercept = mu_hat1.3), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(aes(xintercept = lambda_hat1.3), linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", 
                                   limits = c(-83.95599, -80.96099)) +
  scale_x_continuous(limits = c(0, 6.05), breaks = seq(0, 10, 1)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 1, y = 0.32, label = "mu==0.24", parse = TRUE, color = "red") +
  annotate("text", x = 4, y = 0.7, label = "lambda==5.01", parse = TRUE, color = "red")

heatmap_sim1.6 <- ggplot(df_heatmap1.6) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 850) +
  geom_hline(aes(yintercept = 0.4), color = "black") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  geom_hline(aes(yintercept = mu_hat1.6), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(aes(xintercept = lambda_hat1.6), linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", 
                                   limits = c(-83.95599, -80.96099)) +
  scale_x_continuous(limits = c(0, 6.05), breaks = seq(0, 10, 1)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 5.1, y = 0.3, label = "mu==0.39", parse = TRUE, color = "red") +
  annotate("text", x = 3.4, y = 1.1, label = "lambda==2.3", parse = TRUE, color = "red")

heatmap_sim3 <- ggplot(df_heatmap3) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 1200) +
  geom_hline(aes(yintercept = 0.4), color = "black") +
  geom_vline(aes(xintercept = 2.1), color = "black") +
  geom_hline(aes(yintercept = mu_hat3), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(aes(xintercept = 0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", 
                                   limits = c(-83.95599, -80.96099)) +
  scale_x_continuous(limits = c(0, 6.05), breaks = seq(0, 10, 1)) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 4, y = 0.74, label = "mu==0.81", parse = TRUE, color = "red") +
  annotate("text", x = 1, y = 0.3, label = "lambda %->% 0", parse = TRUE, color = "red")

# 7. Put graphs together ====
figureS1 <- lambda_profile + R0_profile + mu_profile +
  heatmap_sim1.3 + heatmap_sim1.6 + heatmap_sim3 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

ggsave(figureS1, file = "figureS1.pdf", path = "02_plots", width = 10, height = 5)
