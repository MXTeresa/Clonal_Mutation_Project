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
mu_values <- seq(0.01, 3.01, by = 0.02)
prop_values <- seq(0.92, 1, by = 0.01)
k <- 0.5

# 5. Generate data frame ====
listClonal_flu <- list_clonal(n_values, R0_flu, mu_values, maxMuGen, maxFS, maxClonal)
listClonal_sars <- list_clonal(n_values, R0_sars, mu_values, maxMuGen, maxFS, maxClonal)

df_heatmap_flu_nb <- LL_meanNb_nb.df(flu, listClonal_flu, seq(0.001, 0.501, by = 0.01), R0_flu, mu_values, maxMuGen, maxFS, maxIni, k)
#df_heatmap_flu_prop <- LL_prop.df(flu, listClonal_flu, prop_values, R0_flu, mu_values, maxMuGen, maxFS, maxIni)

df_heatmap_flu_nb <- read.csv("03_data/df_heatmap_flu_nb.csv")[, -1]
df_heatmap_flu_prop <- read.csv("03_data/df_heatmap_flu_prop.csv")[, -1]

df_heatmap_sars_nb <- LL_meanNb_nb.df(sars, listClonal_sars, seq(0.001, 6.001, by = 0.1), R0_sars, mu_values, maxMuGen, maxFS, maxIni, k)
#df_heatmap_sars_prop <- LL_prop.df(sars, listClonal_sars, prop_values, R0_sars, mu_values, maxMuGen, maxFS, maxIni)

df_heatmap_sars_nb <- read.csv("03_data/df_heatmap_sars_nb.csv")[, -1]
df_heatmap_sars_prop <- read.csv("03_data/df_heatmap_sars_prop.csv")[, -1]

lambda_flu_hat <-  df_heatmap_flu_pois$lambda[df_heatmap_flu_pois$prob == max(df_heatmap_flu_pois$prob)]
lambda_sars_hat <-  df_heatmap_sars_pois$lambda[df_heatmap_sars_pois$prob == max(df_heatmap_sars_pois$prob)]

mu_flu_hat <-  df_heatmap_flu_pois$mu[df_heatmap_flu_pois$prob == max(df_heatmap_flu_pois$prob)]
mu_flu_hat_nb <-  df_heatmap_flu_nb$mu[df_heatmap_flu_nb$prob == max(df_heatmap_flu_nb$prob)]
mu_flu_hat_prop <-  df_heatmap_flu_prop$mu[df_heatmap_flu_prop$prob == max(df_heatmap_flu_prop$prob)]

mu_sars_hat <-  df_heatmap_sars_pois$mu[df_heatmap_sars_pois$prob == max(df_heatmap_sars_pois$prob)]
mu_sars_hat_nb <-  df_heatmap_sars_nb$mu[df_heatmap_sars_nb$prob == max(df_heatmap_sars_nb$prob)]
mu_sars_hat_prop <-  df_heatmap_sars_prop$mu[df_heatmap_sars_prop$prob == max(df_heatmap_sars_prop$prob)]

melt_distri_flu <- calc_distri_lambda0(flu, R0_flu, mu_flu_hat, maxMuGen, maxFS, clonal)
levels(melt_distri_flu$type) <- c("simulated", "poisson")
melt_distri_sars <- calc_distri_lambda0(sars, R0_sars, mu_sars_hat, maxMuGen, maxFS, clonal)
levels(melt_distri_sars$type) <- c("simulated", "poisson")

melt_distri_flu_nb <- calc_distri_lambda0(flu, R0_flu, mu_flu_hat, maxMuGen, maxFS, clonal)
melt_distri_flu_nb <- melt_distri_flu_nb[melt_distri_flu_nb$type == "calculated",]
melt_distri_flu_nb$type <- "nb"

melt_distri_sars_nb <- calc_distri_lambda0(sars, R0_sars, mu_sars_hat, maxMuGen, maxFS, clonal)
melt_distri_sars_nb <- melt_distri_sars_nb[melt_distri_sars_nb$type == "calculated",]
melt_distri_sars_nb$type <- "nb"

melt_distri_flu_prop <- pmf.P_clonal(1, R0_flu, mu_flu_hat, maxMuGen, maxFS, clonal)
melt_distri_flu_prop$type <- "prop"

melt_distri_sars_prop <- pmf.P_clonal(1, R0_sars, mu_sars_hat, maxMuGen, maxFS, clonal)
melt_distri_sars_prop$type <- "prop"

distri_flu <- rbind(melt_distri_flu, melt_distri_flu_nb, melt_distri_flu_prop)
distri_sars <- rbind(melt_distri_sars, melt_distri_sars_nb, melt_distri_sars_prop)

# 6. Graphs ====
## 6.1 Flu ====
clonal_distri_flu_figure5 <- ggplot(distri_flu, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#255668","#008C80", "#BBDF5F", "#FFC000"),
                    labels = list("observed", "est. (Poisson)", "est. (Neg. Bin.)", "est. (Bimodal)")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  ylim(0, 0.954) +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.8,"line"))

nb_heatmap_flu <- ggplot(df_heatmap_flu_nb) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  geom_hline(aes(yintercept = mu_flu_hat_nb), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.016), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda[NB]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(max(df_heatmap_flu_nb$prob)-2.995, 
                                                                        max(df_heatmap_flu_nb$prob))) +
  xlim(0, 2) +
  ylim(0.01, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 1.3, y = 1.8, label = "mu==1.55", parse = TRUE, color = "red") +
  annotate("text", x = 0.5, y = 0.5, label = "lambda %->% 0", parse = TRUE, color = "red")

prop_heatmap_flu <- ggplot(df_heatmap_flu_prop) +
  geom_contour_fill(aes(x = prop, y = mu, z = prob), na.fill = TRUE, bins = 500) +
  geom_hline(aes(yintercept = mu_flu_hat_prop), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.999), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(proportion~of~Nu==1), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(max(df_heatmap_flu_prop$prob)-2.995, 
                                                                        max(df_heatmap_flu_prop$prob))) +
  xlim(0.8, 1.001) +
  ylim(0.01, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 0.88, y = 1.8, label = "mu==1.55", parse = TRUE, color = "red") +
  annotate("text", x = 0.95, y = 0.5, label = "prop==1", parse = TRUE, color = "red")

## 6.2 Sars ====
clonal_distri_sars_figure5 <- ggplot(distri_sars, aes(x = ClonalMu, y = prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#255668","#008C80", "#BBDF5F", "#FFC000"),
                    labels = list("observed", "est. (Poisson)", "est. (Neg. Bin.)", "est. (Bimodal)")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  ylim(0, 0.954) +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.8,"line"))

nb_heatmap_sars <- ggplot(df_heatmap_sars_nb) +
  geom_contour_fill(aes(x = lambda, y = mu, z = prob), na.fill = TRUE, bins = 200) +
  geom_hline(aes(yintercept = mu_sars_hat_nb), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.04), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(lambda[NB]), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(max(df_heatmap_sars_nb$prob)-2.995,
                                                                        max(df_heatmap_sars_nb$prob))) +
  xlim(0, 6) +
  ylim(0.01, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 4.8, y = 0.8, label = "mu==0.53", parse = TRUE, color = "red") +
  annotate("text", x = 1.3, y = 2.4, label = "lambda %->% 0", parse = TRUE, color = "red")

prop_heatmap_sars <- ggplot(df_heatmap_sars_prop) +
  geom_contour_fill(aes(x = prop, y = mu, z = prob), na.fill = TRUE, bins = 120) +
  geom_hline(aes(yintercept = mu_sars_hat_prop), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(aes(xintercept = 0.999), linetype = "dashed", color = "red", linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  labs(x = expression(proportion~of~Nu==1), y = expression(mu), fill = "") +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = FALSE, 
                                   na.value = "transparent", limits = c(max(df_heatmap_sars_prop$prob)-2.995,
                                                                        max(df_heatmap_sars_prop$prob))) +
  xlim(0.8, 1.001) +
  ylim(0.01, 3.5) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm")) +
  annotate("text", x = 0.88, y = 0.8, label = "mu==0.51", parse = TRUE, color = "red") +
  annotate("text", x = 0.95, y = 2.4, label = "prop==1", parse = TRUE, color = "red")

# 7. Put graphs together ====
figure5 <- clonal_distri_flu_figure5 + nb_heatmap_flu + prop_heatmap_flu + 
  clonal_distri_sars_figure5 + nb_heatmap_sars + prop_heatmap_sars + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

ggsave(figure5, file = "figure5.pdf", path = "02_plots", width = 8, height = 5)
