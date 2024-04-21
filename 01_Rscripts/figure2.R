# 1. Source functions ====
source("01_Rscripts/ClonalMuFunctions.R")

# 2. Load packages ====
library(tidyverse)
library(plotly)
library(reshape)
library(magick)
library(ggplot2)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 17),
                  legend.text = element_text(size = 17)))

# 3. Input parameters ====
n <- 1
N <- 2
R0 <- 1.2
mu <- 0.2
R0_values <- seq(1.1, 7.1, by = 0.1)
mu_values <- seq(0.001, 0.101, by = 0.002)
maxMuGen <- 100
maxFS <- 100
maxClonal <- 5
sample_data <- data.frame(ClonalMu = 0:14, freq = c(2257, 850, 422, 239, 116, 57, 37, 7, 8, 3, 2, 1, 0, 1, 0))
# N = 2, R0 = 1.2, mu = 0.2

# 4. Generate data frame ====
df_mcBozic <- meanClonalBozic.df(R0_values, mu_values)
df_mc <- meanClonal.df(n, R0_values, mu_values, maxMuGen, maxFS, maxClonal)
df_cl <- pmf.P_clonal(N, R0, mu, maxMuGen, maxFS, 14)
df_sim_calc <- df_cl
df_sim_calc$sim_prop <- sample_data$freq/sum(sample_data$freq)
df_sim_calc$ClonalMu <- as.numeric(paste(df_sim_calc$ClonalMu))
df_sim_calc <- melt(df_sim_calc, id = "ClonalMu")
names(df_sim_calc) <- c("ClonalMu", "type", "prob")
df_sim_calc$type <- factor(df_sim_calc$type, levels = c("sim_prop", "prob"))
df1 <- df_sim_calc[1:11, ]
df2 <- df_sim_calc[16:26, ]
df_cl_plot <- rbind(df1, df2)

df_mc$log10mc <- log10(df_mc$mc)

# 5. Graphs
mcBozic_heat <- plot_ly(df_mcBozic, x = ~R0, y = ~mu, z = ~logmc, 
                        type = 'mesh3d', 
                        intensity = ~logmc,
                        colors = c("#255668", "#255668", "#34AF7C", "#BBDF5F", "#EDEF5C", "#EDEF5C")) %>% 
  layout(scene = list(xaxis = list(tickvals = list(1, 3, 5, 7), title = "R₀", range = c(1,7)), 
                      yaxis = list(nticks = 4, title = "µ", range = c(0,0.1)), 
                      zaxis = list(nticks = 8, range = c(-4.5,0), 
                                   title = "mean # of clonal variants(log₁₀)"),
                      camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))) %>%
  hide_colorbar()

mc_heat <- plot_ly(df_mc, x = ~R0, y = ~mu, z = ~log10mc, 
                   type = 'mesh3d', 
                   intensity = ~logmc,
                   colors = c("#255668", "#255668", "#34AF7C", "#BBDF5F", "#EDEF5C", "#EDEF5C")) %>% 
  layout(scene = list(xaxis = list(tickvals = list(1, 3, 5, 7), title = "R₀", range = c(1,7)), 
                      yaxis = list(nticks = 4, title = "µ", range = c(0,0.1)), 
                      zaxis = list(nticks = 8, range = c(-4.5,0), 
                                   title = "mean # of clonal variants(log₁₀)"),
                      camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))) %>%
  hide_colorbar()

sim_calc <- ggplot(df_cl_plot, aes(ClonalMu, prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#008C80", "#BBDF5F"), 
                    labels = list("simulated", "derived")) +
  labs(x = "# clonal variants", y = "proportion", fill = "") +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  theme(legend.position = c(0.78, 0.93))

ggsave(sim_calc, file = "figure2c.pdf", path = "02_plots", width = 4, height = 5)
