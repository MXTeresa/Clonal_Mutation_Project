# 1. Source functions ====
source("01_Rscripts/ClonalMuFunctions.R")

# 2. Load packages ====
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(showtext)
library(metR)
library(scales)
theme_set(theme_minimal() +
            theme(axis.ticks = element_line(colour = "grey50"),
                  axis.line = element_line(colour = "grey50"),
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 14)))
showtext_auto()

# 3. Input parameters and data ====
R0_flu <- 11.1
sites_flu <- 12714
mu_flu <- 1.8*10^(-4)*sites_flu

R0_sars <- 7.4
sites_sars <- 29903
mu_sars <- 1.3*10^(-6)*sites_sars

maxMuGen <- 100
maxFS <- 100
maxIni <- 10
maxClonal <- 10

# 4. Generate data frame ====
## 4.1 Flu ====
df_P_trans1_flu <- P_transAndFix(1, R0_flu, maxIni)
df_P_trans2_flu <- P_transAndFix(2, R0_flu, maxIni)
df_P_trans4_flu <- P_transAndFix(4, R0_flu, maxIni)
df_P_trans6_flu <- P_transAndFix(6, R0_flu, maxIni)
df_P_trans10_flu <- P_transAndFix(10, R0_flu, maxIni)
df_P_trans1_flu$donor_freq <- 1
df_P_trans2_flu$donor_freq <- 2
df_P_trans4_flu$donor_freq <- 4
df_P_trans6_flu$donor_freq <- 6
df_P_trans10_flu$donor_freq <- 10
df_P_trans_flu_combined <- rbind(df_P_trans1_flu, df_P_trans2_flu, df_P_trans4_flu, df_P_trans6_flu, df_P_trans10_flu)

## 4.2 Sars ====
df_P_trans1_sars <- P_transAndFix(1, R0_sars, maxIni)
df_P_trans2_sars <- P_transAndFix(2, R0_sars, maxIni)
df_P_trans4_sars <- P_transAndFix(4, R0_sars, maxIni)
df_P_trans6_sars <- P_transAndFix(6, R0_sars, maxIni)
df_P_trans10_sars <- P_transAndFix(10, R0_sars, maxIni)
df_P_trans1_sars$donor_freq <- 1
df_P_trans2_sars$donor_freq <- 2
df_P_trans4_sars$donor_freq <- 4
df_P_trans6_sars$donor_freq <- 6
df_P_trans10_sars$donor_freq <- 10
df_P_trans_sars_combined <- rbind(df_P_trans1_sars, df_P_trans2_sars, df_P_trans4_sars, df_P_trans6_sars, df_P_trans10_sars)

# 5. Graphs ====
## 5.1 Flu ====
flu_trans_fix <- ggplot(df_P_trans_flu_combined) +
  geom_line(aes(x = N, y = prob, color = factor(donor_freq))) +
  geom_point(aes(x = N, y = prob, color = factor(donor_freq))) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_log10(limits = c(1e-12, 1e0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                name = "probability") +
  scale_color_manual(values = c("#E5815E", "#FFC000", "#BBDF5F", "#008C80", "#255668"))

sars_trans_fix <- ggplot(df_P_trans_sars_combined) +
  geom_line(aes(x = N, y = prob, color = factor(donor_freq))) +
  geom_point(aes(x = N, y = prob, color = factor(donor_freq))) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_log10(limits = c(1e-12, 1e0),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                name = "probability") +
  scale_color_manual(values = c("#E5815E", "#FFC000", "#BBDF5F", "#008C80", "#255668"),
                     labels = list("1%", "2%", "4%", "6%", "10%"),
                     name = "iSNV frequency\nin donor")

# 6. Put graphs together ====
figureS10 <- flu_trans_fix + sars_trans_fix + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave(figureS10, file = "figureS10.pdf", path = "02_plots", width = 9, height = 4)
