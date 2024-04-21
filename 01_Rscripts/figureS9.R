# 1. Source functions ====
source("01_Rscripts/ClonalMuFunctions.R")

# 2. Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(colorspace)
library(showtext)
library(forcats)
theme_set(theme_minimal() +
          theme(axis.ticks = element_line(colour = "grey50"),
                axis.line = element_line(colour = "grey50"),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 15),
                legend.text=element_text(size = 12)))
showtext_auto()

# 3. Input parameters ====
n <- 2
R0 <- 1.2
mu <- 0.2
maxOff <- 8
maxFS <- 100
maxMuGen <- 100
maxClonal <- 5

# 4. Generate colors needed ====
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col2 <- sample(col_vector, 100, replace = TRUE)
col4 <- c(col2, "#008C80", "lightgrey")

# 5. Generate data frame ====
three_pmf <- pmf.alloff(R0,mu,maxOff)
df_FSnormal <- pmf.FSn(n,R0,mu,maxFS)[1:50,]
df_pmfMu_n <- pmf.P_mu_n(n,R0,mu,maxMuGen,maxFS)[1:31,]
df_pmfMu_est <- pmf.P_mu_est(n,R0,mu,maxMuGen,maxFS)[1:11,]
df_S_data <- pmf.mu_est_clonal(n,R0,mu,maxMuGen,maxFS)
df_clonal <- pmf.P_clonal(n,R0,mu,maxMuGen,maxFS,maxClonal)

Px_ref <- Px_analytical(n,R0)

# 6. Graphs ====
offsp <- ggplot(three_pmf, aes(size, prob)) +
  geom_col(aes(fill = type), position = "dodge") +
  scale_fill_manual(values = c("#255668", "#34AF7C", "#BBDF5F"),
                    labels = list("overall", "wild-type", "mutant")) +
  labs(x = "offspring number", y = "probability", fill = "") +
  theme(legend.position = c(0.8, 0.9)) +
  scale_x_continuous(breaks = seq(0, 10, 2))

FinalSize <- ggplot(df_FSnormal, aes(FinalSize, prob_normal)) +
  geom_col(fill = "#008C80") +
  labs(y = "probability", x = "final size")

MuGen <- ggplot(df_pmfMu_n, aes(MuLinGen, prob)) +
  geom_col(fill = "#008C80") +
  labs(y = "probability", x = "# mutant lineages generated")

MuEst <- ggplot(df_pmfMu_est, aes(MuLinEst, prob)) +
  geom_col(fill = "#008C80") +
  labs(y = "probability", x = "# mutant lineages established") +
  scale_x_continuous(breaks = seq(0, 10, 2))

MuEstClonal <- ggplot(df_S_data, aes(fill = fct_rev(fct_inorder(MuLinEst)), 
                                     color = "black",
                                     y = prob_norm, x = Clonal)) + 
  geom_bar(position = "stack", stat = "identity") +
  geom_segment(aes(x = 0, y = Px_ref, xend = 1.45, yend = Px_ref, linetype = "dashed")) +
  scale_linetype_manual(values = "dashed") +
  scale_fill_manual(values = col4) +
  scale_color_manual(values = "black") +
  theme(legend.position = "none") +
  labs(x = "outcome", y = "probability") +
  scale_x_discrete(labels = c(expression(P[x]), expression(P[0]), expression(P[1^{"+"}])))

ClonalNorm <- ggplot(df_clonal, aes(fct_inorder(ClonalMu), prob)) +
  geom_col(fill = "#008C80") +
  labs(y = "probability", x = "# clonal variants")


# 7. Put graphs together ====
figureS9 <- offsp + FinalSize + MuGen + MuEst + MuEstClonal + ClonalNorm +
  plot_annotation(tag_levels = "A")

ggsave(figureS9, file = "figureS9.pdf", path = "02_plots", width = 10, height = 7)
