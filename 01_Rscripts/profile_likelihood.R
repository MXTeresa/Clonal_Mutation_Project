# 1. Source functions ====
source("ClonalMuFunctions.R")

# 3. Input data sets ====
sample <- data.frame(ClonalMu = 0:8, freq = c(78, 12, 7, 1, 0, 1, 0, 1, 0)) # lambda = 2.1, R0 = 1.6, mu = 0.4

# 4. Input parameters ====
lambda <- 2.1
R0 <- 1.6
mu <- 0.4
maxFS <- 100
maxMuGen <- 100
maxClonal <- 10
maxIni <- 10
n_values <- 1:maxIni
R0_values <- seq(1.01, 3.01, by = 0.01)
mu_values <- seq(0.01, 0.9, by = 0.01)
lambda_values <- seq(0.01, 20.01, by = 0.1)

# 5. Generate and modify data frame ====
df_heatmap_profile <- data.frame()

for (r0 in R0_values) {
  listClonal <- list_clonal(n_values, r0, mu_values, maxMuGen, maxFS, maxClonal)
  df_heatmap_temp <- LL_meanNb.df(sample, listClonal, lambda_values, r0, mu_values, maxMuGen, maxFS, maxIni)
  df_heatmap_temp$R0 <- r0
  df_heatmap_profile <- rbind(df_heatmap_profile, df_heatmap_temp)
}

#df_heatmap_profile <- read.csv("03_data/df_heatmap_profile.csv")[,-1]
write.csv(df_heatmap_profile, "df_heatmap_profile.csv")
