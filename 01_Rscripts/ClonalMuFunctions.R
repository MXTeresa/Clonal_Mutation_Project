# 1. Load packages ====
library(parallel)
library(forcats)
library(tidyverse)
library(reshape)

# 2. Offspring distribution pmf ====
pmf.alloff <- function(R0, mu, maxOff) {
  k_overall <- 1
  k_wt <- exp(-mu)
  k_mu <- (1 - exp(-mu))
  p <- 1 / (R0 + 1)
  pmf <- data.frame(
    x = 0:maxOff,
    overall = dnbinom(0:maxOff, k_overall, p),
    wildtype = dnbinom(0:maxOff, k_wt, p),
    mutant = dnbinom(0:maxOff, k_mu, p)
  )
  pmf <- melt(pmf, id.vars = "x")
  names(pmf) <- c("size", "type", "prob")
  return(pmf)
}


# 3. Probability of wild-type establishment ====
P_wt_est <- function(n, R0, mu) {
  lhs <- function(p) {
    p
  }
  rhs <- function(p, R0, mu) {
    rhs <- 1 / (1 + (R0 * exp(-mu) * (1 - p)) / (exp(-mu)))^exp(-mu)
    return(rhs)
  }
  uni <- try(uniroot(function(p) rhs(p, R0, mu) - lhs(p), c(0, 0.9), extendInt = "yes"), silent = TRUE)
  uni <- if (inherits(uni, "try-error")) 1 else uni$root
  if (uni < 1 & uni > 0) {
    ext_1 <- uni
  } else {
    ext_1 <- 1
  }
  est_n <- 1 - ext_1^n
  return(est_n)
}


# 4. Final size ====
FS <- function(y, R0, mu) {
  k <- exp(-mu)
  meanR0 <- R0 * k
  if (y == 1) {
    prob <- 1 / (1 + meanR0 / k)^k
  } else if (y > 1) {
    prod <- 1
    for (j in 0:(y - 2)) {
      prod <- prod * (j / k + y)
    }
    prob <- prod / factorial(y) * (k / (meanR0 + k))^(k * y) * (meanR0 * k / (meanR0 + k))^(y - 1)
  } else {
    prob <- 0
  }
  return(prob)
}


pmf.FS <- function(R0, mu, maxFS) {
  prob <- rep(0, maxFS)
  for (i in 1:maxFS) {
    prob[i] <- FS(i, R0, mu)
  }
  names(prob) <- c(1:maxFS)
  df <- data.frame(FinalSize = names(prob), prob)
  df$FinalSize <- as.numeric(df$FinalSize)
  return(df)
}


FSn <- function(n, y, R0, mu) {
  prob <- 0
  if (y >= n) {
    if (n == 1) {
      prob <- FS(y, R0, mu)
    } else if (n > 1) {
      for (i in 1:(y - n + 1)) {
        prob <- prob + FS(i, R0, mu) * FSn(n - 1, y - i, R0, mu)
      }
    }
  }
  return(prob)
}


pmf.FSn <- function(n, R0, mu, maxFS) {
  df_n1 <- pmf.FS(R0, mu, maxFS)
  prob <- rep(0, maxFS)
  if (n == 1) {
    df <- df_n1
  } else if (n > 1) {
    df_n_one_less <- pmf.FSn(n - 1, R0, mu, maxFS)
    for (i in (1:(maxFS - n))) {
      for (j in ((n - 1):(maxFS - i))) {
        prob[i + j] <- prob[i + j] + df_n1[i, 2] * df_n_one_less[j, 2]
      }
    }
    names(prob) <- seq(1, maxFS)
    df <- data.frame(FinalSize = names(prob), prob)
  }
  pi <- max((1 - P_wt_est(n, R0, mu)), sum(df$prob))
  df$prob_normal <- df$prob / pi
  df$FinalSize <- as.numeric(df$FinalSize)
  return(df)
}


# 5. Mutant offspring generated pmf ====
P_mu_n <- function(n, m, R0, mu, maxFS) {
  prob <- 0
  df <- pmf.FSn(n, R0, mu, maxFS)
  for (i in 1:maxFS) {
    k <- (1 - exp(-mu)) * i
    p <- 1 / (R0 + 1)
    NBmu <- dnbinom(m, k, p)
    prob <- prob + df[i, 3] * NBmu
  }
  return(prob)
}


pmf.P_mu_n <- function(n, R0, mu, maxMuGen, maxFS) {
  prob <- rep(0, maxMuGen + 1)
  df_FS <- pmf.FSn(n, R0, mu, maxFS)
  for (i in 0:maxMuGen) {
    prob_mugen <- 0
    for (j in 1:maxFS) {
      k <- (1 - exp(-mu)) * j
      p <- 1 / (R0 + 1)
      NBmu <- dnbinom(i, k, p)
      prob_mugen <- prob_mugen + df_FS[j, 3] * NBmu
    }
    prob[i + 1] <- prob[i + 1] + prob_mugen
  }
  names(prob) <- c(0:maxMuGen)
  df <- data.frame(MuLinGen = names(prob), prob)
  df$MuLinGen <- as.numeric(df$MuLinGen)
  return(df)
}


# 6. Mutant lineages established pmf ====
P_mu_est <- function(n, m, R0, mu, maxMuGen, maxFS) {
  p_est <- 1 - 1 / R0
  prob <- rep(0, 1)
  pmf_mu_gen <- pmf.P_mu_n(n, R0, mu, maxMuGen, maxFS)
  for (i in m:maxMuGen) {
    p_m <- pmf_mu_gen[i + 1, 2]
    prob <- prob + dbinom(m, i, p_est) * p_m
  }
  return(prob)
}


pmf.P_mu_est <- function(n, R0, mu, maxMuGen, maxFS) {
  p_est <- 1 - 1 / R0
  df_mu_gen <- pmf.P_mu_n(n, R0, mu, maxMuGen, maxFS)
  prob <- rep(0, maxMuGen + 1)
  for (j in 0:maxMuGen) {
    p_m <- df_mu_gen[j + 1, 2]
    for (i in 0:j) {
      prob[i + 1] <- prob[i + 1] + dbinom(i, j, p_est) * p_m
    }
  }
  names(prob) <- c(0:maxMuGen)
  df <- data.frame(MuLinEst = names(prob), prob)
  df$MuLinEst <- as.numeric(df$MuLinEst)
  return(df)
}


# 7. Combine mutant lineages and clonal mutations ====
pmf.mu_est_clonal <- function(n, R0, mu, maxMuGen, maxFS) {
  df <- pmf.P_mu_est(n, R0, mu, maxMuGen, maxFS)
  df$prob_norm <- df$prob * (1 - P_wt_est(n, R0, mu))
  df$Clonal <- "P0"
  df$Clonal[df$MuLinEst == "0"] <- "Px"
  df$Clonal[df$MuLinEst == "1"] <- "P1+"
  df[nrow(df) + 1, ] <- c("S_inf", 0, P_wt_est(n, R0, mu), "P0")
  df$prob <- as.numeric(df$prob)
  df$prob_norm <- as.numeric(df$prob_norm)
  df <- df[order(df$prob_norm, decreasing = TRUE), ]
  df$Clonal <- factor(df$Clonal, order = TRUE, levels = c("Px", "P0", "P1+"))
  return(df)
}


# 8. Clonal mutation pmf ====
Px_analytical <- function(n, R0) {
  px <- 1 / R0^n
  return(px)
}


Px <- function(n, R0, mu, maxMuGen, maxFS) {
  p_wt_ext <- 1 - P_wt_est(n, R0, mu)
  p_mu_est_0 <- P_mu_est(n, 0, R0, mu, maxMuGen, maxFS)
  prob <- p_mu_est_0 * p_wt_ext
  return(prob)
}


P0 <- function(n, R0, mu, maxMuGen, maxFS) {
  S_inf <- P_wt_est(n, R0, mu)
  S_2plus <- (1 - P_mu_est(n, 0, R0, mu, maxMuGen, maxFS) - P_mu_est(n, 1, R0, mu, maxMuGen, maxFS)) * (1 - S_inf)
  prob <- S_inf + S_2plus
  return(prob)
}


P1plus <- function(n, R0, mu, maxMuGen, maxFS) {
  p_wt_ext <- 1 - P_wt_est(n, R0, mu)
  prob <- P_mu_est(n, 1, R0, mu, maxMuGen, maxFS) * p_wt_ext
  return(prob)
}


P_clonal_round <- function(n, m, R0, mu, maxMuGen, maxFS) {
  prob <- P1plus(n, R0, mu, maxMuGen, maxFS) * dpois(m, mu) / (1 - dpois(0, mu))
  return(prob)
}


pmf.P_clonal_round <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  p1plus <- P1plus(n, R0, mu, maxMuGen, maxFS)
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- P0(n, R0, mu, maxMuGen, maxFS)
  if (maxClonal != 0) {
    for (i in 1:maxClonal) {
      prob[i + 2] <- p1plus * dpois(i, mu) / (1 - dpois(0, mu))
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  return(df)
}


P_clonal <- function(n, m, R0, mu, maxMuGen, maxFS) {
  df_clonal_round <- pmf.P_clonal_round(n, R0, mu, maxMuGen, maxFS, m)[-1, ]
  px <- 1 / R0
  prob <- 0
  if (m == 0) {
    prob <- df_clonal_round[1, 2]
  } else {
    for (i in 1:m) {
      prob <- prob + df_clonal_round[i + 1, 2] * P_clonal(1, m - i, R0, mu, maxMuGen, maxFS) / (1 - px)
    }
  }
  return(prob)
}


pmf.P_clonal_withPx <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df_clonal_round <- pmf.P_clonal_round(n,R0,mu,maxMuGen,maxFS,maxClonal)[-1, ]
  px <- 1 / R0
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- df_clonal_round[1, 2]
  if (maxClonal != 0) {
    df_clonal_max_one_less <- pmf.P_clonal_withPx(1, R0, mu, maxMuGen, maxFS, maxClonal - 1)
    for (m in 1:maxClonal) {
      for (i in 1:m) {
        prob[m+2] <- prob[m+2] + df_clonal_round[i+1,2] * df_clonal_max_one_less[m-i+2,2] / (1-px)
      }
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  return(df)
}


pmf.P_clonal <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df_clonal_round <- pmf.P_clonal_round(n, R0, mu, maxMuGen, maxFS, maxClonal)[-1, ]
  px <- 1 / R0
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- df_clonal_round[1, 2]
  if (maxClonal != 0) {
    df_clonal_max_one_less <- pmf.P_clonal_withPx(1, R0, mu, maxMuGen, maxFS, maxClonal - 1)
    for (m in 1:maxClonal) {
      for (i in 1:m) {
        prob[m+2] <- prob[m+2] + df_clonal_round[i+1,2] * df_clonal_max_one_less[m-i+2,2] / (1-px)
      }
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  fold <- max(1 - df$prob[1], sum(df[-1,]$prob))
  df$prob <- df$prob / fold
  df <- df[-1, ]
  return(df)
}


list_clonal <- function(n_values, R0, mu_values, maxMuGen, maxFS, clonal) {
  listClonal <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(n_values) %>%
        lapply(function(y) {
          pmf.P_clonal(n_values[y], R0, mu_values[x], maxMuGen, maxFS, clonal)
        })
    }, mc.cores = detectCores())
  names(listClonal) <- mu_values
  for (i in 1:length(mu_values)) {
    names(listClonal[[i]]) <- n_values
  }
  return(listClonal)
}


list_clonal_mu <- function(n_values, R0_values, mu, maxMuGen, maxFS, clonal) {
  listClonal <- 1:length(R0_values) %>%
    mclapply(function(x) {
      1:length(n_values) %>%
        lapply(function(y) {
          pmf.P_clonal(n_values[y], R0_values[x], mu, maxMuGen, maxFS, clonal)
        })
    }, mc.cores = detectCores())
  names(listClonal) <- R0_values
  for (i in 1:length(R0_values)) {
    names(listClonal[[i]]) <- n_values
  }
  return(listClonal)
}


# 9. Log likelihood ====
pmf.clonal_params <- function(clonal, listClonal, n_values, R0, mu_values, maxMuGen, maxFS) {
  df <- expand.grid(n = paste0("", n_values), mu = paste0("", mu_values))
  df$prob <- 0
  
  for (N in n_values) {
    for (MU in mu_values) {
      muloc <- which(names(listClonal) == MU)
      nloc <- which(names(listClonal[[muloc]]) == N)
      prob <- listClonal[[muloc]][[nloc]]$prob[clonal+1]
      df$prob[which(df$n == N & df$mu == MU)] <- prob
    }
  }
  
  df$n <- as.numeric(df$n)
  df$mu <- as.numeric(paste(df$mu))
  return(df)
}


pmf.clonal_params_mu <- function(clonal, listClonal, n_values, R0_values, mu, maxMuGen, maxFS) {
  df <- expand.grid(n = paste0("", n_values), R0 = paste0("", R0_values))
  df$prob <- 0
  
  for (N in n_values) {
    for (Rnaught in R0_values) {
      R0loc <- which(names(listClonal) == Rnaught)
      nloc <- which(names(listClonal[[R0loc]]) == N)
      prob <- listClonal[[R0loc]][[nloc]]$prob[clonal+1]
      df$prob[which(df$n == N & df$R0 == Rnaught)] <- prob
    }
  }
  
  df$n <- as.numeric(df$n)
  df$R0 <- as.numeric(paste(df$R0))
  return(df)
}


pmf.ini_size <- function(lambda, R0, maxIni) {
  df <- data.frame(N = 0:maxIni, pois = dpois(0:maxIni, lambda))
  df$adj_prob <- 0
  for (n in 0:maxIni) {
    Pest <- 1 - 1 / R0^n
    df$adj_prob[n + 1] <- df$pois[n + 1] * Pest
  }
  df$adj_prob <- df$adj_prob / sum(df$adj_prob)
  df <- melt(df, id.vars = "N")
  names(df) <- c("N", "type", "prob")
  df$type <- factor(df$type)
  levels(df$type) <- c("poisson", "adjusted")
  return(df)
}


pmf.ini_size_nb <- function(lambda, R0, maxIni, k) {
  df <- data.frame(N = 0:maxIni, nbin = dnbinom(0:maxIni, size = k, mu = lambda))
  df$adj_prob <- 0
  for (n in 0:maxIni) {
    Pest <- 1-1/R0^n
    df$adj_prob[n+1] <- df$nbin[n+1] * Pest
  }
  df$adj_prob <- df$adj_prob/sum(df$adj_prob)
  df <- melt(df, id.vars = "N")
  names(df) <- c("N", "type", "prob")
  df$type <- factor(df$type)
  levels(df$type) <- c("nbinom", "adjusted")
  return(df)
}


LL_lambda <- function(df, listClonal, lambda, R0, mu, maxMuGen, maxFS, maxIni) {
  df$prob <- df$freq / sum(df$freq)
  dfPn <- pmf.ini_size(lambda, R0, maxIni)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  prob <- 1
  
  list <- mclapply(df$ClonalMu, 
                   pmf.clonal_params, 
                   listClonal = listClonal, 
                   n_values = 1:maxIni, 
                   R0 = R0, 
                   mu_values = mu, 
                   maxMuGen = maxMuGen, 
                   maxFS = maxFS,
                   mc.cores = detectCores())
  names(list) <- df$ClonalMu
  
  for (i in df$ClonalMu) {
    loc <- which(names(list) == i)
    dfClonalN <- list[[loc]]
    PclonalLam <- 0
    for (n in 1:maxIni) {
      Pn <- dfPn$prob[dfPn$N == n]
      PclonalN <- dfClonalN$prob[dfClonalN$n == n]
      PclonalLam <- PclonalLam + Pn*PclonalN
    }
    k <- df$freq[df$ClonalMu == i]
    prob <- prob * PclonalLam^k
  }
  prob <- log(prob)
  return(prob)
}


LL_lambda_nb <- function(df, listClonal, lambda, R0, mu, maxMuGen, maxFS, maxIni, k) {
  df$prob <- df$freq / sum(df$freq)
  dfPn <- pmf.ini_size_nb(lambda, R0, maxIni, k)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  prob <- 1
  
  list <- mclapply(df$ClonalMu, 
                   pmf.clonal_params, 
                   listClonal = listClonal, 
                   n_values = 1:maxIni, 
                   R0 = R0, 
                   mu_values = mu, 
                   maxMuGen = maxMuGen, 
                   maxFS = maxFS,
                   mc.cores = detectCores())
  names(list) <- df$ClonalMu
  
  for (i in df$ClonalMu) {
    loc <- which(names(list) == i)
    dfClonalN <- list[[loc]]
    PclonalLam <- 0
    for (n in 1:maxIni) {
      Pn <- dfPn$prob[dfPn$N == n]
      PclonalN <- dfClonalN$prob[dfClonalN$n == n]
      PclonalLam <- PclonalLam + Pn*PclonalN
    }
    m <- df$freq[df$ClonalMu == i]
    prob <- prob * PclonalLam^m
  }
  prob <- log(prob)
  return(prob)
}


LL_lambda_mu <- function(df, listClonal, lambda, R0, mu, maxMuGen, maxFS, maxIni) {
  df$prob <- df$freq / sum(df$freq)
  dfPn <- pmf.ini_size(lambda, R0, maxIni)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  prob <- 1
  
  list <- mclapply(df$ClonalMu, 
                   pmf.clonal_params_mu, 
                   listClonal = listClonal, 
                   n_values = 1:maxIni, 
                   R0_values = R0, 
                   mu = mu, 
                   maxMuGen = maxMuGen, 
                   maxFS = maxFS,
                   mc.cores = detectCores())
  names(list) <- df$ClonalMu
  
  for (i in df$ClonalMu) {
    loc <- which(names(list) == i)
    dfClonalN <- list[[loc]]
    PclonalLam <- 0
    for (n in 1:maxIni) {
      Pn <- dfPn$prob[dfPn$N == n]
      PclonalN <- dfClonalN$prob[dfClonalN$n == n]
      PclonalLam <- PclonalLam + Pn*PclonalN
    }
    k <- df$freq[df$ClonalMu == i]
    prob <- prob * PclonalLam^k
  }
  prob <- log(prob)
  return(prob)
}


LL_meanNb.df <- function(df, listClonal, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni) {
  prob <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(lambda_values) %>%
        lapply(function(y) {
          LL_lambda(df, listClonal, lambda_values[y], 
                    R0, mu_values[x], maxMuGen, maxFS, maxIni)
        })
    }, mc.cores = detectCores())
  names(prob) <- mu_values
  for (i in 1:length(mu_values)) {
    names(prob[[i]]) <- lambda_values
  }
  
  LLdf <- expand.grid(lambda = paste0("", lambda_values), mu = paste0("", mu_values))
  LLdf$prob <- unlist(prob)
  LLdf$lambda <- as.numeric(paste(LLdf$lambda))
  LLdf$mu <- as.numeric(paste(LLdf$mu))
  LLdf$meanN <- 0
  LLdf$meanNb <- 0
  
  for (j in lambda_values) {
    # calculate mean N
    dfPn <- pmf.ini_size(j, R0, maxIni)
    dfPn <- dfPn[dfPn$type == "adjusted",]
    dfPn$prob_N <- dfPn$N * dfPn$prob
    mean_n <- sum(dfPn$prob_N)
    LLdf[which(near(LLdf$lambda, j)),]$meanN <- mean_n
    
    # calculate mean Nb
    mean_nb <- 0
    for (k in 1:maxIni) {
      PNb <- 0
      for (n in k:maxIni) {
        PNb <- PNb + dfPn$prob[n+1] * dbinom(k, n, 1-1/R0) / (1-dbinom(0, n, 1-1/R0))
      }
      mean_nb <- mean_nb + k * PNb
    }
    LLdf[which(near(LLdf$lambda, j)),]$meanNb <- mean_nb
  }
  
  return(LLdf)
}


LL_meanNb_nb.df <- function(df, listClonal, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni, k) {
  prob <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(lambda_values) %>%
        lapply(function(y) {
          LL_lambda_nb(df, listClonal, lambda_values[y],
                       R0, mu_values[x], maxMuGen, maxFS, maxIni, k)
        })
    }, mc.cores = detectCores())
  names(prob) <- mu_values
  for (i in 1:length(mu_values)) {
    names(prob[[i]]) <- lambda_values
  }
  
  LLdf <- expand.grid(lambda = paste0("", lambda_values), mu = paste0("", mu_values))
  LLdf$prob <- unlist(prob)
  LLdf$lambda <- as.numeric(paste(LLdf$lambda))
  LLdf$mu <- as.numeric(paste(LLdf$mu))
  LLdf$meanN <- 0
  LLdf$meanNb <- 0
  
  for (j in lambda_values) {
    # calculate mean N
    dfPn <- pmf.ini_size_nb(j, R0, maxIni, k)
    dfPn <- dfPn[dfPn$type == "adjusted",]
    dfPn$prob_N <- dfPn$N * dfPn$prob
    mean_n <- sum(dfPn$prob_N)
    LLdf[which(near(LLdf$lambda, j)),]$meanN <- mean_n
    
    # calculate mean Nb
    mean_nb <- 0
    for (m in 1:maxIni) {
      PNb <- 0
      for (n in m:maxIni) {
        PNb <- PNb + dfPn$prob[n+1] * dbinom(m, n, 1-1/R0) / (1-dbinom(0, n, 1-1/R0))
      }
      mean_nb <- mean_nb + m * PNb
    }
    LLdf[which(near(LLdf$lambda, j)),]$meanNb <- mean_nb
  }
  
  return(LLdf)
}


LL_meanNb_mu.df <- function(df, listClonal, lambda_values, R0_values, mu, maxMuGen, maxFS, maxIni) {
  prob <- 1:length(R0_values) %>%
    mclapply(function(x) {
      1:length(lambda_values) %>%
        lapply(function(y) {
          LL_lambda_mu(df, listClonal, lambda_values[y], 
                       R0_values[x], mu, maxMuGen, maxFS, maxIni)
        })
    }, mc.cores = detectCores())
  names(prob) <- R0_values
  for (i in 1:length(R0_values)) {
    names(prob[[i]]) <- lambda_values
  }
  
  LLdf <- expand.grid(lambda = paste0("", lambda_values), R0 = paste0("", R0_values))
  LLdf$prob <- unlist(prob)
  LLdf$lambda <- as.numeric(paste(LLdf$lambda))
  LLdf$R0 <- as.numeric(paste(LLdf$R0))
  LLdf$meanN <- 0
  LLdf$meanNb <- 0
  
  for (j in 1:nrow(LLdf)) {
    # calculate mean N
    dfPn <- pmf.ini_size(LLdf$lambda[j], LLdf$R0[j], maxIni)
    dfPn <- dfPn[dfPn$type == "adjusted",]
    dfPn$prob_N <- dfPn$N * dfPn$prob
    mean_n <- sum(dfPn$prob_N)
    LLdf$meanN[j] <- mean_n
    
    # calculate mean Nb
    mean_nb <- 0
    for (k in 1:maxIni) {
      PNb <- 0
      for (n in k:maxIni) {
        PNb <- PNb + dfPn$prob[n+1] * dbinom(k, n, 1-1/R0) / (1-dbinom(0, n, 1-1/R0))
      }
      mean_nb <- mean_nb + k * PNb
    }
    LLdf$meanNb[j] <- mean_nb
  }
  
  return(LLdf)
}


LL_prop <- function(df, listClonal, prop, R0, mu, maxMuGen, maxFS) {
  prob <- 1
  list <- mclapply(df$ClonalMu, 
                   pmf.clonal_params, 
                   listClonal = listClonal, 
                   n_values = 1, 
                   R0 = R0, 
                   mu_values = mu, 
                   maxMuGen = maxMuGen, 
                   maxFS = maxFS,
                   mc.cores = detectCores())
  names(list) <- df$ClonalMu
  
  for (i in df$ClonalMu) {
    loc <- which(names(list) == i)
    dfClonalN <- list[[loc]]
    PclonalLam <- prop*dfClonalN$prob
    k <- df$freq[df$ClonalMu == i]
    prob <- prob * PclonalLam^k
  }
  prob <- log(prob)
  return(prob)
}


LL_prop.df <- function(df, listClonal, prop_values, R0,
                       mu_values, maxMuGen, maxFS, maxIni) {
  prob <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(prop_values) %>%
        lapply(function(y) {
          LL_prop(df, listClonal, prop_values[y],R0, mu_values[x], maxMuGen, maxFS)
        })
    }, mc.cores = detectCores())
  names(prob) <- mu_values
  for (i in 1:length(mu_values)) {
    names(prob[[i]]) <- prop_values
  }
  
  LLdf <- expand.grid(prop = paste0("", prop_values), mu = paste0("", mu_values))
  LLdf$prob <- unlist(prob)
  LLdf$prop <- as.numeric(paste(LLdf$prop))
  LLdf$mu <- as.numeric(paste(LLdf$mu))
  LLdf$meanN <- LLdf$prop
  return(LLdf)
}


calc_distri <- function(sample, lambda, R0, mu, maxMuGen, maxFS, maxClonal, maxIni) {
  dfPn <- pmf.ini_size(lambda, R0, maxIni)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  dfcalc <- data.frame(ClonalMu = 0:maxClonal, prob = 0)
  for (n in 1:maxIni) {
    dfClonal <- pmf.P_clonal(n,R0,mu,maxMuGen,maxFS,maxClonal)
    for (clonal in 0:maxClonal) {
      dfcalc$prob[clonal+1] <- dfcalc$prob[clonal+1] + dfPn$prob[n+1] * dfClonal$prob[clonal+1]
    }
  }
  dfcalc$prob <- dfcalc$prob / sum(dfcalc$prob)
  
  sample$prob_meta <- sample$freq / sum(sample$freq)
  comb_distri <- full_join(sample[, c(1,3)], dfcalc, by = "ClonalMu")
  melt_distri <- melt(comb_distri, id = "ClonalMu")
  names(melt_distri) <- c("ClonalMu", "type", "prob")
  melt_distri$type <- factor(melt_distri$type)
  levels(melt_distri$type) <- c("simulated", "calculated")
  
  return(melt_distri)
}


calc_distri_nb <- function(sample, lambda, R0, mu, maxMuGen, maxFS, maxClonal, maxIni, k) {
  dfPn <- pmf.ini_size_nb(lambda, R0, maxIni, k)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  dfcalc <- data.frame(ClonalMu = 0:maxClonal, prob = 0)
  for (n in 1:maxIni) {
    dfClonal <- pmf.P_clonal(n,R0,mu,maxMuGen,maxFS,maxClonal)
    for (clonal in 0:maxClonal) {
      dfcalc$prob[clonal+1] <- dfcalc$prob[clonal+1] + dfPn$prob[n+1] * dfClonal$prob[clonal+1]
    }
  }
  dfcalc$prob <- dfcalc$prob / sum(dfcalc$prob)
  
  sample$prob_meta <- sample$freq / sum(sample$freq)
  comb_distri <- full_join(sample[, c(1,3)], dfcalc, by = "ClonalMu")
  melt_distri <- melt(comb_distri, id = "ClonalMu")
  names(melt_distri) <- c("ClonalMu", "type", "prob")
  melt_distri$type <- factor(melt_distri$type)
  levels(melt_distri$type) <- c("simulated", "calculated")
  
  return(melt_distri)
}


calc_distri_lambda0 <- function(sample, R0, mu, maxMuGen, maxFS, maxClonal) {
  dfcalc <- data.frame(ClonalMu = 0:maxClonal, prob = 0)
  dfClonal <- pmf.P_clonal(1,R0,mu,maxMuGen,maxFS,maxClonal)
  for (clonal in 0:maxClonal) {
    dfcalc$prob[clonal+1] <- dfcalc$prob[clonal+1] + dfClonal$prob[clonal+1]
  }
  dfcalc$prob <- dfcalc$prob / sum(dfcalc$prob)
  
  sample$prob_meta <- sample$freq / sum(sample$freq)
  comb_distri <- full_join(sample[, c(1,3)], dfcalc, by = "ClonalMu")
  melt_distri <- melt(comb_distri, id = "ClonalMu")
  names(melt_distri) <- c("ClonalMu", "type", "prob")
  melt_distri$type <- factor(melt_distri$type)
  levels(melt_distri$type) <- c("simulated", "calculated")
  
  return(melt_distri)
}


# 10. Expected clonal mutations ====
meanClonalBozic <- function(R0, mu) {
  delta <- 1/R0
  mc <- (delta*mu)/(1-delta)
  return(mc)
}

meanClonalBozic.df <- function(R0_values, mu_values) {
  mean <- 1:length(mu_values) %>%
    mclapply(function(x){
      1:length(R0_values) %>%
        lapply(function(y){
          meanClonalBozic(R0_values[y], mu_values[x])
        })
    }, mc.cores = detectCores()) 
  
  mcdf <- expand.grid(R0 = paste0("", R0_values), mu = paste0("", mu_values))
  vs_slope <- c()
  
  mcdf$mc <- unlist(mean)
  mcdf$logmc <- log10(mcdf$mc)
  return(mcdf)
}

meanClonal <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df <- pmf.P_clonal(n, R0, mu, maxMuGen, maxFS, maxClonal)
  sum <- 0
  for (i in 1:maxClonal) {
    sum <- sum + i * df[i + 1, 2]
  }
  return(sum)
}

meanClonal.df <- function(n, R0_values, mu_values, maxMuGen, maxFS, maxClonal) {
  mean <- 1:length(R0_values) %>%
    mclapply(function(x){
      1:length(mu_values) %>%
        lapply(function(y){
          meanClonal(n,R0_values[x], mu_values[y], maxMuGen, maxFS, maxClonal)
        })
    }, mc.cores = detectCores()) 
  
  mcdf <- expand.grid(mu = paste0("", mu_values), R0 = paste0("", R0_values))
  vs_slope = c()
  
  mcdf$mc <- unlist(mean)
  mcdf <- mcdf[, c("R0", "mu", "mc")]
  mcdf$logmc <- log10(mcdf$mc)
  return(mcdf)
}

# 11. Source of mutation ====
P_transAndFix <- function(q, R0, maxIni) {
  freq <- q/100
  df <- data.frame(N = 1:maxIni, prob = 0)
  for (N in 1:maxIni) {
    P_trans_fix <- 0
    for (k in 1:N) {
      P_trans_fix <- P_trans_fix + dbinom(k, N, freq) * (1/R0)^(N-k) * (1-dbinom(0, k, 1-1/R0))
    }
    P_trans_fix <- P_trans_fix/(1-Px_analytical(N, R0))
    df$prob[N] <- P_trans_fix
  }
  return(df)
}


P_denovoAndFix <- function(R0, mu, maxMuGen, maxFS, maxIni, maxClonal, sites) {
  df <- data.frame(N = 1:maxIni, prob = 0)
  for (N in 1:maxIni) {
    df_clonal <- pmf.P_clonal(N, R0, mu, maxMuGen, maxFS, maxClonal)
    sum <- 0
    for (k in 0:maxClonal) {
      p <- df_clonal$prob[k + 1]
      sum <- sum + p*k/sites
    }
    df$prob[N] <- sum
  }
  return(df)
}