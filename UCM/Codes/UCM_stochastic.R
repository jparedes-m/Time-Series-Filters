# Unobserved Components Model - Stochastic Drift

# Written by Jorge Paredes on Oct, 2022.

# Preamble ----
library(FKF)

# Function ----
ucm_stoch <- function(data, p_estimates = FALSE){
  yt <- select_if(data, is.numeric)[[1]]
  date <- select_if(data, is.Date)[[1]]
  n <- length(yt)
  
  
  # Get the guesses for p1, p2, sigma_e, sigma_w, sigma_u, ----
  
  source("https://raw.githubusercontent.com/jparedes-m/Time-Series-Filters/main/Hodrick-Prescott/Codes/Hodrick%20Prescott%20Filter.R")
  hp_est <- hp_filter(data = data, lambda = 1600)
  yt_cycle_hp <- as.matrix(hp_est$data$cycle)
  yt_trend_hp <- as.matrix(hp_est$data$trend)
  
  # Cycle equation
  ce_mle <- arima(ts(yt_cycle_hp), order = c(2, 0, 0), method = "CSS", include.mean = F)
  phi0 <- as.matrix(ce_mle$coef[1:2])
  sigma_u_hp <- var(ce_mle$residuals)
  #cat("\n Guess created by using Hodrick Prescott & MLE estimation of UCM \n")
  
  # Trend equation
  tr_ols <- arima(ts(yt_trend_hp), order = c(1,0,0), method = "CSS",
                  include.mean = FALSE )
  res_trend <- as.matrix(tr_ols$residuals)
  
  drift0 <- stats::filter(res_trend, rep(1 / 8, 8), sides = 1)
  drift0[1 : 8] <- res_trend[1 : 8]
  e0 <- (res_trend - drift0)[9 : (length(yt)-2)]
  w0 <- (drift0 - stats::lag(drift0))[9 : (length(yt)-2)]
  
  # Construct a0 = (tau[2], drift[2], cycle[2], cicle[1]) ----
  a0 <- as.numeric(c(yt_trend_hp[2], drift0[2], yt_cycle_hp[2], yt_cycle_hp[1]))
  
  # Get P0 ----
  trendVar <- diag(c(var(yt_trend_hp), var(drift0)))
  cycleCov <- cov(cbind(lag(yt_cycle_hp)[3 : n], lag(yt_cycle_hp, n = 2)[3 : n]))
  zeros <- matrix(0, ncol = 2, nrow = 2)
  
  P0 <- rbind(cbind(trendVar, zeros), cbind(zeros, cycleCov))
  
  # Variances of innovations
  sigma_e_hp <- var(e0)
  sigma_w_hp <- var(w0)
  
  # Setup initial guess
  par0 <- as.matrix(c(phi0, 
                      sigma_e_hp,
                      sigma_w_hp,
                      sigma_u_hp), ncol = 1, nrow = 5, byrow = TRUE)
  
  g_hp <- c(phi = phi0, 
            sigma_e = sigma_e_hp,
            sigma_w = sigma_w_hp,
            sigma_u = sigma_u_hp)
  
  # Get the constant matrices in the SS form
  ct <- matrix(0)
  dt <- matrix(0, nrow = 4, ncol = 1)
  GGt <- matrix(0)
  Zt <- matrix(c(1, 0, 1, 0), nrow = 1, ncol = 4, byrow =TRUE)
  
  
  # Set Up sample ----
  yt_s <- t(yt[3:n])
  
  
  # Kalman Function ----
  Kalman <- function(par, yt, a0, P0, ct, Zt, GGt, dt){
    
    # Get the Tt and HHt matrices
    
    Tt <- matrix(c(1,1,0,0,
                   0,1,0,0,
                   0,0,par[1], par[2],
                   0,0,1,0), nrow = 4, ncol = 4, byrow = TRUE)
    
    HHt <- diag(c(par[3], par[4], par[5], 0))
    
    results <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, 
                   Tt = Tt, Zt = Zt, HHt = HHt, 
                   GGt = GGt, yt=yt)
    
    return(results)
  }
  
  # Set up optimizing function ----
  fn <- function(x0, ...){
    -Kalman(x0, yt_s,a0, P0, ct, Zt, GGt, dt)$logLik
  }
  
  # Run global optimizer
  fit.fkf <- optim(par0, fn, method = "SANN")
  
  # Run local optimizer over global min
  fit.fkf <- optim(fit.fkf$par, fn)
  fit.obj <- Kalman(fit.fkf$par, yt = yt_s, a0 = a0, 
                    P0 = P0, ct = ct, 
                    Zt = Zt, GGt = GGt, dt = dt)
  
  # Print estimates ----
  
  kalman_est <- c(phi = c(fit.fkf$par[1],fit.fkf$par[2]),
                  sigma_e = fit.fkf$par[3],
                  sigma_w = fit.fkf$par[4],
                  sigma_u = fit.fkf$par[5])
  
  if(p_estimates == TRUE){
    print(cbind(g_hp, kalman_est))
  }
  
  # Smoother ----
  smooth <- fks(fit.obj)
  smooth <- as.data.frame(t(smooth$ahatt)) %>% 
    rename(trend = V1,
           drift = V2,
           cycle = V3) %>% 
    select(-V4)
  
  row_12 <- data.frame(trend = c(NA, NA), drift = c(NA, NA), cycle = c(NA, NA))
  
  smooth <- rbind(row_12, smooth)
  
  df_s <- data.frame(date = seq(head(date, 1), length = n, by = "quarters"),
                     serie = yt,
                     trend = smooth$trend,
                     drift = smooth$drift,
                     cycle = smooth$cycle)
  
  df_s <- df_s[-(1:2),]
  
  rownames(df_s) <- NULL
  
  # Data frame with trend and cycle estimates ----
  df_ <- as.data.frame(t(fit.obj$at)) %>% 
    rename(trend = V1, drift = V2, 
           cycle = V3, L1_cycle = V4) %>% 
    select(trend, drift, cycle)
  
  row1 <- c(trend = NA, drift = NA, cycle = NA)
  
  df_ <- rbind(row1, df_)
  
  df <- data.frame(date = seq(head(date,1), length = n, by = "quarters"),
                   serie = yt, 
                   trend = df_$trend,
                   drift = df_$drift,
                   cycle = df_$cycle)
  
  df <- df[-1,]
  rownames(df) <- NULL
  
  result <- list(estimates = kalman_est, 
                 data_filter = df,
                 smoother = df_s)
  
  return(result)
  
}
