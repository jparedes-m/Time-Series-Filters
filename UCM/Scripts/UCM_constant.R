# Unobserved Components Model - Constant Drift

# Written by Jorge Paredes on Oct, 2022.

# Preamble ----

# install.packages(c("tidyverse","foreach", "doParallel","lubridate","zoo", "FKF", "mFilter", 
#"latex2exp", "readxl", "progress", "doSNOW", "parallel"))
invisible(lapply(c("tidyverse","foreach", "doParallel","lubridate","zoo", "FKF", "mFilter", 
                   "latex2exp", "readxl", "progress", "doSNOW", "parallel", "fredr"), library, character.only = T))


# Working Directory
home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
laptop <- "C:/Users/jpare/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)
wd <- getwd()

# FRED Setting
fredr_set_key("#####")

rm(laptop, home)

# Data ----
gdp_us <- fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

# Function ----

ucm_const <- function(data, p_estimates = FALSE){
  
  # Set the data ----
  yt <- select_if(data, is.numeric)[[1]]
  date <- select_if(data, is.Date)[[1]]
  n <- length(yt)
  
  # Estimate the useful parameters  ----
  hp_est <- hpfilter(yt, freq = 1600,type = "lambda")
  
  trend_hp <- matrix(hp_est$trend, ncol = 1)
  cycle_hp <- matrix(hp_est$cycle, ncol = 1)
  
  ## Estimate parameters for the trend equation
  D_trend_hp <-  trend_hp - lag(trend_hp)
  delta_hp <-  mean(D_trend_hp, na.rm = TRUE)
  sigma_e_hp <-  var((D_trend_hp - delta_hp), na.rm = TRUE)
  
  ## Estimate parameters for the cycle equation
  cycle_eq <- arima(ts(cycle_hp), order = c(2,0,0), method = "CSS", include.mean = FALSE)
  phi0 <- as.matrix(cycle_eq$coef[1:2])
  sigma_w_hp <- var(cycle_eq$residuals)
  
  # Construct a0 and P0 ----
  a0 <- as.numeric(c(trend_hp[2], cycle_hp[2], cycle_hp[1]))
  
  gamma0 <- var(trend_hp)
  vcv_c <- cov(cbind(lag(cycle_hp)[3 : n], lag(cycle_hp, n = 2)[3 : n]))
  zeros <- matrix(0, nrow = 1, ncol = 2)
  
  P0 <- rbind(cbind(gamma0, zeros), cbind(t(zeros), vcv_c))
  rm(gamma0, vcv_c, zeros)
  # Set up initial guess ----
  par0 <- as.matrix(c(phi0, delta_hp, sigma_e_hp, sigma_w_hp), nrow = 1, ncol = 5, byrow = TRUE)
  g_hp <- c(phi = phi0, delta = delta_hp, sigma_e = sigma_e_hp, sigma_w = sigma_w_hp)
  
  # State space constant matrices ----
  ct <- matrix(0)
  GGt <- matrix(0)
  Zt <- matrix(c(1,1,0), nrow = 1, ncol = 3)
  
  # Set up sample ----
  yt_s <- t(yt[3:n])
  
  # Optimal parameters ----
  Kalman <- function(par, yt, a0, P0, ct, Zt, GGt, dt){
    # dt vector
    dt <- matrix(c(par[3],0,0), nrow = 3, ncol =1)
    
    # Tt matrix
    Tt <- matrix(c(1,0,0,
                   0,par[1], par[2],
                   0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
    
    # HT matrix 
    HHt <- diag(c(par[4], par[5], 0))
    
    results <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, 
                   Tt = Tt, Zt = Zt, HHt = HHt, 
                   GGt = GGt, yt=yt)
    return(results)
  }
  # Optimizing function ----
  fn <- function(x0, ...){
    -Kalman(x0, yt_s,a0, P0, ct, Zt, GGt, dt)$logLik
  }
  # Run optimizers -----
  ## Global optimizer
  fit.fkf <- optim(par0, fn, method = "SANN")
  ## Local optimizer
  fit.fkf <- optim(fit.fkf$par, fn)
  ## Fitted object
  fit.obj <- Kalman(fit.fkf$par, yt = yt_s, a0 = a0, 
                    P0 = P0, ct = ct, 
                    Zt = Zt, GGt = GGt, dt = dt)
  # Estimates ----
  kalman_est <- c(phi = c(fit.fkf$par[1],fit.fkf$par[2]),
                  delta = fit.fkf$par[3],
                  sigma_e = fit.fkf$par[4],
                  sigma_w = fit.fkf$par[5])
  
  if(p_estimates == TRUE){
    print(cbind(g_hp, kalman_est))
  }
  
  # Kalman Smoother ----
  smooth <- fks(fit.obj)
  smooth <- as.data.frame(t(smooth$ahatt)) %>% 
    rename(trend = V1,
           cycle = V2) %>% 
    select(-V3) %>% 
    mutate(drift = kalman_est["delta"])
  row_12 <- data.frame(trend = c(NA, NA), drift = c(NA, NA), cycle = c(NA, NA))
  
  smooth <- rbind(row_12, smooth)
  
  df_s <- data.frame(date = seq(head(date, 1), length = n, by = "quarters"),
                     serie = yt,
                     trend = smooth$trend,
                     drift = smooth$drift,
                     cycle = smooth$cycle)
  
  df_s <- df_s[-(1:2),]
  rownames(df_s) <- NULL
  
  
  # Filter ----
  df_ <- as.data.frame(t(fit.obj$at)) %>% 
    rename(trend = V1, 
           cycle = V2, 
           L1_cycle = V3) %>% 
    select(trend, cycle) %>% 
    mutate(drift = kalman_est["delta"])
  
  row1 <- c(trend = NA, drift = NA, cycle = NA)
  
  df_ <- rbind(row1, df_)
  
  df <- data.frame(date = seq(head(date,1), length = n, by = "quarters"),
                   serie = yt, 
                   trend = df_$trend,
                   drift = df_$drift,
                   cycle = df_$cycle)
  
  df <- df[-1,]
  rownames(df) <- NULL
  # Exportation ----
  result <- list(estimates = kalman_est, 
                 data_filter = df,
                 smoother = df_s)
  
  return(result)
}
us_ucm <- ucm_const(data = gdp_us, p_estimates =TRUE)

# Graph ----
par(mfrow = c(1,2))
plot(us_ucm$data_filter$date, exp(us_ucm$data_filter$serie), type = "l", lwd = 1.75,
     main = "Unobserved Components Model \n Constant Drift \n United States",
     ylab = "Billions of US Dollars (2012, chained)",
     xlab = "Time")
lines(us_ucm$data_filter$date, exp(us_ucm$data_filter$trend), col = "red", lwd = 1.75)

plot(us_ucm$data_filter$date, 100*us_ucm$data_filter$cycle, type = "l", lwd = 1.75,
     main = "Unobserved Components Model \n Constant Drift \n United States",
     ylab = "% Deviation from trend",
     xlab = "Time")
abline(h = 0, col = "red", lwd = 2, lty = 2)
par(mfrow  = c(1,1))
#EOF
