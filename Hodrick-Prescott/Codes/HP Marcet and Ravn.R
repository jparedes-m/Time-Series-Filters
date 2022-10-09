# Hodrick - Prescott decomposition for univariate time series
# The Marcet and Ravn modification

# Written by: Jorge Paredes on June, 2022.

# Preamble ----
invisible(lapply(c("readr", "tidyverse", "lubridate", "fredr"), library, character.only =T))

home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
laptop <- "C:/Users/jpare/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)
wd <- getwd()

rm(home, laptop)

# FRED Setting
fredr_set_key("#######")

# Data ----

gdp_us <- fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

gdp_ec <- read_delim("gdp_ec.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(date=as.Date(date), 
         gdp=log(gdp))


gdp_uk <- fredr(series_id = "CLVMNACSCAB1GQUK", 
                observation_start = as.Date("1975-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

# Marcet and Ravn ----
hp_filter_MR <- function(data1, data2, lambda = 1600, rule = "rule 1"){
  
  if(rule == "rule 1" | rule == 1){
    # Rule 1 ----
    ## Find V 
    find_v <- function(data, L){
      
      hp1 <- hp_filter(data, lambda = L)
      
      y <- hp1$data$serie
      trend <- hp1$data$trend
      n <- length(y)
      
      numerator <- c()
      for(t in 2:(n-1)){
        val <- (trend[t+1] - 2*trend[t] + trend[t-1])^2
        numerator[t] <- val
      }
      numerator <- sum(numerator, na.rm = T)
      
      denominator <- c()
      for(t in 1:n){
        val <- (y[t] - trend[t])^2
        denominator[t] <- val
      }
      denominator <- sum(denominator)
      
      V <- numerator/denominator
      
      return(V)
      
    }
    V1 <- find_v(data1, L = 1600)
    ## Find the lambda for data2 
    G <- function(L, V){
      F_ <- find_v(data2, L)
      g=F_-V
      return(g)
    }
    
    G_root <- uniroot(G, tol = 1e-64, interval = c(0, 1e6), V = V1, extendInt = "yes")
    
    cat(paste(" Root finding: F(λ) - V = 0\n ----------------- \n V:", 
              V1,"\n λ data2 under rule 1 is:", round(G_root$root,4),
              "\n (F-V) =", G_root$f.root, "\n")) 
    opt_lambda <- G_root$root
    ## Compute HP for data2 with optimal lambda 
    HP <- hp_filter(data = data2, lambda = opt_lambda)
    
  }else if(rule== "rule 2" | rule == 2){
      # Rule 2 ----
      ## Find W
      find_w <- function(data, L){
        
        hp1 <- hp_filter(data, lambda = L)
        
        y <- hp1$data$serie
        trend <- hp1$data$trend
        n <- length(y)
        
        numerator <- c()
        for(t in 2:(n-1)){
          val <- (trend[t+1] - 2*trend[t] + trend[t-1])^2
          numerator[t] <- val
        }
        numerator <- sum(numerator, na.rm = T)
        W <- numerator /(n-2)
      }
      W1 <- find_w(data1, L=lambda)
      
      ## Find lambda for data2
      G <- function(L, W){
        F_ <- find_w(data2, L)
        g=F_ - W
        return(g)
      }
      G_root <- uniroot(G, tol = 1e-64, interval = c(0, 1e6), W = W1, extendInt = "yes")
      
      cat(paste(" Root finding: F(λ) - W = 0\n ----------------- \n W:", 
                W1,"\n λ for data2 under rule 2 is:", round(G_root$root,4),
                "\n (F-V) =", G_root$f.root, "\n"))
      
      opt_lambda <- G_root$root
      ## Compute HP for data2 with optimal lambda 
      HP <- hp_filter(data = data2, lambda = opt_lambda)
      
    }
  
  # Return ----
  return(HP)
}

# Examples and timing ----
start <- Sys.time()
HP_ec <- hp_filter_MR(data1 = gdp_us, data2 = gdp_ec, lambda = 1600, rule = "rule 1")
end <- Sys.time()
end-start

start <- Sys.time()
HP_uk <- hp_filter_MR(data1 = gdp_us, data2 = gdp_uk, lambda = 1600, rule = "rule 2")
end <- Sys.time()
end-start

rm(end, start, gdp_ec, gdp_uk, gdp_us)
# Graphs ----

par(mfrow=c(1,2))
plot(HP_ec$data$date, exp(HP_ec$data$serie), main = paste0("Serie and trend \n Hodrick - Prescott (λ = ",round(HP_ec$lambda,0),") \n Ecuador"),
     ylab="Billions of US dollars (2007, chained.)", 
     xlab="Time", type="l", lwd = 1.5)
lines(HP_ec$data$date,exp(HP_ec$data$trend), col="red", lwd = 1.5)
legend("topleft", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1)

plot(HP_ec$data$date, 100*HP_ec$data$cycle, main = paste0("Cycle \n Hodrick - Prescott (λ = ",round(HP_ec$lambda,0),") \n Ecuador"),
     ylab="% Deviation from Trend", xlab = "Time", type="l", lwd = 1.5, yaxp = c(-15, 5, 20))
abline(h=0, col="red", lty=2, lwd = 2)
par(mfrow=c(1,1))


par(mfrow=c(1,2))
plot(HP_uk$data$date, exp(HP_uk$data$serie), main = paste0("Serie and trend \n Hodrick - Prescott (λ = ",round(HP_uk$lambda,0),") \n United Kingdom"),
     ylab="Billions of US dollars (2012, chained.)", 
     xlab="Time", type="l", lwd = 1.5)
lines(HP_uk$data$date,exp(HP_uk$data$trend), col="red", lwd = 1.5)
legend("topleft", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1)

plot(HP_uk$data$date, 100*HP_uk$data$cycle, main = paste0("Cycle \n Hodrick - Prescott (λ = ",round(HP_uk$lambda,0),") \n United Kingdom"),
     ylab="% Deviation from Trend", xlab = "Time", type="l", lwd = 1.5, yaxp = c(-30, 5, 35))
abline(h=0, col="red", lty=2, lwd = 2)
par(mfrow=c(1,1))

#EOF
