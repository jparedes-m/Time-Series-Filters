# One Sided Hodrick Prescott Filter

# Written by Jorge Paredes on Oct, 2022

# Based upon the description given by Wolf et. al (2020) 
## https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3536248

# Preamble -----
invisible(lapply(c("readr", "tidyverse", "lubridate", "fredr"), library, character.only =T))

home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
laptop <- "C:/Users/jpare/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)
wd <- getwd()

rm(home, laptop)

# FRED Setting
fredr_set_key("dedfd4de61dbed2e5bb2e5e13e1d972a")

# Data ----

gdp_us <- fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

hp_filter1 <- function(data, lambda){
  y <- select_if(data, is.numeric)[[1]]
  N <- length(y)
  L = lambda
  date <- select_if(data, is.Date)[[1]]
  tau <- numeric(N)
  
  hp_fill <- hp_filter(data = data, lambda = L)
  hp_fill <- hp_fill$data$trend[1:4]
  tau[1:4] <- hp_fill
  
  for(t in 5:N){
    data_s <- data %>% filter(row_number()<=t)
    
    hp2s<- hp_filter(data = data_s, lambda = L)
    
    trend <- tail(hp2s$data$trend,1)
    tau[[t]] <- trend
  }
  
  c = y - tau
  
  results <- data.frame(date = seq(head(date,1), length.out = N, by = "quarters"),
                        serie = y,
                        trend = tau,
                        cycle = c)
  
  return(results)
  
}

us_hp <- hp_filter1(data = gdp_us , lambda = 1600)

# Graphics ----
par(mfrow=c(1,2))
plot(us_hp$date, exp(us_hp$serie), main = "Serie and trend \n Hodrick - Prescott [1s] (λ = 1600)",
     ylab="Billions of US dollars (2012, chained.)", 
     xlab="Time", type="l", lwd = 1.5)
lines(us_hp$date,exp(us_hp$trend), col="red", lwd = 1.5)
legend("topleft", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1)

plot(us_hp$date, 100*us_hp$cycle, main ="Cycle \n Hodrick - Prescott [1s] (λ = 1600)",
     ylab="% Deviation from Trend", xlab = "Time", type="l", lwd = 1.5, yaxp = c(-10, 5, 15))
abline(h=0, col="red", lty=2, lwd = 2)
par(mfrow=c(1,1))

#EOF