# For this example I will assume you have all your packages loaded I will use the FRED data on the US Real GDP.
# Furthermore, I assume you have the functions `ucm_const` and `ucm_stoch` in your global environment.

library(fredr)

## Put your api key here
fredr_set_key("####")

# Data ----
US_GDP <- fredr::fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2023-06-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

# Unobserved Components Model: Constant Drift ----
us_ucm <- ucm_const(data = US_GDP, p_estimates =TRUE)

par(mfrow = c(1,2))
plot(us_ucm$data_filter$date, exp(us_ucm$data_filter$serie), type = "l", lwd = 1.5,
     main = "Unobserved Components Model \n Constant Drift \n United States",
     ylab = "Billions of US Dollars (2012, chained)",
     xlab = "Time")
lines(us_ucm$data_filter$date, exp(us_ucm$data_filter$trend), col = "red", lwd = 1.5)

plot(us_ucm$data_filter$date, 100*us_ucm$data_filter$cycle, type = "l", lwd = 1.5,
     main = "Unobserved Components Model \n Constant Drift \n United States",
     ylab = "% Deviation from trend",
     xlab = "Time",
     yaxp = c(min(round(100*us_ucm$data_filter$cycle,0)), max(round(100*us_ucm$data_filter$cycle, 0)), 10),
     las =2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
par(mfrow  = c(1,1))

# Unobserved Components Model: Stochastic Drift ----
us_ucm <- ucm_stoch(data = US_GDP, p_estimates = TRUE)

par(mfrow = c(1,2))
plot(us_ucm$data_filter$date, exp(us_ucm$data_filter$serie), type = "l", lwd = 1.5,
     main = "Unobserved Components Model \n Stochastic Drift \n United States",
     ylab = "Billions of US Dollars (2012, chained)",
     xlab = "Time")
lines(us_ucm$data_filter$date, exp(us_ucm$data_filter$trend), col = "red", lwd = 1.5)

plot(us_ucm$data_filter$date, 100*us_ucm$data_filter$cycle, type = "l", lwd = 1.5,
     main = "Unobserved Components Model \n Stochastic Drift \n United States",
     ylab = "% Deviation from trend",
     xlab = "Time",
     yaxp = c(min(round(100*us_ucm$data_filter$cycle,0)), max(round(100*us_ucm$data_filter$cycle, 0)), 10),
     las =2)
abline(h = 0, col = "red", lwd = 2, lty = 2)
par(mfrow  = c(1,1))

# EOF
