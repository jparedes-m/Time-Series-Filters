# How to use the linear filter function:

## Data ----
# You could just download data from FRED and get it into R using readr or readxl packages but I just wanna show off.
library(fredr); library(dplyr)
## instead of the hashtags put your fredr key
fredr_set_key("####")
gdp_us <- fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))

## Filtering ----
Linear <- linear(data = gdp_us, p_estimates = TRUE, drift = TRUE)
# Graphs ----
par(mfrow=c(1,2))
plot(Linear$date, exp(Linear$serie), main = "Serie and trend \n Linear filter with drift.",
     ylab="Billions of US dollars (2012, chained.)", 
     xlab="Time", type="l", lwd = 1.5)
lines(Linear$date,exp(Linear$trend), col="red", lwd = 1.5)
legend("topleft", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1)

plot(Linear$date, 100*Linear$cycle, main ="Cycle \n Linear filter with drift.",
     ylab="% Deviation from Trend", xlab = "Time", type="l", lwd = 1.5)
abline(h=0, col="red", lty=2, lwd = 2)
par(mfrow=c(1,1))
