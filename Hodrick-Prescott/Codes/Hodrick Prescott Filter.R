# Hodrick - Prescott decomposition for univariate time series

# Written by: Jorge Paredes on May, 2022.

# Inspired by the first problem set of intermediate macroeconomics lectured at USFQ by Carlos Uribe

# Preamble ----
invisible(lapply(c("readr","tidyverse", "lubridate", "data.table"), library, character.only=T))

home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
laptop <- "C:/Users/jpare/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)
wd <- getwd()

rm(home, laptop)

# Data ----
fredr_set_key("#######")

gdp_us <- fredr(series_id = "GDPC1", 
                observation_start = as.Date("1970-01-01"),
                observation_end = as.Date("2022-04-01")) %>% 
  select(date, value) %>% 
  rename(gdp = value) %>% 
  mutate(gdp = log(gdp))
# HP Filter ----
hp_filter <- function(data, lambda){
  # Set the data ----
  y <- select_if(data, is.numeric)[[1]]
  date <- select_if(data, is.Date)[[1]]
  N <- length(y); L <- lambda
  
  # Build the F matrix ----
  Row1 <- matrix(c(1,-2,1, rep(0, times = (N-3))), nrow =1)
  Row2 <- matrix(c(-2, 5, -4, 1, rep(0, times = (N-4))), nrow = 1)
  
  mid_mat <- function(dim){
    vec <- c(1,-4,6,-4,1, rep(0, times = (N-5)))
    mat <- toeplitz(vec)
    mat[lower.tri(mat)] <- 0
    mat[seq(dim - 4), ]
  }
  Mid_mat <- mid_mat(N)
  
  Row_n_2 <- matrix(c(rep(0, times = (N-4)), 1,-4, 5,-2), nrow = 1)
  Row_n_1 <- matrix(c(rep(0, times = (N-3)), 1,-2,1), nrow = 1)
  
  F_mat <- rbind(Row1, Row2, Mid_mat, Row_n_2, Row_n_1)
  
  # Create the matrix for filtering (I+LF)^-1 ----
  H <- solve(diag(1, nrow = N, ncol = N)+L*F_mat)
  
  # Export the results as a data frame ----
  tau <- H %*% y; c = y-tau
  
  results <- data.frame(date = seq(head(date,1), length.out = N, by = "quarters"),
                        serie = y,
                        trend = tau,
                        cycle = c)
  
  export <- list(data = results, H_matrix = H, lambda = L)
  
  return(export)
}

us_hp <- hp_filter(data=gdp, lambda = 1600)


# Graphs ----
par(mfrow=c(1,2))
plot(us_hp$data$date, exp(us_hp$data$serie), main = "Serie and trend \n Hodrick - Prescott (λ = 1600)",
     ylab="Billions of US dollars (2012, chained.)", 
     xlab="Time", type="l", lwd = 1.5)
lines(us_hp$data$date,exp(us_hp$data$trend), col="red", lwd = 1.5)
legend("topleft", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1)

plot(us_hp$data$date, 100*us_hp$data$cycle, main ="Cycle \n Hodrick - Prescott (λ = 1600)",
     ylab="% Deviation from Trend", xlab = "Time", type="l", lwd = 1.5, yaxp = c(-10, 5, 15))
abline(h=0, col="red", lty=2, lwd = 2)
par(mfrow=c(1,1))


# EOF
