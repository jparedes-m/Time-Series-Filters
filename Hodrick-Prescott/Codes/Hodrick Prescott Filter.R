# Hodrick - Prescott decomposition for univariate time series

# Written by: Jorge Paredes on May, 2022.

# Inspired by the first problem set of intermediate macroeconomics lectured at USFQ by Carlos Uribe

# Preamble ----
invisible(lapply(c("readr", "dplyr", "tidyverse", "lubridate", "data.table"), library, character.only=T))

home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)

# Data ----
# I used the Real Gross Domestic Product for the U.S. Billions of Chained 2012 Dollars, 
## Seasonally Adjusted Annual Rate
# Data is retrieved from FRED St. Louis: https://fred.stlouisfed.org/series/GDPC1 

gdp <- read_csv("GDPC1.csv") %>% 
  rename(gdp=GDPC1, date=DATE) %>% 
  mutate(date=as.Date(date),
         gdp=log(gdp))

# HP Filter ----
hp_filter <- function(data, lambda, show_F=FALSE, additional_info=FALSE){
  # Arrange the data ----
  data <- as.data.frame(data)
  date <- select_if(data, is.Date)
  y <- select_if(data, is.numeric)
  y <- y[[1]]
  # Compute the F matrix ----
  n <- nrow(data)
  row_1 <- c(1,-2,1, rep(0, each=(n-3)))
  row_2 <- c(-2,5,-4,1, rep(0, each=(n-4)))
  intermediate_rows <- list()
  for(i in 3:(n-2)){
    before <- abs(3-i)
    after <- n-length(c(rep(0, each=before), c(1,-4,6,-4,1)))
    row <- c(rep(0, each=before), c(1,-4,6,-4,1), rep(0, each=after))
    intermediate_rows[[i]] <- t(row)
  }
  intermediate_rows <- intermediate_rows[!sapply(intermediate_rows, is.null)]
  intermediate_rows <- as.matrix(do.call(rbind, intermediate_rows))
  row_n_1 <- c(rep(0, each=(n-4)), 1,-4, 5,-2)
  row_n <- c(rep(0, each=(n-3)), 1,-2,1)
  
  F_ <- rbind(row_1, row_2, intermediate_rows, row_n_1, row_n)
  rownames(F_) <- NULL
  # Compute the lambda * F matrix ----
  LF_=lambda*F_
  # Compute the I matrix ----
  I <- diag(n)
  # Compute the (I+LF)^-1 matrix ---
  A <- (I+LF_)
  A <- solve(A)
  # Compute the (I+LF)^-1 * y  (this is the trend component)----
  trend <- A %*% y
  # Compute the cycle -----
  cycle <-  y - trend
  # Pack everything into a data.frame ----
  variable <- "variable"
  
  trend_name <- paste("Trend",variable)
  
  cycle_name <- paste("Cycle",variable)
  
  result <- data.frame(date = date,
                       variable =y,
                       trend = trend, 
                       cycle = cycle)
  # Additional info ----
  if(additional_info==T){
  cat(paste(" Hodrick - Prescott Filter \n -------------------------\n", 
            "λ =", lambda, "\n",
            "Series Type: y = trend + cycle \n\n"))
  }
  
  if(show_F==TRUE){
    print(F_)
  }
  
  return(result)
  
}

gdp7 <- filter(gdp, row_number()<=7)

hp_filter(data= gdp7, lambda = 1600, show_F = TRUE, additional_info = T)

HP <- hp_filter(data=gdp, lambda = 1600)

# EOF
