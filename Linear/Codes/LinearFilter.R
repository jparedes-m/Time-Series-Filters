# Linear Filter

linear <- function(data, p_estimates = FALSE, drift = TRUE){
  yt <- select_if(data, is.numeric)
  date <- select_if(data, is.Date)
  n <- length(yt)
  
  t <- 1:n
  
  # Estimate regression ----
  if (drift == TRUE) {
    reg <- lm(yt ~ t)
  }else{
    reg <- lm(yt ~ 0+t)
  }
  
  # Print coefficients ----
  if(p_estimates == TRUE){
    summary(reg)
  }
  
  # Get the fitted values (trend) and the cycle ----
  y_hat <-  fitted(reg); c = yt - y_hat
  
  # Exportation ----
  df <- data.frame(date = date, 
                   serie = yt,
                   trend = y_hat,
                   cycle = c)
  
  return(df)
}
