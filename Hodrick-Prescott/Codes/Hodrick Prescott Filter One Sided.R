# One Sided Hodrick Prescott Filter

# Written by Jorge Paredes on Oct, 2022

# Based upon the description given by Wolf et. al (2020) 
## https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3536248

hp_filter1s <- function(data, lambda){
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

#EOF
