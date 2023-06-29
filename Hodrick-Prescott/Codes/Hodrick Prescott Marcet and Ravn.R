# Hodrick - Prescott decomposition for univariate time series
# Taken from The HP-filter in cross-country comparisons, written by Albert Marcet & Morten O Ravn
# Paper available here: https://repositori.upf.edu/handle/10230/1203?locale-attribute=es

# Code written by: Jorge Paredes on June, 2022.

# Marcet and Ravn HP Filter ----
hp_filter_MR <- function(data1, data2, lambda = 1600, rule = "rule 1"){
  ## Call my function for the HP Filter (I do not know if this should go inside the code)
  source("https://raw.githubusercontent.com/jparedes-m/Time-Series-Filters/main/Hodrick-Prescott/Codes/Hodrick%20Prescott%20Filter.R")
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

#EOF
