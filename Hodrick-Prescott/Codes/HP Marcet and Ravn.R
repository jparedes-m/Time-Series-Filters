# Hodrick - Prescott decomposition for univariate time series
# The Marcet and Ravn modification

# Written by: Jorge Paredes on June, 2022.

# Preamble ----
invisible(lapply(c("readr", "tidyverse", "dplyr", "lubridate", "readxl"), library, character.only =T))

home <- "C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Time Series/Univariate/Hodrick Prescott/Data"
setwd(home)

# Data ----

gdp <- read_csv("gdp_us.csv") %>% 
  mutate(date=as.Date(date),
         gdp=log(gdp)) %>% 
  filter(year(date)>=2000)

gdp_ec <- read_csv("gdp_ec.csv") %>% 
  mutate(date=as.Date(date), 
         gdp=log(gdp/1000000))

gdp_uk <- read_csv("gdp_uk.csv") %>% 
  mutate(gdp = log(gdp)) %>% 
  filter(year(date)>=2000)


# Marcet and Ravn ----
hp_filter_MR <- function(data1, data2, lambda, rule = "rule 1", start=1, end, tol=1e-10){
  
  if(rule == "rule 1" | rule == 1){
    # Rule 1 ----
    
    ## Find V 
    find_v <- function(data, L){
      
      hp1 <- hp_filter(data, lambda = L, show_F = F, additional_info = F)
      
      y <- hp1$variable
      trend <- hp1$trend
      n <- nrow(hp1)
      
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
    V <- find_v(data1, L = lambda)
    ## Find the lambda for data2 
    G <- function(L){
      F_ <- find_v(data2, L)
      g=F_-V
      return(g)
    }
    G_root <- uniroot(G, interval = c(start, end), tol = tol)
    
    cat(paste(" Root finding: F(λ) - V = 0\n ----------------- \n V:", 
              V,"\n λ for data2 under rule 1 is:", round(G_root$root,4),
              "\n (F-V) =", G_root$f.root, "\n")) 
    opt_lambda <- G_root$root
    ## Compute HP for data2 with optimal lambda 
    HP <- hp_filter(data = data2, lambda = opt_lambda, show_F = F, additional_info = F)
    
  }else{
    if(rule== "rule 2" | rule == 2){
      # Rule 2 ----
      
      ## Find W
      find_w <- function(data, L){
        
        hp1 <- hp_filter(data, lambda = L, show_F = F, additional_info = F)
        
        y <- hp1$variable
        trend <- hp1$trend
        n <- nrow(hp1)
        
        numerator <- c()
        for(t in 2:(n-1)){
          val <- (trend[t+1] - 2*trend[t] + trend[t-1])^2
          numerator[t] <- val
        }
        numerator <- sum(numerator, na.rm = T)
        W <- numerator /(n-2)
      }
      W <- find_w(data1, L=lambda)
      
      ## Find lambda for data2
      G <- function(L){
        F_ <- find_w(data2, L)
        g=F_ - W
        return(g)
      }
      G_root <- uniroot(G, interval = c(start, end), tol = tol, extendInt = "yes")
      
      cat(paste(" Root finding: F(λ) - W = 0\n ----------------- \n W:", 
                W,"\n λ for data2 under rule 2 is:", round(G_root$root,4),
                "\n (F-V) =", G_root$f.root, "\n"))
      
      opt_lambda <- G_root$root
      ## Compute HP for data2 with optimal lambda 
      HP <- hp_filter(data = data2, lambda = opt_lambda, show_F = F, additional_info = F)
      
    }
  }
  
  # Return ----
  return(HP)
}

# Examples and timing ----
start <- Sys.time()
HP_ec <- hp_filter_MR(data1 = gdp, data2 = gdp_ec, lambda = 1600, rule = "rule 1", start = 1200, end=3800)
end <- Sys.time()
end-start

start <- Sys.time()
HP_uk <- hp_filter_MR(data1 = gdp, data2 = gdp_uk, lambda = 1600, rule = "rule 2", end=3800)
end <- Sys.time()
end-start

# Graphs ----

par(mfrow=c(1,2))
plot(HP_ec$date, HP_ec$variable, main = "Variable - Trend",
     ylab="Ln(GDP)", xlab="Time", col="black", type="l")
lines(HP_ec$date,HP_ec$trend, col="red", lty=2)
legend("bottomright", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1:2)

plot(HP_ec$date, HP_ec$cycle, main ="Cycle",
     ylab="% Deviations from Trend", xlab = "Time", col="blue", lty=1, type="l")
abline(h=0, col="black", lty=2)
par(mfrow=c(1,1))


par(mfrow=c(1,2))
plot(HP_uk$date, HP_uk$variable, main = "Variable - Trend",
     ylab="Ln(GDP)", xlab="Time", col="black", type="l")
lines(HP_uk$date,HP_uk$trend, col="red", lty=2)
legend("bottomright", legend = c("Variable", "Trend"), col=c("black", "red"), lty=1:2)

plot(HP_uk$date, HP_uk$cycle, main ="Cycle",
     ylab="% Deviations from Trend", xlab = "Time", col="blue", lty=1, type="l")
abline(h=0, col="black", lty=2)
par(mfrow=c(1,1))

#EOF

HP_ec %>%
  filter(date >= "2012-12-31" & date <= "2021-10-01") %>%
  mutate(quarter = as.yearqtr(date)) %>% 
  ggplot() +
  aes(x = quarter, y = cycle) +
  geom_line(size = 0.75, colour = "#B22222") +
  geom_hline(yintercept = 0,size = 0.75, linetype="dashed")+
  labs(x = "Quarter", y = "% Deviation from Trend", 
       title = "Ecuador's Business Cycle", 
       subtitle = "Obtained by Hodrick Prescott Filter with Marcet and Ravn Modification | λ = 3276.746", 
       caption = "Data obtained by National Quarterly Accounts | BCE") +
  theme_light() +
  theme(plot.title = element_text(size = 20L, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 15L, hjust = 0.5))+
  scale_x_yearqtr(format = '%Yq%q',n=30) 
