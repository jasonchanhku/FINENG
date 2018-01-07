########################1. loading the necessary packages######################## 
library(quantmod)
library(Quandl)
library(xts)
library(plotly)
library(GenSA)

#########################2. getting underlying data##############################
toDate <- function(x) as.Date(x, origin = "1899-12-30")
USDMYR <- read.zoo("data/USDMYR.csv", header = TRUE, sep = ",", FUN = toDate) #read .csv file of underlying data
USDMYR <- as.xts(USDMYR) #convert to xts format 
USDMYR.r <- na.omit(log(USDMYR / lag(USDMYR))) #returns of underlying and omit 'NA'
USDMYR.vol <- sd(USDMYR) #volatility of underlying

#########################3. getting crude oil data###############################
Brent <- read.zoo("data/Brent.csv", header = TRUE, sep = ",", FUN = toDate) #read .csv file of oil price
Brent <- as.xts(Brent) #convert to xts format


###########4. plotting graphs of underyling and oil price to see correlation#################
#subsetting underlying data for plotting
USDMYRdate <- index(USDMYR)
USDMYRval <- data.frame(USDMYR)
USDMYRval <- USDMYRval[ , 1]

#subsetting oil price data for plotting
Brentdate <- index(Brent)
Brentval <- data.frame(Brent)
Brentval <- Brentval[ , 1]

#below are the necessary commands to plot the graph on plotly
ay <- list( tickfont = list(color = "red"), overlaying = "y", side = "right" )
#run this command manually as qqplot overwrites the plotly graph
print(plot_ly(x = USDMYRdate, y = USDMYRval, name = "USDMYR") %>% add_trace(x = Brentdate, y = Brentval, name = "Brent", yaxis = "y2") %>% layout(title = "USDMYR & Brent", yaxis2 = ay))

#########################5. Kappa and Xi Estimation###############################
#Subsetting the qq data of historical data
temp <- qqnorm(USDMYR.r)
temp1 <- as.data.frame(temp$y)
qqUSDMYR <- temp1[ , 1]

r <- 0.02                 #risk-free rate
n <- 497                  #number of points
theta <- USDMYR.vol^2      # Long-run mean.
nu <- s <- rep(0, n + 1) # Initialize vol and stock price.
s.r <- rep(0, n)         #simulated returns vector
difference <- rep(0, n)   #vector to store difference
nu[1] <- theta
s[1] <- USDMYR[length(USDMYR)]

Gen <- function(q){  #function to be optimized
  
  for(i in 1:n){
    dW1 <- sqrt(1/252)  * rnorm(1)
    dW2 <- sqrt(1/252)  * rnorm(1)
    ds <- r*s[i]*(1/252) + sqrt(nu[i])*s[i]*dW1
    dnu <- q[1]*(theta-nu[i])*(1/252) + q[2]*sqrt(nu[i])*dW2
    s[i + 1] <- s[i] + ds
    nu[i + 1] <- max(nu[i] + dnu, 0) # Ensure non-negative 'nu'.
    s.r[i] <- (s[i+1] - s[i]) / s[i]
    difference[i] <- abs(qqUSDMYR[i] - s.r[i])  #calculating the difference of each component
  }
  return(difference) #return the difference numeric vector
}

Parest <- GenSA(fn = Gen, lower = c(0,0), upper = c(10,10))   #optimize function
kappa <- Parest$par[1]   #store value of kappa
xi <- Parest$par[2]      #store value of xi

############################6. Options Pricing#################################
d <- 100000                                  # Number of trials.
T <- 1                                    # 1 year
m <- T*252                                #time steps
I <- 10000                                # principle
K <- USDMYR[length(USDMYR)]               #Strike Price
dt <- 1/m                                 #change in time
r <- 0.02                                 #risk free rate
theta <- USDMYR.vol^2                     # Long-run mean.
nu <- st <- rep(0, m + 1)                 # Initialize vol and stock price.
nu[1] <- theta
st[1] <- USDMYR[length(USDMYR)]           # Initial value of simulated price
f <- rep(0, d)                            # vector for put prices.
g <- rep(0, d)                            # vector for Lookback prices.


for (j in 1:d) {                 # 'j' cycles through trials.
  for (i in 1:m) {               # 'i' cycles through time.
    dW1 <- sqrt(dt) * rnorm(1)
    dW2 <- sqrt(dt) * rnorm(1)
    ds <- r*st[i]*dt + sqrt(nu[i])*st[i]*dW1
    dnu <- kappa*(theta-nu[i])*dt + xi*sqrt(nu[i])*dW2
    st[i + 1] <- st[i] + ds
    nu[i + 1] <- max(nu[i] + dnu, 0) # Ensure non-negative 'nu'.
  }
  f[j] <- exp(-r * T) * max(K - st[m/2], 0)
  g[j] <- exp(-r * T) * max(st[m/2] - min(st[m/2:m]))
}
put <- mean(f)                           #put option price 
putpart <- I*(1-exp(-r*T))/put           #put participation
putpayoff <- I + putpart*put*exp(r * T)*exp(r * T/2)   #put payoff

LB <- mean(g)                            #lookback option price
LBpart <- putpart*mean(f)*exp(r * T)/LB  #lookback option participation
LBpayoff <- I + LBpart*LB*exp(r * T)     #lookback option payoff

cat("Put Option Price Estimate:", round(put, 4), "\n")
cat("Put Option Participation:", round(putpart, 4), "\n")
cat("Put Option Participation Rate:", round(putpart/I, 4), "\n")
cat("Put Payoff:", round(putpayoff, 4), "\n")
cat("Standard Error:", round(sd(f) / sqrt(d), 4), "\n")
cat("###################################################", "\n")
cat("Lookback Option Price Estimate:", round(LB, 4), "\n")
cat("Lookback Option Participation:", round(LBpart, 4), "\n")
cat("Lookback Option Participation Rate:", round(LBpart/I, 4), "\n")
cat("Lookback Option Payoff:", round(LBpayoff, 4), "\n")
cat("Standard Error:", round(sd(g) / sqrt(d), 4), "\n")

# ############################7.Plotting Options Histogram#################################
# plot histogram of put option price. Run this manually.
print(plot_ly(x = f, type = "histogram")%>% add_trace(x = c(put, put), y = c(0, d/1.5)))

# plot histogram of lookback option price. Run this manually.
print(plot_ly(x = g, type = "histogram")%>% add_trace(x = c(LB, LB), y = c(0, d/1.5)))

