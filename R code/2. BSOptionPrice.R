# 2a) 
BSOptionPrice <- function(S, K, r, T_t, sigma, isput) {
  
  d1 <- (log(S / K) + (r + (sigma^2) / 2) * T_t) / (sigma * sqrt(T_t))
  d2 <- d1 - sigma * sqrt(T_t)
  
  N_d1 <- pnorm(d1)
  N_d2 <- pnorm(d2)
  
  # according to the type of option, calculate the price using theoretical formula
  if (isput) {
    price <- K * exp(-r * T_t) * pnorm(-d2) - S * pnorm(-d1)  # if put option
  } else {
    price <- S * N_d1 - K * exp(-r * T_t) * N_d2  # if call option
  }
  return(price)
}


# 2b)

call_price <- BSOptionPrice(S=100, K=105, r=0.02, T_t=0.5, sigma=0.2, isput=FALSE)
cat("The call price is :", call_price, "\n")

put_price <- BSOptionPrice(S=100, K=105, r=0.02, T_t=0.5, sigma=0.2, isput=TRUE)
cat("The put price is :", put_price, "\n")


# 2c) 
BSImplicitVol <- function(OptionPrice, S, K, r, T_t, isput) {
  result_delta <- function(sigma) {
    BSOptionPrice(S, K, r, T_t, sigma, isput) - OptionPrice
  }

  # uniroot allows us to find the implicit volatility
  result <- uniroot(result_delta, lower = 1e-6, upper = 5)
  return(result$root) 
}


# 2d) 
vol_implicite_call <- BSImplicitVol(OptionPrice = 2.7852, S = 100, K = 105, r = 0.02, T_t = 0.25, isput = FALSE)
cat("The implicit volatility for a call option is :", vol_implicite_call, "\n")


# 2e) 
implicit_vol_put <- BSImplicitVol(OptionPrice = 6.8249, S = 100, K = 105, r = 0.02, T_t = 0.75, isput = TRUE)
cat("The implicit volatility for a put option is :", implicit_vol_put, "\n")
