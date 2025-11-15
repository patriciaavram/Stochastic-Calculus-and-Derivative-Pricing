# From previous exercise, functions needed for this to run: 
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



# 3a) 
BinOptionPrice <- function(S, K, r, T_t, mu, sigma, n, isput) {
  # parameters of binomial tree
  h <- T_t / n  # h=1/n is annualized and we have T_t years 
  u <- exp(mu *h + sigma * sqrt(h))
  d <- exp(mu *h - sigma * sqrt(h))
  q <- (exp(r * h) - d) / (u - d)
  
  # initialization of terminal prices 
  terminal_values <- numeric(n + 1)  # empty vector of size n+1 to stock terminal values
  for (i in 0:n) {
    St <- S * u^i * d^(n - i)  # calculate underlying prices for a given trajectory at the last level of tree
    if (isput) {
      terminal_values[i + 1] <- max(K - St, 0) # payoff for put option
    } else {
      terminal_values[i + 1] <- max(St - K, 0) # payoff for call option
    }
  }
  
  # backwards propagation of tree values
  for (j in (n - 1):0) {
    current_values <- numeric(j + 1)  # temporary vector to store values at current level
    for (i in 0:j) {
      current_values[i + 1] <- exp(-r * h) * (q * terminal_values[i + 2] + (1 - q) * terminal_values[i + 1])
    }
    terminal_values <- current_values  # updating values for next level
  }
  
  return(terminal_values[1])  # returns price at time 0
}


# 3b)

call_price_binom <- BinOptionPrice(S = 100, K = 105, r = 0.02, T_t = 0.5, mu = 0.02 - (0.2^2)/2, sigma = 0.20, n = 20, isput = FALSE)
cat("Price of the call option under binomial tree is :", call_price_binom, "\n")

put_price_binom <- BinOptionPrice(S = 100, K = 105, r = 0.02, T_t = 0.5, mu = 0.02 - (0.2^2)/2, sigma = 0.20, n = 20, isput = TRUE)
cat("Price of the put option under binomial tree is :", put_price_binom, "\n")


# 3c) 
r <- 0.02
sigma <- 0.20

mu_crr <- 0  # for Cox-Ross-Rubinstein
mu_lognormal <- r - (sigma^2)/2  # for lognormal

# We first calculate the theoretical price using Black-Scholes
bs_call <- BSOptionPrice(S = 100, K = 105, r, T_t = 0.5, sigma, isput = FALSE)
bs_put <- BSOptionPrice(S = 100, K = 105, r, T_t = 0.5, sigma, isput = TRUE)

cat("Black-Scholes price (Call) :", bs_call, "\n")
cat("Black-Scholes price (Put) :", bs_put, "\n")

# We test convergence for a different number of steps n, trying with different mu values
mu_lognormal_2 <- 0.1
mu_lognormal_3 <- 0.5

n_values <- c(5, 10, 20, 50, 100, 200)
results <- data.frame(n = n_values, CRR_Call = NA, L_Call_mu_1 = NA, L_Call_mu_2 = NA, L_Call_mu_3 = NA, CRR_Put = NA, L_Put_mu_1 = NA,  L_Put_mu_2 = NA,  L_Put_mu_3 = NA) 

# Using BinOptionPrice function to calculate option prices for each n step value and mu (CRR et lognormal).
for (index in 1:length(n_values)) {
  n_value <- n_values[index]
  
  # Call
  results$CRR_Call[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_crr, n_value, isput = FALSE)
  results$L_Call_mu_1[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal, n_value, isput = FALSE)
  results$L_Call_mu_2[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal_2, n_value, isput = FALSE)
  results$L_Call_mu_3[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal_3, n_value, isput = FALSE)
  
  # Put
  results$CRR_Put[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_crr, n_value, isput = TRUE)
  results$L_Put_mu_1[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal, n_value, isput = TRUE)
  results$L_Put_mu_2[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal_2, n_value, isput = TRUE)
  results$L_Put_mu_3[index] <- BinOptionPrice(S = 100, K = 105, r=0.02, T_t = 0.5, sigma = 0.20, mu = mu_lognormal_3, n_value, isput = TRUE)
}

print(results)


# Price comparison for each step n and plotting convergence graphs
library(ggplot2)
results_long <- reshape2::melt(results, id.vars = "n", variable.name = "Model", value.name = "Price")

ggplot(results_long, aes(x = n, y = Price, color = Model)) +
  geom_line() +
  geom_hline(yintercept = bs_call, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = bs_put, linetype = "dashed", color = "red") +
  labs(title = "Convergence of option prices", x = "Number of steps (n)", y = "Option price") + theme_minimal() +
  scale_color_discrete(name = "Modèles", labels = c("CRR Call", "CRR Put", "LN Call mu 1", "LN Put mu 1", "LN Call mu 2", "LN Put mu 2", "LN Call mu 3", "LN Put mu 3"))


# 3e) 

S0 <- 100 
K <- 105
r <- 0.02
T <- 0.5
sigma <- 0.20
n <- 200  # number of steps for the tree (large n)

h <- T / n  # time interval for each step in tree
u <- exp(sigma * sqrt(h))  #since mu=0 under the hypotheses 
d <- exp(-sigma * sqrt(h)) 
q <- (exp(r * h) - d) / (u - d)  # risk neutral probability

levels_terminal <- 0:n  # define levels for which to generate terminal values of underlying
ST_binomial <- S0 * u^levels_terminal * d^(n - levels_terminal)
probas <- dbinom(levels_terminal, size = n, prob = q)  # probability associates with each terminal price

# define density function of lognormal
mu_lognormal <- log(S0) + (r - sigma^2 / 2) * T  # param mu
sigma_lognormal <- sqrt(sigma^2 * T)  # param sigma
x <- seq(min(ST_binomial) * 0.9, max(ST_binomial) * 1.1, length.out = 200)
densite_lognormal <- dlnorm(x, mu_lognormal, sigma_lognormal)

# graph of underlying density in binomial tree, and graph of random lognormal variable 
hist(ST_binomial, weights = probas, breaks = 100, freq = FALSE, col = rgb(0.2, 0.4, 0.6, 0.6), 
     xlim = range(x), ylim = range(0,0.03), main = "Density functions comparisons", xlab = "Final Price ST", ylab = "Density")
lines(x, densite_lognormal, col = "red", lwd = 2)
#legend("topright", legend = c("Arbre Binomial", "Lognormale"), col = c(rgb(0.2, 0.4, 0.6, 0.6), "red"), lwd = c(4, 2))


