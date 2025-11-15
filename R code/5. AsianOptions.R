# Functions needed to run this file: From question 4

SimulateStockPaths <- function(S, mu, r, sigma, nsteps, sizesteps, npaths, theseed, isPmeasure) {
  set.seed(theseed)
  
  dt <- sizesteps # delta t
  N <- nsteps # total number of time steps
  M <- npaths # total number of trajectories to simulate(repeat the process M times)
  drift <- if (isPmeasure) mu else r - 0.5 * sigma^2
  
  paths <- matrix(0, nrow = M, ncol = N + 1)  # initialize the trajectory matrix
  
  paths[, 1] <- S # toutes les trajectoires commencent avec le prix initial S à t0
  
  # Trajectory simulations
  for (m in 1:M) {
    Z <- rnorm(N, mean = 0, sd = 1)  # generate N gaussian random variables
    
    for (i in 1:N) {
      # Calculating price at ti S(m)t_i with Black-Scholes
      paths[m, i + 1] <- paths[m, i] * exp(drift * dt + sigma * sqrt(dt) * Z[i])
    }
  }
  return(paths) 
}


SimulateStockPathsAntithetic <- function(S, mu, r, sigma, nsteps, sizesteps, npaths, theseed, isPmeasure) {
  if (npaths %% 2 != 0) {
    stop("The number of trajectories (npaths) is not a multiple of 2 (cannot use antithetic variables).")
  }
  
  set.seed(theseed)
  
  dt <- sizesteps
  N <- nsteps
  M <- npaths
  
  paths <- matrix(0, nrow = M, ncol = N + 1)
  paths[, 1] <- S
  
  drift <- if (isPmeasure) mu else r - 0.5 * sigma^2
  
  # generating gaussian random variables for the first half of the trajectories
  Z <- matrix(rnorm((M / 2) * N, mean = 0, sd = 1), nrow = M / 2, ncol = N)
  
  # simulating trajectories for the first half until index M/2
  for (m in 1:(M / 2)) {
    for (i in 1:N) {
      paths[m, i + 1] <- paths[m, i] * exp(drift * dt + sigma * sqrt(dt) * Z[m, i])
    }
  }
  # generate the antithetical trajectories for the second half
  for (m in 1:(M / 2)) {
    for (i in 1:N) {
      paths[M / 2 + m, i + 1] <- paths[M / 2 + m, i] * exp(drift * dt - sigma * sqrt(dt) * Z[m, i])
    }
  }
  return(paths)
}

StockPathsQ <- SimulateStockPaths(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                  npaths = M, theseed = seed, isPmeasure = FALSE)
StockPathsQantithetic <- SimulateStockPathsAntithetic(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                                      npaths = M, theseed = seed, isPmeasure = FALSE)

#5a) 

S0 <- 100
K <- 100
r <- 0.02
T <- 0.5
M <- 10000
nsteps <- 26  # number of time steps (weekly)
dt <- 1 / 52
discount_factor <- exp(-r * T)

# calculate arithmetic mean for each trajectory
calculate_arithmetic_mean <- function(stock_paths) {
  rowMeans(stock_paths[, 2:(nsteps + 1)])  # excluding column S0 (t=0)
}

# arithmetic means of prices for StockPathsQ
A_T_Q <- calculate_arithmetic_mean(StockPathsQ)
A_T_Q_antithetic <- calculate_arithmetic_mean(StockPathsQantithetic)

# payoff of asian option
H_Q <- pmax(0, A_T_Q - K)
H_Q_antithetic <- pmax(0, A_T_Q_antithetic - K)

# Monte Carlo price
price_Q <- discount_factor * mean(H_Q)
price_Q_antithetic <- discount_factor * mean(H_Q_antithetic)

# Standard deviation of payoffs
std_Q <- sd(H_Q)
std_Q_antithetic <- sd(H_Q_antithetic)

# IC 95%
IC_Q <- c(price_Q - 1.96 * (std_Q / sqrt(M)), price_Q + 1.96 * (std_Q / sqrt(M)))
IC_Q_antithetic <- c(price_Q_antithetic - 1.96 * (std_Q_antithetic / sqrt(M)),
                     price_Q_antithetic + 1.96 * (std_Q_antithetic / sqrt(M)))

cat("Monte Carlo price (StockPathsQ) :", price_Q, "\n")
cat("Confidence interval 95% (StockPathsQ) :", IC_Q, "\n")
cat("Monte Carlo price (StockPathsQantithetic) :", price_Q_antithetic, "\n")
cat("Confidence interval 95% (StockPathsQantithetic) :", IC_Q_antithetic, "\n")


# 5c) 

n <- 26
Delta_t <- 1 / 52
E_A_ari_Q <- S0 * exp(r * Delta_t) * (1 - exp(r * n * Delta_t)) / (n * (1 - exp(r * Delta_t)))
cat("E[A_T^(ari)] sous Q :", E_A_ari_Q, "\n")

# calculate the value of c*
# extract the simulated trajectories
A_T_ari <- rowMeans(StockPathsQ[, 2:(n + 1)])  # arithmetic mean for each trajectory
H_ari <- exp(-r * T) * pmax(0, A_T_ari - K)    # payoff of asian option

# calculating moments
cov_H_A <- cov(H_ari, A_T_ari)                # Covariance between H and A_T^(ari)
var_A <- var(A_T_ari)                         # Variance of A_T^(ari)

# calculating c*
c_star <- -cov_H_A / var_A
cat("Value of c* :", c_star, "\n")

# estimation of the asian option corrected price
Pi_0_ari_control <- mean(H_ari + c_star * (A_T_ari - E_A_ari_Q))

# calculation of the confidence interval 95% 
std_error <- sd(H_ari + c_star * (A_T_ari - E_A_ari_Q)) / sqrt(length(H_ari))
confidence_interval <- c(Pi_0_ari_control - 1.96 * std_error, Pi_0_ari_control + 1.96 * std_error)

cat("Estimated price of asian option (controlled) :", Pi_0_ari_control, "\n")
cat("Confidence interval at 95% :", confidence_interval, "\n")

