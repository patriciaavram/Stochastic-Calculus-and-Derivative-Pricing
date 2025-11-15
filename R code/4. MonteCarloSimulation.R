
# 4a) 
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



# 4b)
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


# 4c) 

S0 <- 100
mu <- 0.07 # under P
r <- 0.02 # under Q
sigma <- 0.2
T <- 0.5
dt <- 1 / 52
nsteps <- T / dt
M <- 10000 # number of trajectories
seed <- 123

# simulations under real measure 
StockPathsP <- SimulateStockPaths(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                  npaths = M, theseed = seed, isPmeasure = TRUE)
StockPathsPantithetic <- SimulateStockPathsAntithetic(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                                      npaths = M, theseed = seed, isPmeasure = TRUE)

# simulations under risk neutral measure
StockPathsQ <- SimulateStockPaths(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                  npaths = M, theseed = seed, isPmeasure = FALSE)
StockPathsQantithetic <- SimulateStockPathsAntithetic(S = S0, mu, r, sigma, nsteps, sizesteps = dt, 
                                                      npaths = M, theseed = seed, isPmeasure = FALSE)

# checking the dimension of the matrices
cat("Dimensions StockPathsP :", dim(StockPathsP), "\n")
cat("Dimensions StockPathsPantithetic :", dim(StockPathsPantithetic), "\n")
cat("Dimensions StockPathsQ :", dim(StockPathsQ), "\n")
cat("Dimensions StockPathsQantithetic :", dim(StockPathsQantithetic), "\n")



# 4d) 

# exact parameters of the expected value
T <- 0.5
mu <- 0.07
r <- 0.02
S0 <- 100

# exact expected values under P and Q
E_P_exact <- S0 * exp(mu * T)
E_Q_exact <- S0 * exp(r * T)

cat("Exact expected value under P :", E_P_exact, "\n")
cat("Exact expected value under Q :", E_Q_exact, "\n")

# terminal average for P and Q for monte carlo
E_P_hat <- mean(StockPathsP[, nsteps + 1])  # ST in StockPathsP
E_Q_hat <- mean(StockPathsQ[, nsteps + 1])  # ST in StockPathsQ

# terminal standard deviation for P and Q for monte carlo
SD_P_hat <- sd(StockPathsP[, nsteps + 1])  # STD of ST in StockPathsP
SD_Q_hat <- sd(StockPathsQ[, nsteps + 1])  # STD of ST in StockPathsQ

M <- nrow(StockPathsP) # M trajectories

# confidence interval at 95%
IC_P <- c(E_P_hat - 1.96 * (SD_P_hat / sqrt(M)), E_P_hat + 1.96 * (SD_P_hat / sqrt(M)))
IC_Q <- c(E_Q_hat - 1.96 * (SD_Q_hat / sqrt(M)), E_Q_hat + 1.96 * (SD_Q_hat / sqrt(M)))

cat("Exact expected value of Monte Carlo under P :", E_P_hat, "\n")
cat("Confidence interval at 95% confidence under P :", IC_P, "\n")

cat("Exact expected value of Monte Carlo under Q :", E_Q_hat, "\n")
cat("Confidence interval at 95% confidence under Q :", IC_Q, "\n")

# comparison
cat("Relative error under P :", abs(E_P_hat - E_P_exact) / E_P_exact * 100, "%\n")
cat("Relative error under Q :", abs(E_Q_hat - E_Q_exact) / E_Q_exact * 100, "%\n")

