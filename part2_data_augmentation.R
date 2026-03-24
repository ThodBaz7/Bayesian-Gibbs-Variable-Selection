# Data Augmentation for f(x) ∝ x^3 * exp(-5x) / (x+3)^2

# Gibbs sampler with latent variable w
set.seed(2025)

# MCMC parameters
M <- 10000
burnin <- 1000
total_iter <- M + burnin

# Initialize chains
x_samples <- rep(NA, total_iter)
w_samples <- rep(NA, total_iter)
x_samples[1] <- 1  
w_samples[1] <- 1

# Gibbs Sampler
for (i in 2:total_iter) {
  # Sample w | x : Exponential(x + 3)
  w_samples[i] <- rexp(1, rate = x_samples[i-1] + 3)
  
  # Sample x | w : Gamma(shape = 4, rate = 5 + w)
  x_samples[i] <- rgamma(1, shape = 4, rate = 5 + w_samples[i])
}

# Discard burn-in
x_post <- x_samples[(burnin+1):total_iter]
w_post <- w_samples[(burnin+1):total_iter]

# Estimate E[(X+3)/(X+2)]
quantity <- (x_post + 3) / (x_post + 2)
estimated_mean <- mean(quantity)


cat("Estimated E[(X+3)/(X+2)] =", round(estimated_mean, 4), "\n")

par(mfrow = c(1, 1))

# Plot convergence of cumulative mean
cumulative_mean <- cumsum(quantity) / (1:M)
plot(1:M, cumulative_mean, type = "l", lwd = 2,
     xlab = "Iteration", ylab = "Cumulative mean",
     main = expression("Convergence of E" * group("[", (X+3)/(X+2), "]")),
     col = "darkblue")
abline(h = estimated_mean, col = "red", lty = 2, lwd = 1.5)
grid()

# Plot autocorrelation
acf(quantity, main = "Autocorrelation of (X+3)/(X+2)", 
    lag.max = 50, col = "darkred")

