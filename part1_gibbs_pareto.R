# Gibbs sampling for Pareto(a, b)

# 1. Load data
pop_data <- read.table("population_data.txt", header = FALSE, quote = "\"", fill = TRUE)

head(pop_data)

pop_data <- pop_data[-1, ]

colnames(pop_data) <- c("index", "state", "population")

pop_data$population <- as.numeric(as.character(pop_data$population))
pop_data$index <- as.numeric(as.character(pop_data$index))

head(pop_data)

x <- pop_data$population  # population data
n <- length(x)
x_min <- min(x)

# MCMC parameters
M <- 10000
burnin <- 1000
total_iter <- M + burnin

# Create vectors for storage
a_chain <- numeric(total_iter)
b_chain <- numeric(total_iter)

b_chain[1] <- x_min * 0.9  
a_chain[1] <- 2  

sum_log_x <- sum(log(x))
log_x_min <- log(x_min)

set.seed(2025)

# Gibbs Sampler
for (i in 2:total_iter) {
  # Conditional of a | b
  rate_a <- sum_log_x - n * log(b_chain[i-1])
  if (rate_a <= 0) {
    rate_a <- 0.001
  }
  a_chain[i] <- rgamma(1, shape = n + 1, rate = rate_a)
  
  # Conditional of b | a
  y <- rbeta(1, shape1 = a_chain[i] * n + 1, shape2 = 1)
  b_chain[i] <- x_min * y
  if (b_chain[i] >= x_min) {
    b_chain[i] <- x_min * 0.999  # fix boundary
  }
}

# Discard burn-in
a_post <- tail(a_chain, M)
b_post <- tail(b_chain, M)

# Compute estimates (posterior means)
a_hat <- mean(a_post)
b_hat <- mean(b_post)

# 95% credible intervals (as posterior quantiles)
a_ci <- quantile(a_post, probs = c(0.025, 0.975))
b_ci <- quantile(b_post, probs = c(0.025, 0.975))

# Create table
results_table <- data.frame(
  Parameter = c("a", "b"),
  Estimate = c(a_hat, b_hat),
  Lower_95 = c(a_ci[1], b_ci[1]),
  Upper_95 = c(a_ci[2], b_ci[2])
)

print(results_table)
