# --- Load data ---
data = read.table("real_estate_data.txt", header=TRUE)

# Dimensions
n = nrow(data)
p = ncol(data)

# Prepare matrices
y = data[,1]
X = as.matrix(cbind(1, data[,2:4]))

# -------------------------------------------------------------------
# 1. Classical multiple linear regression
# -------------------------------------------------------------------
regression = lm(y ~ X[, -1])
summary_ols = summary(regression)

print(coef(summary_ols))

significant = which(coef(summary_ols)[, 4] < 0.05)
cat("\nStatistically significant variables (p<0.05):\n")
print(names(significant))

# -------------------------------------------------------------------
# 2. Bayesian regression with all variables (non‑informative priors)
# -------------------------------------------------------------------
set.seed(2025)
# Hyperparameters
b0 = rep(0, p)
tau02 = 1e-5
gam = 0.005
del = 0.005

# Pre‑compute quantities
tXy = t(X) %*% y
tXX = t(X) %*% X
inv.tXX = solve(tXX)
b.hat = unname(lm(y ~ X[, -1])$coef)

# MCMC parameters
M = 10000
burnin = 1000
total = M + burnin

# Storage
b = array(NA, c(total, p))
tau2 = rep(NA, total)

# Initial values
b[1,] = b.hat
tau2[1] = 1/n

library(MASS)

for (m in 2:total) {
  Sigma = inv.tXX / (tau2[m-1] + tau02)
  mu = (tau2[m-1] * b.hat + tau02 * b0) / (tau2[m-1] + tau02)
  b[m,] = mvrnorm(1, mu, Sigma)
  tau2[m] = rgamma(1, (gam + n)/2, (del + sum((y - X %*% b[m,])^2))/2)
}

# Discard burn‑in
b = tail(b, M)
tau2 = tail(tau2, M)

# Results
post_means = colMeans(b)
names(post_means) = c("intercept", "sqm", "age", "dist_center")

post_ci = t(apply(b, 2, quantile, probs = c(0.025, 0.975)))
colnames(post_ci) = c("2.5%", "97.5%")
rownames(post_ci) = c("intercept", "sqm", "age", "dist_center")

cat("\nPosterior means:\n")
print(round(post_means, 4))

cat("\n95% credible intervals:\n")
print(round(post_ci, 4))

# -------------------------------------------------------------------
# 3. Hierarchical model with variable selection (binary indicators w)
# -------------------------------------------------------------------
set.seed(2025)
# Hyperparameters (γ0 = δ0 = ε = ζ = 0.01)
b0 = rep(0, p)
gam = 0.01
gam0 = 0.01
del0 = 0.01
eps = 0.01
zet = 0.01
q = rep(0.5, p)           # prior inclusion probabilities

# Storage
b = array(NA, c(total, p))
tau2 = rep(NA, total)
w = array(1, c(total, p))
w[,1] = 1                 # intercept always included
tau02 = rep(NA, total)
del = rep(NA, total)      # hyperparameter for tau2 (delta in the PDF)

# Initial values
b[1,] = b.hat
tau2[1] = 1/n
tau02[1] = gam0 / del0
del[1] = eps / zet

for (m in 2:total) {
  w[m,] = w[m-1,]
  for (j in 2:p) {
    # Try w_j = 1
    w[m,j] = 1
    WW = diag(w[m,])
    a1 = q[j] * exp(-tau2[m-1] * sum((y - X %*% WW %*% b[m-1,])^2) / 2)
    # Try w_j = 0
    w[m,j] = 0
    WW = diag(w[m,])
    a0 = (1 - q[j]) * exp(-tau2[m-1] * sum((y - X %*% WW %*% b[m-1,])^2) / 2)
    prob = a1 / (a1 + a0)
    w[m,j] = rbinom(1, 1, prob)
  }
  
  # Update beta | w, tau2, tau02
  WW = diag(w[m,])
  Sigma = solve(tau2[m-1] * WW %*% tXX %*% WW + tau02[m-1] * tXX)
  mu = Sigma %*% (tau2[m-1] * WW %*% tXy + tau02[m-1] * tXX %*% b0)
  b[m,] = mvrnorm(1, mu, Sigma)
  
  # Update tau2 | beta, w, delta
  WW = diag(w[m,])
  SSR = sum((y - X %*% WW %*% b[m,])^2)
  tau2[m] = rgamma(1, (gam + n)/2, (del[m-1] + SSR)/2)
  
  # Update tau02 | beta
  tau02[m] = rgamma(1, (gam0 + p)/2, (del0 + t(b[m,] - b0) %*% tXX %*% (b[m,] - b0))/2)
  
  # Update delta | tau2
  del[m] = rgamma(1, (eps + gam)/2, (zet + tau2[m])/2)
}

# Discard burn‑in
b = tail(b, M)
tau2 = tail(tau2, M)
w = tail(w, M)
tau02 = tail(tau02, M)
del = tail(del, M)

# Posterior inclusion probabilities
post_prob_w = colMeans(w)
names(post_prob_w) = c("intercept", "sqm", "age", "dist_center")

cat("\nPosterior inclusion probabilities:\n")
print(round(post_prob_w, 4))

# -------------------------------------------------------------------
# 4. Convert w vectors to model numbers and find the most probable model
# -------------------------------------------------------------------
w.to.model = function(w) {
  1 + sum(w[-1] * 2^(0:(length(w)-2)))
}

model = as.vector(apply(w, 1, w.to.model))
model_freq = table(model) / M

cat("\nModel frequencies:\n")
print(round(model_freq, 4))

best_model = as.numeric(names(which.max(model_freq)))
best_model_prob = max(model_freq)

cat("\nMost probable model:", best_model)
cat("\nPosterior probability:", round(best_model_prob, 4), "\n")

model.to.w = function(z) {
  z = z - 1
  res = c()
  for (j in (p-2):0) {
    res = c(res, floor(z / 2^j))
    z = z - (z >= 2^j) * 2^j
  }
  c(1, rev(res))
}

# -------------------------------------------------------------------
# 5. Inference for the most probable model
# -------------------------------------------------------------------
best_indices = which(model == best_model)
b_best = b[best_indices, ]
w_best = model.to.w(best_model)
included = which(w_best == 1)
var_names = c("intercept", "sqm", "age", "dist_center")[included]

post_means_best = colMeans(b_best)[included]
names(post_means_best) = var_names

cat("\nPosterior means of active coefficients:\n")
print(round(post_means_best, 4))

# Convergence: cumulative means of standardized coefficients
library(coda)

b_best_subset = b_best[, included]
b_best_scaled = scale(b_best_subset)
M_best = length(best_indices)
cum_means = apply(b_best_scaled, 2, cumsum) / (1:M_best)

par(mfrow = c(1, 1))
matplot(1:M_best, cum_means, type = "l", 
        col = c("red", "blue", "green")[1:length(included)],
        lty = 1, lwd = 1.5,
        xlab = "Iteration", ylab = "Standardized mean",
        main = paste("Cumulative means - Model", best_model),
        ylim = c(-0.3, 0.3))
abline(h = 0, lty = 2, col = "gray")
legend("topright", legend = var_names, 
       col = c("red", "blue", "green")[1:length(included)], lty = 1)


# -------------------------------------------------------------------
# 6. Predictions for new apartments using the most probable model
# -------------------------------------------------------------------
new_data = data.frame(
  type = c("Small apartment (city center)", "New apartment (suburb)", 
           "Large apartment (semi‑urban)"),
  sqm = c(55, 85, 120),
  age = c(25, 5, 15),
  dist_center = c(1.5, 8.0, 12)
)

Psi = as.matrix(cbind(1, new_data[, c("sqm", "age", "dist_center")]))
k = nrow(Psi)

tau2_best = tau2[best_indices]
z_best = array(NA, c(M_best, k))
set.seed(2025)
for (m in 1:M_best) {
  beta_w = b_best[m, ] * w_best
  mean_z = Psi %*% beta_w
  z_best[m, ] = rnorm(k, mean = mean_z, sd = 1/sqrt(tau2_best[m]))
}

point_best = colMeans(z_best)
ci_best = t(apply(z_best, 2, quantile, probs = c(0.025, 0.975)))

results_best = data.frame(
  Apartment = new_data$type,
  Prediction = round(point_best, 2),
  CI_2.5 = round(ci_best[, 1], 2),
  CI_97.5 = round(ci_best[, 2], 2)
)
print(results_best)

# -------------------------------------------------------------------
# 7. Model averaging predictions
# -------------------------------------------------------------------
z_all = array(NA, c(M, k))
set.seed(2025)
for (m in 1:M) {
  WW = diag(w[m, ])
  mean_z = Psi %*% (WW %*% b[m, ])
  z_all[m, ] = rnorm(k, mean = mean_z, sd = 1/sqrt(tau2[m]))
}

point_avg = colMeans(z_all)
ci_avg = t(apply(z_all, 2, quantile, probs = c(0.025, 0.975)))

results_avg = data.frame(
  Apartment = new_data$type,
  Prediction = round(point_avg, 2),
  CI_2.5 = round(ci_avg[, 1], 2),
  CI_97.5 = round(ci_avg[, 2], 2)
)
print(results_avg)

cat("\n--- Comparison of methods ---\n")
comparison = data.frame(
  Apartment = new_data$type,
  Best_Model = round(point_best, 2),
  Best_CI = paste0("[", round(ci_best[,1], 2), ",", round(ci_best[,2], 2), "]"),
  Averaging = round(point_avg, 2),
  Avg_CI = paste0("[", round(ci_avg[,1], 2), ",", round(ci_avg[,2], 2), "]")
)
print(comparison)
