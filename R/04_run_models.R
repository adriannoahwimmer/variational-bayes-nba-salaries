# ============================================================
# 04_run_models.R
# Set priors, run OLS, MFVB/CAVI and MCMC/Gibbs
# ============================================================

# Requires:
# source("R/01_load_clean_data.R")
# source("R/02_mfvb_function.R")
# source("R/03_gibbs_sampler.R")


# ============================================================
# 1. Set priors
# ============================================================

p <- ncol(X_data)
n <- nrow(X_data)

Mu_beta_prior <- matrix(0, nrow = p, ncol = 1)
Sigma_beta_prior <- diag(100, p)

A_prior <- 1
B_prior <- 1


# ============================================================
# 2. Run OLS
# ============================================================

time_ols <- system.time({
  ols_fit <- lm(Y_data ~ X_data - 1)
  
})

ols_beta_scaled <- coef(ols_fit)
names(ols_beta_scaled) <- colnames(X_data)
summary(ols_fit)


# ============================================================
# 3. Run CAVI / MFVB
# ============================================================

time_cavi <- system.time({
  cavi_result <- mfvb_algo(
    Mu_beta_prior = Mu_beta_prior,
    Sigma_beta_prior = Sigma_beta_prior,
    A_prior = A_prior,
    B_prior = B_prior,
    X_data = X_data,
    Y_data = Y_data,
    tol = 1e-8,
    max_iter = 25,
    verbose = TRUE
  )
})

cavi_beta_scaled <- as.numeric(cavi_result$Mu_q_beta_post_approx)
names(cavi_beta_scaled) <- colnames(X_data)

cavi_sigma2_mean <- as.numeric(
  cavi_result$Mu_q_sigma2_post_approx
)


# ============================================================
# 4. Run MCMC / Gibbs
# ============================================================

set.seed(123)

n_mcmc_iter <- 10000
burn_in <- 2000

time_mcmc <- system.time({
  mcmc_result <- gibbs_bayes_lm(
    X_data = X_data,
    Y_data = Y_data,
    n_iter = n_mcmc_iter,
    burn_in = burn_in,
    A_prior = A_prior,
    B_prior = B_prior,
    Mu_beta_prior = Mu_beta_prior,
    Sigma_beta_prior = Sigma_beta_prior,
    verbose = TRUE
  )
})

beta_samples <- mcmc_result$beta_samples
sigma2_samples <- mcmc_result$sigma2_samples

mcmc_beta_scaled <- colMeans(beta_samples)
mcmc_sigma2_mean <- mean(sigma2_samples)

mcmc_beta_quantiles <- apply(
 beta_samples,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)


# ============================================================
# 5. Compare betas on the standardised scale
# ============================================================

beta_comparison_scaled <- data.frame(
  variable = colnames(X_data),
  OLS = as.numeric(ols_beta_scaled),
  CAVI_mean = as.numeric(cavi_beta_scaled),
  MCMC_mean = as.numeric(mcmc_beta_scaled)
)

beta_comparison_scaled

print(
  cbind(
    beta_comparison_scaled["variable"],
    round(beta_comparison_scaled[, -1], 4)
  ),
  row.names = FALSE
)


# ============================================================
# 6. Runtime comparison
# ============================================================

runtime_results <- data.frame(
  Method = c("OLS", "MFVB / CAVI", "MCMC / Gibbs"),
  Runtime_seconds = c(
    as.numeric(time_ols["elapsed"]),
    as.numeric(time_cavi["elapsed"]),
    as.numeric(time_mcmc["elapsed"])
  ),
  Iterations = c(
    1,
    cavi_result$n_iter,
    n_mcmc_iter
  ),
  n = c(n, n, n),
  p = c(p, p, p)
)

runtime_results


# ============================================================
# 7. Sigma² comparison
# ============================================================

sigma2_comparison <- data.frame(
  Method = c("CAVI", "MCMC"),
  Sigma2_mean = c(
    cavi_sigma2_mean,
    mcmc_sigma2_mean
  )
)

sigma2_comparison