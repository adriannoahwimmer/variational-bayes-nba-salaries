# ============================================================
# 06_simulation_runtime.R
# Runtime simulation: CAVI vs MCMC for different n and p
# ============================================================

# Requires:
# source("R/02_mfvb_function.R")
# source("R/03_gibbs_sampler.R")

library(dplyr)

# ============================================================
# 1. Simulate Bayesian linear regression data
# ============================================================

simulate_blr_data <- function(n, p, sigma = 1, seed = 123) {
  
  set.seed(seed)
  
  X <- matrix(
    rnorm(n * (p - 1)),
    nrow = n,
    ncol = p - 1
  )
  
  X <- cbind("(Intercept)" = 1, X)
  
  beta_true <- matrix(rnorm(p), ncol = 1)
  
  y <- X %*% beta_true +
    rnorm(n, 0, sigma)
  
  return(list(
    X_data = X,
    Y_data = y,
    beta_true = beta_true
  ))
}


# ============================================================
# 2. Runtime function for one scenario
# ============================================================

run_runtime_scenario <- function(
    n,
    p,
    n_mcmc_iter = 1000,
    burn_in = 200,
    cavi_max_iter = 25,
    seed = 123
) {
  
  sim_data <- simulate_blr_data(
    n = n,
    p = p,
    seed = seed
  )
  
  X_data <- sim_data$X_data
  Y_data <- sim_data$Y_data
  
  Mu_beta_prior <- matrix(0, nrow = p, ncol = 1)
  Sigma_beta_prior <- diag(100, p)
  
  A_prior <- 1
  B_prior <- 1
  
  # ------------------------------------------------------------
  # CAVI runtime
  # ------------------------------------------------------------
  
  time_cavi <- system.time({
    cavi_result <- mfvb_algo(
      Mu_beta_prior = Mu_beta_prior,
      Sigma_beta_prior = Sigma_beta_prior,
      A_prior = A_prior,
      B_prior = B_prior,
      X_data = X_data,
      Y_data = Y_data,
      tol = 1e-8,
      max_iter = cavi_max_iter,
      verbose = FALSE
    )
  })
  
  # ------------------------------------------------------------
  # MCMC runtime
  # ------------------------------------------------------------
  
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
      verbose = FALSE
    )
  })
  
  result <- data.frame(
    n = n,
    p = p,
    CAVI_time_seconds = as.numeric(time_cavi["elapsed"]),
    MCMC_time_seconds = as.numeric(time_mcmc["elapsed"]),
    CAVI_iterations = cavi_result$n_iter,
    MCMC_iterations = n_mcmc_iter
  )
  
  return(result)
}


# ============================================================
# 3. Run scenarios
# ============================================================

simulation_scenarios <- data.frame(
  n = c(100, 1000, 5000, 10000),
  p = c(10, 100, 500, 1000) 
)

runtime_simulation_results <- bind_rows(
  lapply(1:nrow(simulation_scenarios), function(i) {
    
    cat("\nRunning scenario:",
        "n =", simulation_scenarios$n[i],
        ", p =", simulation_scenarios$p[i], "\n")
    
    run_runtime_scenario(
      n = simulation_scenarios$n[i],
      p = simulation_scenarios$p[i],
      n_mcmc_iter = 1000,
      burn_in = 200,
      cavi_max_iter = 25,
      seed = 123 + i
    )
  })
)

runtime_simulation_results

# Optional: nicer rounded table
runtime_simulation_results_rounded <- runtime_simulation_results |>
  mutate(
    CAVI_time_seconds = round(CAVI_time_seconds, 4),
    MCMC_time_seconds = round(MCMC_time_seconds, 4),
    Speedup_MCMC_over_CAVI = round(MCMC_time_seconds / CAVI_time_seconds, 2)
  )

runtime_simulation_results_rounded
