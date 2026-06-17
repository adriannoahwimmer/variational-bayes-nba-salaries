# ============================================================
# 02_mfvb_function.R
# Mean Field Variational Bayes Algorithm
# Bayesian Linear Regression:
#
# y = X beta + epsilon, epsilon ~ N(0, sigma^2 I)
#
# Prior:
# beta    ~ N(mu_beta_prior, Sigma_beta_prior)
# sigma^2 ~ IG(a_sigma2_prior, b_sigma2_prior)
#
# Variational approximation:
# q(beta, sigma^2) = q(beta) q(sigma^2)
# q(beta)          = N(mu_q_beta, Sigma_q_beta)
# q(sigma^2)       = IG(a_q_sigma2, b_q_sigma2)
# ============================================================

library(psych)

# ============================================================
# ELBO for Bayesian Linear Regression
# ============================================================

compute_elbo_blr <- function(
    X_data,
    Y_data,
    Mu_q_beta,
    Sigma_q_beta,
    A_q,
    B_q,
    Mu_beta_prior,
    Sigma_beta_prior,
    A_prior,
    B_prior
) {
  
  X_data <- as.matrix(X_data)
  Y_data <- as.matrix(Y_data)
  
  n <- nrow(X_data)
  p <- ncol(X_data)
  
  A_q <- as.numeric(A_q)
  B_q <- as.numeric(B_q)
  
  Sigma_beta_prior_inv <- solve(Sigma_beta_prior)
  
  logdet <- function(M) {
    as.numeric(determinant(M, logarithm = TRUE)$modulus)
  }
  
  # ------------------------------------------------------------
  # Useful expectations under q(sigma^2)
  # ------------------------------------------------------------
  
  E_log_sigma2 <- log(B_q) - digamma(A_q)
  E_inv_sigma2 <- A_q / B_q
  
  # ------------------------------------------------------------
  # E_q(beta)[(y - X beta)'(y - X beta)]
  # =
  # (y - X mu_q)'(y - X mu_q)
  # + tr(X'X Sigma_q)
  # ------------------------------------------------------------
  
  residual_mean <- Y_data - X_data %*% Mu_q_beta
  
  expected_ssr <- as.numeric(
    t(residual_mean) %*% residual_mean +
      tr(t(X_data) %*% X_data %*% Sigma_q_beta)
  )
  
  # ------------------------------------------------------------
  # 1. Expected log likelihood
  # E_q[log p(y | beta, sigma^2)]
  # ------------------------------------------------------------
  
  elbo_likelihood <-
    -n / 2 * log(2 * pi) -
    n / 2 * E_log_sigma2 -
    0.5 * E_inv_sigma2 * expected_ssr
  
  # ------------------------------------------------------------
  # 2. Expected log prior for beta
  # E_q[log p(beta)]
  # ------------------------------------------------------------
  
  beta_diff <- Mu_q_beta - Mu_beta_prior
  
  elbo_beta_prior <-
    -p / 2 * log(2 * pi) -
    0.5 * logdet(Sigma_beta_prior) -
    0.5 * (
      tr(Sigma_beta_prior_inv %*% Sigma_q_beta) +
        as.numeric(t(beta_diff) %*% Sigma_beta_prior_inv %*% beta_diff)
    )
  
  # ------------------------------------------------------------
  # 3. Expected log prior for sigma^2
  # E_q[log p(sigma^2)]
  # ------------------------------------------------------------
  
  elbo_sigma2_prior <-
    A_prior * log(B_prior) -
    lgamma(A_prior) -
    (A_prior + 1) * E_log_sigma2 -
    B_prior * E_inv_sigma2
  
  # ------------------------------------------------------------
  # 4. Entropy of q(beta)
  # -E_q[log q(beta)]
  # ------------------------------------------------------------
  
  entropy_beta <-
    0.5 * logdet(Sigma_q_beta) +
    p / 2 * (1 + log(2 * pi))
  
  # ------------------------------------------------------------
  # 5. Entropy of q(sigma^2)
  # -E_q[log q(sigma^2)]
  # ------------------------------------------------------------
  
  entropy_sigma2 <-
    A_q +
    log(B_q) +
    lgamma(A_q) -
    (A_q + 1) * digamma(A_q)
  
  # ------------------------------------------------------------
  # Total ELBO
  # ------------------------------------------------------------
  
  elbo <-
    elbo_likelihood +
    elbo_beta_prior +
    elbo_sigma2_prior +
    entropy_beta +
    entropy_sigma2
  
  return(as.numeric(elbo))
}



mfvb_algo <- function(
    Mu_beta_prior,
    Sigma_beta_prior,
    A_prior,
    B_prior,
    X_data,
    Y_data,
    tol = 1e-8,
    max_iter = 25,
    verbose = TRUE
) {
  
  # ------------------------------------------------------------
  # Prepare data
  # ------------------------------------------------------------
  
  X_data <- as.matrix(X_data)
  Y_data <- as.matrix(Y_data)
  
  n <- nrow(X_data)
  
  Sigma_beta_prior_inv <- solve(Sigma_beta_prior)
  
  XtX <- t(X_data) %*% X_data
  XtY <- t(X_data) %*% Y_data
  
  
  # ------------------------------------------------------------
  # Equation 1:
  # a_q(sigma^2) = a_prior + n / 2
  #
  # This quantity is constant because a_prior and n are constant.
  # ------------------------------------------------------------
  
  A_q_sigma2_post_approx <- A_prior + n / 2
  
  
  # ------------------------------------------------------------
  # Initialize variational parameters for q(beta)
  #
  # These are only starting values for the CAVI algorithm.
  # The prior itself is not updated.
  # ------------------------------------------------------------
  
  Mu_q_beta_post_approx <- Mu_beta_prior
  Sigma_q_beta_post_approx <- Sigma_beta_prior
  
  
  # ------------------------------------------------------------
  # CAVI loop
  # ------------------------------------------------------------
  
  diff <- Inf
  iter <- 0
  elbo_values <- numeric(0)
  
  while (diff > tol && iter < max_iter) {
    
    mu_old <- Mu_q_beta_post_approx
    
    
    # ------------------------------------------------------------
    # Equation 2:
    # B_q(sigma^2)
    #
    # B_q(sigma^2) =
    # B_prior
    # + 1/2 * E_q(beta)[(y - X beta)'(y - X beta)]
    #
    # The expectation is:
    #
    # (y - X mu_q_beta)'(y - X mu_q_beta)
    # + tr(X'X Sigma_q_beta)
    # ------------------------------------------------------------
    
    residual_mean <- Y_data - X_data %*% Mu_q_beta_post_approx
    
    B_q_sigma2_post_approx <- B_prior +
      0.5 * (
        as.numeric(t(residual_mean) %*% residual_mean) +
          tr(XtX %*% Sigma_q_beta_post_approx)
      )
    
    
    # ------------------------------------------------------------
    # Expected inverse variance:
    #
    # E_q(1 / sigma^2) = a_q(sigma^2) / B_q(sigma^2)
    # ------------------------------------------------------------
    
    E_inv_sigma2 <- as.numeric(
      A_q_sigma2_post_approx / B_q_sigma2_post_approx
    )
    
    
    # ------------------------------------------------------------
    # Equation 3:
    # Sigma_q_beta
    #
    # Sigma_q_beta =
    # ( E_q(1/sigma^2) X'X + Sigma_beta_prior^(-1) )^(-1)
    # ------------------------------------------------------------
    
    Sigma_q_beta_post_approx <- solve(
      E_inv_sigma2 * XtX + Sigma_beta_prior_inv
    )
    
    
    # ------------------------------------------------------------
    # Equation 4:
    # mu_q(beta)
    #
    # mu_q_beta =
    # Sigma_q_beta *
    # ( E_q(1/sigma^2) X'y
    #   + Sigma_beta_prior^(-1) mu_beta_prior )
    # ------------------------------------------------------------
    
    Mu_q_beta_post_approx <- Sigma_q_beta_post_approx %*% (
      E_inv_sigma2 * XtY +
        Sigma_beta_prior_inv %*% Mu_beta_prior
    )
    
    
    # ------------------------------------------------------------
    # Convergence check
    # ------------------------------------------------------------
    
    diff <- max(abs(Mu_q_beta_post_approx - mu_old))
    
    # ------------------------------------------------------------
    # Compute ELBO after current CAVI update
    # ------------------------------------------------------------
    
    elbo_current <- compute_elbo_blr(
      X_data = X_data,
      Y_data = Y_data,
      Mu_q_beta = Mu_q_beta_post_approx,
      Sigma_q_beta = Sigma_q_beta_post_approx,
      A_q = A_q_sigma2_post_approx,
      B_q = B_q_sigma2_post_approx,
      Mu_beta_prior = Mu_beta_prior,
      Sigma_beta_prior = Sigma_beta_prior,
      A_prior = A_prior,
      B_prior = B_prior
    )
    
    elbo_values <- c(elbo_values, elbo_current)
    
    
    # ------------------------------------------------------------
    # Optional output
    # ------------------------------------------------------------
    
    if (verbose) {
      cat("\n----------------------------------\n")
      cat("Iteration:", iter, "\n")
      cat("----------------------------------\n")
      
      cat("A_q_sigma2_post_approx:\n")
      print(A_q_sigma2_post_approx)
      
      cat("B_q_sigma2_post_approx:\n")
      print(B_q_sigma2_post_approx)
      
      cat("Mu_q_sigma2_post_approx:\n")
      print(B_q_sigma2_post_approx / (A_q_sigma2_post_approx - 1))
      
      cat("Mu_q_beta_post_approx:\n")
      print(Mu_q_beta_post_approx)
      
      cat("Sigma_q_beta_post_approx:\n")
      print(Sigma_q_beta_post_approx)
      
      cat("Current diff:\n")
      print(diff)
      
      cat("ELBO:\n")
      print(elbo_current)
    }
    
    iter <- iter + 1
  }
  
 
  
  # ------------------------------------------------------------
  # Store result
  # ------------------------------------------------------------
  
  result <- list(
    A_q_sigma2_post_approx = A_q_sigma2_post_approx,
    B_q_sigma2_post_approx = B_q_sigma2_post_approx,
    Mu_q_sigma2_post_approx = B_q_sigma2_post_approx /
      (A_q_sigma2_post_approx - 1),
    Mu_q_beta_post_approx = Mu_q_beta_post_approx,
    Sigma_q_beta_post_approx = Sigma_q_beta_post_approx,
    elbo_values = elbo_values,
    elbo_final = tail(elbo_values, 1),
    n_iter = iter,
    final_diff = diff
  )
  
  return(result)
}