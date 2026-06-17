# ============================================================
# 03_gibbs_sampler.R
# MCMC with Gibbs Sampler for Bayesian Linear Regression
# ============================================================

gibbs_bayes_lm <- function(
    X_data,
    Y_data,
    n_iter = 10000,
    burn_in = 2000,
    A_prior = 1,
    B_prior = 1,
    Mu_beta_prior = NULL,
    Sigma_beta_prior = NULL,
    verbose = FALSE
) {
  
  # ------------------------------------------------------------
  # Prepare data
  # ------------------------------------------------------------
  
  X <- as.matrix(X_data)
  y <- as.matrix(Y_data)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }
  
  
  # ------------------------------------------------------------
  # Default priors
  # ------------------------------------------------------------
  
  if (is.null(Mu_beta_prior)) {
    Mu_beta_prior <- matrix(0, nrow = p, ncol = 1)
  }
  
  if (is.null(Sigma_beta_prior)) {
    Sigma_beta_prior <- diag(100, p)
  }
  
  Sigma_beta_prior_inv <- solve(Sigma_beta_prior)
  
  # Useful quantities
  XtX <- t(X) %*% X
  XtY <- t(X) %*% y
  
  
  # ------------------------------------------------------------
  # Storage for samples
  # ------------------------------------------------------------
  
  beta_samples <- matrix(NA, nrow = n_iter, ncol = p)
  sigma2_samples <- numeric(n_iter)
  
  colnames(beta_samples) <- colnames(X)
  
  
  # ------------------------------------------------------------
  # Starting values
  # ------------------------------------------------------------
  
  beta_current <- Mu_beta_prior
  sigma2_current <- 1
  
  
  # ------------------------------------------------------------
  # Gibbs Sampler:
  #
  # beta^(t+1)    ~ p(beta | sigma²^(t), y, X)
  # sigma²^(t+1) ~ p(sigma² | beta^(t+1), y, X)
  # ------------------------------------------------------------
  
  for (iter in 1:n_iter) {
    
    
    # ------------------------------------------------------------
    # Sample from:
    #
    # beta | sigma², y, X ~ N(mu_beta_post, Sigma_beta_post)
    #
    # Sigma_beta_post =
    # (Sigma_beta_prior^(-1) + X'X / sigma²)^(-1)
    #
    # Mu_beta_post =
    # Sigma_beta_post *
    # (Sigma_beta_prior^(-1) Mu_beta_prior + X'y / sigma²)
    # ------------------------------------------------------------
    
    Sigma_beta_post <- solve(
      Sigma_beta_prior_inv + XtX / sigma2_current
    )
    
    Mu_beta_post <- Sigma_beta_post %*% (
      Sigma_beta_prior_inv %*% Mu_beta_prior +
        XtY / sigma2_current
    )
    
    beta_current <- Mu_beta_post +
      t(chol(Sigma_beta_post)) %*% rnorm(p)
    
    
    # ------------------------------------------------------------
    # Sample from:
    #
    # sigma² | beta, y, X ~ IG(A_post, B_post)
    #
    # A_post = A_prior + n/2
    #
    # B_post =
    # B_prior +
    # 0.5 * (y - X beta)'(y - X beta)
    # ------------------------------------------------------------
    
    residuals <- y - X %*% beta_current
    
    A_post <- A_prior + n / 2
    
    B_post <- B_prior +
      0.5 * as.numeric(t(residuals) %*% residuals)
    
    sigma2_current <- 1 / rgamma(
      1,
      shape = A_post,
      rate = B_post
    )
    
    
    # ------------------------------------------------------------
    # Store samples
    # ------------------------------------------------------------
    
    beta_samples[iter, ] <- beta_current
    sigma2_samples[iter] <- sigma2_current
    
    
    # ------------------------------------------------------------
    # Optional progress output
    # ------------------------------------------------------------
    
    if (verbose && iter %% 1000 == 0) {
      cat("Gibbs iteration:", iter, "\n")
    }
  }
  
  
  # ------------------------------------------------------------
  # Remove burn-in
  # ------------------------------------------------------------
  
  beta_samples <- beta_samples[(burn_in + 1):n_iter, , drop = FALSE]
  sigma2_samples <- sigma2_samples[(burn_in + 1):n_iter]
  
  
  # ------------------------------------------------------------
  # Return result
  # ------------------------------------------------------------
  
  return(list(
    beta_samples = beta_samples,
    sigma2_samples = sigma2_samples,
    n_iter = n_iter,
    burn_in = burn_in,
    n_samples_used = n_iter - burn_in
  ))
}