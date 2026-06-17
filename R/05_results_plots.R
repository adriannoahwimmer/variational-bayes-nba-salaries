# ============================================================
# 05_results_plots.R
# Results, plots, back-transformation, overpaid/underpaid
# ============================================================

library(dplyr)

# Requires:
# source("R/01_load_clean_data.R")
# source("R/02_mfvb_function.R")
# source("R/03_gibbs_sampler.R")
# source("R/04_run_models.R")


# ============================================================
# 1. Helper function: back-transform betas
# ============================================================
# The model was estimated on:
#
# Y_log_scaled = beta0_scaled + beta1_scaled * X1_scaled + ...
#
# Back-transformation to:
#
# log(Salary) = beta0_original + beta1_original * X1_raw + ...
#
# Slopes:
# beta_j_original = beta_j_scaled * y_sd / x_j_sd
#
# Intercept:
# beta0_original =
# y_mean + y_sd * beta0_scaled
# - sum(beta_j_original * x_j_mean)
# ============================================================

backtransform_betas <- function(beta_scaled, x_means, x_sds, y_mean, y_sd) {
  
  beta_original <- beta_scaled
  
  intercept_name <- "(Intercept)"
  slope_names <- setdiff(names(beta_scaled), intercept_name)
  
  beta_original[slope_names] <-
    beta_scaled[slope_names] * as.numeric(y_sd) / x_sds[slope_names]
  
  beta_original[intercept_name] <-
    as.numeric(y_mean) +
    as.numeric(y_sd) * beta_scaled[intercept_name] -
    sum(beta_original[slope_names] * x_means[slope_names])
  
  return(beta_original)
}


# ============================================================
# 2. Back-transform betas
# ============================================================

ols_beta_original <- backtransform_betas(
  beta_scaled = ols_beta_scaled,
  x_means = x_means,
  x_sds = x_sds,
  y_mean = y_log_mean,
  y_sd = y_log_sd
)

cavi_beta_original <- backtransform_betas(
  beta_scaled = cavi_beta_scaled,
  x_means = x_means,
  x_sds = x_sds,
  y_mean = y_log_mean,
  y_sd = y_log_sd
)

mcmc_beta_original <- backtransform_betas(
  beta_scaled = mcmc_beta_scaled,
  x_means = x_means,
  x_sds = x_sds,
  y_mean = y_log_mean,
  y_sd = y_log_sd
)


# ============================================================
# 3. Comparison table: scaled betas
# ============================================================

beta_comparison_scaled <- data.frame(
  variable = colnames(X_data),
  OLS_scaled = as.numeric(ols_beta_scaled),
  CAVI_scaled = as.numeric(cavi_beta_scaled),
  MCMC_scaled = as.numeric(mcmc_beta_scaled)
)

beta_comparison_scaled


# ============================================================
# 4. Comparison table: back-transformed betas
# ============================================================
# These betas are interpretable on the log(Salary) scale.
# Example:
# exp(beta) - 1 ≈ relative change in salary
# ============================================================

beta_comparison_original <- data.frame(
  variable = names(ols_beta_original),
  OLS_original = as.numeric(ols_beta_original),
  CAVI_original = as.numeric(cavi_beta_original),
  MCMC_original = as.numeric(mcmc_beta_original)
)

beta_comparison_original
print(cbind(
  beta_comparison_original["variable"],
  round(beta_comparison_original[, -1], 4)
))

# ============================================================
# 5. Percentage interpretation of the CAVI betas
# ============================================================
# Only meaningful for non-intercept variables.
# ============================================================

cavi_effects_percent <- data.frame(
  variable = names(cavi_beta_original),
  beta_log_salary = as.numeric(cavi_beta_original),
  percent_effect = (exp(as.numeric(cavi_beta_original)) - 1) * 100
) |>
  filter(variable != "(Intercept)") |>
  arrange(desc(abs(percent_effect)))

cavi_effects_percent

# ============================================================
# Percentage interpretation of the original-scale betas
# ============================================================
# Model on the original scale:
#
# log(Salary) = beta0 + beta1 * X1_raw + ...
#
# For regular variables:
# Effect of +1 unit:
# 100 * (exp(beta_j) - 1)
#
# For percentage variables such as eFG. and FT.,
# if they are stored as decimals, e.g. 0.55:
# Effect of +1 percentage point:
# 100 * (exp(0.01 * beta_j) - 1)
# ============================================================

percent_point_vars <- c("eFG.", "FT.")

# If eFG. and FT. are stored as 0.55, 0.80, etc.,
# then 0.01 corresponds to one percentage point.
# If they were stored as 55, 80, then 1 would be one percentage point.
decimal_percent_vars <- percent_point_vars[
  percent_point_vars %in% names(model_data) &
    sapply(model_data[percent_point_vars], function(x) max(x, na.rm = TRUE) <= 1.5)
]

effect_interpretation <- beta_comparison_original |>
  filter(variable != "(Intercept)") |>
  mutate(
    unit_change = ifelse(variable %in% decimal_percent_vars, 0.01, 1),
    
    interpretation = ifelse(
      variable %in% decimal_percent_vars,
      "+1 percentage point",
      "+1 raw unit"
    ),
    
    OLS_percent_effect =
      100 * (exp(OLS_original * unit_change) - 1),
    
    CAVI_percent_effect =
      100 * (exp(CAVI_original * unit_change) - 1),
    
    MCMC_percent_effect =
      100 * (exp(MCMC_original * unit_change) - 1)
  ) |>
  select(
    variable,
    interpretation,
    unit_change,
    OLS_percent_effect,
    CAVI_percent_effect,
    MCMC_percent_effect
  )

effect_interpretation

print(
  cbind(
    effect_interpretation[, c("variable", "interpretation", "unit_change")],
    round(effect_interpretation[, c(
      "OLS_percent_effect",
      "CAVI_percent_effect",
      "MCMC_percent_effect"
    )], 4)
  ),
  row.names = FALSE
)

# ============================================================
# 6. Table: scaled betas with 95% intervals
# CAVI: mean +/- 1.96 * sd from q(beta)
# Gibbs: empirical 2.5% and 97.5% posterior quantiles
# ============================================================

options(scipen = 999)

# CAVI posterior standard deviations
cavi_sd_scaled <- sqrt(diag(cavi_result$Sigma_q_beta_post_approx))
names(cavi_sd_scaled) <- colnames(X_data)

# CAVI 95% intervals
cavi_lower_scaled <- cavi_beta_scaled - 1.96 * cavi_sd_scaled
cavi_upper_scaled <- cavi_beta_scaled + 1.96 * cavi_sd_scaled

# Gibbs / MCMC 95% credible intervals from samples
mcmc_beta_quantiles_scaled <- apply(
  beta_samples,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

mcmc_lower_scaled <- mcmc_beta_quantiles_scaled[1, ]
mcmc_upper_scaled <- mcmc_beta_quantiles_scaled[2, ]

# Combined table
beta_intervals_scaled <- data.frame(
  variable = colnames(X_data),
  
  CAVI_mean = as.numeric(cavi_beta_scaled),
  CAVI_lower = as.numeric(cavi_lower_scaled),
  CAVI_upper = as.numeric(cavi_upper_scaled),
  
  Gibbs_mean = as.numeric(mcmc_beta_scaled),
  Gibbs_lower = as.numeric(mcmc_lower_scaled),
  Gibbs_upper = as.numeric(mcmc_upper_scaled)
)

beta_intervals_scaled_pretty <- beta_intervals_scaled

beta_intervals_scaled_pretty[, -1] <- lapply(
  beta_intervals_scaled_pretty[, -1],
  function(x) sprintf("%.4f", x)
)

print(beta_intervals_scaled_pretty, row.names = FALSE)



# ============================================================
# Plot: scaled coefficients with 95% intervals
# ============================================================


plot_data <- beta_intervals_scaled |>
  arrange(CAVI_mean)

y_pos <- seq_len(nrow(plot_data))

x_min <- min(plot_data$CAVI_lower, plot_data$Gibbs_lower)
x_max <- max(plot_data$CAVI_upper, plot_data$Gibbs_upper)

plot(
  plot_data$CAVI_mean,
  y_pos + 0.12,
  xlim = c(x_min, x_max),
  ylim = c(0.5, nrow(plot_data) + 0.5),
  yaxt = "n",
  xlab = "Scaled beta coefficient",
  ylab = "",
  main = "Scaled coefficients with 95% posterior intervals",
  pch = 19
)

axis(
  2,
  at = y_pos,
  labels = plot_data$variable,
  las = 2,
  cex.axis = 0.7
)

abline(v = 0, lty = 2)

# CAVI intervals
segments(
  x0 = plot_data$CAVI_lower,
  y0 = y_pos + 0.12,
  x1 = plot_data$CAVI_upper,
  y1 = y_pos + 0.12,
  lwd = 2
)

points(
  plot_data$CAVI_mean,
  y_pos + 0.12,
  pch = 19
)

# Gibbs intervals
segments(
  x0 = plot_data$Gibbs_lower,
  y0 = y_pos - 0.12,
  x1 = plot_data$Gibbs_upper,
  y1 = y_pos - 0.12,
  lwd = 2,
  lty = 2
)

points(
  plot_data$Gibbs_mean,
  y_pos - 0.12,
  pch = 17
)

legend(
  "bottomright",
  legend = c("CAVI", "Gibbs"),
  pch = c(19, 17),
  lty = c(1, 2),
  lwd = c(2, 2)
)


# ============================================================
# Table: original betas with 95% intervals
# Scale: log(Salary)
# CAVI: mean +/- 1.96 * sd after covariance back-transformation
# Gibbs: empirical 2.5% and 97.5% posterior quantiles after sample back-transformation
# ============================================================

options(scipen = 999)

# ------------------------------------------------------------
# Helper: Back-transform covariance matrix for CAVI
# ------------------------------------------------------------

make_backtransform_matrix <- function(beta_names, x_means, x_sds, y_sd) {
  
  p <- length(beta_names)
  
  T_mat <- matrix(
    0,
    nrow = p,
    ncol = p,
    dimnames = list(beta_names, beta_names)
  )
  
  intercept_name <- "(Intercept)"
  slope_names <- setdiff(beta_names, intercept_name)
  
  # Intercept transformation:
  # beta0_orig = y_mean + y_sd * beta0_scaled
  #              - sum_j beta_j_scaled * y_sd / x_sd_j * x_mean_j
  T_mat[intercept_name, intercept_name] <- as.numeric(y_sd)
  
  for (var in slope_names) {
    T_mat[var, var] <- as.numeric(y_sd) / x_sds[var]
    T_mat[intercept_name, var] <- -x_means[var] * as.numeric(y_sd) / x_sds[var]
  }
  
  return(T_mat)
}

# Transformation matrix
T_beta <- make_backtransform_matrix(
  beta_names = colnames(X_data),
  x_means = x_means,
  x_sds = x_sds,
  y_sd = y_log_sd
)

# CAVI covariance on original log(Salary) scale
cavi_cov_original <- T_beta %*%
  cavi_result$Sigma_q_beta_post_approx %*%
  t(T_beta)

cavi_sd_original <- sqrt(diag(cavi_cov_original))

# CAVI 95% intervals on original scale
cavi_lower_original <- cavi_beta_original - 1.96 * cavi_sd_original
cavi_upper_original <- cavi_beta_original + 1.96 * cavi_sd_original

# ------------------------------------------------------------
# Helper: Back-transform all Gibbs beta samples
# ------------------------------------------------------------

backtransform_beta_samples <- function(beta_samples, x_means, x_sds, y_mean, y_sd) {
  
  beta_samples_original <- beta_samples
  
  intercept_name <- "(Intercept)"
  slope_names <- setdiff(colnames(beta_samples), intercept_name)
  
  # Slopes
  for (var in slope_names) {
    beta_samples_original[, var] <-
      beta_samples[, var] * as.numeric(y_sd) / x_sds[var]
  }
  
  # Intercept
  beta_samples_original[, intercept_name] <-
    as.numeric(y_mean) +
    as.numeric(y_sd) * beta_samples[, intercept_name] -
    as.numeric(
      as.matrix(beta_samples_original[, slope_names, drop = FALSE]) %*%
        x_means[slope_names]
    )
  
  return(beta_samples_original)
}

mcmc_beta_samples_original <- backtransform_beta_samples(
  beta_samples = beta_samples,
  x_means = x_means,
  x_sds = x_sds,
  y_mean = y_log_mean,
  y_sd = y_log_sd
)

# Gibbs posterior means and intervals on original scale
mcmc_beta_original_from_samples <- colMeans(mcmc_beta_samples_original)

mcmc_beta_quantiles_original <- apply(
  mcmc_beta_samples_original,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

mcmc_lower_original <- mcmc_beta_quantiles_original[1, ]
mcmc_upper_original <- mcmc_beta_quantiles_original[2, ]



# ============================================================
# Combined table: original betas with 95% posterior intervals
# without intercept
# ============================================================

beta_intervals_original <- data.frame(
  variable = colnames(X_data),
  
  CAVI_mean = as.numeric(cavi_beta_original),
  CAVI_lower = as.numeric(cavi_lower_original),
  CAVI_upper = as.numeric(cavi_upper_original),
  
  Gibbs_mean = as.numeric(mcmc_beta_original_from_samples),
  Gibbs_lower = as.numeric(mcmc_lower_original),
  Gibbs_upper = as.numeric(mcmc_upper_original)
) |>
  filter(variable != "(Intercept)")

# Pretty print with 4 decimals
beta_intervals_original_pretty <- beta_intervals_original

beta_intervals_original_pretty[, -1] <- lapply(
  beta_intervals_original_pretty[, -1],
  function(x) sprintf("%.4f", x)
)

print(beta_intervals_original_pretty, row.names = FALSE)

# ============================================================
# Plot: original coefficients with 95% intervals
# Scale: log(Salary)
# ============================================================

plot_data <- beta_intervals_original |>
  arrange(CAVI_mean)

y_pos <- seq_len(nrow(plot_data))

x_min <- min(plot_data$CAVI_lower, plot_data$Gibbs_lower)
x_max <- max(plot_data$CAVI_upper, plot_data$Gibbs_upper)

plot(
  plot_data$CAVI_mean,
  y_pos + 0.12,
  xlim = c(x_min, x_max),
  ylim = c(0.5, nrow(plot_data) + 0.5),
  yaxt = "n",
  xlab = "Beta coefficient on log(Salary) scale",
  ylab = "",
  main = "Original-scale coefficients with 95% posterior intervals",
  pch = 19
)

axis(
  2,
  at = y_pos,
  labels = plot_data$variable,
  las = 2,
  cex.axis = 0.7
)

abline(v = 0, lty = 2)

# CAVI intervals
segments(
  x0 = plot_data$CAVI_lower,
  y0 = y_pos + 0.12,
  x1 = plot_data$CAVI_upper,
  y1 = y_pos + 0.12,
  lwd = 2
)

points(
  plot_data$CAVI_mean,
  y_pos + 0.12,
  pch = 19
)

# Gibbs intervals
segments(
  x0 = plot_data$Gibbs_lower,
  y0 = y_pos - 0.12,
  x1 = plot_data$Gibbs_upper,
  y1 = y_pos - 0.12,
  lwd = 2,
  lty = 2
)

points(
  plot_data$Gibbs_mean,
  y_pos - 0.12,
  pch = 17
)

legend(
  "bottomright",
  legend = c("CAVI", "Gibbs"),
  pch = c(19, 17),
  lty = c(1, 2),
  lwd = c(2, 2)
)





# ============================================================
# 8. Predictions with CAVI
# ============================================================
# We use the back-transformed CAVI betas:
#
# pred_log_salary = X_raw_intercept %*% beta_original
# pred_salary     = exp(pred_log_salary)
# ============================================================

model_data$pred_log_salary <- as.numeric(
  X_raw_intercept %*% cavi_beta_original
)

model_data$pred_salary <- exp(model_data$pred_log_salary)

model_data$residual_log <- model_data$log_salary - model_data$pred_log_salary

model_data$salary_ratio <- model_data$Salary / model_data$pred_salary


# ============================================================
# 9. Overpaid players
# ============================================================
# Large positive residual_log:
# actual salary > predicted salary
# ============================================================

# overpaid <- model_data |>
#   arrange(desc(residual_log)) |>
#   select(
#     Player,
#     Year,
#     Salary,
#     pred_salary,
#     salary_ratio,
#     residual_log,
#     Age,
#     G,
#     GS,
#     MP,
#     PTS,
#     AST,
#     STL,
#     BLK,
#     TOV,
#     PF
#   ) |>
#   head(20)
# 
# overpaid


# ============================================================
# 10. Underpaid players
# ============================================================
# Strongly negative residual_log:
# actual salary < predicted salary
# ============================================================

# underpaid <- model_data |>
#   arrange(residual_log) |>
#   select(
#     Player,
#     Year,
#     Salary,
#     pred_salary,
#     salary_ratio,
#     residual_log,
#     Age,
#     G,
#     GS,
#     MP,
#     PTS,
#     AST,
#     STL,
#     BLK,
#     TOV,
#     PF
#   ) |>
#   head(20)
# 
# underpaid


# ============================================================
# 11. Show runtime results
# ============================================================

runtime_results


# ============================================================
# 12. Show sigma² comparison
# ============================================================

sigma2_comparison


# ============================================================
# Print results
# ============================================================

cat("\nBeta comparison scaled:\n")
print(beta_comparison_scaled)

cat("\nBeta comparison original log salary scale:\n")
print(beta_comparison_original)

# cat("\nOverpaid players:\n")
# print(overpaid)
# 
# cat("\nUnderpaid players:\n")
# print(underpaid)

cat("\nRuntime results:\n")
print(runtime_results)

cat("\nSigma2 comparison:\n")
print(sigma2_comparison)

cat("\nFinal CAVI ELBO:\n")
print(cavi_result$elbo_final)


# ============================================================
# Plots
# ============================================================

# 1. ELBO convergence plot
plot(
  cavi_result$elbo_values,
  type = "b",
  main = "CAVI convergence: ELBO",
  xlab = "Iteration",
  ylab = "ELBO / lower bound",
  pch = 19
)



# 4. Posterior density for a single beta
var_name <- "PTS"

var_idx <- match(var_name, colnames(X_data))

if (is.na(var_idx)) {
  stop(paste("Variable", var_name, "not found in X_data."))
}

mcmc_draws <- beta_samples[, var_idx]

cavi_mean <- cavi_beta_scaled[var_idx]

cavi_sd <- sqrt(
  cavi_result$Sigma_q_beta_post_approx[var_idx, var_idx]
)

x_grid <- seq(
  min(mcmc_draws, na.rm = TRUE),
  max(mcmc_draws, na.rm = TRUE),
  length.out = 200
)

plot(
  density(mcmc_draws),
  main = paste("Posterior comparison for", var_name),
  xlab = "Beta",
  ylab = "Density"
)

lines(
  x_grid,
  dnorm(x_grid, mean = cavi_mean, sd = cavi_sd),
  lty = 2,
  lwd = 2
)

prior_beta_mean <- 0
mu_beta <- as.numeric(cavi_mean)

abline(v = prior_beta_mean, lty = 3, lwd = 2)
abline(v = mu_beta, lty = 4, lwd = 2)

points(prior_beta_mean, 0, pch = 16, cex = 1.2)
points(mu_beta, 0, pch = 17, cex = 1.2)

legend(
  "topright",
  legend = c("MCMC", "CAVI", "Prior mean = 0", expression(mu[beta])),
  lty = c(1, 2, 3, 4),
  lwd = c(1, 2, 2, 2),
  pch = c(NA, NA, 16, 17),
  bty = "n"
)
# -----------------------------
# Data for LaTeX / pgfplots
# -----------------------------

mcmc_density <- density(mcmc_draws, n = 400)

x_grid <- mcmc_density$x

plot_data <- data.frame(
  x = x_grid,
  MCMC = mcmc_density$y,
  CAVI = dnorm(x_grid, mean = cavi_mean, sd = cavi_sd)
)

prior_beta_mean <- 0
mu_beta <- as.numeric(cavi_mean)

y_max <- max(plot_data$MCMC, plot_data$CAVI) * 1.05

dir.create("figures", showWarnings = FALSE)

write.table(
  plot_data,
  file = "figures/posterior_pts_density.dat",
  row.names = FALSE,
  quote = FALSE
)

cat(
  paste0("\\newcommand{\\posteriorYmax}{", y_max, "}\n"),
  paste0("\\newcommand{\\priorBetaMean}{", prior_beta_mean, "}\n"),
  paste0("\\newcommand{\\muBetaPTS}{", mu_beta, "}\n"),
  file = "figures/posterior_pts_constants.tex"
)


# 5. Coefficient plot on the original scale
plot(
  cavi_beta_original[-1],
  xaxt = "n",
  main = "CAVI coefficients on log(Salary) scale",
  xlab = "Variable",
  ylab = "Beta",
  pch = 19
)

axis(
  1,
  at = seq_along(cavi_beta_original[-1]),
  labels = names(cavi_beta_original[-1]),
  las = 2,
  cex.axis = 0.7
)

abline(h = 0)

# 6. Actual vs predicted log salary
plot(
  model_data$log_salary,
  model_data$pred_log_salary,
  main = "Actual vs predicted log salary",
  xlab = "Actual log salary",
  ylab = "Predicted log salary",
  pch = 19
)

abline(0, 1)

# 7. Residuals
hist(
  model_data$residual_log,
  breaks = 40,
  main = "Distribution of salary residuals",
  xlab = "log(Salary) - predicted log(Salary)"
)

abline(v = 0)

