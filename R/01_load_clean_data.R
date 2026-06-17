# ============================================================
# 01_load_clean_data.R
# Load, clean, preprocess, scale
# ============================================================

library(readr)
library(dplyr)

# ============================================================
# 1. Settings: predictors
# ============================================================

# Predictors used in the final model (p = 19 incl. intercept)
predictors <- c(
  "Year",
  "Age",
  "G",
  "GS",
  "MP",
  "FG",
  "FGA",   
  "eFG.",
  "FT",
  "FT.",
  "FTA",
  "TRB",
  "AST",
  "STL",
  "BLK",
  "TOV",
  "PF",
  "PTS"
)

# Extended predictor set (alternative specification, not used
# in the final model)
predictors_ad <- c(
  "Year",
  "Age",
  "G",
  "GS",
  "MP",
  "FG",
  "FGA",   
  "FG.",
  "X3P",
  "X3PA",
  "X3P.",
  "X2P",
  "X2PA",
  "X2P.",
  "eFG.",
  "FT",
  "FT.",
  "FTA",
  "ORB",
  "DRB",
  "TRB",
  "AST",
  "STL",
  "BLK",
  "TOV",
  "PF",
  "PTS"
)

# ============================================================
# 2. Load CSV
# ============================================================

nba <- read_csv(
  "data/nba_player_stats_salaries_2010-2025.csv",
  show_col_types = FALSE
)

names(nba) <- make.names(names(nba))

# ============================================================
# 3. Clean data
# ============================================================

model_data <- nba |>
  filter(
    !is.na(Salary),
    Salary > 0
  ) |>
  mutate(
    log_salary = log(Salary),
    # Shooting percentages are NA when a player had no attempts;
    # these are imputed with the column mean
    FG.  = ifelse(is.na(FG.),  mean(FG.,  na.rm = TRUE), FG.),
    X3P. = ifelse(is.na(X3P.), mean(X3P., na.rm = TRUE), X3P.),
    FT.  = ifelse(is.na(FT.),  mean(FT.,  na.rm = TRUE), FT.),
    eFG. = ifelse(is.na(eFG.), mean(eFG., na.rm = TRUE), eFG.),
    X2P. = ifelse(is.na(X2P.), mean(X2P., na.rm = TRUE), X2P.)
  ) |>
  select(
    Player,
    Salary,
    log_salary,
    all_of(predictors)
  ) |>
  na.omit()

# ============================================================
# 4. Build Y: log salary, raw and standardised
# ============================================================

Y_log_raw <- matrix(model_data$log_salary, ncol = 1)

Y_log_scaled <- scale(Y_log_raw)

y_log_mean <- attr(Y_log_scaled, "scaled:center")
y_log_sd   <- attr(Y_log_scaled, "scaled:scale")


# ============================================================
# 5. Build X: standardised, with intercept
# ============================================================

formula_no_intercept <- as.formula(
  paste("~", paste(predictors, collapse = " + "), "- 1")
)

X_raw_no_intercept <- model.matrix(
  formula_no_intercept,
  data = model_data
)

X_scaled_no_intercept <- scale(X_raw_no_intercept)

x_means <- attr(X_scaled_no_intercept, "scaled:center")
x_sds   <- attr(X_scaled_no_intercept, "scaled:scale")

X_scaled_intercept <- cbind(
  "(Intercept)" = 1,
  X_scaled_no_intercept
)

X_raw_intercept <- cbind(
  "(Intercept)" = 1,
  X_raw_no_intercept
)


# ============================================================
# 6. Default objects used by the models
# ============================================================

X_data <- X_scaled_intercept
Y_data <- Y_log_scaled

n <- nrow(X_data)
p <- ncol(X_data)

# ============================================================
# 7. Checks
# ============================================================

cat("\nData preparation finished.\n")

cat("\nDimensions:\n")
cat("model_data:", dim(model_data), "\n")
#cat("X_raw_no_intercept:", dim(X_raw_no_intercept), "\n")
cat("X_scaled_intercept:", dim(X_scaled_intercept), "\n")
cat("Y_log_raw:", dim(Y_log_raw), "\n")
cat("Y_log_scaled:", dim(Y_log_scaled), "\n")

cat("\nMissing values:\n")
cat("X_data NA:", sum(is.na(X_data)), "\n")
cat("X_data NaN:", sum(is.nan(X_data)), "\n")
cat("Y_data NA:", sum(is.na(Y_data)), "\n")
cat("Y_data NaN:", sum(is.nan(Y_data)), "\n")

cat("\nn:", n, "\n")
cat("p:", p, "\n")