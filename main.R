# ============================================================
# main.R
# Run full NBA Salary Analysis
#
# Run this script from the repository root, e.g. by opening
# variational-bayes-nba-salaries.Rproj in RStudio or via
#   Rscript main.R
# ============================================================

# Clean workspace
rm(list = ls())

# Optional: set seed for reproducibility
set.seed(123)

# ============================================================
# 1. Load and prepare data
# ============================================================

source("R/01_load_clean_data.R")


# ============================================================
# 2. Load functions
# ============================================================

source("R/02_mfvb_function.R")
source("R/03_gibbs_sampler.R")


# ============================================================
# 3. Run models
# ============================================================

source("R/04_run_models.R")


# ============================================================
# 4. Results and plots
# ============================================================

source("R/05_results_plots.R")


# ============================================================
# 5. Runtime simulation (optional)
# ============================================================
# Not run by default: the largest scenario (n = 10000, p = 1000)
# takes roughly 25 minutes. Uncomment to reproduce the runtime table.

# source("R/06_simulation_runtime.R")
