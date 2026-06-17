# Variational Inference for Bayesian Computation: An Application to NBA Salary Data

Research seminar project for *Selected Topics in Statistical Learning and Data Science* at the Karlsruhe Institute of Technology (KIT), June 2026.

The project implements **mean-field variational Bayes (MFVB)** for Bayesian linear regression from scratch in R, using the **coordinate ascent variational inference (CAVI)** algorithm. It then compares MFVB with a **Gibbs sampler (MCMC)** and OLS, in two settings:

1. **A simulation study** that measures runtime as the problem size grows.
2. **An application to 7,298 NBA player-seasons (2010–2025)** that asks which box-score statistics are associated with player salaries.

📄 **[Slides](docs/presentation.pdf)** · 📝 **[Handout](docs/handout.pdf)**

---

## Main findings

- MFVB gives a deterministic alternative to MCMC that is much cheaper to compute.
- In the simulation study, CAVI is between about **9× and more than 19× faster** than Gibbs sampling.
- On the NBA data, MFVB and MCMC give **almost identical posterior means and 95% intervals**.

---

## Model

Bayesian linear regression with independent conjugate priors:

```math
\mathbf{y} = \mathbf{X}\boldsymbol\beta + \boldsymbol\varepsilon, \qquad
\boldsymbol\varepsilon \sim \mathcal N(\mathbf 0, \sigma^2 I_n), \qquad
\boldsymbol\beta \sim \mathcal N_p(\boldsymbol\mu_\beta, \boldsymbol\Sigma_\beta), \qquad
\sigma^2 \sim \mathrm{IG}(A, B).
```

### Mean-field approximation

The posterior is approximated by $q(\boldsymbol\beta, \sigma^2) = q(\boldsymbol\beta)\,q(\sigma^2)$, where

```math
q(\boldsymbol\beta) = \mathcal N_p\big(\boldsymbol\mu_{q(\beta)}, \boldsymbol\Sigma_{q(\beta)}\big), \qquad
q(\sigma^2) = \mathrm{IG}\big(A_{q(\sigma^2)}, B_{q(\sigma^2)}\big).
```

### CAVI updates

The following updates are repeated until the variational mean converges. Each update increases the evidence lower bound (ELBO).

```math
\begin{aligned}
A_{q(\sigma^2)} &= A + \tfrac{n}{2} \quad\text{(constant)}\\
B_{q(\sigma^2)} &= B + \tfrac12\Big[(\mathbf y - \mathbf X\boldsymbol\mu_{q(\beta)})^\top(\mathbf y - \mathbf X\boldsymbol\mu_{q(\beta)}) + \operatorname{tr}\big(\mathbf X^\top\mathbf X\,\boldsymbol\Sigma_{q(\beta)}\big)\Big]\\
\boldsymbol\Sigma_{q(\beta)} &= \Big(\tfrac{A_{q(\sigma^2)}}{B_{q(\sigma^2)}}\mathbf X^\top\mathbf X + \boldsymbol\Sigma_\beta^{-1}\Big)^{-1}\\
\boldsymbol\mu_{q(\beta)} &= \boldsymbol\Sigma_{q(\beta)}\Big(\tfrac{A_{q(\sigma^2)}}{B_{q(\sigma^2)}}\mathbf X^\top\mathbf y + \boldsymbol\Sigma_\beta^{-1}\boldsymbol\mu_\beta\Big)
\end{aligned}
```

The closed-form ELBO is computed after every iteration to monitor convergence. Its derivation is in the [handout](docs/handout.pdf).

The **Gibbs sampler** alternates between the two full conditionals, $\boldsymbol\beta \mid \sigma^2, \mathbf y$ (normal) and $\sigma^2 \mid \boldsymbol\beta, \mathbf y$ (inverse gamma). It serves as the benchmark for accuracy and runtime.

---

## Results

### Simulation study: runtime

The data are simulated with $p/n$ held fixed, and $p$ includes the intercept. CAVI is capped at 25 iterations. The Gibbs sampler runs 1,000 iterations.

|      n |     p | CAVI (s) | Gibbs (s) | Speed-up |
|-------:|------:|---------:|----------:|---------:|
|    100 |    10 |   < 0.01 |      0.19 |    > 19× |
|  1,000 |   100 |     0.14 |      1.70 |    12.1× |
|  5,000 |   500 |    16.05 |    151.15 |     9.4× |
| 10,000 | 1,000 |   128.53 |  1,414.59 |    11.0× |

Absolute times depend on the machine.

### NBA salary application

- **Response:** log salary, standardised.
- **Predictors:** Year, Age, G, GS, MP, FG, FGA, eFG%, FT, FT%, FTA, TRB, AST, STL, BLK, TOV, PF and PTS, all standardised. Together with the intercept this gives $p = 19$.
- **Priors:** $\boldsymbol\beta \sim \mathcal N(\mathbf 0, 100 I)$ and $\sigma^2 \sim \mathrm{IG}(1, 1)$.
- **Convergence:** CAVI converges after 5 iterations.
- **Error variance:** the posterior mean of $\sigma^2$ is 0.4788 under both methods.

Coefficients are transformed back to the original log-salary scale via $\beta_j^{\text{orig}} = \beta_j^{\text{scaled}} \, s_y / s_{x_j}$. They are then reported as percentage effects on salary, $100\,(\exp(\beta_j^{\text{orig}}\Delta x_j) - 1)$, holding all other predictors constant.

The table shows the predictors whose 95% posterior intervals exclude zero:

| Variable | Unit change          | MFVB (%) | MCMC (%) |
|----------|----------------------|---------:|---------:|
| Year     | +1 season            |     5.74 |     5.74 |
| Age      | +1 year              |     7.74 |     7.74 |
| G        | +1 game played       |     1.52 |     1.52 |
| GS       | +1 game started      |     0.22 |     0.22 |
| FGA      | +1 field goal attempt per game |     7.23 |     7.22 |
| eFG%     | +1 percentage point  |     0.39 |     0.39 |
| TRB      | +1 rebound per game  |     7.18 |     7.18 |
| AST      | +1 assist per game   |     3.22 |     3.20 |
| BLK      | +1 block per game    |    24.72 |    24.59 |
| PF       | +1 foul per game     |    -8.20 |    -8.14 |

The full coefficient plot comparing MFVB and MCMC is on slides 16–17 of the [presentation](docs/presentation.pdf).

---

## Repository structure

```
.
├── main.R                       # Runs the full NBA analysis (scripts 01–05)
├── R/
│   ├── 01_load_clean_data.R     # Load CSV, clean, log-transform and standardise
│   ├── 02_mfvb_function.R       # CAVI algorithm + closed-form ELBO
│   ├── 03_gibbs_sampler.R       # Gibbs sampler (MCMC benchmark)
│   ├── 04_run_models.R          # Set priors; fit OLS, CAVI and Gibbs
│   ├── 05_results_plots.R       # Back-transformation, tables, plots
│   └── 06_simulation_runtime.R  # Runtime study CAVI vs. Gibbs (optional, slow)
├── data/
│   └── nba_player_stats_salaries_2010-2025.csv
└── docs/
    ├── presentation.pdf         # Seminar slides
    └── handout.pdf              # Two-page handout with derivations
```

---

## How to run

**Requirements:** R ≥ 4.1, because the code uses the native pipe `|>`, plus three packages:

```r
install.packages(c("readr", "dplyr", "psych"))
```

**Run the NBA analysis.** Run this from the repository root, either by opening `variational-bayes-nba-salaries.Rproj` in RStudio and sourcing `main.R`, or from a terminal:

```bash
Rscript main.R
```

This takes a few seconds. The 10,000 Gibbs iterations account for most of that time. The script prints the comparison tables to the console and produces these plots:

- the coefficients with 95% intervals, for MFVB and Gibbs;
- the ELBO convergence;
- the posterior density of the PTS coefficient;
- actual vs. predicted salaries;
- the residuals.

When the script runs via `Rscript`, the plots are written to `Rplots.pdf`. The script also writes `figures/posterior_pts_density.dat` and `figures/posterior_pts_constants.tex`, which were used for the pgfplots figures in the slides.

**Reproduce the runtime study.** Uncomment the last line in `main.R`, or run:

```r
source("R/02_mfvb_function.R")
source("R/03_gibbs_sampler.R")
source("R/06_simulation_runtime.R")
```

⚠️ The largest scenario ($n = 10{,}000$, $p = 1{,}000$) takes about 25 minutes.

---

## Data

The data come from [**NBA Player Stats and Salaries 2010–2025**](https://www.kaggle.com/datasets/ratin21/nba-player-stats-and-salaries-2010-2025) by Ratin21 on Kaggle, released under **CC0: Public Domain**. The original sources are HoopsHype (salaries) and Basketball-Reference (statistics).

Preprocessing in `01_load_clean_data.R`:

- Keep player-seasons with a positive salary and take the log.
- Shooting percentages are missing when a player had zero attempts. These missing values are imputed with the column mean.
- Standardise the response and all predictors, then add an intercept column.

---

## Possible extensions

- Identify potentially underpaid players from the model residuals (the code is already sketched in `05_results_plots.R`) and check whether their next contracts adjust upward.
- Use richer variational families or stochastic variational inference.
- Add advanced basketball metrics to improve salary prediction.

---

## References

1. Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and Hall/CRC.
2. Ormerod, J. T., & Wand, M. P. (2010). Explaining variational approximations. *The American Statistician*, 64(2), 140–153.
3. Casella, G., & George, E. I. (1992). Explaining the Gibbs sampler. *The American Statistician*, 46(3), 167–174.
4. Ratin21 (n.d.). *NBA Player Stats and Salaries 2010–2025*. Kaggle.

---

## License

The code is released under the [MIT License](LICENSE). The dataset is CC0 (public domain); see [Data](#data).

**Author:** Adrian Wimmer
