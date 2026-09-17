# Bayesian Portfolio Optimization

**Bayesian Analysis course project** — Master in Data Science for Economics

*Beatrice Clavarino*

Bayesian asset allocation across 4 ETFs (SPY, GLD, TLT, EEM) representing US equities, gold, long-term US bonds, and emerging market equities. Monthly returns are modeled as multivariate Normal with independent priors (Normal on µ, Inverse-Wishart on Σ), and inference is carried out via a Gibbs sampler in R.

An informative prior on µ, built from Research Affiliates' long-run return forecasts, is compared against a diffuse prior.

### Repository structure:

- **bayes_report.pdf**: The project report, written in LaTeX, provides a comprehensive overview of the methodology, experiments, and results.
- **bayes_code.R**: The R script containing the code used for the project, implementing the Gibbs sampler, MCMC diagnostics, and portfolio optimization.

## Key results

- Informative prior regularizes posterior means, especially for TLT (2.74% vs. −7.42% historical)
- Bayesian portfolio (informative prior) is more diversified: SPY 43.9%, GLD 41.2%, TLT 10.3%, EEM 4.6%
- MCMC diagnostics: well-mixed chains, ESS ≈ 8,000, Gelman-Rubin R̂ = 1.00

## Reproduce

```r
install.packages(c("quantmod", "MASS", "coda", "ggplot2", "quadprog"))
source("code/BayesianCode.R")
```
