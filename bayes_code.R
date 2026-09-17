# ============================================
# Project: Bayesian Portfolio Optimization
# Course: Bayesian Analysis
# Author: Beatrice Clavarino
# ============================================

# ============================================
# LIBRARIES
# ============================================
# install.packages(c("quantmod", "MASS", "coda", "ggplot2"))
# install.packages("quadprog", repos = "https://cloud.r-project.org")
library(quantmod)
library(MASS)
library(coda)
library(quadprog)

# ============================================
# SETUP GENERALE
# ============================================
set.seed(123)
options(scipen = 999)

my_palette <- c("#674846", "#676767", "#DBD7D2", "#98817B",
                "#960018", "#BC8F8F", "#FFDCE2", "#58111A")

etf_colors <- c(
  SPY = "#960018",
  GLD = "#FFDCE2",
  TLT = "#676767",
  EEM = "#98817B"
)

etf_names <- c("SPY", "GLD", "TLT", "EEM")

# ============================================
# LOAD DATA
# ============================================
tickers <- c("SPY", "GLD", "TLT", "EEM")
getSymbols(tickers, from = "2010-01-01", to = "2024-12-31")

prices <- merge(
  monthlyReturn(Cl(SPY)),
  monthlyReturn(Cl(GLD)),
  monthlyReturn(Cl(TLT)),
  monthlyReturn(Cl(EEM))
)
colnames(prices) <- etf_names

# ============================================
# CLEANING AND PROCESSING
# ============================================
prices <- na.omit(prices)

# Converti in matrice e dividi i periodi
Y <- as.matrix(prices)
Y_prior <- Y[index(prices) < "2020-01-01", ]   # Periodo A: 2010-2019
Y_data  <- Y[index(prices) >= "2020-01-01", ]  # Periodo B: 2020-2024

cat("Osservazioni periodo A (prior):", nrow(Y_prior), "mesi\n")
cat("Osservazioni periodo B (dati):", nrow(Y_data), "mesi\n")

# ============================================
# STATISTICHE DESCRITTIVE
# ============================================

cat("\n=== Statistiche Periodo B (2020-2024) ===\n")
cat("\nRendimenti medi mensili:\n")
print(round(colMeans(Y_data), 4))
cat("\nRendimenti medi annualizzati:\n")
print(round(colMeans(Y_data) * 12, 4))
cat("\nDeviazione standard mensile:\n")
print(round(apply(Y_data, 2, sd), 4))
cat("\nMatrice di correlazione:\n")
print(round(cor(Y_data), 3))
cat("\nMatrice di covarianza:\n")
print(round(cov(Y_data), 6))

cat("\n=== Statistiche Periodo A (2010-2019) ===\n")
cat("\nRendimenti medi annualizzati:\n")
print(round(colMeans(Y_prior) * 12, 4))
cat("\nMatrice di correlazione:\n")
print(round(cor(Y_prior), 3))

# ============================================
# DATI PER LA STIMA
# ============================================
Y     <- as.matrix(Y_data)
T_obs <- nrow(Y)
p     <- ncol(Y)

# ============================================
# PRIOR INFORMATIVA
# ============================================

# Fonte: Research Affiliates Asset Allocation Interactive
# Expected 10Y Nominal Returns, as of 03/31/2026
# https://interactive.researchaffiliates.com/asset-allocation
mu_annual <- c(
  SPY = 0.035,   # US Large
  GLD = 0.050,   # Commodities (proxy for Gold)
  TLT = 0.049,   # US Treasury Long
  EEM = 0.073    # Emerging Markets
)
mu0_inf     <- (1 + mu_annual)^(1/12) - 1
Lambda0_inf <- diag(p) * (0.002)^2       # sd = 0.2% mensile
nu0_inf     <- p + 2 + nrow(Y_prior)     # = 126
S0_inf      <- cov(Y_prior) * (nu0_inf - p - 1)

cat("\nPrior mu mensile (informativa):\n")
print(round(mu0_inf, 5))
cat("Gradi di liberta prior informativa (nu0):", nu0_inf, "\n")

# ============================================
# PRIOR DIFFUSA
# ============================================
mu0_dif     <- rep(0, p)
Lambda0_dif <- diag(p) * (0.10)^2   # sd = 10% mensile
nu0_dif     <- p + 2                 # = 6, minimo possibile
S0_dif      <- diag(p) * 0.001

cat("Gradi di liberta prior diffusa (nu0):", nu0_dif, "\n")

# ============================================
# GIBBS SAMPLER
# ============================================
gibbs_sampler <- function(Y, mu0, Lambda0, nu0, S0, n_iter = 10000) {

  T_obs <- nrow(Y)
  p     <- ncol(Y)

  mu_samples    <- matrix(NA, nrow = n_iter, ncol = p)
  Sigma_samples <- array(NA, dim = c(p, p, n_iter))

  Sigma_current <- S0 / (nu0 - p - 1)
  y_bar         <- colMeans(Y)
  Lambda0_inv   <- solve(Lambda0)
  nu_n          <- nu0 + T_obs

  cat("Avvio Gibbs sampler...\n")

  for (s in 1:n_iter) {

    # STEP 1: campiona mu | Sigma, Y
    # Full conditional: N(mu_n, Lambda_n)
    Sigma_inv    <- solve(Sigma_current)
    Lambda_n_inv <- Lambda0_inv + T_obs * Sigma_inv
    Lambda_n     <- solve(Lambda_n_inv)
    mu_n         <- Lambda_n %*% (Lambda0_inv %*% mu0 + T_obs * Sigma_inv %*% y_bar)
    mu_current   <- mvrnorm(1, mu = mu_n, Sigma = Lambda_n)

    # STEP 2: campiona Sigma | mu, Y
    # Full conditional: IW(nu_n, S_n)
    S_mu <- matrix(0, p, p)
    for (t in 1:T_obs) {
      e_t  <- Y[t, ] - mu_current
      S_mu <- S_mu + outer(e_t, e_t)
    }
    S_n           <- S0 + S_mu
    W             <- rWishart(1, df = nu_n, Sigma = solve(S_n))[,,1]
    Sigma_current <- solve(W)

    mu_samples[s, ]      <- mu_current
    Sigma_samples[, , s] <- Sigma_current

    if (s %% 2000 == 0) cat("Iterazione", s, "di", n_iter, "\n")
  }

  cat("Gibbs sampler completato.\n")
  return(list(mu = mu_samples, Sigma = Sigma_samples))
}

# ============================================
# LANCIA IL GIBBS SAMPLER
# ============================================
cat("\n=== PRIOR INFORMATIVA ===\n")
results_inf <- gibbs_sampler(Y, mu0_inf, Lambda0_inf, nu0_inf, S0_inf, n_iter = 10000)

cat("\n=== PRIOR DIFFUSA ===\n")
results_dif <- gibbs_sampler(Y, mu0_dif, Lambda0_dif, nu0_dif, S0_dif, n_iter = 10000)

# ============================================
# BURN-IN
# ============================================
burnin <- 2000

mu_inf    <- results_inf$mu[(burnin+1):10000, ]
mu_dif    <- results_dif$mu[(burnin+1):10000, ]
Sigma_inf <- results_inf$Sigma[,,(burnin+1):10000]
Sigma_dif <- results_dif$Sigma[,,(burnin+1):10000]

colnames(mu_inf) <- etf_names
colnames(mu_dif) <- etf_names

cat("Campioni utili dopo burn-in:", nrow(mu_inf), "\n")

# ============================================
# CODA DIAGNOSTICS
# ============================================
mu_inf_mcmc <- mcmc(mu_inf)
mu_dif_mcmc <- mcmc(mu_dif)

cat("\n=== Effective Sample Size - Prior informativa ===\n")
print(round(effectiveSize(mu_inf_mcmc)))

cat("\n=== Effective Sample Size - Prior diffusa ===\n")
print(round(effectiveSize(mu_dif_mcmc)))

# Gelman-Rubin: seconda catena
set.seed(456)
results_inf_2 <- gibbs_sampler(Y, mu0_inf, Lambda0_inf, nu0_inf, S0_inf, n_iter = 10000)
mu_inf_2 <- results_inf_2$mu[(burnin+1):10000, ]
colnames(mu_inf_2) <- etf_names

chain_list <- mcmc.list(mcmc(mu_inf), mcmc(mu_inf_2))

cat("\n=== Gelman-Rubin - Prior informativa ===\n")
print(gelman.diag(chain_list))

# ============================================
# SOMMARIO POSTERIOR
# ============================================
cat("\n=== POSTERIOR - Prior informativa ===\n")
cat("Media posteriore annualizzata:\n")
print(round(colMeans(mu_inf) * 12, 4))
cat("Intervalli di credibilita 95%:\n")
print(round(apply(mu_inf, 2, quantile, c(0.025, 0.975)) * 12, 4))

cat("\n=== POSTERIOR - Prior diffusa ===\n")
cat("Media posteriore annualizzata:\n")
print(round(colMeans(mu_dif) * 12, 4))
cat("Intervalli di credibilita 95%:\n")
print(round(apply(mu_dif, 2, quantile, c(0.025, 0.975)) * 12, 4))

cat("\n=== DATI STORICI (confronto) ===\n")
print(round(colMeans(Y) * 12, 4))

# ============================================
# PESI OTTIMALI DEL PORTAFOGLIO
# ============================================
portfolio_weights <- function(mu, Sigma, gamma = 3) {
  p <- length(mu)
  
  # Programmazione quadratica — soluzione esatta
  # Massimizza: w'mu - (gamma/2) w'Sigma w
  # Equivale a minimizzare: (gamma/2) w'Sigma w - w'mu
  
  # Matrice Dmat e vettore dvec per solve.QP
  Dmat <- gamma * Sigma          # matrice quadratica
  dvec <- mu                     # vettore lineare
  
  # Vincoli:
  # 1. sum(w) = 1  → vincolo di uguaglianza
  # 2. w_i >= 0    → vincoli di non-negatività
  
  # Matrice dei vincoli (p+1 vincoli: 1 uguaglianza + p disuguaglianze)
  Amat <- cbind(
    rep(1, p),   # vincolo sum(w) = 1
    diag(p)      # vincoli w_i >= 0
  )
  bvec <- c(1, rep(0, p))  # RHS dei vincoli
  
  # meq = 1 indica che il primo vincolo è di uguaglianza
  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) NULL
  )
  
  # Se solve.QP fallisce, usa l'euristica come fallback
  if (is.null(result)) {
    Sigma_inv <- solve(Sigma)
    ones      <- rep(1, p)
    A         <- as.numeric(t(ones) %*% Sigma_inv %*% ones)
    B         <- as.numeric(t(ones) %*% Sigma_inv %*% mu)
    w         <- (1/gamma) * Sigma_inv %*% mu +
                 ((1 - B/gamma) / A) * (Sigma_inv %*% ones)
    w         <- pmax(as.vector(w), 0)
    w         <- w / sum(w)
    return(w)
  }
  
  return(result$solution)
}

n_samples   <- nrow(mu_inf)
weights_inf <- matrix(NA, nrow = n_samples, ncol = 4)
weights_dif <- matrix(NA, nrow = n_samples, ncol = 4)

for (s in 1:n_samples) {
  weights_inf[s, ] <- portfolio_weights(mu_inf[s, ], Sigma_inf[,, s])
  weights_dif[s, ] <- portfolio_weights(mu_dif[s, ], Sigma_dif[,, s])
}

colnames(weights_inf) <- etf_names
colnames(weights_dif) <- etf_names

cat("\n=== Pesi ottimali - Prior informativa ===\n")
print(round(colMeans(weights_inf), 3))
cat("\n=== Pesi ottimali - Prior diffusa ===\n")
print(round(colMeans(weights_dif), 3))
cat("\n=== Markowitz classico (confronto) ===\n")
w_markowitz <- portfolio_weights(colMeans(Y), cov(Y))
names(w_markowitz) <- etf_names
print(round(w_markowitz, 3))

# ============================================
# DISTRIBUZIONE PREDITTIVA
# ============================================
n_pred  <- nrow(mu_inf)
pred_inf <- matrix(NA, nrow = n_pred, ncol = 4)
pred_dif <- matrix(NA, nrow = n_pred, ncol = 4)

for (s in 1:n_pred) {
  pred_inf[s, ] <- mvrnorm(1, mu = mu_inf[s, ], Sigma = Sigma_inf[,, s])
  pred_dif[s, ] <- mvrnorm(1, mu = mu_dif[s, ], Sigma = Sigma_dif[,, s])
}

colnames(pred_inf) <- etf_names
colnames(pred_dif) <- etf_names

cat("\n=== Intervalli predittivi 95% - Prior informativa ===\n")
print(round(apply(pred_inf, 2, quantile, c(0.025, 0.5, 0.975)) * 12, 4))
cat("\n=== Intervalli predittivi 95% - Prior diffusa ===\n")
print(round(apply(pred_dif, 2, quantile, c(0.025, 0.5, 0.975)) * 12, 4))

# ============================================
# GRAFICI
# ============================================

my_palette <- c("#674846", "#676767", "#DBD7D2", "#98817B",
                "#960018", "#BC8F8F", "#FFDCE2", "#58111A")

etf_colors <- c(
  SPY = "#960018",
  GLD = "#FFDCE2",
  TLT = "#676767",
  EEM = "#98817B"
)

# --- Grafico 1: Traceplots ---
pdf("plot1_traceplots.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot(mu_inf[, i],
       type = "l",
       col  = etf_colors[i],
       main = paste("Traceplot -", etf_names[i]),
       xlab = "Iterazione",
       ylab = "mu",
       lwd  = 0.8)
}
dev.off()

# --- Grafico 2: Posterior per ETF ---
for (i in 1:4) {
  d_inf <- density(mu_inf[, i] * 12)
  d_dif <- density(mu_dif[, i] * 12)
  xlim  <- range(c(d_inf$x, d_dif$x))
  ylim  <- range(c(d_inf$y, d_dif$y))

  pdf(paste0("plot2_posterior_", etf_names[i], ".pdf"), width = 8, height = 6)
  plot(d_dif,
       col  = "#676767", lwd = 2,
       main = paste("Posterior -", etf_names[i]),
       xlab = "Rendimento annuale",
       ylab = "Densita",
       xlim = xlim, ylim = ylim)
  lines(d_inf, col = etf_colors[i], lwd = 2)
  abline(v = mu_annual[i],        col = "#960018", lty = 2, lwd = 1.5)
  abline(v = colMeans(Y)[i] * 12, col = "black",   lty = 3, lwd = 1.5)
  legend("topright",
         legend = c("Post. diffusa", "Post. informativa",
                    "Prior analisti", "Dato storico"),
         col    = c("#676767", etf_colors[i], "#960018", "black"),
         lty    = c(1, 1, 2, 3), lwd = 2, cex = 0.8, bty = "n")
  dev.off()
}

# --- Grafico 3: Pesi portafoglio ---
pdf("plot3_portfolio_weights.pdf", width = 9, height = 6)
weights_matrix <- rbind(
  Markowitz           = w_markowitz,
  "Post. diffusa"     = colMeans(weights_dif),
  "Post. informativa" = colMeans(weights_inf)
)
barplot(weights_matrix,
        beside      = TRUE,
        col         = c("#676767", "#BC8F8F", "#960018"),
        names.arg   = etf_names,
        ylab        = "Peso ottimale",
        ylim        = c(0, 0.7),
        legend.text = rownames(weights_matrix),
        args.legend = list(x = "topright", bty = "n", cex = 0.8))
title("Confronto pesi portafoglio ottimale")
dev.off()

# --- Grafico 4: Intervalli di credibilita sui pesi ---
pdf("plot4_weights_credible_intervals.pdf", width = 9, height = 6)
par(mfrow = c(1, 2))
for (label in c("informativa", "diffusa")) {
  w      <- if (label == "informativa") weights_inf else weights_dif
  means  <- colMeans(w)
  lower  <- apply(w, 2, quantile, 0.025)
  upper  <- apply(w, 2, quantile, 0.975)
  plot(1:4, means,
       ylim = c(0, 1), xlim = c(0.5, 4.5),
       pch  = 19, col = unname(etf_colors),
       xaxt = "n", xlab = "",
       ylab = "Peso ottimale",
       main = paste("Prior", label))
  axis(1, at = 1:4, labels = etf_names)
  for (i in 1:4) {
    lines(c(i, i), c(lower[i], upper[i]),
          col = unname(etf_colors)[i], lwd = 2)
  }
  abline(h = 0.25, lty = 2, col = "#676767")
}
dev.off()

# --- Grafico 5: Distribuzione predittiva ---
pdf("plot5_predictive_distribution.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (i in 1:4) {
  d_inf <- density(pred_inf[, i] * 12)
  d_dif <- density(pred_dif[, i] * 12)
  xlim  <- range(c(d_inf$x, d_dif$x))
  ylim  <- range(c(0, max(d_inf$y, d_dif$y)))
  plot(d_dif,
       col  = "#DBD7D2", lwd = 2,
       main = paste("Predictive -", etf_names[i]),
       xlab = "Rendimento annuale",
       ylab = "Densita",
       xlim = xlim, ylim = ylim)
  lines(d_inf, col = etf_colors[i], lwd = 2)
  legend("topright",
         legend = c("Diffusa", "Informativa"),
         col    = c("#DBD7D2", etf_colors[i]),
         lty = 1, lwd = 2, cex = 0.8, bty = "n")
}
dev.off()

cat("\nTutti i grafici salvati.\n")
getwd()
