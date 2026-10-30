## p4: two numerical choices the engine has to make, measured without any fit.
##
## a. SCORING a predictive distribution that is a mixture over parameter and
##    random-effect draws: Poisson given lambda, log lambda ~ N(m, s^2).
##      RB      average the family's CDF / pmf over the M draws, then score
##              ("Rao-Blackwellised": no response is simulated)
##      sample  simulate one response per draw and score the sample as an
##              empirical distribution (the CRPS in its energy form)
##    Truth by 200-node Gauss-Hermite quadrature over log lambda.
## b. INTEGRATING a population-averaged mean E[plogis(eta + sigma Z)]:
##    Gauss-Hermite with K nodes against Monte Carlo with n draws, and the
##    same for a finite-difference slope of that mean.
##
## A prototype, kept as it was run on 24 September 2026 (base R; no model is
## fitted); the findings are in p4_scoring_integration.md. Run from the
## package root: Rscript studies/prototypes/p4_scoring_integration.R
set.seed(20260924)
res_dir <- file.path("studies", "prototypes", "results")

gh <- function(K) {                      # nodes/weights for E[f(Z)], Z ~ N(0,1)
  i <- seq_len(K - 1); b <- sqrt(i / 2)
  J <- matrix(0, K, K); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values * sqrt(2), w = e$vectors[1, ]^2)
}
G200 <- gh(200)

## ---- a. scoring -------------------------------------------------------------
crps_disc <- function(Fk, y) {           # CRPS of a count distribution, CDF on 0..K
  k <- seq_along(Fk) - 1L
  sum((Fk - (k >= y))^2)
}
truth <- function(m, s, y, K) {
  lam <- exp(m + s * G200$x)
  Fk <- colSums(G200$w * outer(lam, 0:K, function(l, k) ppois(k, l)))
  c(crps = crps_disc(Fk, y),
    logs = -log(sum(G200$w * dpois(y, lam))))
}
cfg <- expand.grid(m = c(0.5, 1.5), s = c(0.3, 0.8), q = c(0.5, 0.95))
Ms <- c(100L, 400L, 1600L); nrep <- 400L
rows <- list()
for (i in seq_len(nrow(cfg))) {
  m <- cfg$m[i]; s <- cfg$s[i]
  ## an observation at the q-th quantile of the predictive distribution
  K <- qpois(1 - 1e-12, exp(m + 6 * s))
  lamT <- exp(m + s * G200$x)
  FkT <- colSums(G200$w * outer(lamT, 0:K, function(l, k) ppois(k, l)))
  y <- which(FkT >= cfg$q[i])[1] - 1L
  tr <- truth(m, s, y, K)
  for (M in Ms) {
    est <- replicate(nrep, {
      lam <- exp(m + s * rnorm(M))
      ## the CDF is 1 to working precision beyond this, so the sum can stop
      Kr <- max(y, qpois(1 - 1e-12, max(lam)))
      Fk <- rowMeans(outer(0:Kr, lam, ppois))
      Y <- rpois(M, lam)
      ## energy form: E|Y - y| - E|Y - Y'| / 2, the pairwise mean by sorting
      Ys <- sort(Y); pw <- 2 * sum((2 * seq_len(M) - M - 1) * Ys) / M^2
      c(rb_crps = crps_disc(Fk, y), rb_logs = -log(mean(dpois(y, lam))),
        smp_crps = mean(abs(Y - y)) - pw / 2,
        smp_logs = { p <- mean(Y == y); if (p > 0) -log(p) else Inf })
    })
    rows[[length(rows) + 1L]] <- data.frame(
      m = m, s = s, q = cfg$q[i], y = y, M = M,
      true_crps = tr[["crps"]], true_logs = tr[["logs"]],
      rb_crps_rmse = sqrt(mean((est["rb_crps", ] - tr[["crps"]])^2)),
      smp_crps_rmse = sqrt(mean((est["smp_crps", ] - tr[["crps"]])^2)),
      smp_crps_bias = mean(est["smp_crps", ]) - tr[["crps"]],
      rb_logs_rmse = sqrt(mean((est["rb_logs", ] - tr[["logs"]])^2)),
      smp_logs_rmse = { z <- est["smp_logs", ]; sqrt(mean((z[is.finite(z)] - tr[["logs"]])^2)) },
      smp_logs_inf = mean(!is.finite(est["smp_logs", ])))
  }
}
sc <- do.call(rbind, rows)
sc$crps_rmse_ratio <- sc$smp_crps_rmse / sc$rb_crps_rmse
sc$logs_rmse_ratio <- sc$smp_logs_rmse / sc$rb_logs_rmse
cat("== a. scoring estimators: RMSE, sample-based over Rao-Blackwellised ==\n")
print(format(sc[, c("m", "s", "q", "y", "M", "true_crps", "rb_crps_rmse", "smp_crps_rmse",
                    "crps_rmse_ratio", "rb_logs_rmse", "smp_logs_rmse", "smp_logs_inf")],
             digits = 3), row.names = FALSE)
cat("\nmedian RMSE ratio, CRPS:", round(median(sc$crps_rmse_ratio), 2),
    " log score:", round(median(sc$logs_rmse_ratio), 2), "\n")

## ---- b. integration ----------------------------------------------------------
pa_exact <- function(eta, s) sum(G200$w * plogis(eta + s * G200$x))
rowsb <- list()
for (s in c(0.5, 1.5, 3)) for (eta in c(-2, 0, 1.5)) {
  ex <- pa_exact(eta, s)
  for (K in c(5L, 10L, 20L, 40L)) {
    g <- gh(K)
    rowsb[[length(rowsb) + 1L]] <- data.frame(sigma = s, eta = eta, method = sprintf("GH%02d", K),
      abs_err = abs(sum(g$w * plogis(eta + s * g$x)) - ex))
  }
  for (n in c(100L, 200L, 1000L)) {
    e <- replicate(500, mean(plogis(eta + s * rnorm(n))))
    rowsb[[length(rowsb) + 1L]] <- data.frame(sigma = s, eta = eta, method = sprintf("MC%04d", n),
      abs_err = sqrt(mean((e - ex)^2)))
  }
}
ib <- do.call(rbind, rowsb)
cat("\n== b. population-averaged mean: max abs error over eta (MC: RMSE over 500 seeds) ==\n")
print(round(xtabs(abs_err ~ method + sigma, aggregate(abs_err ~ method + sigma, ib, max)), 6))

## the slope: d/dx of E[plogis(b0 + b1 x + sigma Z)], at sigma = 1.5, b1 = 1
b0 <- -0.5; b1 <- 1; s <- 1.5; xs <- seq(-2, 2, length.out = 41); h <- 1e-4
slope_exact <- mean(vapply(xs, function(x) sum(G200$w * dlogis(b0 + b1 * x + s * G200$x)) * b1, 0))
g20 <- gh(20)
slope_gh20 <- mean(vapply(xs, function(x)
  (sum(g20$w * plogis(b0 + b1 * (x + h) + s * g20$x)) -
   sum(g20$w * plogis(b0 + b1 * (x - h) + s * g20$x))) / (2 * h), 0))
slope_mc <- function(n) replicate(500, {
  z <- rnorm(n)                                   # common random numbers across x +- h
  mean(vapply(xs, function(x)
    (mean(plogis(b0 + b1 * (x + h) + s * z)) - mean(plogis(b0 + b1 * (x - h) + s * z))) / (2 * h), 0))
})
m100 <- slope_mc(100L); m1000 <- slope_mc(1000L)
sl <- c(exact = slope_exact, gh20 = slope_gh20,
        mc100_rmse = sqrt(mean((m100 - slope_exact)^2)),
        mc100_worst_5pct = unname(quantile(abs(m100 - slope_exact), 0.95)),
        mc1000_rmse = sqrt(mean((m1000 - slope_exact)^2)))
cat("\n== b2. average marginal slope, sigma = 1.5 ==\n"); print(signif(sl, 4))
saveRDS(list(scoring = sc, integration = ib, slope = sl), file.path(res_dir, "p4_scoring_integration.rds"))
