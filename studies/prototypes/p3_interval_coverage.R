## p3: which sources of uncertainty must a prediction interval carry?
##
## Coverage of central 90% and 50% predictive intervals for NEW observations,
## judged by the randomised PIT of the realised value (exact for a count, no
## discreteness artefact), under four uncertainty models:
##   M0 plug-in        beta, theta fixed; a seen group at its fitted effect
##   M1 RE posterior   M0, but the seen group's effect drawn from its
##                     conditional posterior given the estimates
##   M2 + fixed eff.   beta and the group effects drawn jointly given theta,
##                     with theta and the dispersion held at their estimates
##   M3 joint          everything, theta included, from the joint precision
## Targets: a new row in a SEEN group, a row in a NEW group, and an AR(1)
## forecast h = 1..4 steps ahead. Few groups (8) is where they should differ;
## 40 groups is the control where they should not.
##
## M0 to M3 became ilm_predict()'s uncertainty = "none", "latent", "fixed" and
## "joint". Its "fixed" draws the dispersion, which M2 holds.
##
## A prototype, kept as it was run on 24 September 2026 with illume 0.0.7.9000
## (300 replicates per iid cell and 120 AR replicates); the findings are in
## p3_interval_coverage.md. It reads the fit's joint precision and parameter
## vectors by name, as that version laid them out, so a later illume may need
## changes. Run from the package root:
##   Rscript studies/prototypes/p3_interval_coverage.R [R_iid] [R_ar] [cores]
args <- commandArgs(TRUE)
R_iid <- if (length(args) >= 1) as.integer(args[1]) else 300L
R_ar  <- if (length(args) >= 2) as.integer(args[2]) else 150L
ncore <- if (length(args) >= 3) as.integer(args[3]) else 2L
res_dir <- file.path("studies", "prototypes", "results")

worker_setup <- function() {
  suppressPackageStartupMessages({ library(illume); library(Matrix) })
  TRUE
}

## ---- shared pieces (shipped to workers) -------------------------------------
mu_joint <- function(f) {
  rn <- rownames(f$jointPrecision); fx <- f$opt$par; rd <- f$sdr$par.random
  mu <- numeric(length(rn))
  for (nm in unique(rn))
    mu[rn == nm] <- if (nm %in% names(rd)) rd[names(rd) == nm] else fx[names(fx) == nm]
  mu
}
## N(0, Q^-1) draws, holding directions with no curvature (a boundary) fixed
rmvn_prec <- function(Q, S) {
  e <- eigen(as.matrix(Q), symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values)
  Z <- matrix(rnorm(sum(keep) * S), sum(keep), S)
  list(D = e$vectors[, keep, drop = FALSE] %*% (Z / sqrt(e$values[keep])),
       held = sum(!keep))
}
gh <- function(K) {
  i <- seq_len(K - 1); b <- sqrt(i / 2)
  J <- matrix(0, K, K); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE); list(x = e$values * sqrt(2), w = e$vectors[1, ]^2)
}
## randomised PIT of y under a mixture: eta_s (S draws), sd_s for gaussian,
## re_sd_s > 0 for a new group (integrated by quadrature)
pit <- function(y, fam, eta_s, sd_s = NULL, re_sd_s = NULL, G20) {
  if (fam == "gaussian") {
    tot <- if (is.null(re_sd_s)) sd_s else sqrt(sd_s^2 + re_sd_s^2)
    return(mean(pnorm(y, eta_s, tot)))
  }
  cdf <- function(q) {
    if (q < 0) return(0)
    if (is.null(re_sd_s)) return(mean(ppois(q, exp(eta_s))))
    mean(vapply(seq_along(eta_s), function(s)
      sum(G20$w * ppois(q, exp(eta_s[s] + re_sd_s[s] * G20$x))), 0))
  }
  lo <- cdf(y - 1); hi <- cdf(y)
  lo + runif(1) * (hi - lo)
}

rep_iid <- function(rep, fam, G, n, sb, S = 1000L) {
  set.seed(700000 + rep + (if (fam == "poisson") 50000 else 0) + G * 1000)
  b0 <- if (fam == "gaussian") 1 else log(3); b1 <- 0.5
  d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
  b <- rnorm(G, 0, sb)
  gen <- function(eta) if (fam == "gaussian") rnorm(length(eta), eta, 1) else rpois(length(eta), exp(eta))
  d$y <- gen(b0 + b1 * d$x + b[d$g])
  ts <- data.frame(g = factor(seq_len(G), levels = levels(d$g)), x = rnorm(G))
  ts$y <- gen(b0 + b1 * ts$x + b)
  nn <- 5L; tn <- data.frame(x = rnorm(nn)); tn$y <- gen(b0 + b1 * tn$x + rnorm(nn, 0, sb))
  f <- tryCatch(suppressWarnings(suppressMessages(
         ilm_model(y ~ x + (1 | g), data = d, family = fam, joint = TRUE, verbose = FALSE))),
       error = function(e) NULL)
  if (is.null(f) || is.null(f$jointPrecision)) return(NULL)
  Q <- f$jointPrecision; rn <- rownames(Q); mu <- mu_joint(f)
  ib <- which(rn == "beta"); it <- which(rn == "theta"); iu <- which(rn == "bvec")
  idp <- which(rn == "logdisp")
  G20 <- gh(20)
  ## M3: everything
  jd <- rmvn_prec(Q, S); D3 <- jd$D + mu
  ## M2: (beta, b) given theta -- the precision's (beta, b) block
  ic <- c(ib, iu); d2 <- rmvn_prec(Q[ic, ic], S)
  D2 <- matrix(mu, length(mu), S); D2[ic, ] <- D2[ic, ] + d2$D
  ## M1: b given (beta, theta)
  d1 <- rmvn_prec(Q[iu, iu], S)
  D1 <- matrix(mu, length(mu), S); D1[iu, ] <- D1[iu, ] + d1$D
  D0 <- matrix(mu, length(mu), 1L)
  sd_of <- function(D) if (length(idp)) exp(D[idp, ]) else NULL
  re_sd <- function(D) exp(D[it, ])
  out <- list()
  for (M in c("M0", "M1", "M2", "M3")) {
    D <- switch(M, M0 = D0, M1 = D1, M2 = D2, M3 = D3)
    u_seen <- vapply(seq_len(G), function(i) {
      eta <- as.numeric(c(1, ts$x[i]) %*% D[ib, , drop = FALSE]) + D[iu[i], ]
      pit(ts$y[i], fam, eta, sd_of(D), NULL, G20)
    }, 0)
    u_new <- vapply(seq_len(nn), function(i) {
      eta <- as.numeric(c(1, tn$x[i]) %*% D[ib, , drop = FALSE])
      pit(tn$y[i], fam, eta, sd_of(D), re_sd(D), G20)
    }, 0)
    out[[M]] <- data.frame(rep = rep, fam = fam, G = G, method = M,
                           target = c(rep("seen", G), rep("new", nn)),
                           u = c(u_seen, u_new))
  }
  r <- do.call(rbind, out)
  r$held <- jd$held; r$sb_hat <- sqrt(f$Sigma$g[1, 1])
  r
}

rep_ar <- function(rep, S = 1000L) {
  set.seed(900000 + rep)
  G <- 20; Tn <- 16; H <- 4; rho <- 0.7; sa <- 0.8; sb <- 0.5
  full <- expand.grid(time = 1:(Tn + H), g = factor(seq_len(G)))
  full$x <- rnorm(nrow(full))
  a <- unlist(lapply(seq_len(G), function(g) {
    z <- numeric(Tn + H); z[1] <- rnorm(1, 0, sa)
    for (t in 2:(Tn + H)) z[t] <- rho * z[t - 1] + rnorm(1, 0, sa * sqrt(1 - rho^2)); z }))
  full$y <- rpois(nrow(full), exp(0.3 + 0.3 * full$x + rnorm(G, 0, sb)[full$g] + a))
  obs <- full[full$time <= Tn, ]; fut <- full[full$time > Tn, ]
  f <- tryCatch(suppressWarnings(suppressMessages(
         ilm_model(y ~ x + (1 | g), data = obs, family = "poisson", joint = TRUE,
                   ar = ilm_ar1(obs$time, obs$g, verbose = FALSE), verbose = FALSE))),
       error = function(e) NULL)
  if (is.null(f) || is.null(f$jointPrecision)) return(NULL)
  Q <- f$jointPrecision; rn <- rownames(Q); mu <- mu_joint(f)
  ib <- which(rn == "beta"); iu <- which(rn == "bvec"); ia <- which(rn == "B_ar")
  last <- ia[(seq_len(G) - 1L) * f$ar$Tt + f$ar$Tt]
  jd <- rmvn_prec(Q, S); D3 <- jd$D + mu
  ir <- c(iu, ia); d1 <- rmvn_prec(Q[ir, ir], S)
  D1 <- matrix(mu, length(mu), S); D1[ir, ] <- D1[ir, ] + d1$D
  D0 <- matrix(mu, length(mu), 1L)
  h <- fut$time - Tn; gi <- as.integer(fut$g)
  out <- list()
  for (M in c("M0", "M1", "M3")) {
    D <- switch(M, M0 = D0, M1 = D1, M3 = D3)
    rho_s <- tanh(D[rn == "rho_raw", ]); sa2_s <- exp(2 * D[rn == "lchol_ar", ])
    u <- vapply(seq_len(nrow(fut)), function(i) {
      ph <- rho_s^h[i]
      ## the innovation part is integrated exactly by quadrature, so only
      ## parameter and state draws are simulated
      eta <- as.numeric(c(1, fut$x[i]) %*% D[ib, , drop = FALSE]) + D[iu[gi[i]], ] +
        ph * D[last[gi[i]], ]
      pit(fut$y[i], "poisson", eta, NULL, sqrt(sa2_s * (1 - ph^2)), gh(20))
    }, 0)
    out[[M]] <- data.frame(rep = rep, method = M, h = h, u = u)
  }
  r <- do.call(rbind, out); r$held <- jd$held; r$rho_hat <- f$rho; r
}

## ---- run --------------------------------------------------------------------
library(parallel)
cl <- makeCluster(ncore)
invisible(clusterCall(cl, worker_setup))
clusterExport(cl, c("mu_joint", "rmvn_prec", "gh", "pit", "rep_iid", "rep_ar"))
t0 <- proc.time()[3]
jobs <- expand.grid(rep = seq_len(R_iid), fam = c("gaussian", "poisson"), G = c(8L, 40L),
                    stringsAsFactors = FALSE)
iid <- parLapplyLB(cl, seq_len(nrow(jobs)), function(j, jobs)
  tryCatch(rep_iid(jobs$rep[j], jobs$fam[j], jobs$G[j], n = if (jobs$G[j] == 8L) 6L else 6L, sb = 0.8),
           error = function(e) NULL), jobs = jobs)
t1 <- proc.time()[3]
ar <- parLapplyLB(cl, seq_len(R_ar), function(r) tryCatch(rep_ar(r), error = function(e) NULL))
t2 <- proc.time()[3]
stopCluster(cl)
iid_ok <- vapply(iid, Negate(is.null), TRUE); ar_ok <- vapply(ar, Negate(is.null), TRUE)
iid <- do.call(rbind, iid[iid_ok]); ar <- do.call(rbind, ar[ar_ok])

cov_tab <- function(z, by) {
  z$c90 <- z$u > 0.05 & z$u < 0.95; z$c50 <- z$u > 0.25 & z$u < 0.75
  ## MC error from replicate-level means, since PITs within a replicate share a fit
  agg <- aggregate(cbind(c90, c50) ~ ., data = z[c(by, "rep", "c90", "c50")], FUN = mean)
  f <- function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
  o <- aggregate(cbind(c90, c50) ~ ., data = agg[c(setdiff(by, "rep"), "c90", "c50")],
                 FUN = f)
  do.call(data.frame, o)
}
ti <- cov_tab(iid, c("fam", "G", "target", "method"))
ta <- cov_tab(ar, c("h", "method"))
cat(sprintf("iid: %d of %d replicates usable (%.1f min); ar: %d of %d (%.1f min)\n",
            sum(iid_ok), length(iid_ok), (t1 - t0) / 60, sum(ar_ok), length(ar_ok), (t2 - t1) / 60))
cat("share of iid fits with a held (flat) direction:",
    round(mean(tapply(iid$held > 0, iid$rep + 1e4 * (iid$fam == "poisson") + 1e6 * iid$G, any)), 3), "\n")
cat("\n== new observations: coverage of central 90% and 50% intervals ==\n")
print(format(ti, digits = 3), row.names = FALSE)
cat("\n== Poisson AR(1) forecasts, 20 groups x 16 times ==\n")
print(format(ta, digits = 3), row.names = FALSE)
saveRDS(list(iid = iid, ar = ar, tab_iid = ti, tab_ar = ta,
             minutes = c(iid = (t1 - t0) / 60, ar = (t2 - t1) / 60)),
        file.path(res_dir, "p3_interval_coverage.rds"))
