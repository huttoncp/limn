## p11: four measurements taken before conformal recalibration was built in
## limn: within-group coverage (a), few groups (b) and forecasts (c, and d,
## added after c's result). Each part is run by name.
##
##   a. WITHIN-GROUP coverage under misspecification: a Poisson model fitted to
##      negative binomial counts, 30 groups of 10, seen-group rows. Conformal
##      calibration over all rows targets coverage over all rows; this
##      measures coverage by the group's true effect (terciles), and tries
##      conformal by bins ("Mondrian"), with the level chosen separately in
##      each tercile of the predicted mean.
##   b. FEW GROUPS, new-group rows, with a CORRECT Poisson model (P3's design:
##      groups of 6, group SD 0.8). ML, REML (P10's remedy), and
##      conformal on the ML fit from leave-one-group-out scores: one row per
##      group (group-level exchangeability; unbounded below 9 groups at 90%)
##      and all rows (more scores, but rows within a group are not exchangeable
##      with other groups' rows).
##   c. FORECASTING. P3's AR(1) design (20 series x 16 times, horizons 1-4),
##      Poisson AR(1) fitted to Poisson (correct) or NB (overdispersed) counts.
##      Conformal scores come from rolling origins inside the training data
##      (origins 8-12, each refitted), per horizon and pooled.
##   d. FORECASTING again, with three counts per series-time cell (see rep_d).
## Coverage is judged by the randomised PIT, with one uniform per test row
## shared by every method.
##
## A prototype, kept as it was run on 24 September 2026 with the illume of that
## day (the script did not record its version), 200 replicates per cell; the
## findings are in p11_conformal_followups.md. It reads the fit's joint
## precision and parameter vectors by name, as that version laid them out, so
## a later illume may need changes. Run from the package root:
##   Rscript studies/prototypes/p11_conformal_followups.R <a|b|c|d> [R] [cores]
args <- commandArgs(TRUE)
part <- if (length(args) >= 1) args[1] else "a"
R <- if (length(args) >= 2) as.integer(args[2]) else 100L
ncore <- if (length(args) >= 3) as.integer(args[3]) else 2L
res_dir <- file.path("studies", "prototypes", "results")

worker_setup <- function() {
  suppressPackageStartupMessages({ library(illume); library(Matrix) }); TRUE
}

## ---- shared pieces -------------------------------------------------------------
mu_joint <- function(f) {
  rn <- rownames(f$jointPrecision); fx <- f$opt$par; rd <- f$sdr$par.random
  mu <- numeric(length(rn))
  for (nm in unique(rn)) mu[rn == nm] <- if (nm %in% names(rd)) rd[names(rd) == nm] else fx[names(fx) == nm]
  mu
}
rmvn_prec <- function(Q, S) {
  e <- eigen(as.matrix(Q), symmetric = TRUE); keep <- e$values > 1e-8 * max(e$values)
  e$vectors[, keep, drop = FALSE] %*% (matrix(rnorm(sum(keep) * S), sum(keep), S) / sqrt(e$values[keep]))
}
gh <- local({ K <- 20; i <- seq_len(K - 1); b <- sqrt(i / 2)
  J <- matrix(0, K, K); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE); list(x = e$values * sqrt(2), w = e$vectors[1, ]^2) })
conf_q <- function(s, alpha = 0.1) {
  n <- length(s); k <- ceiling((n + 1) * (1 - alpha))
  if (k > n) Inf else sort(s)[k]
}
fit_glmm <- function(dd, fam, reml = FALSE) suppressWarnings(suppressMessages(
  ilm_model(y ~ x + (1 | g), data = dd, family = fam, joint = TRUE, reml = reml,
            verbose = FALSE)))

## randomised PITs (and plug-in conditional means) of the rows of `test` under
## the GLMM fit f: "seen" uses the row's group effect, "new" integrates a fresh
## one by Gauss-Hermite. `fam` is the FITTED family.
glmm_pits <- function(f, test, target, v, fam, S = 300L) {
  Q <- f$jointPrecision; rn <- rownames(Q); mu <- mu_joint(f)
  D <- rmvn_prec(Q, S) + mu
  ib <- which(rn == "beta"); it <- which(rn == "theta"); iu <- which(rn == "bvec")
  k_s <- if (fam == "nbinom") exp(D[rn == "logdisp", ]) else NULL
  lev <- levels(f$model$g)
  cdf <- function(q, eta, sd_re) {
    if (q < 0) return(0)
    if (is.null(sd_re)) {
      lam <- exp(eta)
      return(if (fam == "poisson") mean(ppois(q, lam)) else mean(pnbinom(q, size = k_s, mu = lam)))
    }
    lam <- exp(eta + outer(sd_re, gh$x))
    P <- if (fam == "poisson") ppois(q, lam) else pnbinom(q, size = k_s, mu = lam)
    mean(P %*% gh$w)
  }
  u <- numeric(nrow(test)); m <- numeric(nrow(test))
  for (i in seq_len(nrow(test))) {
    eta <- as.numeric(c(1, test$x[i]) %*% D[ib, , drop = FALSE])
    eta0 <- sum(c(1, test$x[i]) * mu[ib])
    if (target == "seen") {
      j <- iu[match(as.character(test$g[i]), lev)]
      eta <- eta + D[j, ]; eta0 <- eta0 + mu[j]; sd_re <- NULL
    } else sd_re <- exp(D[it, ])
    lo <- cdf(test$y[i] - 1, eta, sd_re); hi <- cdf(test$y[i], eta, sd_re)
    u[i] <- lo + v[i] * (hi - lo); m[i] <- exp(eta0)
  }
  list(u = u, mean = m)
}

## ---- a. within-group coverage ------------------------------------------------
rep_a <- function(rep, fam) {
  set.seed(50000 + rep)
  G <- 30; n <- 10; nt <- 5
  gen <- function(x, b) rnbinom(length(x), size = 1.5, mu = exp(log(3) + 0.5 * x + b))
  b <- rnorm(G, 0, 0.5)
  d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n)); d$y <- gen(d$x, b[d$g])
  ts <- data.frame(g = factor(rep(seq_len(G), each = nt), levels = levels(d$g)),
                   x = rnorm(G * nt)); ts$y <- gen(ts$x, b[as.integer(ts$g)])
  v <- runif(nrow(ts))
  f <- fit_glmm(d, fam)
  pt <- glmm_pits(f, ts, "seen", v, fam)
  kr <- sample(rep(1:5, length.out = nrow(d)))
  cal <- do.call(rbind, lapply(1:5, function(j) {
    fj <- fit_glmm(d[kr != j, ], fam); te <- d[kr == j, ]
    pj <- glmm_pits(fj, te, "seen", runif(nrow(te)), fam)
    data.frame(s = abs(pj$u - 0.5), m = pj$mean) }))
  q <- conf_q(cal$s)
  cuts <- stats::quantile(cal$m, c(1 / 3, 2 / 3))
  qb <- vapply(1:3, function(k) conf_q(cal$s[findInterval(cal$m, cuts) + 1L == k]), 0)
  bin <- findInterval(pt$mean, cuts) + 1L
  true_t <- findInterval(b[as.integer(ts$g)], stats::qnorm(c(1 / 3, 2 / 3), 0, 0.5)) + 1L
  data.frame(rep = rep, fam = fam, true_tercile = true_t, pred_bin = bin,
             model = abs(pt$u - 0.5) <= 0.45, conformal = abs(pt$u - 0.5) <= q,
             mondrian = abs(pt$u - 0.5) <= qb[bin])
}

## ---- b. few groups ---------------------------------------------------------------
rep_b <- function(rep, G) {
  set.seed(60000 + rep + 1000 * G)
  n <- 6; sb <- 0.8
  gen <- function(x, b) rpois(length(x), exp(log(3) + 0.5 * x + b))
  b <- rnorm(G, 0, sb)
  d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n)); d$y <- gen(d$x, b[d$g])
  tn <- data.frame(x = rnorm(5)); tn$y <- gen(tn$x, rnorm(5, 0, sb))
  v <- runif(5)
  f_ml <- fit_glmm(d, "poisson"); f_re <- fit_glmm(d, "poisson", reml = TRUE)
  u_ml <- glmm_pits(f_ml, tn, "new", v, "poisson")$u
  u_re <- glmm_pits(f_re, tn, "new", v, "poisson")$u
  one <- numeric(G); all <- list()
  for (g in seq_len(G)) {
    fg <- fit_glmm(droplevels(d[as.integer(d$g) != g, ]), "poisson")
    te <- d[as.integer(d$g) == g, ]
    s <- abs(glmm_pits(fg, te, "new", runif(nrow(te)), "poisson")$u - 0.5)
    one[g] <- s[sample.int(length(s), 1L)]; all[[g]] <- s
  }
  q1 <- conf_q(one); qa <- conf_q(unlist(all))
  data.frame(rep = rep, G = G, model_ml = abs(u_ml - 0.5) <= 0.45,
             model_reml = abs(u_re - 0.5) <= 0.45,
             conf_one = abs(u_ml - 0.5) <= q1, conf_all = abs(u_ml - 0.5) <= qa,
             unbounded_one = is.infinite(q1), unbounded_all = is.infinite(qa),
             level_one = min(2 * q1, 1), level_all = min(2 * qa, 1))
}

## ---- c. forecasting -----------------------------------------------------------------
fit_ar <- function(dd) suppressWarnings(suppressMessages(
  ilm_model(y ~ x + (1 | g), data = dd, family = "poisson", joint = TRUE,
            ar = ilm_ar1(dd$time, dd$g, verbose = FALSE), verbose = FALSE)))
## randomised PITs of rows `fut`, forecast from fit f whose last time is To
ar_pits <- function(f, fut, To, v, S = 300L) {
  Q <- f$jointPrecision; rn <- rownames(Q); mu <- mu_joint(f)
  D <- rmvn_prec(Q, S) + mu
  ib <- which(rn == "beta"); iu <- which(rn == "bvec"); ia <- which(rn == "B_ar")
  G <- nlevels(f$model$g); lev <- levels(f$model$g)
  last <- ia[(seq_len(G) - 1L) * f$ar$Tt + f$ar$Tt]
  rho_s <- tanh(D[rn == "rho_raw", ]); sa2_s <- exp(2 * D[rn == "lchol_ar", ])
  vapply(seq_len(nrow(fut)), function(i) {
    h <- fut$time[i] - To; gi <- match(as.character(fut$g[i]), lev)
    ph <- rho_s^h
    eta <- as.numeric(c(1, fut$x[i]) %*% D[ib, , drop = FALSE]) + D[iu[gi], ] + ph * D[last[gi], ]
    sd_in <- sqrt(sa2_s * (1 - ph^2))
    cdf <- function(q) if (q < 0) 0 else mean(ppois(q, exp(eta + outer(sd_in, gh$x))) %*% gh$w)
    lo <- cdf(fut$y[i] - 1); hi <- cdf(fut$y[i])
    lo + v[i] * (hi - lo)
  }, 0)
}
rep_c <- function(rep, fam) {
  set.seed(70000 + rep + if (fam == "nbinom") 5000 else 0)
  G <- 20; Tn <- 16; H <- 4; rho <- 0.7; sa <- 0.8; sb <- 0.5
  full <- expand.grid(time = 1:(Tn + H), g = factor(seq_len(G)))
  full$x <- rnorm(nrow(full))
  a <- unlist(lapply(seq_len(G), function(g) {
    z <- numeric(Tn + H); z[1] <- rnorm(1, 0, sa)
    for (t in 2:(Tn + H)) z[t] <- rho * z[t - 1] + rnorm(1, 0, sa * sqrt(1 - rho^2)); z }))
  mu <- exp(0.3 + 0.3 * full$x + rnorm(G, 0, sb)[full$g] + a)
  full$y <- if (fam == "poisson") rpois(nrow(full), mu) else rnbinom(nrow(full), size = 1.5, mu = mu)
  obs <- full[full$time <= Tn, ]; fut <- full[full$time > Tn, ]
  v <- runif(nrow(fut))
  u <- ar_pits(fit_ar(obs), fut, Tn, v)
  cal <- do.call(rbind, lapply(8:12, function(o) {
    tr <- obs[obs$time <= o, ]; te <- obs[obs$time > o & obs$time <= o + H, ]
    data.frame(h = te$time - o, s = abs(ar_pits(fit_ar(tr), te, o, runif(nrow(te))) - 0.5)) }))
  qh <- vapply(1:H, function(k) conf_q(cal$s[cal$h == k]), 0)
  qp <- conf_q(cal$s)
  hf <- fut$time - Tn
  data.frame(rep = rep, fam = fam, h = hf, model = abs(u - 0.5) <= 0.45,
             conf_h = abs(u - 0.5) <= qh[hf], conf_pooled = abs(u - 0.5) <= qp,
             level_h = pmin(2 * qh[hf], 1), level_pooled = min(2 * qp, 1))
}

## ---- d. forecasting, misspecified in a way the latent cannot absorb ------------
## Part c's NB counts turned out not to miscover: with one count per series-time
## cell, the AR latent soaks up the overdispersion. Here each cell has three
## counts sharing its latent value, so the overdispersion WITHIN a cell is left
## to the Poisson family, which cannot carry it.
rep_d <- function(rep, fam) {
  set.seed(80000 + rep + if (fam == "nbinom") 5000 else 0)
  G <- 20; Tn <- 16; H <- 4; rho <- 0.7; sa <- 0.8; sb <- 0.5; k <- 3L
  cells <- expand.grid(time = 1:(Tn + H), g = factor(seq_len(G)))
  a <- unlist(lapply(seq_len(G), function(g) {
    z <- numeric(Tn + H); z[1] <- rnorm(1, 0, sa)
    for (t in 2:(Tn + H)) z[t] <- rho * z[t - 1] + rnorm(1, 0, sa * sqrt(1 - rho^2)); z }))
  bg <- rnorm(G, 0, sb)
  full <- cells[rep(seq_len(nrow(cells)), each = k), ]
  full$x <- rnorm(nrow(full))
  mu <- exp(0.3 + 0.3 * full$x + bg[full$g] + rep(a, each = k))
  full$y <- if (fam == "poisson") rpois(nrow(full), mu) else rnbinom(nrow(full), size = 1.5, mu = mu)
  obs <- full[full$time <= Tn, ]; fut <- full[full$time > Tn, ]
  v <- runif(nrow(fut))
  u <- ar_pits(fit_ar(obs), fut, Tn, v)
  cal <- do.call(rbind, lapply(8:12, function(o) {
    tr <- obs[obs$time <= o, ]; te <- obs[obs$time > o & obs$time <= o + H, ]
    data.frame(h = te$time - o, s = abs(ar_pits(fit_ar(tr), te, o, runif(nrow(te))) - 0.5)) }))
  qh <- vapply(1:H, function(k) conf_q(cal$s[cal$h == k]), 0)
  qp <- conf_q(cal$s)
  hf <- fut$time - Tn
  data.frame(rep = rep, fam = fam, h = hf, model = abs(u - 0.5) <= 0.45,
             conf_h = abs(u - 0.5) <= qh[hf], conf_pooled = abs(u - 0.5) <= qp,
             level_h = pmin(2 * qh[hf], 1), level_pooled = min(2 * qp, 1))
}

## ---- run --------------------------------------------------------------------------
library(parallel)
cl <- makeCluster(ncore); invisible(clusterCall(cl, worker_setup))
clusterExport(cl, c("mu_joint", "rmvn_prec", "gh", "conf_q", "fit_glmm", "glmm_pits",
                    "rep_a", "rep_b", "fit_ar", "ar_pits", "rep_c", "rep_d"))
jobs <- switch(part,
  a = expand.grid(rep = seq_len(R), arg = c("poisson", "nbinom"), stringsAsFactors = FALSE),
  b = expand.grid(rep = seq_len(R), arg = c(6L, 8L, 10L, 15L, 20L)),
  c = expand.grid(rep = seq_len(R), arg = c("poisson", "nbinom"), stringsAsFactors = FALSE),
  d = expand.grid(rep = seq_len(R), arg = c("poisson", "nbinom"), stringsAsFactors = FALSE))
fn <- switch(part, a = "rep_a", b = "rep_b", c = "rep_c", d = "rep_d")
t0 <- proc.time()[3]
out <- parLapplyLB(cl, seq_len(nrow(jobs)), function(j, jobs, fn)
  tryCatch(get(fn)(jobs$rep[j], jobs$arg[j]), error = function(e) NULL), jobs = jobs, fn = fn)
stopCluster(cl)
ok <- !vapply(out, is.null, TRUE); z <- do.call(rbind, out[ok])
mins <- (proc.time()[3] - t0) / 60
cat(sprintf("part %s: %d of %d replicates usable, %.1f min\n", part, sum(ok), length(ok), mins))
f2 <- function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
if (part == "a") {
  ## replicate-level means first, then across replicates
  for (by in c("true_tercile", "pred_bin")) {
    z$bin <- z[[by]]
    agg <- aggregate(cbind(model, conformal, mondrian) ~ rep + fam + bin, z, mean)
    tab <- do.call(data.frame, aggregate(cbind(model, conformal, mondrian) ~ fam + bin, agg, f2))
    cat("\n== coverage of nominal 90% by", by, "(1 = lowest) ==\n"); print(format(tab, digits = 3), row.names = FALSE)
  }
  agg <- aggregate(cbind(model, conformal, mondrian) ~ rep + fam, z, mean)
  cat("\n== marginal ==\n"); print(format(do.call(data.frame, aggregate(cbind(model, conformal, mondrian) ~ fam, agg, f2)), digits = 3), row.names = FALSE)
} else if (part == "b") {
  agg <- aggregate(cbind(model_ml, model_reml, conf_one, conf_all) ~ rep + G, z, mean)
  tab <- do.call(data.frame, aggregate(cbind(model_ml, model_reml, conf_one, conf_all) ~ G, agg, f2))
  cat("\n== new-group coverage of nominal 90% ==\n"); print(format(tab, digits = 3), row.names = FALSE)
  lv <- aggregate(cbind(unbounded_one, unbounded_all, level_one, level_all) ~ G,
                  unique(z[c("rep", "G", "unbounded_one", "unbounded_all", "level_one", "level_all")]),
                  function(v) mean(v))
  cat("\n== conformal: share unbounded and mean level chosen ==\n"); print(format(lv, digits = 3), row.names = FALSE)
} else if (part %in% c("c", "d")) {
  agg <- aggregate(cbind(model, conf_h, conf_pooled) ~ rep + fam + h, z, mean)
  tab <- do.call(data.frame, aggregate(cbind(model, conf_h, conf_pooled) ~ fam + h, agg, f2))
  cat("\n== forecast coverage of nominal 90%, by horizon ==\n"); print(format(tab, digits = 3), row.names = FALSE)
  lv <- aggregate(cbind(level_h, level_pooled) ~ fam + h, z, mean)
  cat("\n== mean level chosen ==\n"); print(format(lv, digits = 3), row.names = FALSE)
}
saveRDS(list(rows = z, minutes = mins), file.path(res_dir, paste0("p11", part, "_conformal.rds")))
