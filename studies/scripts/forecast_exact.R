## forecast_exact: are ilm_forecast()'s forecasts the exact conditional
## distribution, and do they agree with an independent implementation?
##
## Gaussian panels of 8 series with a covariate, for AR(1), CAR(1) and a
## random walk, each plus observation noise. New times: past each series'
## last time (horizons 1 to 4), and between observed times for CAR(1) and the
## walk.
##
##   1. Exact multivariate-normal conditioning at illume's estimates: the
##      latent's covariance over the observed and new times, generalised least
##      squares for beta, and the conditional mean and variance of the linear
##      predictor at the new times (universal kriging). No limn or illume code
##      is used for it. Compared with uncertainty = "fixed", which holds the
##      latent's parameters at their estimates (it still draws the residual SD,
##      which the reference holds, so the reference should be narrower by
##      about that SD's share of the uncertainty).
##   2. glmmTMB, with the new rows in the data at weight zero: ar1() for
##      AR(1), ou() for CAR(1); it has no random walk. Its predictions at those
##      rows are its conditional modes, and its se.fit the generalised delta
##      method.
##      a. With its variance parameters mapped to illume's estimates and beta
##         free, it computes the same thing as the exact reference, which it
##         should reproduce to the optimiser's tolerance.
##      b. Fitted freely, its delta-method SEs are compared with limn's joint
##         draws by horizon, from the last observed time (h = 0) on. At h = 0
##         the two are the same quantity. Beyond it, limn draws each new time
##         from its exact conditional given each draw's parameters, while the
##         delta method linearises in them, so they may part.
##
## Run from the package root: Rscript studies/scripts/forecast_exact.R
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(glmmTMB)
})
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "forecast_exact")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
S <- 40000L

## one panel: G series of Tn times (or a subset of them), with latent a and noise
panel <- function(type, seed, G = 8L, Tn = 20L) {
  set.seed(seed)
  tt <- if (type == "car1") lapply(seq_len(G), function(g) sort(sample(1:30, 16)))
        else rep(list(1:Tn), G)
  d <- do.call(rbind, lapply(seq_len(G), function(g) {
    t <- tt[[g]]
    a <- switch(type,
      ar1 = as.numeric(stats::arima.sim(list(ar = 0.7), length(t), sd = 0.6)),
      car1 = { z <- numeric(length(t)); z[1] <- stats::rnorm(1, 0, 0.8)
               for (i in seq_along(t)[-1]) { r <- 0.6^(t[i] - t[i - 1])
                 z[i] <- r * z[i - 1] + stats::rnorm(1, 0, 0.8 * sqrt(1 - r^2)) }; z },
      rw1 = c(0, cumsum(stats::rnorm(length(t) - 1, 0, 0.35))))
    data.frame(id = sprintf("s%d", g), t = t, x = stats::rnorm(length(t)), a = a)
  }))
  d$id <- factor(d$id)
  d$y <- 1 + 0.5 * d$x + d$a + stats::rnorm(nrow(d), 0, 0.4)
  d
}

## the new rows: horizons 1 to 4 past each series' last time, and for CAR(1) and
## the walk two times between observed ones in every series
new_rows <- function(d, type) {
  last <- tapply(d$t, d$id, max)
  fut <- do.call(rbind, lapply(names(last), function(g)
    data.frame(id = g, t = last[[g]] + 1:4, kind = "forecast")))
  if (type != "ar1") {                              # CAR(1) in its gaps, the walk at half steps
    mid <- do.call(rbind, lapply(names(last), function(g) {
      t <- d$t[d$id == g]
      k <- if (type == "car1") which(diff(t) > 1)[1:2] else c(3L, 10L)
      data.frame(id = g, t = t[k] + 0.5, kind = "between") }))
    fut <- rbind(fut, mid[!is.na(mid$t), ])
  }
  set.seed(99); fut$x <- stats::rnorm(nrow(fut))
  fut
}

## exact conditional mean and variance of eta at new rows, at the given parameters
exact_kriging <- function(d, nd, type, s2, rho, q, sig_e) {
  X <- cbind(1, d$x); Xn <- cbind(1, nd$x)
  n <- nrow(d); m <- nrow(nd)
  id_all <- c(as.character(d$id), as.character(nd$id)); t_all <- c(d$t, nd$t)
  C <- matrix(0, n + m, n + m)
  for (g in unique(id_all)) {
    k <- which(id_all == g); tk <- t_all[k]
    C[k, k] <- switch(type,
      ar1 = , car1 = s2 * rho^abs(outer(tk, tk, "-")),
      rw1 = { t0 <- min(d$t[d$id == g]); q * outer(tk - t0, tk - t0, pmin) })
  }
  io <- seq_len(n); iu <- n + seq_len(m)
  V <- C[io, io] + diag(sig_e^2, n)
  Vi <- solve(V)
  A <- solve(t(X) %*% Vi %*% X)
  b <- A %*% t(X) %*% Vi %*% d$y
  K <- C[iu, io] %*% Vi
  mean <- drop(Xn %*% b + K %*% (d$y - X %*% b))
  R <- Xn - K %*% X
  var <- diag(C[iu, iu] - K %*% C[io, iu] + R %*% A %*% t(R))
  list(mean = mean, sd = sqrt(var))
}

rows <- list(); hz <- list()
for (type in c("ar1", "car1", "rw1")) {
  d <- panel(type, seed = match(type, c("ar1", "car1", "rw1")))
  nd <- new_rows(d, type)
  ar <- switch(type, ar1 = ilm_ar1(~ t | id), car1 = ilm_car1(~ t | id), rw1 = ilm_rw1(~ t | id))
  f <- suppressWarnings(illume::ilm_model(y ~ x, data = d, family = "gaussian", ar = ar,
                                          verbose = FALSE))
  vc <- illume::ilm_varcorr(f)
  ex <- exact_kriging(d, nd, type, s2 = as.numeric(vc$latent$Sigma), rho = vc$latent$rho,
                      q = as.numeric(vc$latent$var_per_time), sig_e = vc$dispersion$value[1])
  pf <- ilm_forecast(f, nd, target = "link", uncertainty = "fixed", nsim = S, seed = 11L)
  m_e <- rowMeans(pf$draws); s_e <- apply(pf$draws, 1L, stats::sd)
  for (k in unique(nd$kind)) {
    on <- nd$kind == k
    rows[[paste(type, k, "exact")]] <- data.frame(
      structure = type, rows = k, against = "exact MVN at the estimates (fixed)",
      n = sum(on), mean_diff_in_mc_se = max(abs(m_e[on] - ex$mean[on]) / (ex$sd[on] / sqrt(S))),
      sd_ratio_min = min(s_e[on] / ex$sd[on]), sd_ratio_max = max(s_e[on] / ex$sd[on]),
      mean_diff_in_sd = max(abs(m_e[on] - ex$mean[on]) / ex$sd[on]))
  }
  if (type != "rw1") {
    ## glmmTMB with the new rows at weight zero; the last observed time too, for h = 0
    last_t <- tapply(d$t, d$id, max)
    lastrow <- do.call(rbind, lapply(split(d, d$id), function(z) z[which.max(z$t), c("id", "t", "x")]))
    nw <- rbind(data.frame(nd[c("id", "t", "x")],
                           h = ifelse(nd$kind == "forecast", nd$t - last_t[as.character(nd$id)], NA)),
                data.frame(lastrow, h = 0))
    dd <- rbind(data.frame(d[c("id", "t", "x", "y")], w = 1),
                data.frame(id = nd$id, t = nd$t, x = nd$x, y = 0, w = 0))
    dd$id <- factor(dd$id)
    tm_fit <- function(start = NULL, map = NULL) {
      if (type == "ar1") {
        dd$tf <- factor(dd$t, levels = sort(unique(dd$t)))
        glmmTMB(y ~ x + ar1(tf + 0 | id), data = dd, weights = w, start = start, map = map)
      } else {
        dd$tn <- numFactor(dd$t)
        glmmTMB(y ~ x + ou(tn + 0 | id), data = dd, weights = w, start = start, map = map)
      }
    }
    tm_pred <- function(m) {
      z <- nw; z$id <- factor(z$id, levels = levels(dd$id)); z$w <- 1; z$y <- 0
      if (type == "ar1") z$tf <- factor(z$t, levels = sort(unique(dd$t))) else z$tn <- numFactor(z$t)
      predict(m, z, re.form = NULL, se.fit = TRUE)
    }
    ## a. variances mapped to illume's estimates
    s_lat <- sqrt(as.numeric(vc$latent$Sigma)); rho <- vc$latent$rho
    th2 <- if (type == "ar1") rho / sqrt(1 - rho^2) else log(-log(rho))
    m_map <- tm_fit(start = list(theta = c(log(s_lat), th2), betadisp = log(vc$dispersion$value[1])),
                    map = list(theta = factor(c(NA, NA)), betadisp = factor(NA)))
    pm <- tm_pred(m_map); k_new <- seq_len(nrow(nd))
    rows[[paste(type, "glmmTMB mapped")]] <- data.frame(
      structure = type, rows = paste(unique(nd$kind), collapse = " and "),
      against = "glmmTMB, variances mapped, against the exact reference", n = nrow(nd),
      mean_diff_in_mc_se = NA, sd_ratio_min = min(pm$se.fit[k_new] / ex$sd),
      sd_ratio_max = max(pm$se.fit[k_new] / ex$sd),
      mean_diff_in_sd = max(abs(pm$fit[k_new] - ex$mean) / ex$sd))
    ## b. fitted freely: the delta method against joint draws, by horizon
    po <- tm_pred(tm_fit())
    pj <- ilm_forecast(f, nw[c("id", "t", "x")], target = "link", nsim = S, seed = 12L)
    mj <- rowMeans(pj$draws); sj <- apply(pj$draws, 1L, stats::sd)
    for (h in sort(unique(stats::na.omit(nw$h)))) {
      on <- which(nw$h %in% h)
      hz[[paste(type, h)]] <- data.frame(structure = type, h = h, n = length(on),
        mean_diff_in_sd = max(abs(mj[on] - po$fit[on]) / po$se.fit[on]),
        sd_ratio_min = min(sj[on] / po$se.fit[on]), sd_ratio_max = max(sj[on] / po$se.fit[on]))
    }
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
hz <- do.call(rbind, hz); rownames(hz) <- NULL
print(res); print(hz)
write.csv(res, file.path(out_dir, "forecast_exact.csv"), row.names = FALSE)
write.csv(hz, file.path(out_dir, "forecast_delta_by_horizon.csv"), row.names = FALSE)

f2 <- function(v) sprintf("%.2f", v); f3 <- function(v) sprintf("%.3f", v)
md <- c(
  "# Forecasts against exact conditioning and glmmTMB", "",
  sprintf("## limn %s (%s; illume %s; glmmTMB %s)", version, format(Sys.Date()),
          utils::packageVersion("illume"), utils::packageVersion("glmmTMB")), "",
  "Generated by `studies/scripts/forecast_exact.R`. Do not edit by hand.", "",
  sprintf(paste("Gaussian panels of 8 series with a covariate and observation noise. New times: horizons 1 to 4",
                "past each series' last time, and for CAR(1) and the walk two times between observed ones. %d",
                "draws per forecast, so a draw-mean's Monte Carlo SE is its SD / %.0f, and an SD's relative",
                "error is about %.1f%%."), S, sqrt(S), 100 / sqrt(2 * S)), "",
  "**Exact.** limn with the latent's parameters held (`uncertainty = \"fixed\"`) against exact conditioning at illume's estimates, and glmmTMB with its variances mapped to the same estimates against that exact reference:", "",
  "| structure | new times | comparison | rows | largest mean difference | SD ratio |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %d | %s | %s to %s |", res$structure, res$rows, res$against, res$n,
          ifelse(is.na(res$mean_diff_in_mc_se), sprintf("%.1e SD", res$mean_diff_in_sd),
                 sprintf("%s MC SEs", f2(res$mean_diff_in_mc_se))),
          ifelse(is.na(res$mean_diff_in_mc_se), sprintf("%.6f", res$sd_ratio_min), f3(res$sd_ratio_min)),
          ifelse(is.na(res$mean_diff_in_mc_se), sprintf("%.6f", res$sd_ratio_max), f3(res$sd_ratio_max))), "",
  "**Joint draws against the delta method.** limn's joint draws against glmmTMB's own fit and se.fit, by horizon; h = 0 is each series' last observed time, where the two compute the same thing:", "",
  "| structure | horizon | rows | largest mean difference (SDs) | SD ratio (limn / glmmTMB) |",
  "|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %s | %s to %s |", hz$structure, hz$h, hz$n, f3(hz$mean_diff_in_sd),
          f3(hz$sd_ratio_min), f3(hz$sd_ratio_max)), "",
  sprintf("A mean's Monte Carlo SE is %.3f SDs.", 1 / sqrt(S)))
writeLines(md, file.path("studies", "findings", "forecast_exact.md"))
cat(md, sep = "\n")
