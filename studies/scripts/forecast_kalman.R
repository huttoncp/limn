## forecast_kalman: with every parameter held (uncertainty = "latent"), are
## ilm_forecast()'s forecasts a Kalman filter's?
##
## A gaussian model with a latent AR(1), CAR(1) or random walk plus
## observation noise is a linear state-space model. Given its parameters, the
## forecast of each series is what the Kalman filter gives: filter the series
## to its last time, then carry the state forward (stats::KalmanRun() and
## stats::KalmanForecast(), from base R, sharing no code with limn or illume).
## uncertainty = "latent" holds every parameter at illume's estimates and
## draws only the latent cells, so it should reproduce that forecast.
##
## Panels of 4 series with a covariate: AR(1) and the walk at times 1 to 40,
## CAR(1) at 30 irregular integer times out of 1 to 50, which the filter sees as
## a grid with the unobserved times missing. The filter's state model:
##   AR(1), CAR(1): transition rho per time unit, innovation variance
##     Sigma (1 - rho^2), started from the stationary variance Sigma;
##   random walk: transition 1, innovation variance var_per_time, started at
##     zero at the series' first time (illume anchors each walk there).
## Its observation is y - X beta, with variance sigma^2.
##
##   1. The point forecast: limn's estimate (the fitted cells carried forward
##      with no new shocks) against the filter's mean, horizons 1 to 6.
##   2. The distribution: limn's draws of the linear predictor against the
##      filter's state mean and variance, and limn's predictive quantiles of a
##      new observation against the filter's normal forecast.
##   3. For contrast, the spread of the other settings against the filter's:
##      "none" leaves out the state's own uncertainty, "fixed" and "joint" add
##      the parameters'.
##
## Run from the package root: Rscript studies/scripts/forecast_kalman.R
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "forecast_kalman")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
S <- 40000L
H <- 6L
probs <- c(0.05, 0.5, 0.95)
S_other <- 10000L

panel <- function(type, seed, G = 4L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(G), function(g) {
    t <- if (type == "car1") sort(sample(1:50, 30)) else 1:40
    a <- switch(type,
      ar1 = as.numeric(stats::arima.sim(list(ar = 0.7), length(t), sd = 0.6)),
      car1 = { z <- numeric(length(t)); z[1] <- stats::rnorm(1, 0, 0.8)
               for (i in seq_along(t)[-1]) { r <- 0.8^(t[i] - t[i - 1])
                 z[i] <- r * z[i - 1] + stats::rnorm(1, 0, 0.8 * sqrt(1 - r^2)) }; z },
      rw1 = c(0, cumsum(stats::rnorm(length(t) - 1, 0, 0.35))))
    x <- stats::rnorm(length(t))
    data.frame(id = sprintf("s%d", g), t = t, x = x,
               y = 1 + 0.5 * x + a + stats::rnorm(length(t), 0, 0.4))
  }))
}

## the Kalman filter's forecast of one series, H steps past its last time
kalman <- function(ds, beta, sigma, lt) {
  grid <- seq(min(ds$t), max(ds$t))
  z <- rep(NA_real_, length(grid))
  z[match(ds$t, grid)] <- ds$y - beta[1] - beta[2] * ds$x
  if (lt$type == "rw1") { phi <- 1; Q <- c(lt$var_per_time); P0 <- 0 }
  else { phi <- lt$rho; Q <- c(lt$Sigma) * (1 - phi^2); P0 <- c(lt$Sigma) }
  mod <- list(T = matrix(phi), Z = 1, h = sigma^2, V = matrix(Q), a = 0,
              P = matrix(P0), Pn = matrix(P0))
  run <- stats::KalmanRun(z, mod, update = TRUE)
  fc <- stats::KalmanForecast(H, attr(run, "mod"))
  data.frame(id = ds$id[1], t = max(ds$t) + seq_len(H), h = seq_len(H),
             state_mean = fc$pred, state_var = fc$var - sigma^2, obs_var = fc$var)
}

one <- function(type) {
  d <- panel(type, seed = switch(type, ar1 = 11, car1 = 12, rw1 = 13))
  d$id <- factor(d$id)
  ar <- switch(type, ar1 = illume::ilm_ar1(~ t | id), car1 = illume::ilm_car1(~ t | id),
               rw1 = illume::ilm_rw1(~ t | id))
  fit <- suppressWarnings(illume::ilm_model(y ~ x, data = d, family = "gaussian", ar = ar,
                                            verbose = FALSE))
  vc <- illume::ilm_varcorr(fit)
  dr0 <- illume::ilm_draws(fit, nsim = 1L, seed = 1L, given = "parameters", natural = FALSE)
  beta <- dr0$mode[dr0$map$block == "beta"]
  ## the residual SD the fit's likelihood uses, exp(logdisp). For a gaussian fit with
  ## nothing integrated out, ilm_varcorr() reports lm's sqrt(RSS / (n - p)) instead;
  ## here the two should agree, and the findings say whether they do
  sigma <- exp(dr0$mode[dr0$map$block == "logdisp"][[1]])
  kf <- do.call(rbind, lapply(split(d, d$id), kalman, beta = beta, sigma = sigma, lt = vc$latent))
  set.seed(99)
  kf$x <- stats::rnorm(nrow(kf))
  kf$lp_mean <- beta[1] + beta[2] * kf$x + kf$state_mean
  nd <- kf[c("id", "t", "x")]

  fc <- function(target, unc, nsim) ilm_forecast(fit, nd, target = target, uncertainty = unc,
                                                 nsim = nsim, seed = 1L)
  lk <- fc("link", "latent", S)
  est <- lk$estimate
  dm <- rowMeans(lk$draws); dsd <- apply(lk$draws, 1L, stats::sd)
  rs <- fc("response", "latent", S_other)
  q <- vapply(probs, function(pr) lmn_pred_quantile(rs, pr), numeric(nrow(kf)))
  qk <- outer(kf$lp_mean, rep(1, length(probs))) +
        outer(sqrt(kf$obs_var), stats::qnorm(probs))
  h_corr <- vapply(split(seq_len(nrow(kf)), kf$id), function(i)
    stats::cor(lk$draws[i[1], ], lk$draws[i[H], ]), 0)
  ## the filter's correlation between horizons 1 and H: phi^(H-1) Var(1) / sd(1) sd(H)
  phi <- if (type == "rw1") 1 else vc$latent$rho
  k_corr <- phi^(H - 1) * sqrt(kf$state_var[kf$h == 1] / kf$state_var[kf$h == H])

  contrast <- do.call(rbind, lapply(c("none", "latent", "fixed", "joint"), function(u) {
    p <- if (u == "latent") lk else fc("link", u, S_other)
    r <- apply(p$draws, 1L, stats::sd) / sqrt(kf$state_var)
    data.frame(structure = type, uncertainty = u, h = kf$h, ratio = r)
  }))
  list(
    rows = data.frame(structure = type, id = kf$id, h = kf$h,
                      est_diff = est - kf$lp_mean,
                      mean_z = (dm - kf$lp_mean) / (dsd / sqrt(S)),
                      sd_ratio = dsd / sqrt(kf$state_var),
                      q_diff = apply(abs(q - qk), 1L, max) / sqrt(kf$obs_var)),
    corr = data.frame(structure = type, id = names(h_corr), limn = h_corr, kalman = k_corr),
    sigma = data.frame(structure = type, fit = sigma, varcorr = vc$dispersion$value[[1]]),
    contrast = contrast)
}

res <- lapply(c("ar1", "car1", "rw1"), one)
rows <- do.call(rbind, lapply(res, `[[`, "rows"))
corr <- do.call(rbind, lapply(res, `[[`, "corr"))
contrast <- do.call(rbind, lapply(res, `[[`, "contrast"))
sig <- do.call(rbind, lapply(res, `[[`, "sigma"))
write.csv(rows, file.path(out_dir, "rows.csv"), row.names = FALSE)
write.csv(corr, file.path(out_dir, "corr.csv"), row.names = FALSE)
write.csv(contrast, file.path(out_dir, "contrast.csv"), row.names = FALSE)

tab1 <- do.call(rbind, lapply(split(rows, rows$structure), function(s)
  data.frame(structure = s$structure[1], n = nrow(s), est = max(abs(s$est_diff)),
             z = max(abs(s$mean_z)), sd_lo = min(s$sd_ratio), sd_hi = max(s$sd_ratio),
             q = max(s$q_diff))))
tab2 <- do.call(rbind, lapply(split(corr, corr$structure), function(s)
  data.frame(structure = s$structure[1], diff = max(abs(s$limn - s$kalman)),
             lo = min(s$kalman), hi = max(s$kalman))))
tab3 <- stats::aggregate(ratio ~ structure + uncertainty + h, data = contrast, FUN = mean)
tab3 <- tab3[tab3$h %in% c(1L, H), ]
tab3$uncertainty <- factor(tab3$uncertainty, c("none", "latent", "fixed", "joint"))
tab3 <- tab3[order(tab3$structure, tab3$h, tab3$uncertainty), ]
wide <- stats::reshape(tab3, idvar = c("structure", "h"), timevar = "uncertainty",
                       direction = "wide")

md <- c(
  "# Forecasts with the parameters held, against the Kalman filter", "",
  sprintf("## limn %s (%s; illume %s; %s)", version, format(Sys.Date()),
          utils::packageVersion("illume"), R.version.string), "",
  "Generated by `studies/scripts/forecast_kalman.R`. Do not edit by hand.", "",
  paste("Gaussian panels of 4 series with a covariate and observation noise; AR(1) and the random walk",
        "at times 1 to 40, CAR(1) at 30 irregular times out of 1 to 50. Each series is forecast 1 to",
        sprintf("%d steps past its last time, and compared with `stats::KalmanRun()` and", H),
        "`stats::KalmanForecast()` at illume's estimates. limn used `uncertainty = \"latent\"`, which",
        sprintf("holds every parameter and draws only the latent cells, with %d draws: a draw-mean's", S),
        "Monte Carlo SE is its SD / 200, and an SD's relative error about 0.35%."), "",
  "**Point forecast and distribution.** The estimate is limn's point forecast, the fitted cells carried forward with no new shocks. The mean and SD are of the draws of the linear predictor, against the filter's state mean and variance. The quantiles are a new observation's at 5%, 50% and 95% (from 10,000 draws), against the filter's normal forecast:", "",
  "| structure | rows | largest estimate difference | largest mean difference (MC SEs) | SD ratio (limn / filter) | largest quantile difference (SDs) |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %d | %.1e | %.2f | %.4f to %.4f | %.4f |", tab1$structure, tab1$n, tab1$est,
          tab1$z, tab1$sd_lo, tab1$sd_hi, tab1$q), "",
  sprintf("**Joint over horizons.** The correlation of each series' draws at horizons 1 and %d, against the filter's (phi^%d Var(1)^1/2 / Var(%d)^1/2); its Monte Carlo SE is about 0.005 or less:", H, H - 1, H), "",
  "| structure | largest difference | filter's correlation |", "|---|---|---|",
  sprintf("| %s | %.4f | %.3f to %.3f |", tab2$structure, tab2$diff, tab2$lo, tab2$hi), "",
  "**What each setting adds.** The SD of the linear predictor's draws over the filter's, averaged over the series (10,000 draws for the settings other than \"latent\"):", "",
  "| structure | horizon | none | latent | fixed | joint |", "|---|---|---|---|---|---|",
  sprintf("| %s | %d | %.3f | %.3f | %.3f | %.3f |", wide$structure, wide$h, wide$ratio.none,
          wide$ratio.latent, wide$ratio.fixed, wide$ratio.joint), "",
  sprintf(paste("The filter uses the residual SD the fit's likelihood uses, exp(logdisp): %s.",
                "`ilm_varcorr()` reports %s. (For a gaussian model with nothing integrated out it",
                "reports lm's sqrt(RSS / (n - p)) instead, which is not the likelihood's.)"),
          paste(sprintf("%s %.4f", sig$structure, sig$fit), collapse = ", "),
          if (isTRUE(all.equal(sig$fit, sig$varcorr, tolerance = 1e-10))) "the same for these fits"
          else paste(sprintf("%s %.4f", sig$structure, sig$varcorr), collapse = ", ")), "",
  "\"none\" carries the fitted cells forward from their estimates, so it leaves out the uncertainty in the last cell; its ratio is below 1 at short horizons and approaches 1 as the new shocks dominate (for the random walk they never stop adding). \"fixed\" and \"joint\" add the uncertainty in the fixed effects and, for \"joint\", in the variances; the filter has neither.")
writeLines(md, file.path("studies", "findings", "forecast_kalman.md"))
cat(md, sep = "\n")
cat("\ncommit:", system("git rev-parse HEAD", intern = TRUE), "\n")
