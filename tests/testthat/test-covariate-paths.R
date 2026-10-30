ar_series <- function(seed, Tn = 40, G = 1) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(sprintf("s%d", seq_len(G)), each = Tn)), t = rep(seq_len(Tn), G))
  d$x <- 1 + as.numeric(unlist(lapply(seq_len(G), function(g)
    stats::arima.sim(list(ar = 0.7), Tn, sd = 0.6)))) + stats::rnorm(Tn * G, 0, 0.3)
  d
}
fit_x <- function(d, f = x ~ 1)
  suppressMessages(suppressWarnings(illume::ilm_model(f, data = d, family = "gaussian",
                                     ar = illume::ilm_ar1(~ t | id), verbose = FALSE)))

test_that("with the parameters held, a covariate's paths are its Kalman filter forecast", {
  d <- ar_series(21); Tn <- 40
  fx <- fit_x(d)
  vc <- illume::ilm_varcorr(fx)
  dr <- illume::ilm_draws(fx, nsim = 1, given = "parameters", natural = FALSE)
  b0 <- dr$mode[dr$map$block == "beta"][[1]]; s <- exp(dr$mode[dr$map$block == "logdisp"][[1]])
  rho <- vc$latent$rho; v <- c(vc$latent$Sigma)
  mod <- list(T = matrix(rho), Z = 1, h = s^2, V = matrix(v * (1 - rho^2)), a = 0,
              P = matrix(v), Pn = matrix(v))
  kf <- stats::KalmanForecast(4, attr(stats::KalmanRun(d$x - b0, mod, update = TRUE), "mod"))
  fut <- data.frame(id = "s1", t = Tn + 1:4, x = NA)
  S <- 4000
  p <- ilm_covariate_paths(fx, fut, nsim = S, uncertainty = "latent")
  expect_equal(nrow(p), 4 * S)
  expect_equal(p$.draw, rep(seq_len(S), each = 4))
  expect_equal(p$t, rep(Tn + 1:4, S))
  X <- matrix(p$x, 4)                                     # rows: times; columns: paths
  sdk <- sqrt(kf$var)
  expect_true(all(abs(rowMeans(X) - (b0 + kf$pred)) < 4 * sdk / sqrt(S)))
  expect_equal(apply(X, 1L, stats::sd), sdk, tolerance = 0.05)
  ## joint over times: horizons 1 and 2 correlate as the filter says
  k12 <- rho * (kf$var[1] - s^2) / sqrt(kf$var[1] * kf$var[2])
  expect_equal(stats::cor(X[1, ], X[2, ]), k12, tolerance = 0.1)
})

test_that("paths go straight into a forecast, and several covariates are drawn independently", {
  d <- ar_series(22, Tn = 30, G = 2)
  d$z <- ave(stats::rnorm(60), d$id, FUN = cumsum)
  d$y <- 2 + 0.5 * d$x - 0.3 * d$z + stats::rnorm(60, 0, 0.3)
  fx <- fit_x(d); fz <- fit_x(d, z ~ 1)
  fy <- illume::ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  fut <- data.frame(id = rep(c("s1", "s2"), each = 3), t = rep(31:33, 2), x = NA, z = NA)
  p <- ilm_covariate_paths(list(fx, fz), fut, nsim = 400)
  expect_false(anyNA(p$x) || anyNA(p$z))
  X <- matrix(p$x, 6); Z <- matrix(p$z, 6)
  expect_lt(abs(stats::cor(X[1, ], Z[1, ])), 0.15)       # independent models, independent draws
  expect_no_message(fc <- ilm_forecast(fy, p, nsim = 400))
  expect_identical(fc$meta$covariate_draws, 400L)
  expect_setequal(fc$meta$drawn, c("x", "z"))
  expect_equal(nrow(as.data.frame(fc)), 6)
})

test_that("covariate paths refuse what they cannot draw", {
  d <- ar_series(23); d$w <- exp(d$x)
  fx <- fit_x(d)
  fut <- data.frame(id = "s1", t = 41:42, x = NA)
  expect_error(ilm_covariate_paths(fx, cbind(.draw = 1, fut)), "already has a `.draw` column")
  expect_error(ilm_covariate_paths(fit_x(d, log(w) ~ 1), fut), "not an expression [(]log[(]w[)][)]")
  expect_error(ilm_covariate_paths(list(fx, fx), fut), "same covariate, 'x'")
  expect_error(ilm_covariate_paths("x", fut), "fitted ilm_model")
})

test_that("paths from another tool become blocks of draws, row for row", {
  fut <- data.frame(id = rep(c("a", "b"), each = 2), t = rep(31:32, 2), x = NA, season = 1)
  sims <- expand.grid(t = 31:32, id = c("a", "b"), .rep = c(7, 3, 5))
  sims$.sim <- seq_len(nrow(sims)) / 10
  sims <- sims[sample(nrow(sims)), ]                        # any order
  p <- ilm_as_draws(sims, fut, by = c("id", "t"), values = c(x = ".sim"))
  expect_equal(nrow(p), 12)
  expect_equal(p$.draw, rep(1:3, each = 4))                  # paths numbered as they first appear
  ids <- unique(sims$.rep)
  for (k in 1:3) {
    blk <- p[p$.draw == k, ]
    ref <- sims[sims$.rep == ids[k], ]
    expect_equal(blk$x, ref$.sim[match(paste(fut$id, fut$t), paste(ref$id, ref$t))])
    expect_equal(blk$season, fut$season)
  }
  ## drawn columns found by name when they share one
  names(sims)[names(sims) == ".sim"] <- "x"
  expect_equal(ilm_as_draws(sims, fut, by = c("id", "t"))$x, p$x)
  ## every path must cover every row, once
  expect_error(ilm_as_draws(sims[-1, ], fut, by = c("id", "t")), "every path must cover every row")
  expect_error(ilm_as_draws(rbind(sims, sims[1, ]), fut, by = c("id", "t")), "more than one value")
  expect_error(ilm_as_draws(sims, fut, by = "id"), "identify each row")
  expect_error(ilm_as_draws(sims, fut, by = c("id", "t"), path = "draw"), "'draw' is missing")
})
