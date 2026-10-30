test_that("ilm_forecast() is ilm_predict() with a forecast's settings", {
  set.seed(1)
  d <- data.frame(id = factor(rep(c("a", "b"), each = 20)), t = rep(1:20, 2))
  d$y <- 5 + ave(rnorm(40, 0, 0.4), d$id, FUN = cumsum) + rnorm(40, 0, 0.4)
  fit <- suppressWarnings(illume::ilm_model(y ~ 1, data = d, family = "gaussian",
                                            ar = illume::ilm_rw1(~ t | id), verbose = FALSE))
  fut <- ilm_future(d, time = "t", group = "id", h = 3, response = "y")
  fc <- ilm_forecast(fit, fut, nsim = 50, seed = 2)
  expect_equal(fc, ilm_predict(fit, fut, target = "response", groups = "fitted",
                               nsim = 50, seed = 2))
  expect_identical(fc$meta$target, "response")
  ## the forecasts widen with the horizon, as a walk's should
  w <- with(as.data.frame(ilm_forecast(fit, fut, nsim = 2000)), upper_0.9 - lower_0.9)
  expect_true(all(diff(w[1:3]) > 0) && all(diff(w[4:6]) > 0))
})

test_that("a seen series' point forecast carries its last cell forward without shocks", {
  set.seed(5); G <- 4; Tn <- 15
  d <- expand.grid(t = 1:Tn, id = factor(letters[1:G]))
  d$y <- 2 + unlist(lapply(1:G, function(g) as.numeric(arima.sim(list(ar = 0.6), Tn, sd = 0.5)))) +
    rnorm(nrow(d), 0, 0.3)
  fit <- suppressWarnings(illume::ilm_model(y ~ 1, data = d, family = "gaussian",
                                            ar = illume::ilm_ar1(~ t | id), verbose = FALSE))
  nd <- data.frame(id = c("a", "a", "c"), t = c(16, 18, 17))
  p <- ilm_forecast(fit, nd, target = "link", nsim = 200)
  rho <- illume::ilm_varcorr(fit)$latent$rho
  rf <- illume::ilm_ranef(fit); rf <- rf[rf$type == "ar1", ]
  last <- function(g) rf$mode[as.character(rf$level) == g & rf$time == Tn]
  expect_equal(p$estimate, coef(fit)[[1]] + c(rho * last("a"), rho^3 * last("a"), rho^2 * last("c")))
  ## it does not depend on the draws
  expect_equal(ilm_forecast(fit, nd, target = "link", nsim = 50, seed = 9)$estimate, p$estimate)
})

test_that("forecasting a trend in time warns, and names a random walk", {
  set.seed(6); G <- 5; Tn <- 12
  d <- expand.grid(t = 1:Tn, id = factor(letters[1:G]))
  d$y <- 1 + 0.1 * d$t + sin(2 * pi * d$t / 6) + rnorm(nrow(d), 0, 0.3)
  fit_ar <- function(f, ar) suppressMessages(suppressWarnings(illume::ilm_model(
    f, data = d, family = "gaussian", ar = ar, verbose = FALSE)))
  ahead <- data.frame(id = "a", t = Tn + 1:2)
  ## every parameter held: with no walk or AR in these data, their variances and
  ## the residual SD sit at the boundary, which draws cannot use until illume
  ## holds it (with the variances alone held, some draws' SD overflows)
  f1 <- fit_ar(y ~ t + ilm_fourier(t, 6), illume::ilm_ar1(~ t | id))
  expect_warning(ilm_forecast(f1, ahead, uncertainty = "latent", nsim = 20), "trends in `t` [(]t[)].*ilm_rw1[(]~ t [|] id[)]")
  ## a periodic term alone is safe to forecast
  f2 <- fit_ar(y ~ ilm_fourier(t, 6), illume::ilm_ar1(~ t | id))
  expect_no_warning(ilm_forecast(f2, ahead, uncertainty = "latent", nsim = 20))
  ## a random walk already lets the level move
  f3 <- fit_ar(y ~ t, illume::ilm_rw1(~ t | id))
  expect_no_warning(ilm_forecast(f3, ahead, uncertainty = "latent", nsim = 20))
  ## times inside the data are not forecasts
  expect_no_warning(ilm_forecast(f1, data.frame(id = "a", t = 3), uncertainty = "latent", nsim = 20))
  ## with no correlation over time, ilm_future()'s rows say which column is time
  f4 <- suppressMessages(suppressWarnings(illume::ilm_model(
    y ~ poly(t, 2) + (1 | id), data = rbind(d, d), family = "gaussian", verbose = FALSE)))
  fut <- ilm_future(d, time = "t", group = "id", h = 2, response = "y")
  expect_warning(ilm_forecast(f4, fut, groups = "fitted", uncertainty = "latent", nsim = 20),
                 "trends in `t` [(]poly[(]t, 2[)][)]")
})

test_that("with every parameter held, a gaussian AR(1) forecast is the Kalman filter's", {
  ## A well-identified case, so the fit is interior. The earlier 40-point case was
  ## weakly identified: its likelihood basins were under 0.15 log-likelihood units
  ## apart, its ML fit under illume 0.0.8.9001 was a local optimum, and REML is now
  ## illume's default. Weak identification of such fits is an open question for the
  ## package family.
  set.seed(8); Tn <- 300
  d <- data.frame(id = factor("a"), t = 1:Tn)
  d$y <- 1 + as.numeric(arima.sim(list(ar = 0.9), Tn, sd = 0.3)) + rnorm(Tn, 0, 0.6)
  fit <- suppressWarnings(illume::ilm_model(y ~ 1, data = d, family = "gaussian",
                                            ar = illume::ilm_ar1(~ t | id), verbose = FALSE))
  vc <- illume::ilm_varcorr(fit)
  dr <- illume::ilm_draws(fit, nsim = 1, given = "parameters", natural = FALSE)
  b0 <- coef(fit)[[1]]; s <- exp(dr$mode[dr$map$block == "logdisp"][[1]]); rho <- vc$latent$rho; v <- c(vc$latent$Sigma)
  mod <- list(T = matrix(rho), Z = 1, h = s^2, V = matrix(v * (1 - rho^2)), a = 0,
              P = matrix(v), Pn = matrix(v))
  kf <- stats::KalmanForecast(4, attr(stats::KalmanRun(d$y - b0, mod, update = TRUE), "mod"))
  nd <- data.frame(id = "a", t = Tn + 1:4)
  p <- ilm_forecast(fit, nd, target = "link", uncertainty = "latent", nsim = 4000)
  expect_equal(p$estimate, b0 + kf$pred, tolerance = 1e-6)
  expect_equal(unname(apply(p$draws, 1L, stats::sd)), sqrt(kf$var - s^2), tolerance = 0.05)
  ## "none" leaves out the last cell's own uncertainty, so it is narrower
  p0 <- ilm_forecast(fit, nd, target = "link", uncertainty = "none", nsim = 4000)
  expect_true(stats::sd(p0$draws[1, ]) < 0.9 * stats::sd(p$draws[1, ]))
})

test_that("uncertainty = \"latent\" holds every parameter, the dispersion included", {
  set.seed(9)
  d <- data.frame(x = rnorm(60)); d$y <- 1 + 0.5 * d$x + rnorm(60, 0, 0.7)
  fit <- illume::ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  nd <- data.frame(x = c(-1, 2))
  m <- ilm_predict(fit, nd, uncertainty = "latent", nsim = 50)
  expect_equal(m$draws, matrix(m$estimate, 2, 50))          # nothing left to draw
  r <- ilm_predict(fit, nd, target = "response", uncertainty = "latent", nsim = 50)
  dr <- illume::ilm_draws(fit, nsim = 1, given = "parameters", natural = FALSE)
  s <- exp(dr$mode[dr$map$block == "logdisp"][[1]])        # the fit's own residual SD
  expect_equal(lmn_pred_quantile(r, 0.9), stats::qnorm(0.9, m$estimate, s), tolerance = 1e-8)
  ## where "fixed" still draws the fixed effects and the dispersion
  f <- ilm_predict(fit, nd, uncertainty = "fixed", nsim = 50)
  expect_gt(stats::sd(f$draws[1, ]), 0)
})
