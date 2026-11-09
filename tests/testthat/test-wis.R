## The weighted interval score, against Bracher et al. (2021)'s definitions:
## equations (1) and (2) (p. 4), the K = 0 case (p. 4), equation (3)'s
## approximation to the CRPS (p. 4) and equation (4)'s quantile form (p. 5).

test_that("the WIS is equation (1) with equation (2)'s weights", {
  y <- c(-2.3, 0.1, 1.7, 6)
  m <- c(0, 0.5, 1, 2); s <- 1.3
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(s)))
  alpha <- c(0.02, 0.05, 1:9 / 10)
  by_hand <- 0.5 * abs(y - m)                                     # w0 = 1/2, median m
  for (a in alpha) {
    l <- stats::qnorm(a / 2, m, s); u <- stats::qnorm(1 - a / 2, m, s)
    is <- (u - l) + 2 / a * (l - y) * (y < l) + 2 / a * (y - u) * (y > u)
    by_hand <- by_hand + a / 2 * is                               # w_k = alpha_k / 2
  }
  by_hand <- by_hand / (length(alpha) + 1 / 2)
  expect_equal(lmn_pred_scores(p, y, "wis")$wis, by_hand, tolerance = 1e-8)
  ## the default is Bracher et al.'s set: 11 intervals and the median
  expect_equal(lmn_pred_scores(p, y, "wis", alpha = alpha)$wis,
               lmn_pred_scores(p, y, "wis")$wis)
  expect_equal(lmn_wis_alpha, c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  expect_equal(eval(formals(ilm_scores.ilm_pred)$alpha), lmn_wis_alpha)
})

test_that("the WIS is the mean of twice the quantile scores (equation (4))", {
  y <- c(0, 2, 5, 17); lam <- c(0.5, 2, 3, 9)
  pp <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(log(lam))))
  alpha <- c(0.1, 0.5, 0.8)
  tau <- sort(c(alpha / 2, 0.5, 1 - alpha / 2))                   # 2K + 1 quantiles
  pin <- lmn_pred_scores(pp, y, "pinball", probs = tau)
  q_form <- unname(2 * rowSums(pin[paste0("pinball_", vapply(tau, format, ""))]) / length(tau))
  expect_equal(lmn_pred_scores(pp, y, "wis", alpha = alpha)$wis, q_form, tolerance = 1e-12)
  ## the same with the default set, outcomes below, inside and above every interval,
  ## the pinball loss written as Bracher et al. write it, (1{y <= q} - tau)(q - y)
  y <- c(-9, -1.2, -0.3, 0, 0.4, 2.2, 11)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(rep(0, 7))), list(logsd = 0))
  tau <- sort(c(lmn_wis_alpha / 2, 0.5, 1 - lmn_wis_alpha / 2))
  q <- stats::qnorm(tau)
  pin <- vapply(y, function(v) sum(((v <= q) - tau) * (q - v)), 0)
  expect_equal(lmn_pred_scores(p, y, "wis")$wis, pin / (length(lmn_wis_alpha) + 1 / 2),
               tolerance = 1e-12)
})

test_that("the default set is the median and 11 intervals: 23 levels, 0.01 to 0.99", {
  tau <- sort(c(lmn_wis_alpha / 2, 0.5, 1 - lmn_wis_alpha / 2))
  expect_length(unique(tau), 23L)
  expect_equal(range(tau), c(0.01, 0.99))
})

test_that("the worked example of Bracher et al., as their correction gives it", {
  ## two negative binomial forecasts scored at y = 190 (section 3.2, p. 7): 105.3
  ## for F (mean 60, size 4) and 88.9 for G (mean 80, size 10), as corrected in
  ## Bracher et al. (2022, p. 1); the article printed 103.9 and 87.8
  wis <- function(mu, size) {
    p <- lmn_new_pred(test_dist_nbinom(), list(mu = matrix(log(mu))), list(logsize = log(size)))
    tau <- sort(c(lmn_wis_alpha / 2, 0.5, 1 - lmn_wis_alpha / 2))
    expect_equal(vapply(tau, function(t) lmn_pred_quantile(p, t), 0),
                 stats::qnbinom(tau, size = size, mu = mu))
    lmn_pred_scores(p, 190, "wis")$wis
  }
  expect_equal(round(wis(60, 4), 1), 105.3)
  expect_equal(round(wis(80, 10), 1), 88.9)
})

test_that("the default WIS is on the CRPS's scale", {
  ## a standard normal forecast, outcomes 0 to 6 SDs from the median: the
  ## default WIS is 0.85 to 1 times the CRPS (limn measures 0.878 to 0.966)
  y <- seq(0, 6, by = 0.1)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(rep(0, length(y)))), list(logsd = 0))
  r <- lmn_pred_scores(p, y, c("wis", "crps"))
  expect_true(all(r$wis / r$crps > 0.85 & r$wis / r$crps < 1))
})

test_that("with no intervals the WIS is the absolute error", {
  y <- c(-2.3, 0.1, 1.7, 6); m <- c(0, 0.5, 1, 2)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(1.3)))
  sc <- lmn_pred_scores(p, y, c("wis", "absolute"), alpha = numeric(0))
  expect_equal(sc$wis, sc$absolute)
})

test_that("with many evenly spaced intervals the WIS approaches the CRPS (equation (3))", {
  y <- c(-2.3, 0.1, 1.7, 6); m <- c(0, 0.5, 1, 2)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(1.3)))
  crps <- lmn_pred_scores(p, y, "crps")$crps
  w99 <- lmn_pred_scores(p, y, "wis", alpha = 1:99 / 100)$wis
  w9 <- lmn_pred_scores(p, y, "wis", alpha = 1:9 / 10)$wis
  ## limn's own measurement of the gap; Bracher et al. state only the approximation
  expect_lt(max(abs(w99 / crps - 1)), 0.01)
  expect_true(all(abs(w99 - crps) < abs(w9 - crps)))
})

test_that("the WIS goes through ilm_scores() and checks its levels", {
  y <- c(-2.3, 0.1, 1.7, 6); m <- c(0, 0.5, 1, 2)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(1.3)))
  s <- ilm_scores(p, y, rules = c("crps", "wis"), alpha = c(0.1, 0.5))
  expect_equal(attr(s, "rules"), c("crps", "wis"))
  expect_equal(s$wis, lmn_pred_scores(p, y, "wis", alpha = c(0.1, 0.5))$wis)
  expect_equal(lmn_default_rule(data.frame(x = 1, wis = 2), NULL), "wis")
  expect_error(ilm_scores(p, y, rules = "wis", alpha = 1), "strictly between 0 and 1")
  expect_error(ilm_scores(p, y, rules = "wis", alpha = c(0, 0.5)), "strictly between 0 and 1")
  expect_error(ilm_scores(p, y, rules = "wis", alpha = c(0.2, 0.2)), "more than once")
  expect_error(ilm_scores(p, y, rules = "wis", alpha = NA_real_), "strictly between 0 and 1")
  ## other rules don't read alpha
  expect_silent(ilm_scores(p, y, rules = "crps", alpha = 2))
})

test_that("ilm_cv() scores the WIS, with alpha passed through", {
  set.seed(1)
  d <- data.frame(site = factor(rep(1:8, each = 8)), x = rnorm(64))
  d$y <- rpois(64, exp(0.5 + 0.3 * d$x + rnorm(8, 0, 0.3)[d$site]))
  fit <- suppressWarnings(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                                            verbose = FALSE))
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 4)
  cv <- suppressWarnings(ilm_cv(fit, d, folds, rules = c("wis", "absolute"), nsim = 50,
                                alpha = numeric(0)))
  expect_true("wis" %in% attr(cv, "rules"))
  expect_equal(cv$wis, cv$absolute)
})
