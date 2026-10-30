## Scores are checked against scoringRules, an independent implementation, where
## it has the same score, and against direct formulas otherwise.

test_that("scores reproduce Gneiting and Raftery's (2007) formulas", {
  y <- c(-2.3, 0.1, 1.7, 6)
  m <- c(0, 0.5, 1, 2); s <- 1.3
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(s)))
  ## the gaussian CRPS in closed form (section 4.2, p. 367), turned into a loss
  z <- (y - m) / s
  gr <- -s * (1 / sqrt(pi) - 2 * stats::dnorm(z) - z * (2 * stats::pnorm(z) - 1))
  expect_equal(lmn_pred_scores(p, y, "crps")$crps, gr, tolerance = 1e-10)
  ## the interval score, equation (43), for the central 80% interval
  l <- stats::qnorm(0.1, m, s); u <- stats::qnorm(0.9, m, s); a <- 0.2
  expect_equal(lmn_pred_scores(p, y, "interval", level = 0.8)[["interval_0.8"]],
               (u - l) + 2 / a * (l - y) * (y < l) + 2 / a * (y - u) * (y > u),
               tolerance = 1e-8)
  ## the quantile score, equation (41), turned into the pinball loss
  q <- stats::qnorm(0.3, m, s)
  expect_equal(lmn_pred_scores(p, y, "pinball", probs = 0.3)[["pinball_0.3"]],
               -(y - q) * ((y <= q) - 0.3), tolerance = 1e-8)
})

test_that("one-draw scores match scoringRules: gaussian, Poisson, negative binomial", {
  skip_if_not_installed("scoringRules")
  y <- c(-2.3, 0.1, 1.7, 6)
  m <- c(0, 0.5, 1, 2); s <- 1.3
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(s)))
  sc <- lmn_pred_scores(p, y, c("log", "crps"))
  expect_equal(sc$log, scoringRules::logs_norm(y, m, s), tolerance = 1e-12)
  expect_equal(sc$crps, scoringRules::crps_norm(y, m, s), tolerance = 1e-10)

  yc <- c(0, 2, 5, 17)
  lam <- c(0.5, 2, 3, 9)
  pp <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(log(lam))))
  sp <- lmn_pred_scores(pp, yc, c("log", "crps"))
  expect_equal(sp$log, scoringRules::logs_pois(yc, lam), tolerance = 1e-12)
  expect_equal(sp$crps, scoringRules::crps_pois(yc, lam), tolerance = 1e-10)

  size <- 2.5
  pn <- lmn_new_pred(test_dist_nbinom(), list(mu = matrix(log(lam))),
                     list(logsize = log(size)))
  sn <- lmn_pred_scores(pn, yc, c("log", "crps"))
  expect_equal(sn$log, scoringRules::logs_nbinom(yc, size = size, mu = lam),
               tolerance = 1e-12)
  ## scoringRules' negative binomial CRPS is a closed form through 'hypergeo'.
  ## In scoringRules 1.1.3 it returns -Inf for some sizes (2.5, with a mean of
  ## 2 or more), so it is compared at sizes where it is finite; the test below
  ## checks every size against the exact energy form.
  skip_if_not_installed("hypergeo")
  for (sz in c(0.7, 10)) {
    pz <- lmn_new_pred(test_dist_nbinom(), list(mu = matrix(log(lam))),
                       list(logsize = log(sz)))
    expect_equal(lmn_pred_scores(pz, yc, "crps")$crps,
                 scoringRules::crps_nbinom(yc, size = sz, mu = lam), tolerance = 1e-10)
  }
})

test_that("the negative binomial CRPS equals the exact energy form over its pmf", {
  ## E|X - y| - E|X - X2| / 2, both summed over the probability mass function:
  ## a different computation from limn's sum over the CDF
  energy <- function(y, size, mu) {
    k <- 0:qnbinom(1 - 1e-14, size = size, mu = mu)
    p <- dnbinom(k, size = size, mu = mu); p <- p / sum(p)
    sum(abs(k - y) * p) - 0.5 * sum(abs(outer(k, k, `-`)) * outer(p, p))
  }
  yc <- c(0, 2, 5, 17); lam <- c(0.5, 2, 3, 9)
  for (sz in c(0.7, 2.5, 10)) {
    pz <- lmn_new_pred(test_dist_nbinom(), list(mu = matrix(log(lam))),
                       list(logsize = log(sz)))
    expect_equal(lmn_pred_scores(pz, yc, "crps")$crps,
                 mapply(energy, yc, sz, lam), tolerance = 1e-10)
  }
})

test_that("mixture scores match scoringRules' closed forms for normal mixtures", {
  skip_if_not_installed("scoringRules")
  set.seed(2)
  S <- 40
  mu <- rbind(rnorm(S, 0, 0.4), rnorm(S, 3, 2.5))    # a tight and a wide mixture
  sd <- exp(rnorm(S, 0, 0.3))
  y <- c(0.2, -4)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = mu), list(logsd = log(sd)))
  sc <- lmn_pred_scores(p, y, c("log", "crps"))
  w <- matrix(1 / S, 2, S)
  sdm <- matrix(sd, 2, S, byrow = TRUE)
  expect_equal(sc$crps, scoringRules::crps_mixnorm(y, mu, sdm, w), tolerance = 1e-10)
  expect_equal(sc$log, scoringRules::logs_mixnorm(y, mu, sdm, w), tolerance = 1e-12)
})

test_that("a count mixture's CRPS is the sum over its support", {
  p <- test_pred(test_dist_poisson(), m = log(c(0.7, 12)), S = 25, spread = 0.6)
  y <- c(3, 4)
  brute <- vapply(1:2, function(i) {
    k <- 0:2000
    Fk <- colMeans(outer(exp(p$lp$mu[i, ]), k, function(l, kk) ppois(kk, l)))
    sum((Fk - (k >= y[i]))^2)
  }, 0)
  expect_equal(lmn_pred_scores(p, y, "crps")$crps, brute, tolerance = 1e-10)
})

test_that("outcomes far in a tail still score correctly", {
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(0)), list(logsd = 0))
  sc <- lmn_pred_scores(p, 60, c("log", "crps"))
  expect_equal(sc$log, -dnorm(60, log = TRUE))          # finite, on the log scale
  expect_equal(sc$crps, 60 - 1 / sqrt(pi), tolerance = 1e-8)  # CRPS ~ |y| - 1/sqrt(pi)
  pc <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(log(2))))
  expect_equal(lmn_pred_scores(pc, 80, "log")$log, -dpois(80, 2, log = TRUE))
})

test_that("Brier, pinball, interval, squared and absolute follow their definitions", {
  m <- c(0, 1); s <- 2
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(s)))
  y <- c(3.5, 0.8)
  sc <- lmn_pred_scores(p, y, c("pinball", "interval", "squared", "absolute"),
                        probs = c(0.25, 0.9), level = 0.8)
  for (pr in c(0.25, 0.9)) {
    q <- qnorm(pr, m, s)
    expect_equal(sc[[paste0("pinball_", pr)]], (y - q) * (pr - (y < q)), tolerance = 1e-10)
  }
  l <- qnorm(0.1, m, s); u <- qnorm(0.9, m, s)
  wink <- (u - l) + (2 / 0.2) * (l - y) * (y < l) + (2 / 0.2) * (y - u) * (y > u)
  expect_equal(sc$interval_0.8, wink, tolerance = 1e-10)
  expect_equal(sc$squared, (y - m)^2)
  expect_equal(sc$absolute, abs(y - m), tolerance = 1e-10)

  pb <- lmn_new_pred(test_dist_binomial(), list(mu = matrix(qlogis(c(0.2, 0.9)))))
  expect_equal(lmn_pred_scores(pb, c(1, 1), "brier")$brier, (c(0.2, 0.9) - 1)^2)
  pc <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(0)))
  expect_error(lmn_pred_scores(pc, 1, "brier"), "binary")
})

test_that("missing outcomes score NA, and row descriptions are kept", {
  p <- test_pred(test_dist_gaussian(), m = c(0, 1, 2), S = 5, par = list(logsd = 0))
  p$rows <- data.frame(row = 1:3, group = c("a", "a", "b"))
  sc <- lmn_pred_scores(p, c(0.5, NA, 1), c("log", "crps"))
  expect_s3_class(sc, "ilm_pred_scores")
  expect_equal(sc$group, c("a", "a", "b"))
  expect_true(is.na(sc$log[2]) && is.na(sc$crps[2]))
  expect_false(anyNA(sc$log[-2]))
  expect_identical(attr(sc, "orientation"), "lower is better")
})

test_that("bad requests stop with a reason", {
  p <- test_pred(test_dist_poisson(), m = c(0, 1), S = 5)
  expect_error(lmn_pred_scores(p, c(1, 2), "logg"), "unknown rule")
  expect_error(lmn_pred_scores(p, c(1, 2.5), "log"), "whole numbers")
  expect_error(lmn_pred_scores(p, 1, "log"), "1 values for 2 predicted rows")
  expect_error(lmn_pred_scores(p, c(1, 2), "interval", level = 1.2), "between 0 and 1")
})

test_that("a wide discrete support is summed by Euler-Maclaurin, within its error estimate", {
  ## the exact sum over every integer of the support, for reference
  exact <- function(p, y) vapply(seq_along(y), function(i) {
    lo <- min(lmn_draw_quantiles(p, i, 1e-12)); hi <- max(lmn_draw_quantiles(p, i, 1 - 1e-12))
    Fk <- lmn_mix_cdf_run(p, i, lo, hi - lo + 1)[1, ]; k <- lo:hi
    sum((Fk - (k >= y[i]))^2) + max(0, lo - y[i]) + max(0, y[i] - hi - 1)
  }, 0)
  pn <- test_pred(test_dist_nbinom(), m = log(c(3000, 20000)), S = 60, spread = 0.05,
                  par = list(logsize = rep(log(20), 60)))
  y <- c(2500, 21000)
  w <- lmn_crps_discrete(pn, y)
  ex <- exact(pn, y)
  expect_true(all(attr(w, "step") > 1))                       # both rows took the wide path
  expect_true(all(abs(as.numeric(w) - ex) <= 1e-6 * ex))
  expect_true(all(abs(as.numeric(w) - ex) <= attr(w, "error_estimate")))
  expect_false(any(attr(w, "flagged")))
  expect_true(all(attr(w, "bound") > 0))                      # the rigorous bracket, for the record
  ## a row either side of the threshold: summed exactly, and by the wide path
  w1 <- lmn_crps_discrete(pn, y, max_support = 1e6)
  expect_identical(attr(w1, "error_estimate"), c(0, 0))
  expect_equal(as.numeric(w1), ex, tolerance = 1e-12)
  ## narrow rows are summed exactly, with nothing to estimate
  pp <- test_pred(test_dist_poisson(), m = log(c(200, 900)), S = 40, spread = 0.2)
  wn <- lmn_crps_discrete(pp, c(190, 950))
  expect_identical(attr(wn, "error_estimate"), c(0, 0))
  expect_null(attr(lmn_pred_scores(pp, c(190, 950), rules = "crps"), "crps_bound"))
  ## the scores carry the error estimate and the bracket
  s <- lmn_pred_scores(pn, y, rules = "crps")
  expect_equal(attr(s, "crps_bound"), attr(w, "error_estimate"))
  expect_equal(attr(s, "crps_bracket"), attr(w, "bound"))
})

test_that("an average's CRPS error estimate is its total's, scaled like the score", {
  dp <- test_dist_poisson()
  p <- test_pred(dp, m = log(1e6), S = 5, spread = 0.05)
  tot <- lmn_crps_discrete(p, 1000400)
  avg <- p; avg$scale <- 0.5
  s <- lmn_pred_scores(avg, 500200, rules = "crps")
  expect_equal(s$crps, 0.5 * as.numeric(tot))
  expect_equal(attr(s, "crps_bound"), 0.5 * attr(tot, "error_estimate"))
  expect_true(attr(s, "crps_bound") > 0)
})

test_that("the discrete CRPS by running sums equals the sum over the CDF at every point", {
  ## the direct form: the mixture's CDF evaluated at each point of the support
  direct <- function(p, y) vapply(seq_along(y), function(i) {
    lo <- min(lmn_draw_quantiles(p, i, 1e-12)); hi <- max(lmn_draw_quantiles(p, i, 1 - 1e-12))
    k <- lo:hi
    Fk <- lmn_mix(p, i, matrix(k, 1L), "p")[1, ]
    sum((Fk - (k >= y[i]))^2) + max(0, lo - y[i]) + max(0, y[i] - hi - 1)
  }, 0)
  ## rows of very different widths, so that they fall in different width groups,
  ## and outcomes inside, below and above the support
  m <- log(c(0.3, 4, 60, 900, 2500))
  y <- c(0, 3, 400, 10, 2600)
  pp <- test_pred(test_dist_poisson(), m = m, S = 30, spread = 0.2)
  expect_gt(length(lmn_width_groups(vapply(seq_along(y), function(i)
    max(lmn_draw_quantiles(pp, i, 1 - 1e-12)) - min(lmn_draw_quantiles(pp, i, 1e-12)) + 1, 0))), 2L)
  ## the widest row exceeds 2,000 values, so max_support is raised to test the
  ## exact sum (the wide path has its own test)
  expect_equal(as.numeric(lmn_crps_discrete(pp, y, max_support = 1e6)), direct(pp, y), tolerance = 1e-12)
  pn <- test_pred(test_dist_nbinom(), m = m, S = 30, spread = 0.2, par = list(logsize = rep(log(3), 30)))
  ## likewise for the negative binomial rows
  expect_equal(as.numeric(lmn_crps_discrete(pn, y, max_support = 1e6)), direct(pn, y), tolerance = 1e-12)
  ## the narrowest supports: two points (a Poisson mean of 1e-9 reaches 1 at 1 - 1e-12)
  p1 <- test_pred(test_dist_poisson(), m = rep(log(1e-9), 2), S = 5, spread = 0)
  expect_equal(as.numeric(lmn_crps_discrete(p1, c(0, 2))), direct(p1, c(0, 2)), tolerance = 1e-12)
  ## and one point (a mean of 1e-14 stays at 0 even at 1 - 1e-12)
  p0 <- test_pred(test_dist_poisson(), m = rep(log(1e-14), 2), S = 5, spread = 0)
  expect_equal(as.numeric(lmn_crps_discrete(p0, c(0, 3))), direct(p0, c(0, 3)), tolerance = 1e-12)
})

test_that("width groups hold rows within a factor of 1.5 of the group's narrowest", {
  w <- c(1, 3, 6, 7, 100, 140, 200, 5000)
  g <- lmn_width_groups(w)
  expect_setequal(unlist(g), seq_along(w))
  for (k in g) expect_true(max(w[k]) <= 1.5 * min(w[k]) + 5)
})

test_that("rows whose outcome lies outside the support take the tighter bounds, and say so", {
  p <- test_pred(test_dist_poisson(), m = log(c(20, 20, 20)), S = 50, spread = 0.1)
  inside <- lmn_crps_discrete(p, c(18, 20, 25))
  expect_identical(attr(inside, "recomputed"), 0L)
  far <- lmn_crps_discrete(p, c(0, 20, 400))
  expect_identical(attr(far, "recomputed"), 2L)
  ## the recomputed rows are exact against the sum over the tighter support
  tight <- lmn_crps_discrete(p, c(0, 20, 400), eps = 1e-12, eps_outside = 1e-12)
  expect_equal(as.numeric(far), as.numeric(tight), tolerance = 1e-12)
})
