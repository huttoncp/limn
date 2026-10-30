test_that("a continuous PIT is the predictive CDF at the outcome", {
  m <- c(0, 2); s <- 1.5
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = log(s)))
  u <- ilm_pit(p, c(1, 1))
  expect_equal(as.vector(u), pnorm(1, m, s))
  expect_false(attr(u, "randomised"))
})

test_that("a discrete PIT is randomised within its jump, reproducibly", {
  p <- test_pred(test_dist_poisson(), m = log(c(1, 4, 9)), S = 20)
  y <- c(0, 3, 11)
  lo <- lmn_pred_cdf(p, y - 1); hi <- lmn_pred_cdf(p, y)
  set.seed(5); before <- .Random.seed
  u <- ilm_pit(p, y, seed = 7)
  expect_identical(.Random.seed, before)                 # the caller's stream untouched
  expect_true(all(u >= lo & u <= hi))
  expect_true(attr(u, "randomised"))
  expect_identical(u, ilm_pit(p, y, seed = 7))
})

test_that("PITs of outcomes drawn from the predictive distribution are uniform", {
  set.seed(11)
  n <- 3000
  m <- rnorm(n, 1)
  y <- rpois(n, exp(m))
  p <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(m)))
  u <- ilm_pit(p, y, seed = 1)
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
  expect_equal(mean(u), 0.5, tolerance = 0.02)
})

test_that("missing outcomes give a missing PIT", {
  p <- test_pred(test_dist_gaussian(), m = c(0, 1), S = 3, par = list(logsd = 0))
  u <- ilm_pit(p, c(NA, 0.3))
  expect_true(is.na(u[1]))
  expect_false(is.na(u[2]))
  expect_error(ilm_pit(list(), 1), "ilm_pred")
})
