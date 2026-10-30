## The latent's steps against conditioning the process's own covariance, and
## the sequential paths against the joint conditional distribution.

cov_stationary <- function(t, s2, rho) s2 * rho^abs(outer(t, t, "-"))
cov_walk <- function(t, q, t0) q * outer(t - t0, t - t0, pmin)
cond <- function(C, known, a) {                  # mean and covariance given C[known]
  u <- setdiff(seq_len(nrow(C)), known)
  W <- C[u, known, drop = FALSE] %*% solve(C[known, known])
  list(mean = drop(W %*% a), cov = C[u, u, drop = FALSE] - W %*% C[known, u, drop = FALSE])
}

test_that("each step is the conditional distribution of the process", {
  par_c <- list(type = "car1", sig2 = 0.7, rho = 0.55, q = NULL, step = 1)
  par_a <- list(type = "ar1", sig2 = 0.7, rho = -0.4, q = NULL, step = 2)
  par_w <- list(type = "rw1", sig2 = NULL, rho = NULL, q = 0.3, step = 1)
  aL <- 0.8; aR <- -0.5
  ## forward (or backward) from one neighbour: t = 3.7 from 1.2
  for (par in list(par_c, par_a)) {
    tt <- if (par$type == "ar1") c(1, 7) else c(1.2, 3.7)
    ref <- cond(cov_stationary(tt / par$step, par$sig2, par$rho), 1, aL)
    st <- lmn_latent_step("step", matrix(aL), NULL, diff(tt), NA, par)
    expect_equal(c(st$mean, st$var), c(ref$mean, ref$cov), tolerance = 1e-12)
  }
  ## a bridge: t = 2.5 between 1.2 and 4.1
  tt <- c(1.2, 4.1, 2.5)
  ref <- cond(cov_stationary(tt, par_c$sig2, par_c$rho), 1:2, c(aL, aR))
  st <- lmn_latent_step("bridge", matrix(aL), matrix(aR), 2.5 - 1.2, 4.1 - 2.5, par_c)
  expect_equal(c(st$mean, st$var), c(ref$mean, ref$cov), tolerance = 1e-12)
  ## a walk anchored at t0 = 0: forward, and a bridge
  ref <- cond(cov_walk(c(1.2, 3.7), 0.3, 0), 1, aL)
  st <- lmn_latent_step("step", matrix(aL), NULL, 2.5, NA, par_w)
  expect_equal(c(st$mean, st$var), c(ref$mean, ref$cov), tolerance = 1e-12)
  ref <- cond(cov_walk(tt, 0.3, 0), 1:2, c(aL, aR))
  st <- lmn_latent_step("bridge", matrix(aL), matrix(aR), 2.5 - 1.2, 4.1 - 2.5, par_w)
  expect_equal(c(st$mean, st$var), c(ref$mean, ref$cov), tolerance = 1e-12)
  ## the stationary start
  st <- lmn_latent_step("start", matrix(0), NULL, NA, NA, par_c)
  expect_equal(c(st$mean, st$var), c(0, 0.7))
})

## one series with cells at 1, 2, 5 and 6, and new times before, between and
## after them; the paths' draws against the joint conditional distribution
path_case <- function() {
  cells <- data.frame(group = factor("a"), time = c(1, 2, 5, 6))
  a_c <- c(0.4, -0.2, 0.9, 0.6)
  tt <- c(0, 3, 4, 8, 9, 3)                            # 3 twice: one shared value
  ar <- data.frame(group = "a", time = tt, cell = NA_integer_,
                   prev_cell = c(NA, 2, 2, 4, 4, 2), next_cell = c(1, 3, 3, NA, NA, 3))
  list(cells = cells, a_c = a_c, ar = ar, tt = tt)
}

test_that("a seen series' new times are one draw from their joint conditional", {
  S <- 40000
  x <- path_case()
  ## CAR(1): the cells at 1, 2, 5 and 6 are known; the new times are 0, 3, 4, 8, 9
  par <- list(type = "car1", sig2 = rep(0.8, S), rho = rep(0.6, S), step = 1)
  got <- lmn_latent_paths(x$ar, x$cells, matrix(x$a_c, 4, S), par, rep(FALSE, 6), seed = 7)
  expect_equal(got$values[2, ], got$values[6, ])        # the same time shares its value
  ref <- cond(cov_stationary(c(1, 2, 5, 6, 0, 3, 4, 8, 9), 0.8, 0.6), 1:4, x$a_c)
  v <- got$values[1:5, ]
  expect_true(all(abs(rowMeans(v) - ref$mean) < 4 * sqrt(diag(ref$cov) / S)))
  expect_lt(max(abs(stats::cov(t(v)) - ref$cov)), 0.02)
  ## a walk anchored at zero at 1: the cells at 2, 5 and 6 are known, and the
  ## time before the anchor walks back from it independently
  par_w <- list(type = "rw1", q = rep(0.25, S))
  a_w <- matrix(c(0, x$a_c[2:4]), 4, S)
  got <- lmn_latent_paths(x$ar, x$cells, a_w, par_w, rep(FALSE, 6), seed = 7)
  ref <- cond(cov_walk(c(2, 5, 6, 3, 4, 8, 9), 0.25, 1), 1:3, x$a_c[2:4])
  v <- got$values[2:5, ]
  expect_true(all(abs(rowMeans(v) - ref$mean) < 4 * sqrt(diag(ref$cov) / S)))
  expect_lt(max(abs(stats::cov(t(v)) - ref$cov)), 0.02)
  expect_lt(abs(mean(got$values[1, ])), 4 * sqrt(0.25 / S))
  expect_lt(abs(stats::var(got$values[1, ]) - 0.25), 0.02)
  expect_lt(max(abs(stats::cor(got$values[1, ], t(v)))), 0.03)
})

test_that("a new series is a fresh path from the stationary distribution", {
  S <- 40000
  par <- list(type = "car1", sig2 = rep(0.5, S), rho = rep(0.7, S), step = 1)
  cells <- data.frame(group = factor("a"), time = 1:3)
  ar <- data.frame(group = "z", time = c(2, 3.5, 6), cell = NA_integer_,
                   prev_cell = NA_integer_, next_cell = NA_integer_)
  got <- lmn_latent_paths(ar, cells, matrix(0, 3, S), par, rep(TRUE, 3), seed = 3)
  expect_equal(stats::cov(t(got$values)), cov_stationary(c(2, 3.5, 6), 0.5, 0.7), tolerance = 0.05)
  ## a new walk starts from zero at its first time
  par_w <- list(type = "rw1", q = rep(0.2, S))
  got <- lmn_latent_paths(ar, cells, matrix(0, 3, S), par_w, rep(TRUE, 3), seed = 3)
  expect_equal(got$values[1, ], rep(0, S))
  expect_equal(stats::cov(t(got$values[2:3, ])), cov_walk(c(3.5, 6), 0.2, 2), tolerance = 0.05)
})
