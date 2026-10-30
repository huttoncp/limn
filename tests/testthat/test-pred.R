test_that("the constructor checks its pieces against the family", {
  dg <- test_dist_gaussian()
  lp <- list(mu = matrix(0, 3, 4))
  expect_error(lmn_new_pred(dg, lp), "needs parameter")
  expect_error(lmn_new_pred(dg, list(eta = matrix(0, 3, 4)), list(logsd = 0)),
               "linear predictor")
  expect_error(lmn_new_pred(dg, lp, list(logsd = c(0, 0))), "one value per draw")
  expect_error(lmn_new_pred(dg, lp, list(logsd = 0), weights = c(1, -1, 1, 1)),
               "non-negative")
  expect_error(lmn_new_pred(list(d = 1), lp), "must supply functions")
  p <- lmn_new_pred(dg, lp, list(logsd = 0))
  expect_s3_class(p, "ilm_pred")
  expect_output(print(p), "a new observation for 3 rows, from 4 draws")
})

test_that("with one draw, the mixture is the family's own distribution", {
  dg <- test_dist_gaussian()
  m <- c(-1, 0, 2.5)
  p <- lmn_new_pred(dg, list(mu = matrix(m)), list(logsd = log(1.7)))
  expect_equal(lmn_pred_mean(p), m)
  expect_equal(lmn_pred_cdf(p, c(0, 0, 0)), pnorm(0, m, 1.7))
  for (pr in c(0.05, 0.5, 0.9))
    expect_equal(lmn_pred_quantile(p, pr), qnorm(pr, m, 1.7), tolerance = 1e-10)
})

test_that("a mixture's CDF and mean are the weighted averages over draws", {
  dg <- test_dist_gaussian()
  mu <- matrix(c(0, 1, 3, -2, 0.5, 4), 2)      # 2 rows, 3 draws
  w <- c(0.2, 0.3, 0.5)
  sd <- exp(c(0, 0.5, -0.2))
  p <- lmn_new_pred(dg, list(mu = mu), list(logsd = log(sd)), weights = w)
  q <- c(0.7, -1)
  want <- vapply(1:2, function(i) sum(w * pnorm(q[i], mu[i, ], sd)), 0)
  expect_equal(lmn_pred_cdf(p, q), want)
  expect_equal(lmn_pred_mean(p), as.vector(mu %*% w))
})

test_that("a mixture quantile inverts the mixture CDF", {
  dg <- test_dist_gaussian()
  p <- test_pred(dg, m = c(0, 5), S = 40, spread = 1.5, par = list(logsd = 0))
  for (pr in c(0.01, 0.3, 0.97)) {
    q <- lmn_pred_quantile(p, pr)
    expect_equal(lmn_pred_cdf(p, q), rep(pr, 2), tolerance = 1e-10)
  }
})

test_that("the capped root takes at most the bisection's passes plus 2, and its answer", {
  ## Two three-draw mixtures on which the first build's guard, which counted
  ## from the bracket's current tolerance rather than the smallest it could
  ## still reach, overran the cap by one pass (45 against 42, 43 against 40).
  ## The draws are hex doubles so they are these exact mixtures everywhere.
  ## A third case hands the root a wide bracket closing towards a quantile
  ## near zero, where the stopping width shrinks by a factor of 1,000 on the
  ## way: the case the guard and the growing budget must both get right.
  dw <- list(family = "weibull", discrete = FALSE, lp_names = "mu", par_names = "logsd",
             d = function(y, lp, par, log = FALSE) stop("unused"),
             p = function(q, lp, par) stats::pweibull(q, 1 / exp(par$logsd), exp(lp$mu)),
             q = function(prob, lp, par) stats::qweibull(prob, 1 / exp(par$logsd), exp(lp$mu)),
             r = function(lp, par) stop("unused"), mean = function(lp, par) stop("unused"))
  hex <- function(...) as.numeric(c(...))
  cases <- list(
    list(dist = test_dist_gaussian(), p = 0.1,
         m = hex("0x1.0e1ad232faef5p+4", "0x1.828f94a46b504p+5", "-0x1.18c56bf2ac3fbp+4"),
         s = hex("0x1.5c33a000c1c3p-5", "0x1.ef11fa4726febp-6", "0x1.123357775564fp-5"),
         w = hex("0x1.06eee19f57c46p-1", "0x1.ce55657083362p-2", "0x1.1e66ba866a08ap-5")),
    list(dist = dw, p = 0.975,
         m = hex("-0x1.0b876f877e3dep-1", "-0x1.a9eb708d22a0ep-2", "-0x1.cc98de21fd759p-2"),
         s = hex("0x1.999999999999ap-4", "0x1.4eae46599bba6p-1", "0x1.9e24bf5d58a6cp-3"),
         w = hex("0x1.31c7aa3a3bbc9p-6", "0x1.07b1760f64a81p-1", "0x1.dd80993d92f43p-2")),
    list(dist = test_dist_gaussian(), p = 0.5, m = 0.001, s = 1, w = 1, bracket = c(-3, 1000)))
  real_mix <- lmn_mix
  n <- 0L
  local_mocked_bindings(lmn_mix = function(...) { n <<- n + 1L; real_mix(...) })
  for (cs in cases) {
    pr <- lmn_new_pred(cs$dist, list(mu = matrix(cs$m, 1)), par = list(logsd = log(cs$s)), weights = cs$w)
    if (is.null(cs$bracket)) {
      Q <- lmn_draw_quantiles(pr, 1L, cs$p)
      b <- lmn_widen_bracket(pr, 1L, cs$p, min(Q), max(Q), FALSE)
    } else {
      b <- list(lo = cs$bracket[1], hi = cs$bracket[2],
                flo = real_mix(pr, 1L, matrix(cs$bracket[1]), "p")[, 1],
                fhi = real_mix(pr, 1L, matrix(cs$bracket[2]), "p")[, 1])
    }
    ## the bisection's passes and answer, from the same bracket to the same width
    lo <- b$lo; hi <- b$hi; ref <- 0L
    while (hi - lo > 1e-12 * max(1, abs(lo), abs(hi))) {
      mid <- (lo + hi) / 2; ref <- ref + 1L
      if (real_mix(pr, 1L, matrix(mid), "p")[, 1] >= cs$p) hi <- mid else lo <- mid
    }
    n <- 0L
    q <- lmn_quantile_root(pr, 1L, cs$p, b$lo, b$hi, b$flo, b$fhi)
    expect_lte(n, ref + 2L)
    expect_lte(abs(q - (lo + hi) / 2), 2e-12 * max(1, abs(q)))
  }
})

test_that("a discrete mixture quantile is the smallest k with F(k) >= p", {
  dp <- test_dist_poisson()
  p <- test_pred(dp, m = log(c(0.4, 3, 40)), S = 30, spread = 0.5)
  for (pr in c(0.05, 0.5, 0.95)) {
    q <- lmn_pred_quantile(p, pr)
    brute <- vapply(1:3, function(i) {
      k <- 0:400
      Fk <- colMeans(outer(exp(p$lp$mu[i, ]), k, function(l, kk) ppois(kk, l)))
      k[which(Fk >= pr)[1]]
    }, 0)
    expect_equal(q, brute)
  }
})
