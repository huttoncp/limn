scores2 <- function(seed = 1, G = 12, per = 6) {
  set.seed(seed)
  g <- rep(seq_len(G), each = per)
  shared <- rnorm(G)[g]                 # rows of a group share their errors
  a <- data.frame(group = g, log = 1 + shared + rnorm(G * per, 0, 0.5))
  b <- data.frame(group = g, log = a$log - 0.05 + 0.4 * shared + rnorm(G * per, 0, 0.3))
  list(a = a, b = b)
}

test_that("independent rows give the paired t-test", {
  s <- scores2()
  cmp <- ilm_compare(a = s$a, b = s$b, rule = "log")
  tt <- t.test(s$b$log, s$a$log, paired = TRUE)
  expect_equal(cmp$difference, unname(tt$estimate))
  expect_equal(cmp$statistic, unname(tt$statistic))
  expect_equal(cmp$p_value, tt$p.value)
  expect_equal(cmp$df, unname(tt$parameter))
})

test_that("a clustered standard error matches sandwich::vcovCL (HC1)", {
  skip_if_not_installed("sandwich")
  s <- scores2()
  cmp <- ilm_compare(s$a, s$b, rule = "log", by = "group")
  d <- s$b$log - s$a$log
  se <- sqrt(sandwich::vcovCL(lm(d ~ 1), cluster = s$a$group)[1, 1])
  expect_equal(cmp$se, se, tolerance = 1e-12)
  expect_equal(cmp$df, 11)
  expect_equal(cmp$clusters, 12)
  ## rows sharing errors: the clustered SE is larger than the row-wise one
  expect_gt(cmp$se, ilm_compare(s$a, s$b, rule = "log")$se)
})

test_that("Diebold-Mariano at lag 0 is the paired t-test on origin means", {
  set.seed(4)
  o <- rep(1:15, each = 3)
  a <- data.frame(origin = o, horizon = 1L, crps = rexp(45))
  b <- data.frame(origin = o, horizon = 1L, crps = a$crps * runif(45, 0.7, 1.2))
  cmp <- ilm_compare(a, b, rule = "crps", by = "origin", method = "dm")
  mo <- tapply(b$crps - a$crps, o, mean)
  tt <- t.test(mo)
  expect_equal(cmp$statistic, unname(tt$statistic), tolerance = 1e-12)
  expect_equal(cmp$p_value, tt$p.value, tolerance = 1e-12)
  expect_equal(cmp$clusters, 15)
})

test_that("the default lag is the longest horizon minus one", {
  set.seed(6)
  o <- rep(1:20, each = 3)
  a <- data.frame(origin = o, horizon = rep(1:3, 20), log = rnorm(60))
  b <- data.frame(origin = o, horizon = rep(1:3, 20), log = a$log + rnorm(60, 0.1))
  cmp <- suppressWarnings(ilm_compare(a, b, rule = "log", by = "origin", method = "dm"))
  expect_match(cmp$method[1], "lag 2")
  expect_match(suppressWarnings(ilm_compare(a, b, rule = "log", by = "origin", method = "dm",
                                            lag = 0))$method[1], "lag 0")
})

test_that("comparisons need the same rows and the named columns", {
  s <- scores2()
  expect_error(ilm_compare(s$a), "at least two")
  expect_error(ilm_compare(s$a, s$b[-1, ], rule = "log"), "different numbers of rows")
  expect_error(ilm_compare(s$a, s$b, rule = "crps"), "no 'crps' column")
  expect_error(ilm_compare(s$a, s$b, rule = "log", by = "site"), "no 'site' column")
  b2 <- s$b; b2$group <- rev(b2$group)
  expect_error(ilm_compare(s$a, b2, rule = "log", by = "group"), "disagree on 'group'")
  ## a plain table is read by its rules' column names ...
  expect_identical(ilm_compare(s$a, s$b)$rule, "log")
  ## ... and one with none has to be told
  v <- data.frame(value = s$a$log); w <- data.frame(value = s$b$log)
  expect_error(ilm_compare(v, w), "name the score")
  expect_identical(ilm_compare(v, w, rule = "value")$rule, "value")
})

test_that("the printout says which model is better, by how much, and how sure", {
  s <- scores2()
  expect_output(print(ilm_compare(ar = s$a, iid = s$b, rule = "log", by = "group")),
                "iid's log score is (lower|higher) than ar's by .* clustered by group, 12 clusters; p")
  expect_output(print(ilm_compare(ar = s$a, iid = s$b, rule = "log")), "lower is better")
})

test_that("the fixed-smoothing Diebold-Mariano statistics match outside computations", {
  set.seed(4); x <- as.numeric(stats::arima.sim(list(ar = 0.4), 37)) + 0.3
  ## WPE: the Daniell periodogram at the first m Fourier frequencies, as stats::spec.pgram()
  w <- lmn_dm_wpe(x)
  sp <- stats::spec.pgram(x, taper = 0, detrend = FALSE, demean = TRUE, fast = FALSE, plot = FALSE)
  expect_identical(w$m, 3)
  expect_equal(w$statistic, sqrt(37) * mean(x) / sqrt(mean(sp$spec[1:3])))
  expect_equal(w$p, 2 * stats::pt(-abs(w$statistic), 6))
  ## WCE: Bartlett weights to M, as sandwich's Newey-West estimate with lag M - 1
  skip_if_not_installed("sandwich")
  wc <- lmn_dm_wce(x)
  v <- sandwich::lrvar(x, type = "Newey-West", prewhite = FALSE, adjust = FALSE, lag = wc$M - 1) * 37
  expect_equal(wc$statistic, sqrt(37) * mean(x) / sqrt(v))
  ## Coroneo and Iacone's critical values at T = 20: 2.57 and 2.09 (WCE), 2.78 (WPE, t on 4 df)
  expect_equal(unname(round(lmn_dm_wce(rnorm(20))$cv, 2)), c(2.57, 2.09))
  expect_equal(round(stats::qt(0.975, lmn_dm_wpe(rnorm(20))$df), 2), 2.78)
})

test_that("over origins the default is the fixed-smoothing test, pooled and per horizon", {
  set.seed(7)
  o <- rep(1:30, each = 3)
  a <- data.frame(origin = o, horizon = rep(1:3, 30), crps = rexp(90))
  b <- data.frame(origin = o, horizon = rep(1:3, 30), crps = a$crps + rnorm(90, 0.05, 0.3))
  cmp <- ilm_compare(a, b, rule = "crps", by = "origin")
  expect_equal(cmp$horizon, c(NA, 1:3))
  expect_match(cmp$method, "fixed-smoothing")
  ## the pooled row is the WPE statistic on the origin means, with t on 2m df
  w <- lmn_dm_wpe(as.numeric(tapply(b$crps - a$crps, o, mean)))
  expect_equal(cmp$statistic[1], w$statistic)
  expect_equal(cmp$p_value[1], w$p)
  expect_equal(cmp$df[1], 6)                                 # m = floor(30^(1/3)) = 3
  expect_equal(cmp$difference[1], mean(b$crps - a$crps))
  ## a horizon's row uses that horizon's rows only
  k <- a$horizon == 2
  w2 <- lmn_dm_wpe(as.numeric(tapply(b$crps[k] - a$crps[k], o[k], mean)))
  expect_equal(cmp$statistic[3], w2$statistic)
  expect_equal(cmp$n, c(90, 30, 30, 30))
  expect_error(ilm_compare(a, b, rule = "crps", by = "origin", lag = 1), "applies to method")
  expect_error(ilm_compare(a, b, rule = "crps", by = "row", method = "dm"), "apply to by = \"origin\" only")
  expect_error(ilm_compare(a, b, rule = "crps", by = "horizon", lag = 1), "clustered")
  ## a lag you give applies to the pooled row and every horizon's row
  lg <- suppressWarnings(ilm_compare(a, b, rule = "crps", by = "origin", method = "dm", lag = 1))
  expect_match(lg$method, "lag 1")
  lg0 <- suppressWarnings(ilm_compare(a, b, rule = "crps", by = "origin", method = "dm"))
  expect_equal(sub(".*[(]lag ([0-9]+)[)].*", "\\1", lg0$method), c("2", "0", "1", "2"))
  ## the usual test warns with few origins over several horizons, in the registry's words
  expect_warning(ilm_compare(a, b, rule = "crps", by = "origin", method = "dm"),
                 "30 origins over 3 horizons.*method = \"wpe\"", class = "lmn_dm_few_origins")
  expect_no_warning(ilm_compare(a[a$horizon == 1, ], b[b$horizon == 1, ], rule = "crps",
                                by = "origin", method = "dm"))
  expect_output(print(cmp), "horizon 2")
  expect_output(print(cmp), "pooled over horizons")
})

test_that("the test is given as data beside its description", {
  set.seed(8)
  o <- rep(1:30, each = 3)
  a <- data.frame(origin = o, horizon = rep(1:3, 30), group = rep(1:10, 9), crps = rexp(90))
  b <- a; b$crps <- a$crps + rnorm(90, 0.05, 0.3)
  w <- ilm_compare(a, b, rule = "crps", by = "origin")
  expect_equal(w$test, rep("wpe", 4))
  expect_equal(w$bandwidth, rep(3L, 4))                      # floor(30^(1/3))
  expect_equal(w$lag, rep(NA_integer_, 4))
  expect_equal(w$df, 2 * w$bandwidth)
  d <- suppressWarnings(ilm_compare(a, b, rule = "crps", by = "origin", method = "dm"))
  expect_equal(d$test, rep("dm", 4))
  expect_equal(d$lag, c(2L, 0L, 1L, 2L))                     # pooled, then horizons 1 to 3
  expect_equal(d$bandwidth, rep(NA_integer_, 4))
  expect_equal(ilm_compare(a, b, rule = "crps", by = "group")$test, "clustered")
  expect_equal(ilm_compare(a, b, rule = "crps")$test, "independent")
  expect_true(all(c("method", "test", "bandwidth", "lag") %in% names(w)))
})

test_that("a negative truncated variance is taken as 0 and the null rejected (Diebold and Mariano)", {
  set.seed(11)
  Tn <- 20L
  tt <- seq_len(Tn)
  dd <- 0.3 + (-1)^tt + rnorm(Tn, 0, 0.1)      # strongly negatively autocorrelated
  a <- data.frame(origin = tt, horizon = 1L, log = rnorm(Tn))
  b <- data.frame(origin = tt, horizon = 1L, log = a$log + dd)
  ## the truncated variance at lag 1, written afresh: autocovariances with divisor T
  dbar <- mean(dd)
  g <- function(k) sum((dd[(k + 1):Tn] - dbar) * (dd[1:(Tn - k)] - dbar)) / Tn
  expect_lt(g(0) + 2 * g(1), 0)
  expect_warning(cmp <- ilm_compare(a, b, rule = "log", by = "origin", method = "dm", lag = 1),
                 class = "lmn_dm_negative_variance")
  expect_identical(cmp$se, 0)
  expect_identical(cmp$statistic, Inf)                  # b's scores are higher on average
  expect_identical(cmp$p_value, 0)
  expect_identical(cmp$lag, 1L)
  expect_match(cmp$method, "variance taken as 0", fixed = TRUE)
  expect_true("dm_negative_variance" %in% attr(cmp, "checks")$id)     # as the help says
  ## the same series at lag 0, and a positively autocorrelated one at lag 1, are tested as usual
  seen <- NULL
  keep <- function(cnd) { seen <<- c(seen, class(cnd)); invokeRestart("muffleWarning") }
  r0 <- withCallingHandlers(ilm_compare(a, b, rule = "log", by = "origin", method = "dm", lag = 0),
                            warning = keep)
  dp <- 0.3 + sin(tt / 3) + rnorm(Tn, 0, 0.1)
  bp <- data.frame(origin = tt, horizon = 1L, log = a$log + dp)
  r1 <- withCallingHandlers(ilm_compare(a, bp, rule = "log", by = "origin", method = "dm", lag = 1),
                            warning = keep)
  expect_false("lmn_dm_negative_variance" %in% seen)
  expect_true(is.finite(r0$statistic) && is.finite(r1$statistic))
  expect_no_match(c(r0$method, r1$method), "taken as 0")
})

test_that("a variance of exactly 0 is not a negative one: no warning, the ordinary path", {
  ## Diebold and Mariano's rule is for a negative estimate (1995, Section 1.1)
  seen <- NULL
  keep <- function(cnd) { seen <<- c(seen, class(cnd)); invokeRestart("muffleWarning") }
  res <- withCallingHandlers(lmn_se_dm(rep(0.5, 20), 1:20, 1L), warning = keep)
  expect_false("lmn_dm_negative_variance" %in% seen)
  expect_identical(res$se, 0)
  expect_identical(res$method, "Diebold-Mariano over origins (lag 1)")
})

test_that("the modified test reproduces Harvey, Leybourne and Newbold's Table 1", {
  skip_on_cran()
  ## Harvey, Leybourne and Newbold (1997), Table 1 (p. 285), the fourth row of each h
  ## block: the modified statistic (9) with t(n - 1) critical values, nominal 10%,
  ## two-sided. Their design (p. 284): independent N(0, 1) errors, d = e1^2 - e2^2,
  ## horizon h with truncation h - 1, 10,000 replications. The cells are those where a
  ## negative variance estimate is very rare (under 0.5% of replications, so that the rule
  ## for it moves the size by under half a point), and the comparison rests on eq. (9)
  ## and the t critical values. The paper does not say how its simulations treated a
  ## negative estimate.
  cells <- data.frame(h = c(1, 1, 1, 2, 2, 3), n = c(8, 16, 64, 32, 128, 64),
                      table1 = c(8.4, 9.6, 10.1, 12.2, 10.8, 12.2))
  R <- 10000L
  set.seed(19970281)
  for (i in seq_len(nrow(cells))) {
    h <- cells$h[i]; n <- cells$n[i]
    neg <- 0L
    rej <- vapply(seq_len(R), function(r) {
      d <- rnorm(n)^2 - rnorm(n)^2
      res <- withCallingHandlers(lmn_se_dm(d, seq_len(n), h - 1L),
                                 lmn_dm_negative_variance = function(w) {
                                   neg <<- neg + 1L; invokeRestart("muffleWarning") })
      res$p < 0.10
    }, TRUE)
    size <- 100 * mean(rej)
    ## both percentages carry Monte Carlo error from 10,000 replications
    se <- 100 * sqrt(2 * 0.1 * 0.9 / R)
    expect_lt(neg / R, 0.005)
    expect_lt(abs(size - cells$table1[i]), 3 * se + 0.05)    # + 0.05: Table 1 is rounded
  }
})
