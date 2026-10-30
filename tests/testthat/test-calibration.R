test_that("a correct predictive distribution reads as calibrated", {
  set.seed(4)
  n <- 800; m <- rnorm(n, 1, 0.5)
  y <- rnbinom(n, size = 1.5, mu = exp(m))
  p <- lmn_new_pred(test_dist_nbinom(), list(mu = matrix(m)), list(logsize = log(1.5)))
  cal <- lmn_pred_calibration(p, y)
  expect_s3_class(cal, "ilm_pred_calibration")
  expect_match(cal$verdict, "^Consistent with calibration")
  expect_equal(cal$coverage$coverage, cal$coverage$level, tolerance = 0.05)
  expect_gt(cal$uniformity[["p_value"]], 0.01)
  expect_output(print(cal), "coverage by level")
  ## a coverage is a share to three figures; its standard error an estimate in
  ## percentage points
  expect_output(print(cal), "[0-9]+[.][0-9]%  ?[(]SE [0-9.]+ percentage points[)]")
})

test_that("Poisson on overdispersed counts is too narrow, not biased, and names nbinom", {
  set.seed(4)
  n <- 800; m <- rnorm(n, 1, 0.5)
  y <- rnbinom(n, size = 1.5, mu = exp(m))
  p <- lmn_new_pred(test_dist_poisson(), list(mu = matrix(m)))
  cal <- lmn_pred_calibration(p, y)
  expect_true(any(grepl("^Too narrow", cal$verdict)))
  expect_true(any(grepl("family = \"nbinom\"", cal$verdict)))
  ## the mean is right, so it must not be called biased, even though the
  ## mean PIT moves when only the shape is wrong
  expect_false(any(grepl("^Biased", cal$verdict)))
})

test_that("a biased mean is called biased, in the outcome's units", {
  set.seed(8)
  n <- 600; m <- rnorm(n)
  y <- rnorm(n, m + 0.4)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = 0))
  cal <- lmn_pred_calibration(p, y)
  expect_true(any(grepl("^Biased: outcomes are on average 0\\.[34].* above", cal$verdict)))
})

test_that("binary outcomes are judged by the calibration intercept and slope", {
  set.seed(4)
  n <- 800; pb <- runif(n, 0.05, 0.95)
  mk <- function(yb) lmn_pred_calibration(
    lmn_new_pred(test_dist_binomial(), list(mu = matrix(qlogis(pb)))), yb)
  good <- mk(rbinom(n, 1, pb))
  expect_match(good$verdict, "^Consistent with calibration")
  expect_equal(good$binary$measure, c("intercept", "slope"))
  moderate <- mk(rbinom(n, 1, plogis(1.8 * qlogis(pb))))       # truth more extreme
  expect_true(any(grepl("^Underconfident", moderate$verdict)))
  expect_equal(moderate$binary$estimate[2], 1.8, tolerance = 0.2)
  extreme <- mk(rbinom(n, 1, plogis(0.5 * qlogis(pb))))        # truth more moderate
  expect_true(any(grepl("^Overconfident", extreme$verdict)))
  expect_false(any(grepl("^Too (narrow|wide)", c(moderate$verdict, extreme$verdict))))
})

test_that("calibration needs enough outcomes", {
  p <- test_pred(test_dist_gaussian(), m = rnorm(10), S = 3, par = list(logsd = 0))
  expect_error(lmn_pred_calibration(p, rnorm(10)), "at least 20")
})

test_that("the PIT plot returns its bars and band, and checks its input", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(1)
  u <- runif(400)
  d <- ilm_plot_pit(u, bins = 8)
  expect_equal(nrow(d), 8L)
  expect_equal(sum(d$count), 400L)
  expect_equal(mean(d$density), 1)
  expect_true(all(d$band_lower < 1 & d$band_upper > 1))
  expect_error(ilm_plot_pit(c(0.2, 1.3)), "between 0 and 1")
  expect_error(ilm_plot_pit(NA_real_), "no PITs")
})

## ilm_calibration() on a prediction: limn's method for illume's generic
calib_fits <- function() {
  set.seed(1)
  d <- data.frame(x = rnorm(600))
  d$y <- rnbinom(600, size = 1.5, mu = exp(1 + 0.5 * d$x))
  train <- d[1:300, ]; test <- d[301:600, ]
  fit <- function(fam) suppressMessages(suppressWarnings(
    illume::ilm_model(y ~ x, data = train, family = fam, verbose = FALSE)))
  list(test = test, pois = ilm_predict(fit("poisson"), test, target = "response", nsim = 200),
       nb = ilm_predict(fit("nbinom"), test, target = "response", nsim = 200))
}

test_that("illume's generic dispatches to limn's method for a prediction", {
  f <- calib_fits()
  cal <- illume::ilm_calibration(f$pois, f$test$y)
  expect_s3_class(cal, "ilm_pred_calibration")
  expect_identical(cal, lmn_pred_calibration(f$pois, f$test$y))
  expect_identical(illume::ilm_calibration(f$pois, f$test$y, bins = 8, levels = c(0.8, 0.95), seed = 3),
                   lmn_pred_calibration(f$pois, f$test$y, bins = 8, levels = c(0.8, 0.95), seed = 3))
  expect_equal(cal$coverage$level, c(0.5, 0.8, 0.9, 0.95))
  expect_equal(cal$family, "poisson")
})

test_that("on fitted models, Poisson on overdispersed counts is too narrow and names nbinom", {
  f <- calib_fits()
  pois <- illume::ilm_calibration(f$pois, f$test$y)
  expect_true(any(grepl("^Too narrow", pois$verdict)))
  expect_true(any(grepl(paste0("^Remedy: Overdispersion is the usual cause for counts ",
                               "\\(limn's own reasoning\\): refit with family = \"nbinom\""),
                        pois$verdict)))
  expect_identical(pois$checks$id[1], "calib_narrow_counts")
  nb <- illume::ilm_calibration(f$nb, f$test$y)
  expect_false(any(grepl("^Too (narrow|wide)", nb$verdict)))
  expect_output(print(pois), "Remedy: ")
})

test_that("the method checks its arguments and leaves the caller's random state", {
  f <- calib_fits()
  expect_error(illume::ilm_calibration(f$pois, f$test$y, levels = 1.2), "`levels`")
  expect_error(illume::ilm_calibration(f$pois, f$test$y, levels = numeric()), "`levels`")
  set.seed(9); a <- runif(1)
  set.seed(9); invisible(illume::ilm_calibration(f$pois, f$test$y, seed = 2)); b <- runif(1)
  expect_identical(a, b)
})

## a biased mean comes first, and a reversed slope is named
test_that("a biased mean is reported first, and the shortfall beside it waits on the mean", {
  set.seed(8)
  n <- 600; m <- rnorm(n)
  y <- rnorm(n, m + 0.5)                                # sd right, mean shifted
  cal <- lmn_pred_calibration(lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)),
                                           list(logsd = 0)), y)
  expect_identical(cal$checks$id, c("calib_bias", "calib_narrow_bias"))
  expect_match(cal$verdict[1], "^Biased: outcomes are on average 0[.][45].* above")
  expect_match(cal$verdict[3], "^Too narrow: .*A biased mean alone can cause this")
  expect_match(cal$verdict[4], "^Remedy: Fix the mean first")
  ## no remedy that widens intervals is offered while the mean is off
  expect_false(any(grepl("dispformula|nbinom|reml", cal$verdict)))
  ## with the mean right, the same shortfall gets its dispersion remedy
  y2 <- rnorm(n, m, 1.6)
  cal2 <- lmn_pred_calibration(lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)),
                                            list(logsd = 0)), y2)
  expect_identical(cal2$checks$id, "calib_narrow_dispersion")
  expect_true(any(grepl("dispformula", cal2$verdict)))
})

test_that("a calibration slope below 0 is reversed, not only overconfident", {
  set.seed(4)
  n <- 800; pb <- runif(n, 0.05, 0.95)
  mk <- function(yb) lmn_pred_calibration(
    lmn_new_pred(test_dist_binomial(), list(mu = matrix(qlogis(pb)))), yb)
  rev <- mk(rbinom(n, 1, plogis(-0.8 * qlogis(pb))))     # truth runs the other way
  expect_true("binary_reversed" %in% rev$checks$id)
  expect_false("binary_overconfident" %in% rev$checks$id)
  expect_true(any(grepl("^Reversed: the calibration slope is -.*order the outcomes the wrong way",
                        rev$verdict)))
  expect_true(any(grepl("^Remedy: Use fewer terms", rev$verdict)))
  ## a slope between 0 and 1 is still overconfident
  ext <- mk(rbinom(n, 1, plogis(0.5 * qlogis(pb))))
  expect_true("binary_overconfident" %in% ext$checks$id)
  expect_false("binary_reversed" %in% ext$checks$id)
})
