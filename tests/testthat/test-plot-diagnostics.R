horizon_scores <- function() {
  set.seed(1)
  s <- expand.grid(horizon = 1:4, origin = 1:12, site = 1:5)
  shared <- rnorm(12)[s$origin]                       # origins share their surprises
  ar <- data.frame(s, log = 1 + 0.08 * s$horizon + 0.3 * shared + rnorm(nrow(s), 0, 0.3))
  iid <- data.frame(s, log = ar$log + 0.03 * s$horizon + rnorm(nrow(s), 0, 0.1))
  list(ar = ar, iid = iid)
}

test_that("the horizon plot's intervals are clustered by origin", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  s <- horizon_scores()
  d <- ilm_plot_horizon(ar = s$ar, iid = s$iid)
  expect_equal(nrow(d), 8L)
  expect_equal(d$model, rep(c("ar", "iid"), each = 4))
  h2 <- s$ar[s$ar$horizon == 2, ]
  cl <- lmn_se_clustered(h2$log, h2$origin)
  expect_equal(d$se[d$model == "ar" & d$horizon == 2], cl$se)
  expect_equal(d$mean[d$model == "ar" & d$horizon == 2], mean(h2$log))
  expect_equal(d$origins, rep(12L, 8))
  ## without origins the rows are independent
  d1 <- ilm_plot_horizon(s$ar[setdiff(names(s$ar), "origin")])
  expect_equal(d1$se[2], sd(h2$log) / sqrt(nrow(h2)))
  expect_error(ilm_plot_horizon(data.frame(log = 1:3)), "'horizon' column")
})

test_that("the comparison plot's intervals match the comparison's own test", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  s <- horizon_scores()
  cmp <- ilm_compare(ar = s$ar, iid = s$iid, rule = "log", by = "origin")
  d <- ilm_plot_compare(cmp)
  tq <- qt(0.975, cmp$df)
  expect_equal(d$lower, cmp$difference - tq * cmp$se)
  expect_equal(d$label[1], "iid vs ar (log, by origin)")               # pooled over horizons
  expect_equal(d$label[3], "iid vs ar (log, by origin, horizon 2)")
  ## the interval excludes zero exactly when the test rejects at 5%
  expect_identical(d$lower > 0 | d$upper < 0, cmp$p_value < 0.05)
  expect_error(ilm_plot_compare(data.frame()), "ilm_compare")
})

test_that("calibration plots coverage, or reliability for binary outcomes", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  set.seed(4)
  n <- 600; m <- rnorm(n, 1, 0.5)
  cal <- lmn_pred_calibration(lmn_new_pred(test_dist_poisson(), list(mu = matrix(m))),
                              rnbinom(n, size = 1.5, mu = exp(m)))
  d <- ilm_plot_calibration(cal)
  expect_equal(d$level, c(0.5, 0.8, 0.9, 0.95))
  expect_true(all(d$upper < d$level))                 # too narrow at every level
  pb <- runif(n, 0.05, 0.95)
  calb <- lmn_pred_calibration(
    lmn_new_pred(test_dist_binomial(), list(mu = matrix(qlogis(pb)))), rbinom(n, 1, pb))
  r <- ilm_plot_calibration(calb)
  expect_equal(nrow(r), 10L)
  expect_equal(sum(r$n), n)
  expect_true(all(r$lower <= r$observed & r$observed <= r$upper))
  expect_true(mean(r$lower <= r$predicted & r$predicted <= r$upper) >= 0.8)
  expect_error(ilm_plot_calibration(list()), "calibration of a predictive distribution")
})

test_that("reliability bins shrink with few outcomes, and use Wilson intervals", {
  r <- lmn_reliability(seq(0.1, 0.9, length.out = 50), rep(c(0, 1), 25), bins = 10)
  expect_equal(nrow(r), 2L)                          # 50 outcomes: bins of at least 20
  w <- lmn_reliability(rep(0.5, 40), c(rep(1, 10), rep(0, 30)), bins = 1)
  ## Wilson's interval for 10 of 40
  z <- qnorm(0.975); ph <- 0.25; n <- 40
  expect_equal(w$lower, (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n))
})
