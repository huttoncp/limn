## A prediction names its time and series columns and keeps, for each row
## ahead of a series' data, the time it was forecast from.

test_that("a forecast keeps its time and series columns and each row's origin", {
  set.seed(5); G <- 3; Tn <- 12
  d <- expand.grid(t = 1:Tn, id = factor(letters[1:G]))
  d <- d[!(d$id == "c" & d$t > 10), ]                    # series c ends at 10
  d$y <- 2 + rnorm(nrow(d), 0, 0.5)
  d$y[d$id == "b" & d$t == 12] <- NA                  # b last seen at 11
  fit <- suppressWarnings(illume::ilm_model(y ~ 1, data = d, family = "gaussian",
                                            ar = illume::ilm_ar1(~ t | id), verbose = FALSE))
  nd <- data.frame(id = c("a", "a", "c", "b", "b", "new"), t = c(13, 15, 12, 5, 12, 13))
  p <- suppressWarnings(ilm_forecast(fit, nd, nsim = 30, seed = 1))
  expect_identical(p$meta$time, "t")
  expect_identical(p$meta$series, "id")
  ## ahead of a seen series' data: its last observed time, even on the fit's
  ## own grid of cells (c at 12, b at 12); within its data, or a new series: NA
  expect_equal(p$origin, c(12, 12, 10, NA, 11, NA))
  ## a subset keeps its rows' origins
  expect_equal(lmn_subset_pred(p, c(3, 1))$origin, c(10, 12))
})

test_that("a model without a latent process has no forecast metadata", {
  set.seed(2)
  d <- data.frame(g = factor(rep(letters[1:10], each = 6)), x = rnorm(60))
  d$y <- rpois(60, exp(0.3 + 0.5 * d$x))
  fit <- suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE))
  p <- suppressWarnings(ilm_predict(fit, d[1:4, ], groups = "fitted", nsim = 20, seed = 1))
  expect_null(p$meta$time)
  expect_null(p$meta$series)
  expect_equal(p$origin, rep(NA_real_, 4))
})
