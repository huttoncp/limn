fan_pred <- function() {
  rows <- expand.grid(month = 7:10, site = c("a", "b"), stringsAsFactors = FALSE)
  m <- 0.3 * rows$month + (rows$site == "b")
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)),
                    list(logsd = matrix(log(0.2 * (rows$month - 6)))),
                    rows = rows)
  list(pred = p, m = m, sd = 0.2 * (rows$month - 6), rows = rows)
}

test_that("the fan chart's bands are each row's own central intervals", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  x <- fan_pred()
  d <- ilm_plot_fan(x$pred, time = "month", group = "site", levels = c(0.5, 0.9))
  o <- order(x$rows$site, x$rows$month)
  expect_equal(d$median, x$m[o], tolerance = 1e-10)
  expect_equal(d$lower_0.9, qnorm(0.05, x$m, x$sd)[o], tolerance = 1e-10)
  expect_equal(d$upper_0.5, qnorm(0.75, x$m, x$sd)[o], tolerance = 1e-10)
  ## the bands widen as the predictive distributions do
  wa <- (d$upper_0.9 - d$lower_0.9)[d$group == "a"]
  expect_true(all(diff(wa) > 0))
})

test_that("history is drawn before the predictions, series by series", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  x <- fan_pred()
  h <- expand.grid(month = 1:6, site = c("a", "b"), stringsAsFactors = FALSE)
  h$y <- rnorm(nrow(h))
  d <- ilm_plot_fan(x$pred, time = "month", group = "site", history = h, response = "y")
  expect_equal(nrow(d), nrow(h) + nrow(x$rows))
  expect_equal(sum(!is.na(d$observed)), nrow(h))
  expect_true(all(is.na(d$median[!is.na(d$observed)])))
  expect_equal(d$time[d$group == "a"], 1:10)          # history, then the forecast
})

test_that("panels can take their own scales, or a given one, through ...", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  x <- fan_pred()
  h <- expand.grid(month = 1:6, site = c("a", "b"), stringsAsFactors = FALSE)
  h$y <- rnorm(nrow(h))
  d0 <- ilm_plot_fan(x$pred, time = "month", group = "site", history = h, response = "y")
  d1 <- ilm_plot_fan(x$pred, time = "month", group = "site", history = h, response = "y",
                     facet.args = list(free = TRUE))
  d2 <- ilm_plot_fan(x$pred, time = "month", group = "site", history = h, response = "y",
                     ylim = c(-5, 10))                  # once an error: ylim was set twice
  expect_identical(d1, d0)
  expect_identical(d2, d0)
  expect_silent(ilm_plot_fan(x$pred, time = "month"))  # one panel, no history
})

test_that("a fan chart needs its columns", {
  x <- fan_pred()
  expect_error(ilm_plot_fan(x$pred, time = "week"), "not in `data`")
  h <- data.frame(month = 1:3, site = "a", y = 1:3)
  expect_error(ilm_plot_fan(x$pred, time = "month", history = h), "give `response`")
  expect_error(ilm_plot_fan(x$pred, time = "month", levels = 1.5), "strictly between")
})
