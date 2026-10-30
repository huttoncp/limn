scaled_setup <- function() {
  set.seed(9)
  hist <- data.frame(site = rep(c("a", "b"), each = 8), month = rep(1:8, 2),
                     y = c(cumsum(rnorm(8)), 10 + cumsum(rnorm(8, 0, 3))))
  rows <- data.frame(site = rep(c("a", "b"), each = 2), month = rep(9:10, 2))
  m <- c(1, 1.5, 12, 13)
  p <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)), list(logsd = 0),
                    rows = rows)
  list(hist = hist, rows = rows, m = m, pred = p, y = c(0.4, 2.2, 15, 9))
}

test_that("scaled errors divide by each series' in-sample naive error", {
  s <- scaled_setup()
  nv <- lmn_naive_scale(s$hist, "month", "site", "y")
  qa <- vapply(c("a", "b"), function(k) mean(abs(diff(s$hist$y[s$hist$site == k]))), 0)
  q2 <- vapply(c("a", "b"), function(k) mean(diff(s$hist$y[s$hist$site == k])^2), 0)
  expect_equal(nv$abs, unname(qa)); expect_equal(nv$sq, unname(q2))
  sc <- lmn_pred_scores(s$pred, s$y, c("mase", "rmsse"), naive = nv)
  k <- match(s$rows$site, c("a", "b"))
  expect_equal(sc$mase, abs(s$y - s$m) / qa[k], tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(sc$rmsse, (s$y - s$m)^2 / q2[k], tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("a seasonal scale uses the seasonal naive error", {
  s <- scaled_setup()
  nv <- lmn_naive_scale(s$hist, "month", "site", "y", period = 4)
  ya <- s$hist$y[s$hist$site == "a"]
  expect_equal(nv$abs[1], mean(abs(diff(ya, lag = 4))))
  expect_identical(attr(nv, "period"), 4L)
})

test_that("scaled errors need a scale for every series", {
  s <- scaled_setup()
  expect_error(lmn_pred_scores(s$pred, s$y, "mase"), "need `naive`")
  nv <- lmn_naive_scale(s$hist[s$hist$site == "a", ], "month", "site", "y")
  expect_error(lmn_pred_scores(s$pred, s$y, "mase", naive = nv), "no in-sample scale for series 'b'")
  flat <- data.frame(month = 1:4, y = 3)
  expect_error(lmn_naive_scale(flat, "month", NULL, "y"), "never change")
  expect_error(lmn_naive_scale(s$hist, "month", "site", "y", period = 8), "too few")
})
