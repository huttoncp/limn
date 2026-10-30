## ilm_scores() on a prediction, and the two ways to give MASE and RMSSE their
## scale: history = through the prediction's own columns, or
## naive = from ilm_naive_scale(). They give identical scores where both apply.

scores_fit <- function() {
  set.seed(1)
  d <- data.frame(site = factor(rep(c("a", "b", "c"), each = 16)), t = rep(1:16, 3))
  d$y <- rpois(nrow(d), exp(1.5 + ave(rnorm(nrow(d), 0, 0.15), d$site, FUN = cumsum)))
  train <- d[d$t <= 12, ]; test <- d[d$t > 12, ]
  fit <- suppressMessages(suppressWarnings(illume::ilm_model(y ~ 1, data = train, family = "poisson",
                                                             ar = illume::ilm_rw1(~ t | site), verbose = FALSE)))
  list(train = train, test = test, fc = ilm_forecast(fit, test[c("site", "t")], nsim = 100, seed = 1))
}

test_that("a prediction records its response, time and series columns", {
  s <- scores_fit()
  expect_identical(s$fc$meta$response, "y")
  expect_identical(s$fc$meta$time, "t")
  expect_identical(s$fc$meta$series, "site")
})

test_that("history = and naive = give identical scaled errors", {
  s <- scores_fit()
  a <- ilm_scores(s$fc, s$test$y, rules = c("mase", "rmsse"), history = s$train)
  b <- ilm_scores(s$fc, s$test$y, rules = c("mase", "rmsse"),
                  naive = ilm_naive_scale(s$train, time = "t", group = "site", response = "y"))
  expect_identical(a$mase, b$mase)
  expect_identical(a$rmsse, b$rmsse)
  ## with a seasonal lag too
  a4 <- ilm_scores(s$fc, s$test$y, rules = "mase", history = s$train, period = 4L)
  b4 <- ilm_scores(s$fc, s$test$y, rules = "mase",
                   naive = ilm_naive_scale(s$train, time = "t", group = "site", response = "y", period = 4L))
  expect_identical(a4$mase, b4$mase)
  expect_false(identical(a4$mase, a$mase))
  ## the other rules are unchanged by either
  expect_identical(ilm_scores(s$fc, s$test$y, rules = "crps", history = s$train)$crps,
                   ilm_scores(s$fc, s$test$y, rules = "crps")$crps)
  expect_identical(attr(a, "orientation"), "lower is better")
})

test_that("history = says what to do when it cannot work, and one route is asked for", {
  s <- scores_fit()
  expect_error(ilm_scores(s$fc, s$test$y, rules = "mase", history = s$train,
                          naive = ilm_naive_scale(s$train, "t", "site", "y")), "not both")
  ## a prediction with no correlation over time has no time column
  set.seed(2)
  d <- data.frame(g = factor(rep(letters[1:25], each = 4)), x = rnorm(100))
  d$y <- rpois(100, exp(1 + 0.3 * d$x))
  f <- suppressMessages(suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                                                           verbose = FALSE)))
  p <- ilm_predict(f, d[1:4, ], target = "response", groups = "fitted", nsim = 50)
  expect_null(p$meta$time)
  expect_error(ilm_scores(p, d$y[1:4], rules = "mase", history = d),
               "naive = ilm_naive_scale\\(history, time = .*response = \"y\"")
  expect_error(ilm_scores(s$fc, s$test$y, rules = "mase", history = "train"), "data frame")
})
