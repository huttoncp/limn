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

## A model with time only in its formula: a trend, and a random slope in time.
slope_case <- function() {
  set.seed(7)
  d <- expand.grid(t = 1:16, id = factor(letters[1:6]))
  d$y <- 5 + 0.3 * d$t + rnorm(6, 0, 1)[d$id] + rnorm(6, 0, 0.1)[d$id] * d$t +
    rnorm(nrow(d), 0, 0.5)
  train <- d[d$t <= 12, ]
  fit <- suppressMessages(suppressWarnings(illume::ilm_model(
    y ~ t + (1 + t | id), data = train, family = "gaussian", verbose = FALSE)))
  list(d = d, train = train, fit = fit,
       fut = ilm_future(train, time = "t", group = "id", h = 4, response = "y"))
}

test_that("rows from ilm_future() name the time and series without a correlation over time", {
  s <- slope_case()
  p <- suppressWarnings(ilm_forecast(s$fit, s$fut, nsim = 30, seed = 1))
  expect_identical(p$meta$time, "t")
  expect_identical(p$meta$series, "id")
  ## so a scaled error can be read from the history, and it is the one the
  ## scale computed directly gives
  y <- s$d$y[match(paste(s$fut$id, s$fut$t), paste(s$d$id, s$d$t))]
  expect_equal(ilm_scores(p, y, rules = "mase", history = s$train)$mase,
               ilm_scores(p, y, rules = "mase", naive = ilm_naive_scale(
                 s$train, time = "t", group = "id", response = "y"))$mase)
})

test_that("rows changed by merge() or subset() need time and group named", {
  s <- slope_case()
  merged <- merge(s$fut[c("id", "t")], data.frame(id = letters[1:6], z = 1))
  subsetted <- subset(s$fut, t > 13)
  for (nd in list(merged, subsetted)) {
    expect_null(attr(nd, "time"))
    expect_null(suppressWarnings(ilm_forecast(s$fit, nd, nsim = 20))$meta$time)
    p <- suppressWarnings(ilm_forecast(s$fit, nd, nsim = 20, time = "t", group = "id"))
    expect_identical(c(p$meta$time, p$meta$series), c("t", "id"))
  }
  ## named columns take precedence over the ones ilm_future() recorded
  fut2 <- s$fut; fut2$t2 <- fut2$t
  p <- suppressWarnings(ilm_predict(s$fit, fut2, groups = "fitted", nsim = 20, time = "t2"))
  expect_identical(c(p$meta$time, p$meta$series), c("t2", "id"))
})

test_that("time and group name columns of the new rows", {
  s <- slope_case()
  expect_error(ilm_forecast(s$fit, s$fut, nsim = 20, time = s$fut$t), "not the values themselves")
  expect_error(ilm_forecast(s$fit, s$fut, nsim = 20, group = c("id", "t")), "not the values themselves")
  expect_error(ilm_forecast(s$fit, s$fut, nsim = 20, time = "week"), "'week', which is not in `newdata`")
})

## A correlation over time built by name, and the same one built from vectors.
ar_pair <- function() {
  set.seed(3)
  a <- expand.grid(t = 1:12, g = factor(c("a", "b", "c")))
  a$y <- rnorm(nrow(a), 10, 1) + rep(as.numeric(arima.sim(list(ar = 0.5), 12)), 3)
  fit <- function(ar) suppressMessages(suppressWarnings(illume::ilm_model(
    y ~ 1, data = a, family = "gaussian", ar = ar, verbose = FALSE)))
  list(by_name = fit(illume::ilm_ar1(~ t | g)), from_vectors = fit(illume::ilm_ar1(a$t, a$g)),
       nd = data.frame(t = c(13, 14, 13), g = c("a", "a", "b")))
}

test_that("a correlation over time built by name keeps its own columns", {
  f <- ar_pair()
  expect_error(ilm_forecast(f$by_name, cbind(f$nd, tt = f$nd$t), nsim = 20, time = "tt"),
               "ilm_ar1[(]~ t [|] g[)], uses 't'.*Leave `time` out, or give `time = \"t\"`")
  expect_error(ilm_forecast(f$by_name, cbind(f$nd, gg = f$nd$g), nsim = 20, group = "gg"),
               "uses 'g'.*Leave `group` out, or give `group = \"g\"`")
  ## the same names pass, and change nothing
  expect_equal(ilm_forecast(f$by_name, f$nd, target = "link", nsim = 20, seed = 1, time = "t", group = "g"),
               ilm_forecast(f$by_name, f$nd, target = "link", nsim = 20, seed = 1))
})

test_that("rows from ilm_future() built around another time column say so", {
  f <- ar_pair()
  set.seed(3)
  a <- expand.grid(t = 1:12, g = factor(c("a", "b", "c")))
  a$y <- 10; a$week <- a$t + 100
  ## ilm_future() treats the model's own time column as one that varies, and
  ## leaves it missing
  fut <- suppressMessages(ilm_future(a, time = "week", group = "g", h = 2, response = "y"))
  expect_error(ilm_forecast(f$by_name, fut, nsim = 20),
               "record 'week' as their time column.*ilm_ar1[(]~ t [|] g[)], uses 't'.*time = \"t\", or fill in 't'")
  ## filled in, the term's own column is used
  fut$t <- fut$week - 100
  expect_equal(ilm_forecast(f$by_name, fut, target = "link", nsim = 20, seed = 1)$estimate,
               ilm_forecast(f$by_name, fut[c("t", "g")], target = "link", nsim = 20, seed = 1)$estimate)
})

test_that("a correlation over time built from vectors takes its time and series as columns", {
  f <- ar_pair()
  nd <- data.frame(week = f$nd$t, series = f$nd$g)
  expect_error(ilm_forecast(f$from_vectors, nd, nsim = 20), "built from vectors.*give `time` and `group`")
  ref <- ilm_forecast(f$by_name, f$nd, target = "link", nsim = 50, seed = 1)
  p <- ilm_forecast(f$from_vectors, nd, target = "link", nsim = 50, seed = 1,
                    time = "week", group = "series")
  expect_equal(p$estimate, ref$estimate)
  expect_equal(p$origin, ref$origin)
  expect_identical(c(p$meta$time, p$meta$series), c("week", "series"))
  ## the columns follow the rows when `set` repeats them and scenarios stack them
  ps <- ilm_forecast(f$from_vectors, nd, target = "link", nsim = 50, seed = 1,
                     time = "week", group = "series", set = list(z = 1:2))
  expect_equal(ps$estimate, rep(ref$estimate, each = 2))
  pl <- ilm_forecast(f$from_vectors, list(low = nd, high = nd[3:1, ]), target = "link", nsim = 50,
                     seed = 1, time = "week", group = "series")
  expect_equal(pl$estimate, c(ref$estimate, rev(ref$estimate)))
})
