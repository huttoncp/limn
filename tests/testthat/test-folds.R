panel <- function(G = 6, Tn = 10) {
  data.frame(site = rep(sprintf("s%02d", seq_len(G)), each = Tn),
             month = rep(seq_len(Tn), G), y = seq_len(G * Tn))
}

test_that("random folds hold out every row exactly once, in near-equal folds", {
  d <- panel()
  f <- ilm_folds(d, "random", K = 7)
  tests <- unlist(lapply(f, `[[`, "test"))
  expect_setequal(tests, seq_len(nrow(d)))
  expect_length(tests, nrow(d))
  sz <- vapply(f, function(k) length(k$test), 0L)
  expect_lte(max(sz) - min(sz), 1L)
  for (k in f) expect_setequal(k$train, setdiff(seq_len(nrow(d)), k$test))
})

test_that("random folds are reproducible, and leave the caller's RNG alone", {
  d <- panel()
  set.seed(99); before <- .Random.seed
  a <- ilm_folds(d, "random", K = 4, seed = 3)
  expect_identical(.Random.seed, before)
  expect_identical(a[[1]]$test, ilm_folds(d, "random", K = 4, seed = 3)[[1]]$test)
  expect_false(identical(a[[1]]$test, ilm_folds(d, "random", K = 4, seed = 4)[[1]]$test))
})

test_that("grouped folds never split a group between training and test", {
  d <- panel(G = 11)
  f <- ilm_folds(d, "grouped", K = 4, group = "site")
  for (k in f) expect_length(intersect(d$site[k$train], d$site[k$test]), 0L)
  held <- lapply(f, function(k) unique(d$site[k$test]))
  expect_setequal(unlist(held), unique(d$site))       # every group held out once
  expect_length(unlist(held), 11L)
  expect_lte(diff(range(lengths(held))), 1L)          # group counts balanced
})

test_that("blocked folds hold out contiguous time and leave a gap either side", {
  d <- panel(Tn = 12)
  f <- ilm_folds(d, "blocked", K = 3, time = "month", gap = 1)
  expect_equal(lapply(f, function(k) sort(unique(d$month[k$test]))),
               list(1:4, 5:8, 9:12))
  expect_equal(sort(unique(d$month[f[[2]]$train])), c(1:3, 10:12))
  expect_setequal(unlist(lapply(f, `[[`, "test")), seq_len(nrow(d)))
})

test_that("rolling origins train up to the origin and test the next horizon", {
  d <- panel(Tn = 10)
  f <- ilm_folds(d, "rolling", time = "month", initial = 4, horizon = 2)
  expect_length(f, 5L)                                # origins 4..8
  expect_equal(vapply(f, `[[`, 0, "origin"), 4:8)
  for (k in f) {
    o <- k$origin
    expect_equal(sort(unique(d$month[k$train])), seq_len(o))
    expect_equal(sort(unique(d$month[k$test])), o + 1:2)
    expect_equal(k$horizon, d$month[k$test] - o)
  }
})

test_that("a sliding window keeps its size, and a gap pushes the test set out", {
  d <- panel(Tn = 12)
  f <- ilm_folds(d, "rolling", time = "month", initial = 5, horizon = 3,
                 window = "sliding", size = 3, gap = 1, step = 2)
  expect_equal(vapply(f, `[[`, 0, "origin"), c(5, 7))
  for (k in f) {
    o <- k$origin
    expect_equal(sort(unique(d$month[k$train])), (o - 2):o)
    expect_equal(sort(unique(d$month[k$test])), o + 2:4)
    expect_equal(range(k$horizon), c(2, 4))
  }
})

test_that("rolling folds count distinct time points, dates included", {
  d <- data.frame(day = as.Date("2026-01-01") + c(0, 0, 3, 3, 9, 9, 10, 10),
                  y = 1:8)
  f <- ilm_folds(d, "rolling", time = "day", initial = 2, horizon = 1)
  expect_length(f, 2L)
  expect_s3_class(f[[1]]$origin, "Date")
  expect_equal(d$day[f[[1]]$test], rep(as.Date("2026-01-10"), 2))
  expect_equal(f[[2]]$horizon, c(1L, 1L))
})

test_that("impossible requests stop with a reason", {
  d <- panel(G = 3, Tn = 5)
  expect_error(ilm_folds(d, "grouped", K = 4, group = "site"), "more folds than the 3 groups")
  expect_error(ilm_folds(d, "grouped", K = 2), "`group` is needed")
  expect_error(ilm_folds(d, "grouped", K = 2, group = "nope"), "not in `data`")
  expect_error(ilm_folds(d, "rolling", time = "month", initial = 4, horizon = 2),
               "no room")
  expect_error(ilm_folds(d, "rolling", time = "month", initial = 3, window = "sliding",
                         size = 4), "longer than the first training window")
  expect_error(ilm_folds(d, "random", K = 1), "at least 2")
})

test_that("as.data.frame and print describe the folds", {
  d <- panel(Tn = 8)
  f <- ilm_folds(d, "rolling", time = "month", initial = 5, horizon = 2,
                 known = "month")
  df <- as.data.frame(f)
  expect_named(df, c("fold", "row", "set", "origin", "horizon"))
  expect_equal(nrow(df), sum(vapply(f, function(k) length(k$train) + length(k$test), 0L)))
  expect_true(all(is.na(df$horizon[df$set == "train"])))
  expect_output(print(f), "rolling origins")
  expect_output(print(f), "known in advance: month")
  expect_output(print(ilm_folds(d, "rolling", time = "month", initial = 5)), "ex-post")
})
