test_that("the split-conformal quantile is the ceiling((n+1)(1-alpha))-th score", {
  s <- c(0.05, 0.3, 0.12, 0.44, 0.2, 0.38, 0.01, 0.27, 0.41)   # 9 scores
  expect_equal(lmn_conformal_quantile(s, 0.1), max(s))          # k = 9
  expect_equal(lmn_conformal_quantile(s[-1], 0.1), Inf)         # 8 scores: k = 9 > 8
  expect_equal(lmn_conformal_quantile(s, 0.5), sort(s)[5])      # k = ceiling(5) = 5
})

test_that("recalibration gives nominal coverage, and rescues too-narrow intervals", {
  set.seed(3)
  cover <- function(narrow) mean(replicate(1500, {
    ## PITs of outcomes under a model whose intervals are right (narrow = 1)
    ## or too narrow by a factor: the outcomes are more spread than it thinks
    u <- stats::pnorm(stats::rnorm(201, 0, narrow))
    lv <- lmn_conformal_levels(u[-201], 0.9)$level
    abs(u[201] - 0.5) <= lv / 2
  }))
  expect_equal(cover(1), 0.9, tolerance = 0.02)
  set.seed(4)
  model_only <- mean(abs(stats::pnorm(stats::rnorm(1e4, 0, 1.6)) - 0.5) <= 0.45)
  expect_lt(model_only, 0.8)
  expect_equal(cover(1.6), 0.9, tolerance = 0.02)
})

test_that("CDF pooling weighs every group the same, whatever its size", {
  set.seed(5)
  g <- rep(c("big", "s1", "s2", "s3"), c(400, 10, 10, 10))
  s <- c(runif(400, 0, 0.1), runif(30, 0.3, 0.5))    # the big group scores low
  q <- lmn_pooled_quantile(s, g, 0.1)
  ## brute force: the average of the four groups' empirical CDFs at q
  Fbar <- function(t) mean(vapply(split(s, g), function(v) mean(v <= t), 0))
  expect_gte(Fbar(q), 0.9)
  expect_lt(Fbar(max(s[s < q])), 0.9)                 # and nothing smaller reaches it
  ## weighting rows instead would let the big group dominate
  expect_lt(lmn_conformal_quantile(s, 0.1), q)
})

test_that("repeated subsampling needs enough groups, and is reproducible", {
  set.seed(6)
  g <- rep(seq_len(8), each = 5); s <- runif(40, 0, 0.5)
  expect_equal(lmn_subsample_quantile(s, g, 0.1), Inf)          # 8 groups: k <= 9
  g <- rep(seq_len(30), each = 5); s <- runif(150, 0, 0.5)
  set.seed(1); before <- .Random.seed
  q <- lmn_subsample_quantile(s, g, 0.1, seed = 2)
  expect_identical(.Random.seed, before)
  expect_identical(q, lmn_subsample_quantile(s, g, 0.1, seed = 2))
  expect_true(q > stats::quantile(s, 0.7) && q <= max(s))
})

test_that("levels come per bin, and say when a bin is too small to bound", {
  set.seed(7)
  u <- runif(90); by <- rep(c("h1", "h2", "h3"), c(40, 40, 10))
  lv <- lmn_conformal_levels(u[1:80], 0.9, by = by[1:80])
  expect_equal(lv$bin, c("h1", "h2"))
  lv2 <- lmn_conformal_levels(u[c(1:10, 81:88)], 0.9, by = c(rep("a", 10), rep("b", 8)))
  expect_equal(lv2$unbounded, c(FALSE, TRUE))
  expect_equal(lv2$level[2], 1)
  expect_error(lmn_conformal_levels(u, 0.9, method = "pooled"), "needs `group`")
  expect_error(lmn_conformal_levels(u, 1.1), "between 0 and 1")
})

test_that("predicted-mean bins are thirds, or fewer when the scores are few", {
  b <- lmn_mean_bins(1:300)
  expect_length(b$cuts, 2L)
  expect_equal(as.integer(b$bin(c(1, 150, 300))), 1:3)
  expect_length(lmn_mean_bins(1:30)$cuts, 0L)          # 30 scores: one bin, not three
})
