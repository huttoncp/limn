## limn's results keep the diagnostics they raise, with the values behind the
## words, beside the seed and the draws asked for. The warnings are still
## raised as before.

few_fit <- function() {
  set.seed(2); G <- 8; n <- 10                 # 8 groups: the few-groups warning
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:G), each = n)), x = rnorm(G * n))
  d$y <- rpois(nrow(d), exp(0.5 + 0.4 * d$x + rnorm(G, 0, 0.6)[d$g]))
  suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE))
}

test_that("a filled diagnostic keeps the values its words were made from", {
  d <- lmn_diag("dm_few_origins", T = 12, H = 4)
  expect_identical(d$values[[1]], list(T = 12, H = 4))
  expect_match(d$detail, "12 origins over 4 horizons")
  expect_identical(nrow(lmn_no_checks()), 0L)
  expect_true("values" %in% names(lmn_no_checks()))
})

test_that("a prediction keeps its checks, still warns, and keeps its seed and draws", {
  f <- few_fit()
  nd <- data.frame(g = c("g01", "new"), x = c(0, 1))
  expect_warning(p <- ilm_predict(f, nd, groups = "fitted", nsim = 40, seed = 7), class = "lmn_few_groups_means")
  expect_true(all(c("few_groups_means", "unseen_levels") %in% p$checks$id))
  expect_identical(p$checks$values[[match("unseen_levels", p$checks$id)]], list(n = 1L))
  expect_equal(p$meta$nsim, 40)
  expect_equal(p$meta$seed$value, 7)
  expect_match(p$meta$seed$kind, "Mersenne-Twister")
  ## a total keeps the prediction's checks
  tot <- suppressWarnings(ilm_aggregate(p))
  expect_identical(tot$checks, p$checks)
  ## no checks: an empty frame, not NULL
  x <- suppressWarnings(ilm_predict(f, data.frame(g = "g01", x = 0), groups = "fitted", nsim = 20))
  expect_false("unseen_levels" %in% x$checks$id)
})

test_that("a comparison keeps the checks it raises", {
  set.seed(6); o <- rep(1:20, each = 3)
  a <- data.frame(origin = o, horizon = rep(1:3, 20), log = rnorm(60))
  b <- data.frame(origin = o, horizon = rep(1:3, 20), log = a$log + rnorm(60, 0.1))
  expect_warning(cmp <- ilm_compare(a, b, rule = "log", by = "origin", method = "dm"),
                 class = "lmn_dm_few_origins")
  ck <- attr(cmp, "checks")
  expect_identical(ck$id, "dm_few_origins")
  expect_equal(ck$values[[1]], list(T = 20, H = 3))
  quiet_cmp <- ilm_compare(a, b, rule = "log", by = "origin")          # wpe: no warning
  expect_identical(nrow(attr(quiet_cmp, "checks")), 0L)
})

test_that("cross-validation keeps its own checks and its folds' predictions'", {
  set.seed(3)
  d <- data.frame(site = factor(rep(sprintf("s%02d", 1:8), each = 10)), t = rep(1:10, 8), x = rnorm(80))
  d$y <- rpois(80, exp(0.5 + 0.3 * d$x + rnorm(8, 0, 0.4)[d$site]))
  fit <- suppressMessages(suppressWarnings(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                                                             ar = illume::ilm_ar1(~ t | site), verbose = FALSE)))
  folds <- ilm_folds(d, type = "random", K = 3)
  cv <- suppressWarnings(ilm_cv(fit, d, folds, rules = "crps", nsim = 30))
  ck <- attr(cv, "checks")
  expect_true("cv_interpolation" %in% ck$id)                  # its own, for random folds
  ## grouped folds predict held-out sites as new ones, which with 8 sites
  ## warns in each fold's prediction: kept once on the result
  gf <- ilm_folds(d, type = "grouped", group = "site", K = 4)
  cg <- suppressWarnings(ilm_cv(fit, d, gf, rules = "crps", nsim = 30))
  ckg <- attr(cg, "checks")
  expect_true("few_groups_new" %in% ckg$id)
  expect_false(anyDuplicated(paste(ckg$id, ckg$detail)) > 0)
})
