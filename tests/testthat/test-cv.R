cv_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(site = factor(rep(sprintf("s%02d", 1:12), each = 10)), t = rep(1:10, 12),
                  x = stats::rnorm(120), size = rep(stats::runif(12, 1, 2), each = 10))
  d$y <- stats::rpois(120, exp(0.5 + 0.3 * d$x + stats::rnorm(12, 0, 0.4)[d$site]))
  d
}
## fitted with a direct call: illume::ilm_refit() re-evaluates the fit's call, which a
## wrapper passing ... would leave unevaluable
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("grouped cross-validation scores every row once, as a manual refit does", {
  d <- cv_data()
  fit <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", verbose = FALSE))
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 4)
  cv <- ilm_cv(fit, d, folds, rules = c("log", "crps"), nsim = 100)
  expect_s3_class(cv, "ilm_cv")
  expect_setequal(cv$row, seq_len(nrow(d)))
  expect_false(anyDuplicated(cv$row) > 0)
  expect_identical(attr(cv, "groups"), "new")               # held-out sites are new units
  ## fold 1 by hand
  fo <- folds[[1]]
  f1 <- suppressMessages(suppressWarnings(illume::ilm_refit(fit, data = d[fo$train, ])))
  p1 <- suppressWarnings(ilm_predict(f1, d[fo$test, ], target = "response", groups = "new",
                                     nsim = 100, seed = 2))
  s1 <- lmn_pred_scores(p1, d$y[fo$test], rules = c("log", "crps"))
  expect_equal(cv$crps[cv$fold == 1], s1$crps)
  expect_equal(cv$log[cv$fold == 1], s1$log)
  expect_output(print(cv), "Mean scores")
})

test_that("rolling folds carry origin and horizon, and say whether they are ex-ante", {
  d <- cv_data(2)
  fit <- quiet(illume::ilm_model(y ~ x + size + (1 | site), data = d, family = "poisson",
                                 ar = illume::ilm_ar1(~ t | site), verbose = FALSE))
  folds <- ilm_folds(d, type = "rolling", time = "t", group = "site", initial = 7, horizon = 2)
  cv <- ilm_cv(fit, d, folds, rules = "crps", nsim = 100)
  expect_true(all(c("origin", "horizon") %in% names(cv)))
  expect_identical(attr(cv, "groups"), "fitted")
  ## x changes over time and was not named known; size is constant within a site
  expect_identical(attr(cv, "evaluation"), "ex-post")
  expect_identical(attr(cv, "unknown_covariates"), "x")
  expect_output(print(cv), "By horizon")
  expect_output(print(cv), "Ex-post: x")
  folds_k <- ilm_folds(d, type = "rolling", time = "t", group = "site", initial = 7, horizon = 2,
                       known = "x")
  expect_identical(attr(ilm_cv(fit, d, folds_k, rules = "crps", nsim = 50), "evaluation"), "ex-ante")
})

test_that("two models' cross-validations go straight into ilm_compare() and ilm_plot_horizon()", {
  d <- cv_data(3)
  fit_ar <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                                    ar = illume::ilm_ar1(~ t | site), verbose = FALSE))
  fit_iid <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson",
                                     verbose = FALSE))
  folds <- ilm_folds(d, type = "rolling", time = "t", group = "site", initial = 6, horizon = 2,
                     known = "x")
  cv_ar <- ilm_cv(fit_ar, d, folds, rules = "crps", nsim = 50)
  cv_iid <- ilm_cv(fit_iid, d, folds, rules = "crps", nsim = 50)
  cmp <- quiet(ilm_compare(ar = cv_ar, iid = cv_iid, by = "origin"))
  expect_s3_class(cmp, "ilm_compare")
  ## the first model is the reference: the pooled row is iid's mean score minus ar's
  expect_equal(cmp$difference[1], mean(cv_iid$crps - cv_ar$crps), tolerance = 1e-12)
  expect_true(all(cmp$lag >= 0, na.rm = TRUE))
  d_h <- ilm_plot_horizon(ar = cv_ar, iid = cv_iid, rule = "crps")
  expect_setequal(unique(d_h$horizon), 1:2)
})

test_that("random folds on a model with a correlation over time warn that they score interpolation", {
  d <- cv_data(3)
  fit <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", ar = illume::ilm_ar1(~ t | site), verbose = FALSE))
  folds <- ilm_folds(d, type = "random", K = 3)
  expect_warning(ilm_cv(fit, d, folds, rules = "crps", nsim = 50),
                 "folds are random.*interpolation.*rolling", class = "lmn_cv_interpolation")
})

test_that("a fold that cannot be refitted is recorded, not fatal", {
  d <- cv_data(4)
  fit <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", verbose = FALSE))
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 3)
  folds[[2]]$train <- integer(0)                            # nothing to fit
  cv <- ilm_cv(fit, d, folds, rules = "crps", nsim = 50)
  expect_true(all(is.na(cv$crps[cv$fold == 2])) && all(!is.na(cv$failed[cv$fold == 2])))
  expect_true(all(!is.na(cv$crps[cv$fold != 2])))
  expect_output(print(cv), "1 fold\\(s\\) failed")
  expect_error(ilm_cv(fit, d[-1, ], folds), "data frame the folds were built on")
})

test_that("cross-validation stops, in the registry's words, when no fold can be refitted", {
  d <- cv_data(5)
  fit <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", verbose = FALSE))
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 3)
  for (k in seq_along(folds)) folds[[k]]$train <- integer(0)   # nothing can be fitted
  expect_error(ilm_cv(fit, d, folds, rules = "crps", nsim = 50),
               "None of the 3 folds could be refitted", class = "lmn_cv_failed")
})

test_that("the results are the same whatever the number of workers", {
  skip_on_cran()
  d <- cv_data(6)
  fit <- quiet(illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", verbose = FALSE))
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 4)
  one <- ilm_cv(fit, d, folds, rules = c("log", "crps"), nsim = 100, ncores = 1)
  two <- ilm_cv(fit, d, folds, rules = c("log", "crps"), nsim = 100, ncores = 2)
  expect_equal(two, one)
})

test_that("workers are told to load limn from source only when it was loaded from source", {
  ## an installed package also has a DESCRIPTION and an R/ directory; loading
  ## that with load_all() gave workers an empty limn under R CMD check
  tmp <- function() { d <- tempfile(); dir.create(d); d }
  inst <- tmp()
  dir.create(file.path(inst, "R")); dir.create(file.path(inst, "Meta"))
  file.create(file.path(inst, "DESCRIPTION")); file.create(file.path(inst, "Meta", "package.rds"))
  expect_null(lmn_source_path(inst))
  src <- tmp()
  dir.create(file.path(src, "R")); file.create(file.path(src, "DESCRIPTION"))
  expect_identical(lmn_source_path(src), src)
  expect_null(lmn_source_path(tmp()))
  unlink(c(inst, src), recursive = TRUE)
  expect_null(lmn_source_path(find.package("stats")))
})

test_that("workers load the same illume, from the same library, as the parent", {
  skip_on_cran()
  ## a worker quietly loading another library's illume would give different
  ## numbers without an error; this makes it fail loudly
  here <- list(version = as.character(utils::packageVersion("illume")),
               path = normalizePath(find.package("illume")),
               limn = as.character(utils::packageVersion("limn")))
  there <- lmn_pool_lapply(1:2, function(i)
    list(version = as.character(utils::packageVersion("illume")),
         path = normalizePath(find.package("illume")),
         limn = as.character(utils::packageVersion("limn"))), ncores = 2)
  for (w in there) expect_identical(w, here)
})
