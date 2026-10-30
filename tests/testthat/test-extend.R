## The two extension points: a package that fits illume's models another way
## registers methods for ilm_predict_draws() and ilm_cv_refit(), and every verb
## then works on its fits. Here a stub class plays that package: a subclass of
## "ilm_model" whose methods count their calls and hand on to illume's.

ext_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(site = factor(rep(sprintf("s%02d", 1:12), each = 10)), x = stats::rnorm(120))
  d$y <- stats::rpois(120, exp(0.5 + 0.3 * d$x + stats::rnorm(12, 0, 0.4)[d$site]))
  d
}
ext_fit <- function(d) suppressMessages(suppressWarnings(
  illume::ilm_model(y ~ x + (1 | site), data = d, family = "poisson", verbose = FALSE)))

## methods for the stub class, registered on limn's generics; each test
## removes them again with on.exit(drop_stub_methods())
stub_methods <- function(draws = NULL, refit = NULL) {
  calls <- new.env()
  calls$draws <- 0L; calls$refit <- 0L
  if (is.null(draws)) draws <- function(fit, nsim, seed, uncertainty, natural, ...) {
    calls$draws <- calls$draws + 1L
    class(fit) <- setdiff(class(fit), "lmn_stub_fit")
    ilm_predict_draws(fit, nsim = nsim, seed = seed, uncertainty = uncertainty, natural = natural)
  }
  if (is.null(refit)) refit <- function(fit, data, ...) {
    calls$refit <- calls$refit + 1L
    class(fit) <- setdiff(class(fit), "lmn_stub_fit")
    out <- ilm_cv_refit(fit, data)
    class(out) <- c("lmn_stub_fit", class(out))
    out
  }
  ns <- asNamespace("limn")
  registerS3method("ilm_predict_draws", "lmn_stub_fit", draws, envir = ns)
  registerS3method("ilm_cv_refit", "lmn_stub_fit", refit, envir = ns)
  calls
}
drop_stub_methods <- function() {
  tbl <- get(".__S3MethodsTable__.", envir = asNamespace("limn"))
  rm(list = intersect(c("ilm_predict_draws.lmn_stub_fit", "ilm_cv_refit.lmn_stub_fit"),
                      ls(tbl, all.names = TRUE)), envir = tbl)
}
as_stub <- function(fit) { class(fit) <- c("lmn_stub_fit", class(fit)); fit }

test_that("an illume fit's draws through the generic are the engine's own", {
  d <- ext_data(); fit <- ext_fit(d); vc <- illume::ilm_varcorr(fit)
  for (u in c("joint", "fixed", "latent", "none")) {
    a <- ilm_predict_draws(fit, nsim = 20L, seed = 3L, uncertainty = u, natural = TRUE)
    b <- lmn_param_draws(fit, u, 20L, 3L, vc, natural = TRUE)
    expect_identical(a, b)
  }
})

test_that("a fit of another class predicts through its own draws method", {
  d <- ext_data(); fit <- ext_fit(d)
  calls <- stub_methods(); on.exit(drop_stub_methods())
  nd <- d[1:6, ]
  for (u in c("joint", "none")) {
    p0 <- suppressWarnings(ilm_predict(fit, nd, target = "response", groups = "new",
                                       uncertainty = u, nsim = 50))
    p1 <- suppressWarnings(ilm_predict(as_stub(fit), nd, target = "response", groups = "new",
                                       uncertainty = u, nsim = 50))
    expect_identical(as.data.frame(p1), as.data.frame(p0))
  }
  expect_identical(calls$draws, 2L)
  ## and the scoring verbs follow
  p1 <- suppressWarnings(ilm_predict(as_stub(fit), nd, target = "response", groups = "new",
                                     nsim = 50))
  expect_true(all(is.finite(lmn_pred_scores(p1, nd$y, rules = "crps")$crps)))
})

test_that("draws in the wrong shape stop, saying what is wrong", {
  d <- ext_data(); fit <- ext_fit(d)
  short <- function(fit, nsim, seed, uncertainty, natural, ...) {
    class(fit) <- setdiff(class(fit), "lmn_stub_fit")
    dr <- ilm_predict_draws(fit, nsim = nsim, seed = seed, uncertainty = uncertainty,
                            natural = natural)
    dr$draws <- dr$draws[, -1L, drop = FALSE]
    dr
  }
  stub_methods(draws = short); on.exit(drop_stub_methods())
  expect_error(ilm_predict(as_stub(fit), d[1:3, ], target = "mean", groups = "fitted", nsim = 20),
               "`draws` has 19 columns, not 20", class = "lmn_draws_shape")
  stub_methods(draws = function(...) list(draws = "no"))
  expect_error(ilm_predict(as_stub(fit), d[1:3, ], target = "mean", groups = "fitted", nsim = 20),
               "no numeric matrix", class = "lmn_draws_shape")
})

test_that("cross-validation refits a fit of another class through its own method", {
  d <- ext_data(); fit <- ext_fit(d)
  folds <- ilm_folds(d, type = "grouped", group = "site", K = 3)
  cv0 <- ilm_cv(fit, d, folds, rules = "crps", nsim = 50)
  calls <- stub_methods(); on.exit(drop_stub_methods())
  cv1 <- ilm_cv(as_stub(fit), d, folds, rules = "crps", nsim = 50)
  expect_identical(calls$refit, 3L)
  expect_identical(calls$draws, 3L)              # each refitted fold is still the stub class
  expect_identical(cv1$crps, cv0$crps)
})

test_that("workers are told to load the package whose methods a fit dispatches to", {
  d <- ext_data(); fit <- ext_fit(d)
  expect_identical(lmn_method_packages(fit), character(0))   # illume's fit: limn's own methods
  m <- function(fit, data, ...) NULL
  environment(m) <- asNamespace("stats")                    # as if stats had registered it
  stub_methods(refit = m); on.exit(drop_stub_methods())
  expect_identical(lmn_method_packages(as_stub(fit)), "stats")
})
