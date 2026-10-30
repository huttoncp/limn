## The diagnostics registry: one entry per check, from which warnings, errors,
## printed verdicts and the help pages all take their words.

test_that("every entry has the fields, and the help's words fill every placeholder", {
  for (id in names(lmn_diagnostics)) {
    e <- lmn_diagnostics[[id]]
    expect_true(all(c("check", "status", "detail", "cause", "suggestion", "doc") %in% names(e)),
                info = id)
    expect_true(e$status %in% c("OK", "NOTE", "WARN", "FAIL"), info = id)
    expect_true(nzchar(e$detail), info = id)
    if (e$status != "OK")                               # a problem says why, and what to do
      expect_true(nzchar(e$cause) && nzchar(e$suggestion), info = id)
    d <- do.call(lmn_diag, c(list(id), e$doc))
    expect_false(grepl("[{][a-z_]+[}]", lmn_diag_text(d)), info = id)
  }
})

test_that("a warning carries its entry's class and fields", {
  set.seed(3); d <- data.frame(g = factor(rep(letters[1:8], each = 6)), x = rnorm(48))
  d$y <- rpois(48, exp(1 + 0.5 * d$x + rnorm(8, 0, 0.8)[d$g]))
  fit <- suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                                            verbose = FALSE))
  w <- tryCatch(ilm_predict(fit, data.frame(g = "zz", x = 0), target = "response",
                            groups = "fitted", nsim = 20), warning = function(w) w)
  expect_s3_class(w, "lmn_few_groups_new")
  expect_s3_class(w, "lmn_diagnostic")
  expect_identical(w$diagnostic$status, "WARN")
  expect_match(w$diagnostic$detail, "'g' has 8 levels")
  expect_identical(conditionMessage(w), lmn_diag_text(w$diagnostic))
  ## a new level is noted when the prediction is printed, in the registry's words
  p <- suppressWarnings(ilm_predict(fit, data.frame(g = c("a", "zz"), x = 0),
                                    target = "response", groups = "fitted", nsim = 20))
  expect_output(print(p), "note: 1 of the rows name a level the fit has not seen")
})

test_that("the new-unit few-groups check suggests REML only to a fit that didn't use it", {
  set.seed(3); d <- data.frame(g = factor(rep(letters[1:8], each = 6)), x = rnorm(48))
  d$y <- rpois(48, exp(1 + 0.5 * d$x + rnorm(8, 0, 0.8)[d$g]))
  d$z <- d$x + rnorm(8, 0, 0.8)[d$g] + rnorm(48)
  warn_of <- function(fam, reml, target = "response") {
    fit <- suppressWarnings(illume::ilm_model(if (fam == "gaussian") z ~ x + (1 | g) else y ~ x + (1 | g),
                                              data = d, family = fam, reml = reml, verbose = FALSE))
    tryCatch(ilm_predict(fit, data.frame(g = "a", x = 0), target = target, groups = "new", nsim = 20),
             warning = function(w) w)
  }
  for (fam in c("poisson", "gaussian")) {
    w_ml <- warn_of(fam, FALSE); w_reml <- warn_of(fam, TRUE)
    expect_s3_class(w_ml, "lmn_few_groups_new")
    expect_match(conditionMessage(w_ml), "Refit with reml = TRUE", fixed = TRUE)
    expect_s3_class(w_reml, "lmn_few_groups_new_reml")
    expect_no_match(conditionMessage(w_reml), "Refit with reml = TRUE", fixed = TRUE)
    expect_match(conditionMessage(w_reml), "already used restricted maximum likelihood", fixed = TRUE)
    ## the two share what was found and why; a mean's check is the same either way
    expect_identical(w_reml$diagnostic$cause, w_ml$diagnostic$cause)
    expect_s3_class(warn_of(fam, TRUE, target = "mean"), "lmn_few_groups_means")
  }
})

test_that("verdicts and the help are composed from the same entries", {
  d <- lmn_diag("calib_narrow_counts", coverage = "80%", level = 90)
  v <- lmn_verdict_lines(d)
  ## the reference follows what was found and why, not the remedy
  expect_identical(v[1], paste0(paste(d$detail, d$cause), " Reference: ", d$reference, "."))
  expect_match(v[1], "Reference: Gneiting.*pp. 252-253\\.$")
  expect_match(v[2], "^Remedy: Overdispersion is the usual cause for counts \\(limn's own reasoning\\)")
  expect_no_match(v[2], "Reference")
  rd <- lmn_diag_rd(c("trend", "unseen_levels"))
  expect_identical(rd[1], "@section Checks:")
  e <- lmn_diagnostics$trend
  expect_true(any(grepl(lmn_diag_text(do.call(lmn_diag, c(list("trend"), e$doc))), rd, fixed = TRUE)))
})

test_that("every reference names the pages or section that support its entry", {
  ## a check's reference supports what was found and why; anything in an entry that
  ## the source does not give is marked as limn's own reasoning, or is limn's study
  refs <- Filter(Negate(is.null), lapply(lmn_diagnostics, `[[`, "reference"))
  expect_gt(length(refs), 0)
  for (id in names(refs))
    expect(grepl(", (p\\. [0-9]+|pp\\. [0-9]|section [0-9])", refs[[id]]),
           paste0("the reference for '", id, "' names no page or section: ", refs[[id]]))
  ## and the help says how to read them
  expect_match(lmn_diag_rd("trend")[2], "A reference supports what was found and why")
})
