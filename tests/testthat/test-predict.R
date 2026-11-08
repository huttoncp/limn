## The engine is checked against illume's own numbers here, and against an
## independent implementation (glmmTMB) and exact integrals in studies/.

pois_fit <- function() {
  set.seed(1); G <- 20; n <- 10                 # 20 groups: no few-groups warning (B2)
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:G), each = n)), x = rnorm(G * n))
  b0 <- rnorm(G, 0, 0.6); b1 <- rnorm(G, 0, 0.3)
  d$y <- rpois(nrow(d), exp(0.5 + 0.4 * d$x + b0[d$g] + b1[d$g] * d$x))
  list(d = d, fit = suppressWarnings(illume::ilm_model(y ~ x + (1 + x | g), data = d,
                                                        family = "poisson", verbose = FALSE)))
}
nd3 <- data.frame(g = c("g01", "g07", "zz"), x = c(-1, 0.5, 1))

test_that("at the estimates the engine gives illume's own predictions", {
  x <- pois_fit(); fit <- x$fit
  typ <- ilm_predict(fit, nd3, groups = "typical", nsim = 50)
  pop <- ilm_predict(fit, nd3, groups = "population", nsim = 50)
  expect_equal(typ$estimate, as.numeric(predict(fit, nd3, groups = "typical", type = "response")))
  expect_equal(pop$estimate, as.numeric(predict(fit, nd3, groups = "population", type = "response")))
  lnk <- ilm_predict(fit, nd3, target = "link", groups = "typical", nsim = 50)
  expect_equal(exp(lnk$estimate), typ$estimate)
  ## a known group's mean is exp(X beta-hat + Z b-hat), with ilm_ranef()'s modes
  fitd <- ilm_predict(fit, nd3, groups = "fitted", nsim = 50)
  rf <- illume::ilm_ranef(fit); b <- coef(fit)
  man <- vapply(1:2, function(i) {
    r <- rf[rf$level == nd3$g[i], ]
    exp(b[1] + b[2] * nd3$x[i] + r$mode[r$dim == "(Intercept)"] + r$mode[r$dim == "x"] * nd3$x[i])
  }, 0)
  expect_equal(fitd$estimate[1:2], unname(man))
})

test_that("a fit with smooths predicts as illume fits it", {
  ## illume's design leaves a smooth's unpenalised column unnamed
  set.seed(9); d <- expand.grid(t = 1:20, g = factor(sprintf("g%02d", 1:20)))
  d$x <- rnorm(nrow(d)); d$z <- rnorm(nrow(d))
  d$y <- 2 + sin(d$x) + 0.5 * d$z^2 + rnorm(20, 0, 0.3)[d$g] + rnorm(nrow(d), 0, 0.2)
  fit <- suppressWarnings(illume::ilm_model(y ~ s(x) + s(z) + (1 | g), data = d,
                                           family = "gaussian", verbose = FALSE))
  p <- ilm_predict(fit, d[1:6, ], groups = "fitted", uncertainty = "none", target = "link")
  expect_equal(unname(p$estimate), unname(stats::predict(fit)[1:6]), tolerance = 1e-8)
})

test_that("a random intercept alone enters a known unit's prediction", {
  ## its design column comes unnamed, which once dropped the effect silently
  set.seed(2); d <- data.frame(g = factor(rep(sprintf("u%02d", 1:20), each = 6)), x = rnorm(120))
  d$y <- rnorm(120, 1 + 0.5 * d$x + rnorm(20)[d$g], 0.7)
  fit <- suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                                             verbose = FALSE))
  nd <- data.frame(g = c("u02", "u11"), x = c(0.3, -1))
  p <- ilm_predict(fit, nd, groups = "fitted", uncertainty = "none")
  rf <- illume::ilm_ranef(fit); b <- coef(fit)
  man <- b[1] + b[2] * nd$x + rf$mode[match(nd$g, rf$level)]
  expect_equal(p$estimate, unname(man))
  expect_false(any(p$rows$new_level))
})

test_that("an unseen level under \"fitted\" is predicted as a new unit, and marked", {
  x <- pois_fit()
  p <- ilm_predict(x$fit, nd3, groups = "fitted", nsim = 300)
  expect_identical(p$rows$new_level, c(FALSE, FALSE, TRUE))
  w <- function(i) diff(stats::quantile(p$draws[i, ], c(0.05, 0.95)))
  expect_gt(w(3), 2 * w(2))                     # a new unit is far less certain
  pn <- ilm_predict(x$fit, nd3[3, ], groups = "new", nsim = 300)
  expect_equal(unname(stats::quantile(pn$draws[1, ], 0.5)), unname(stats::quantile(p$draws[3, ], 0.5)),
               tolerance = 0.25)
})

test_that("rows naming the same new unit share its effect within a draw", {
  x <- pois_fit()
  nd <- data.frame(g = c("zz", "zz", "yy"), x = c(0, 0, 0))
  p <- ilm_predict(x$fit, nd, target = "link", groups = "new", nsim = 200)
  expect_equal(p$draws[1, ], p$draws[2, ])       # same unit, same draw
  expect_false(isTRUE(all.equal(p$draws[1, ], p$draws[3, ])))
})

test_that("groups must be named, for every term, from the four choices", {
  x <- pois_fit()
  expect_error(ilm_predict(x$fit, nd3), "`groups` is needed")
  expect_error(ilm_predict(x$fit, nd3, groups = "average"), "must be one of")
  expect_error(ilm_predict(x$fit, nd3, groups = list(site = "new")), "no term 'site'")
  expect_s3_class(ilm_predict(x$fit, nd3, groups = list(g = "typical"), nsim = 20), "ilm_pred")
  ## a model with no random terms needs no `groups`; illume's ilm_draws() does
  ## not yet draw for a fit without random effects (reported), so until it
  ## does this part is skipped
  f0 <- illume::ilm_model(y ~ x, data = x$d, family = "poisson", verbose = FALSE)
  skip_if(inherits(try(illume::ilm_draws(f0, 2L), silent = TRUE), "try-error"),
          "illume::ilm_draws() cannot yet draw for a fit without random effects")
  p0 <- ilm_predict(f0, nd3, nsim = 20)
  expect_equal(p0$estimate, as.numeric(predict(f0, nd3, type = "response")))
})

test_that("scenarios stack, and `set` repeats every row per combination", {
  x <- pois_fit()
  sc <- ilm_predict(x$fit, list(low = transform(nd3, x = -1), high = transform(nd3, x = 1)),
                    groups = "typical", nsim = 20)
  expect_equal(sc$rows$.scenario, rep(c("low", "high"), each = 3))
  st <- ilm_predict(x$fit, nd3[1:2, ], groups = "typical", set = list(x = c(-1, 0, 1)), nsim = 20)
  expect_equal(nrow(st$rows), 6L)
  expect_equal(st$rows$x, rep(c(-1, 0, 1), 2))
})

test_that("plug-in predictions carry no parameter uncertainty", {
  x <- pois_fit()
  p <- ilm_predict(x$fit, nd3[1:2, ], groups = "fitted", uncertainty = "none")
  expect_equal(ncol(p$draws), 1L)
  df <- as.data.frame(p)
  expect_equal(df$lower_0.9, df$estimate); expect_equal(df$upper_0.9, df$estimate)
  ## but a new unit's effect is still drawn, over repeated copies of the estimates
  pn <- ilm_predict(x$fit, nd3[3, ], groups = "new", uncertainty = "none", nsim = 200)
  expect_equal(ncol(pn$draws), 200L)
  expect_gt(diff(range(pn$draws)), 0)
  pf <- ilm_predict(x$fit, nd3[1:2, ], groups = "fitted", uncertainty = "fixed", nsim = 100)
  expect_gt(stats::sd(pf$draws[1, ]), 0)
})

test_that("response predictions work with the scoring verbs; means are refused", {
  x <- pois_fit()
  r <- ilm_predict(x$fit, nd3, target = "response", groups = "fitted", nsim = 100)
  expect_identical(r$kind, "response")
  u <- ilm_pit(r, c(1, 3, 2))
  expect_true(all(u >= 0 & u <= 1))
  expect_named(lmn_pred_scores(r, c(1, 3, 2), c("log", "crps")),
               c(names(r$rows), "log", "crps"))
  m <- ilm_predict(x$fit, nd3, groups = "fitted", nsim = 50)
  expect_error(ilm_pit(m, c(1, 3, 2)), "prediction of the mean")
  expect_error(lmn_pred_scores(m, c(1, 3, 2), "log"), "target = \"response\"")
})

test_that("a correlation over time enters at its fitted cells, and a walk starts at zero", {
  set.seed(3); G <- 5; n <- 12
  d <- data.frame(g = factor(rep(sprintf("s%d", 1:G), each = n)), t = rep(1:n, G), x = rnorm(G * n))
  lvl <- unlist(lapply(1:G, function(k) cumsum(rnorm(n, 0, 0.3))))
  d$y <- rnorm(nrow(d), 1 + 0.5 * d$x + lvl, 0.3)
  fit <- suppressWarnings(illume::ilm_model(y ~ x, data = d, family = "gaussian",
                                            ar = illume::ilm_rw1(~ t | g), verbose = FALSE))
  nd <- data.frame(g = c("s1", "s1", "s3"), t = c(1, 6, 9), x = 0)
  p <- ilm_predict(fit, nd, target = "link", groups = "fitted", nsim = 50)
  rf <- illume::ilm_ranef(fit); rf <- rf[rf$type == "rw1", ]
  cell <- function(g, t) rf$mode[as.character(rf$level) == g & rf$time == t]
  expect_equal(p$estimate, coef(fit)[[1]] + c(cell("s1", 1), cell("s1", 6), cell("s3", 9)))
  expect_equal(cell("s1", 1), 0)                                  # the anchor
  ## two steps past the last cell, a walk's forecast is centred on that cell,
  ## with the variance of two steps (the parameters held at their estimates)
  pf <- ilm_predict(fit, data.frame(g = "s1", t = 14, x = 0), target = "link",
                    groups = "fitted", uncertainty = "none", nsim = 4000)
  q <- as.numeric(illume::ilm_varcorr(fit)$latent$var_per_time)
  expect_lt(abs(mean(pf$draws) - coef(fit)[[1]] - cell("s1", 12)), 4 * sqrt(2 * q / 4000))
  expect_equal(stats::sd(pf$draws), sqrt(2 * q), tolerance = 0.05)
  ## a typical series has no latent part; the population mean integrates it
  pt <- ilm_predict(fit, nd, target = "link", groups = "typical", nsim = 20)
  expect_equal(pt$estimate, rep(coef(fit)[[1]], 3))
})

test_that("zero-inflated, dispersion-model and ordinal fits predict", {
  set.seed(4); G <- 20; n <- 12
  d <- data.frame(g = factor(rep(1:G, each = n)), x = rnorm(G * n))
  d$y <- rpois(nrow(d), exp(0.8 + 0.4 * d$x + rnorm(G, 0, 0.4)[d$g]))
  d$y[runif(nrow(d)) < 0.25] <- 0L
  fz <- suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                                           ziformula = ~ 1, verbose = FALSE))
  nd <- data.frame(g = c("1", "2"), x = c(0, 1))
  pz <- ilm_predict(fz, nd, groups = "typical", nsim = 30)
  expect_equal(pz$estimate, as.numeric(predict(fz, nd, groups = "typical", type = "response")))
  rz <- ilm_predict(fz, nd, target = "response", groups = "fitted", nsim = 30)
  expect_gt(lmn_pred_cdf(rz, c(0, 0))[1], stats::ppois(0, pz$estimate[1]))   # the extra zeros
  d$z <- rnorm(nrow(d), 1 + 0.5 * d$x, exp(0.3 * d$x))
  fd <- suppressWarnings(suppressMessages(illume::ilm_model(z ~ x + (1 | g), data = d,
                                           family = "gaussian", dispformula = ~ x, verbose = FALSE)))
  rd <- ilm_predict(fd, data.frame(g = "1", x = c(-1.5, 1.5)), target = "response",
                    groups = "typical", uncertainty = "none")
  wd <- with(as.data.frame(rd), upper_0.9 - lower_0.9)
  expect_gt(wd[2], 2 * wd[1])                                     # the spread grows with x
  d$o <- cut(d$z, 3, labels = c("lo", "mid", "hi"), ordered_result = TRUE)
  fo <- suppressWarnings(suppressMessages(illume::ilm_model(o ~ x + (1 | g), data = d,
                                           family = "ordinal", verbose = FALSE)))
  ro <- ilm_predict(fo, data.frame(g = "1", x = c(-2, 2)), target = "response",
                    groups = "typical", nsim = 30)
  ## an ordered outcome is summarised in its categories, with no mean
  so <- as.data.frame(ro)
  expect_s3_class(so$estimate, "ordered")
  expect_true(so$estimate[1] < so$estimate[2])
  expect_false("mean" %in% names(so))
  expect_error(ilm_predict(fo, nd, groups = "typical"), "probability for each category")
})

test_that("variances held under \"fixed\" are the estimates' in every draw", {
  x <- pois_fit(); vc <- illume::ilm_varcorr(x$fit)
  dr <- lmn_param_draws(x$fit, "fixed", 20L, 1L, vc, natural = TRUE)
  own <- illume::ilm_draws(x$fit, nsim = 20L, seed = 1L, given = "theta")
  expect_equal(dr$natural$re$g, own$natural$re$g, ignore_attr = TRUE)
  expect_equal(dr$draws, own$draws)
  ## and are not computed when no unit's effects are drawn
  expect_null(lmn_param_draws(x$fit, "joint", 20L, 1L, vc, natural = FALSE)$natural$re)
})

test_that("variance draws that are not finite stop, naming the way round them", {
  vc <- list(re = list(g = matrix(1e-9, 1, 1)))
  dr <- list(draws = matrix(0, 2, 3),
             natural = list(re = list(g = array(c(0.1, Inf, 0.2), c(1, 1, 3)))))
  expect_error(lmn_check_natural(dr, vc, "g"), "1 of 3 draws .*uncertainty = \"fixed\"")
  dr$natural$re$g[] <- 0.1
  expect_true(lmn_check_natural(dr, vc, "g"))
})

test_that("terms built from a variable are rebuilt on the new rows", {
  ## illume's grids lose the variable behind log(x) and the like; the engine
  ## builds designs from the new rows, as an lm() fit's own terms do
  set.seed(4); G <- 20; n <- 12
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:G), each = n)), x = rexp(G * n) + 0.2,
                  z = rnorm(G * n), w = runif(G * n, 0, 10), t = rep(1:n, G))
  d$y <- 1 + 0.8 * log(d$x) - 0.3 * d$z^2 + sin(d$w / 2) + sin(2 * pi * d$t / 12) +
    rnorm(G, 0, 0.7)[d$g] + rnorm(G * n, 0, 0.5)
  rhs <- ~ log(x) + poly(z, 2) + splines::ns(w, df = 3) + ilm_fourier(t, 12, K = 2)
  fit <- suppressWarnings(illume::ilm_model(update(rhs, y ~ . + (1 | g)), data = d,
                                            family = "gaussian", verbose = FALSE))
  nd <- data.frame(g = c("g02", "g09", "zz"), x = c(0.5, 2, 4.5), z = c(-1.5, 0.2, 2.1),
                   w = c(0.5, 5, 9.7), t = c(3, 13, 30))     # t = 30 is past the data
  tt <- delete.response(terms(lm(update(rhs, y ~ .), data = d)))
  X <- model.matrix(tt, model.frame(tt, nd))
  p <- ilm_predict(fit, nd, groups = "typical", uncertainty = "none")
  expect_equal(p$estimate, unname(drop(X %*% coef(fit))))
})

test_that("few groups warn: for means always, for a new observation of a new unit", {
  set.seed(3); d <- data.frame(g = factor(rep(letters[1:8], each = 6)), x = rnorm(48))
  d$y <- rpois(48, exp(1 + 0.5 * d$x + rnorm(8, 0, 0.8)[d$g]))
  fit <- suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                                            verbose = FALSE))
  nd <- data.frame(g = c("a", "zz"), x = 0)
  resp <- function(...) ilm_predict(fit, target = "response", nsim = 20, ...)
  expect_warning(resp(nd, groups = "new"), "8 levels.*new unit's observation.*reml = TRUE")
  expect_warning(resp(nd, groups = "fitted"), "8 levels")              # "zz" is a new unit
  expect_warning(resp(nd, groups = "population"), "8 levels")
  expect_no_warning(resp(nd[1, ], groups = "fitted"))
  expect_no_warning(resp(nd, groups = "typical"))
  ## a mean's interval falls short whichever unit it is for
  expect_warning(ilm_predict(fit, nd[1, ], groups = "fitted", nsim = 20), "means.*z-intervals")
  expect_warning(ilm_predict(fit, nd, groups = "typical", target = "link", nsim = 20),
                 "z-intervals")
})

test_that("draws that could not be computed stop, naming the way round them", {
  vals <- matrix(c(1, 2, NaN, 4, Inf, 6), 2)
  lpd <- list(v_int = matrix(c(0, 0, 190, 0, 0, 0), 2))
  expect_error(lmn_check_values(vals, lpd, "mean"),
               "1 of 3 draws of the mean .*variance reaches 190.*uncertainty = \"fixed\"")
  expect_true(lmn_check_values(vals[, -2], lpd, "mean"))      # Inf alone is kept
})

test_that("a mean's summary gives the value at the estimates and the expected value", {
  x <- pois_fit()
  p <- ilm_predict(x$fit, nd3, groups = "fitted", nsim = 200)
  df <- as.data.frame(p)
  expect_equal(df$estimate, p$estimate)
  expect_equal(df$mean, rowMeans(p$draws))
  ## a plug-in prediction for seen units has one draw, the estimates, so the two agree
  p0 <- ilm_predict(x$fit, nd3[1:2, ], groups = "fitted", uncertainty = "none")
  expect_equal(as.data.frame(p0)$mean, p0$estimate)
})

test_that("covariates given as draws mix into the prediction, block by block", {
  set.seed(3)
  d <- data.frame(x = rnorm(80)); d$y <- 1 + 0.8 * d$x + rnorm(80, 0, 0.5)
  fit <- illume::ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  C <- 40; xs <- rnorm(C, 1, 0.6)
  nd <- data.frame(.draw = rep(seq_len(C), each = 2), row = rep(1:2, C),
                   x = c(rbind(xs, 0.5)))                 # row 1 drawn, row 2 known
  expect_message(p <- ilm_predict(fit, nd, target = "response", uncertainty = "latent", nsim = 50),
                      "rounded up from 50 to 80, a multiple of the 40")
  expect_identical(p$meta$covariate_draws, 40L)
  expect_identical(p$meta$drawn, "x")
  expect_equal(p$rows$x, c(NA, 0.5)); expect_equal(p$rows$row, 1:2)
  ## with every parameter held, row 1 is the mixture over the covariate draws, exactly
  dr <- illume::ilm_draws(fit, nsim = 1, given = "parameters", natural = FALSE)
  b <- dr$mode[dr$map$block == "beta"]; s <- exp(dr$mode[dr$map$block == "logdisp"][[1]])
  q <- lmn_pred_quantile(p, 0.9)
  expect_equal(mean(stats::pnorm(q[1], b[1] + b[2] * xs, s)), 0.9, tolerance = 1e-8)
  expect_equal(q[2], stats::qnorm(0.9, b[1] + b[2] * 0.5, s), tolerance = 1e-8)
  ## a drawn row's estimate is its predictive median; a known row's is the plug-in
  m <- ilm_predict(fit, nd, uncertainty = "latent", nsim = 80)
  expect_equal(m$estimate[1], stats::median(b[1] + b[2] * rep(xs, 2)))
  expect_equal(m$estimate[2], unname(b[1] + b[2] * 0.5))
})

test_that("identical covariate blocks give the prediction without them", {
  x <- pois_fit()
  nd <- nd3[1:2, ]
  blocks <- cbind(.draw = rep(1:4, each = 2), nd[rep(1:2, 4), ])
  a <- ilm_predict(x$fit, blocks, groups = "fitted", nsim = 40, seed = 3)
  b <- ilm_predict(x$fit, nd, groups = "fitted", nsim = 40, seed = 3)
  expect_equal(a$draws, b$draws, ignore_attr = TRUE)
  expect_equal(a$estimate, b$estimate)
  expect_length(a$meta$drawn, 0L)
})

test_that("covariate blocks must list the same rows", {
  x <- pois_fit()
  nd <- nd3[1:2, ]
  bad <- cbind(.draw = c(1, 1, 2), nd[c(1, 2, 1), ])
  expect_error(ilm_predict(x$fit, bad, groups = "fitted"), "block 2 has 1 rows and block 1 has 2")
  swapped <- cbind(.draw = rep(1:2, each = 2), nd[c(1, 2, 2, 1), ])
  expect_error(ilm_predict(x$fit, swapped, groups = "fitted"),
               "same rows in the same order.*block 2 differs from block 1 in the units of 'g'")
})

test_that("a new observation whose family overflows in some draws stops in the registry's words", {
  x <- pois_fit()
  p <- ilm_predict(x$fit, nd3, target = "response", groups = "fitted", nsim = 20)
  expect_true(lmn_check_family(p))
  p$lp$mu[2, 3] <- 1e6                                  # exp() overflows: no quantiles
  expect_error(lmn_check_family(p), "1 of 20 draws of a new observation's distribution.*uncertainty = \"latent\"",
               class = "lmn_nan_family")
})

test_that("a family whose quantile function stops on a NaN parameter still stops in the registry's words", {
  set.seed(4); d <- data.frame(g = factor(rep(sprintf("g%02d", 1:24), each = 6)), x = rnorm(144))
  d$y <- ifelse(runif(144) < 0.3, 0L, rpois(144, exp(0.4 + 0.5 * d$x + rnorm(24, 0, 0.4)[d$g])))
  f <- suppressWarnings(suppressMessages(illume::ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                                                           ziformula = ~ 1, verbose = FALSE)))
  p <- suppressWarnings(ilm_predict(f, d[1:4, c("g", "x")], target = "response", groups = "fitted", nsim = 20))
  expect_true(lmn_check_family(p))
  zi <- setdiff(names(p$lp), "mu")[1]
  p$lp[[zi]][, c(3, 11)] <- NaN                         # the zero part's quantile stops, not NaN
  expect_error(lmn_check_family(p), "2 of 20 draws of a new observation's distribution",
               class = "lmn_nan_family")
})

test_that("an ordered outcome's draw that cannot be computed is caught too (its quantile is NA, not NaN)", {
  set.seed(5); d <- data.frame(g = factor(rep(sprintf("g%02d", 1:24), each = 6)), x = rnorm(144))
  d$y <- factor(cut(0.5 * d$x + rnorm(24, 0, 0.4)[d$g] + rnorm(144), c(-Inf, -0.5, 0.3, Inf),
                    labels = c("lo", "mid", "hi")), ordered = TRUE)
  f <- suppressWarnings(suppressMessages(illume::ilm_model(y ~ x + (1 | g), data = d, family = "ordinal",
                                                           verbose = FALSE)))
  p <- suppressWarnings(ilm_predict(f, d[1:4, c("g", "x")], target = "response", groups = "fitted", nsim = 20))
  expect_true(lmn_check_family(p))
  p$lp$mu[, c(3, 11)] <- NaN
  expect_error(lmn_check_family(p), "2 of 20 draws of a new observation's distribution",
               class = "lmn_nan_family")
})
