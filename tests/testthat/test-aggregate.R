## Totals from joint draws: exact where a sum stays in its family, simulated
## otherwise. Checked against sums of each draw's own distribution done by hand.

agg_fit <- function(family, seed = 1) {
  set.seed(seed); G <- 20; n <- 6
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:G), each = n)), x = rnorm(G * n))
  eta <- 0.5 + 0.4 * d$x + rnorm(G, 0, 0.5)[d$g]
  d$y <- switch(family, poisson = rpois(nrow(d), exp(eta)),
                nbinom = rnbinom(nrow(d), mu = exp(eta), size = 2),
                gaussian = rnorm(nrow(d), eta, 0.8))
  suppressMessages(suppressWarnings(illume::ilm_model(y ~ x + (1 | g), data = d, family = family,
                                                      verbose = FALSE)))
}
agg_rows <- data.frame(g = rep(c("g01", "g02"), each = 3), x = c(-1, 0, 1, 0.5, 0, -0.5))
agg_groups <- rep(1:2, each = 3)
agg_out <- data.frame(g = c("g01", "g02"))

test_that("a Poisson total is, draw by draw, the Poisson of the summed means", {
  p <- ilm_predict(agg_fit("poisson"), agg_rows, target = "response", groups = "fitted", nsim = 200)
  a <- lmn_aggregate(p, agg_groups, agg_out)
  expect_identical(a$kind, "response")
  expect_match(a$meta$aggregate, "exact")
  mu <- rowsum(exp(p$lp$mu), agg_groups)                         # 2 x 200
  for (q in c(3, 8, 15))
    expect_equal(lmn_pred_cdf(a, c(q, q)), unname(rowMeans(stats::ppois(q, mu))), tolerance = 1e-12)
})

test_that("a gaussian total adds the means and the variances within each draw", {
  p <- ilm_predict(agg_fit("gaussian"), agg_rows, target = "response", groups = "fitted", nsim = 200)
  a <- lmn_aggregate(p, agg_groups, agg_out)
  expect_match(a$meta$aggregate, "exact")
  m <- rowsum(p$lp$mu, agg_groups)
  sd <- sqrt(3) * matrix(exp(p$par$logdisp), 2, 200, byrow = TRUE)  # three rows, one SD a draw
  for (q in c(0, 2, 5))
    expect_equal(lmn_pred_cdf(a, c(q, q)), unname(rowMeans(stats::pnorm(q, m, sd))), tolerance = 1e-12)
})

test_that("other families' totals are summed from simulated outcomes, and add up", {
  p <- ilm_predict(agg_fit("nbinom"), agg_rows, target = "response", groups = "fitted", nsim = 500)
  a <- lmn_aggregate(p, agg_groups, agg_out, outcomes = 8L, seed = 4)
  expect_identical(a$kind, "sample")
  expect_match(a$meta$aggregate, "simulated")
  expect_true(all(a$draws == round(a$draws)))
  ## the totals' mean is the sum of the rows' means, within Monte Carlo error
  want <- as.numeric(rowsum(lmn_pred_mean(p), agg_groups))
  se <- apply(a$draws, 1L, stats::sd) / sqrt(ncol(a$draws))
  expect_true(all(abs(lmn_pred_mean(a) - want) < 4 * se))
  ## the same outcomes summed at two levels agree sample by sample
  all_one <- lmn_aggregate(p, rep(1L, 6), data.frame(g = "all"), outcomes = 8L, seed = 4)
  expect_equal(colSums(a$draws), all_one$draws[1, ])
  ## draws of means sum exactly, with their estimates
  pm <- ilm_predict(agg_fit("nbinom"), agg_rows, groups = "fitted", nsim = 50)
  am <- lmn_aggregate(pm, agg_groups, agg_out)
  expect_equal(am$draws, unname(rowsum(pm$draws, agg_groups)), ignore_attr = TRUE)
  expect_equal(am$estimate, as.numeric(rowsum(pm$estimate, agg_groups)))
  expect_error(lmn_aggregate(ilm_predict(agg_fit("nbinom"), agg_rows, target = "link",
                                         groups = "fitted", nsim = 20), agg_groups, agg_out),
               "no meaning")
})

test_that("samples are scored by the fair CRPS, and have no log score", {
  set.seed(8); X <- matrix(rpois(2 * 60, 5), 2); y <- c(4, 9)
  sp <- lmn_new_sample_pred(X, discrete = TRUE)
  brute <- vapply(1:2, function(i) mean(abs(X[i, ] - y[i])) -
                    sum(abs(outer(X[i, ], X[i, ], "-"))) / (2 * 60 * 59), 0)
  expect_equal(lmn_crps_sample(sp, y), brute, tolerance = 1e-12)
  skip_if_not_installed("scoringRules")
  ## scoringRules scores the samples' empirical distribution; the fair estimator
  ## differs from it by the pairs term's 1 / (M - 1) in place of 1 / M
  edf <- scoringRules::crps_sample(y, X)
  pairs <- vapply(1:2, function(i) sum(abs(outer(X[i, ], X[i, ], "-"))), 0)
  expect_equal(lmn_crps_sample(sp, y), edf - pairs / (2 * 60^2 * 59), tolerance = 1e-12)
  expect_error(lmn_pred_scores(sp, y, rules = "log"), "no density|has none")
  ## the PIT of a whole-number sample is randomised between F(y - 1) and F(y)
  u <- ilm_pit(sp, y, seed = 1)
  expect_true(all(u >= lmn_pred_cdf(sp, y - 1) & u <= lmn_pred_cdf(sp, y)))
})

test_that("ilm_aggregate() totals within `by`, and never across scenarios", {
  fit <- agg_fit("poisson")
  nd <- list(low = data.frame(g = c("g01", "g01", "g02"), x = -1),
             high = data.frame(g = c("g01", "g01", "g02"), x = 1))
  p <- ilm_predict(fit, nd, target = "response", groups = "fitted", nsim = 100)
  a <- ilm_aggregate(p, by = "g")
  expect_equal(a$rows$n_rows, c(2L, 1L, 2L, 1L))          # two scenarios by two groups
  expect_setequal(names(a$rows), c(".scenario", "g", "n_rows"))
  expect_equal(ilm_aggregate(p)$rows$n_rows, c(3L, 3L))   # one total per scenario
  expect_output(print(a), "a total of new observations.*exact")
  expect_error(ilm_aggregate(p, by = "site"), "do not have: site")
})

test_that("an average is the total over the number of rows, scored through it", {
  p <- ilm_predict(agg_fit("poisson"), agg_rows, target = "response", groups = "fitted", nsim = 200)
  tot <- lmn_aggregate(p, agg_groups, agg_out)
  avg <- lmn_aggregate(p, agg_groups, agg_out, fun = "mean")
  y <- c(7, 9)
  expect_equal(lmn_pred_cdf(avg, y / 3), lmn_pred_cdf(tot, y))
  expect_equal(lmn_pred_quantile(avg, 0.8), lmn_pred_quantile(tot, 0.8) / 3)
  expect_equal(lmn_pred_mean(avg), lmn_pred_mean(tot) / 3)
  s_tot <- lmn_pred_scores(tot, y, rules = c("log", "crps"))
  s_avg <- lmn_pred_scores(avg, y / 3, rules = c("log", "crps"))
  expect_equal(s_avg$log, s_tot$log)                        # a count's average keeps its steps
  expect_equal(s_avg$crps, s_tot$crps / 3)
  u <- ilm_pit(avg, y / 3, seed = 2)
  expect_true(all(u >= lmn_pred_cdf(tot, y - 1) - 1e-12 & u <= lmn_pred_cdf(tot, y) + 1e-12))
  expect_error(lmn_pred_scores(avg, c(7.1, 9) / 3, rules = "crps"), "must be whole")
  ## a gaussian average is itself gaussian: the mean over m, the SD over m
  pg <- ilm_predict(agg_fit("gaussian"), agg_rows, target = "response", groups = "fitted", nsim = 200)
  ag <- lmn_aggregate(pg, agg_groups, agg_out, fun = "mean")
  sd <- sqrt(3) * matrix(exp(pg$par$logdisp), 2, 200, byrow = TRUE) / 3
  expect_equal(lmn_pred_cdf(ag, c(1, 1)),
               unname(rowMeans(stats::pnorm(1, rowsum(pg$lp$mu, agg_groups) / 3, sd))),
               tolerance = 1e-12)
  ## simulated averages are the simulated totals over m
  pn <- ilm_predict(agg_fit("nbinom"), agg_rows, target = "response", groups = "fitted", nsim = 200)
  sn <- lmn_aggregate(pn, agg_groups, agg_out)
  an <- lmn_aggregate(pn, agg_groups, agg_out, fun = "mean")
  expect_equal(lmn_pred_quantile(an, 0.5), lmn_pred_quantile(sn, 0.5) / 3)
  expect_equal(as.data.frame(an)$estimate, as.data.frame(sn)$estimate / 3)
})

test_that("missing simulated totals stop in the registry's words", {
  expect_error(lmn_signal("missing_totals", bad = 3, S = 8, type = "error"),
               "3 of 8 simulated totals are missing.*family = \"poisson\"",
               class = "lmn_missing_totals")
})
