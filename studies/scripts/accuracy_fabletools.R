## accuracy_fabletools: do limn's point and interval accuracy measures agree
## with fabletools::accuracy() on the same forecasts?
##
## Outside implementation: fabletools (not a limn dependency), on normal
## predictive distributions for three monthly series, 36 months of history and
## 6 ahead. Per series:
##   MASE and RMSSE, seasonal (lag 12, fabletools' default for monthly data) and
##     non-seasonal (lag 1);
##   MAE and RMSE; the 80% Winkler score; the quantile score at 0.1, 0.5, 0.9.
## For a normal distribution the median is the mean, so limn's median-based
## absolute errors and fabletools' mean-based ones should agree exactly here.
## (For a skewed distribution they differ by design: limn scores the absolute
## error with the predictive median, the point forecast that minimises it.)
##
## The normal predictive distributions are built with the test suite's
## distribution object (tests/testthat/helper-dist.R), which stands in for
## illume's ilm_dist() until that exists.
##
## Run from the package root: Rscript studies/scripts/accuracy_fabletools.R
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(fabletools); library(tsibble); library(distributional)
})
source(file.path("tests", "testthat", "helper-dist.R"))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "accuracy_fabletools")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260924)
keys <- c("north", "south", "west"); Ttr <- 36L; H <- 6L
full <- expand.grid(i = seq_len(Ttr + H), key = keys, stringsAsFactors = FALSE)
lvl <- c(north = 20, south = 50, west = 5)[full$key]
full$y <- lvl + 3 * sin(2 * pi * full$i / 12) + 0.1 * full$i +
  stats::rnorm(nrow(full), 0, c(north = 1, south = 4, west = 0.5)[full$key])
full$month <- yearmonth("2020 Jan") + (full$i - 1L)
train <- full[full$i <= Ttr, ]; test <- full[full$i > Ttr, ]
m <- test$y + stats::rnorm(nrow(test), 0, 1.5)           # a forecast that is off
s <- 1 + 0.3 * (test$i - Ttr)                            # and widens with horizon

## ---- fabletools ---------------------------------------------------------------
fc <- as_fable(tsibble(data.frame(key = test$key, month = test$month,
                                  y = dist_normal(m, s), .mean = m),
                       key = key, index = month),
               response = "y", distribution = y)
ts_full <- tsibble(full[c("key", "month", "y")], key = key, index = month)
acc <- accuracy(fc, ts_full, measures = list(
  MASE12 = MASE, RMSSE12 = RMSSE,
  MASE1 = function(.resid, .train, ...) MASE(.resid, .train, .period = 1),
  RMSSE1 = function(.resid, .train, ...) RMSSE(.resid, .train, .period = 1),
  MAE = MAE, RMSE = RMSE,
  winkler80 = function(.dist, .actual, ...) winkler_score(.dist, .actual, level = 80),
  qscore = function(.dist, .actual, ...) quantile_score(.dist, .actual, probs = c(0.1, 0.5, 0.9))))
acc <- as.data.frame(acc)[order(acc$key), ]

## ---- limn --------------------------------------------------------------------------
pred <- lmn_new_pred(test_dist_gaussian(), list(mu = matrix(m)),
                     list(logsd = matrix(log(s))), rows = data.frame(key = test$key, i = test$i))
sc12 <- lmn_pred_scores(pred, test$y, c("mase", "rmsse", "absolute", "squared", "interval", "pinball"),
                        probs = c(0.1, 0.5, 0.9), level = 0.8,
                        naive = lmn_naive_scale(train, "i", "key", "y", period = 12))
sc1 <- lmn_pred_scores(pred, test$y, c("mase", "rmsse"),
                       naive = lmn_naive_scale(train, "i", "key", "y", period = 1))
by_key <- function(v, f = mean) as.numeric(tapply(v, test$key, f)[acc$key])
pin <- rowMeans(sc12[c("pinball_0.1", "pinball_0.5", "pinball_0.9")])
lim <- data.frame(
  key = acc$key,
  MASE12 = by_key(sc12$mase), RMSSE12 = sqrt(by_key(sc12$rmsse)),
  MASE1 = by_key(sc1$mase), RMSSE1 = sqrt(by_key(sc1$rmsse)),
  MAE = by_key(sc12$absolute), RMSE = sqrt(by_key(sc12$squared)),
  winkler80 = by_key(sc12$interval_0.8), pinball = by_key(pin))
measures <- c("MASE12", "RMSSE12", "MASE1", "RMSSE1", "MAE", "RMSE", "winkler80")
cmp <- data.frame(measure = measures, max_rel_diff = vapply(measures, function(k)
  max(abs(lim[[k]] - acc[[k]]) / pmax(1, abs(acc[[k]]))), 0))
## fabletools' quantile score is twice the pinball loss (Gneiting & Raftery's
## definition); limn reports the pinball loss itself
ratio <- acc$qscore / lim$pinball
cmp <- rbind(cmp, data.frame(measure = "quantile score / pinball",
                             max_rel_diff = max(abs(ratio - 2) / 2)))
## ---- mixtures: limn's quantiles by bisection against distributional's ----------
## 20 draws per row, so each forecast is a mixture of normals. fabletools reads
## its quantiles and mean from distributional's dist_mixture(); limn from its
## own mixture machinery. The measures read from those agree or not.
S <- 20L
mu <- m + outer(rep(1, length(m)), stats::rnorm(S, 0, 1.2)) +
  matrix(stats::rnorm(length(m) * S, 0, 0.4), length(m), S)
mix <- do.call(dist_mixture, c(lapply(seq_len(S), function(k) dist_normal(mu[, k], s)),
                               list(weights = rep(1 / S, S))))
fcm <- as_fable(tsibble(data.frame(key = test$key, month = test$month, y = mix,
                                   .mean = rowMeans(mu)), key = key, index = month),
                response = "y", distribution = y)
accm <- accuracy(fcm, ts_full, measures = list(
  RMSE = RMSE, RMSSE12 = RMSSE,
  winkler80 = function(.dist, .actual, ...) winkler_score(.dist, .actual, level = 80),
  qscore = function(.dist, .actual, ...) quantile_score(.dist, .actual, probs = c(0.1, 0.5, 0.9))))
accm <- as.data.frame(accm)[order(accm$key), ]
predm <- lmn_new_pred(test_dist_gaussian(), list(mu = mu),
                      list(logsd = matrix(log(s), length(s), S)),
                      rows = data.frame(key = test$key, i = test$i))
scm <- lmn_pred_scores(predm, test$y, c("rmsse", "squared", "interval", "pinball"),
                       probs = c(0.1, 0.5, 0.9), level = 0.8,
                       naive = lmn_naive_scale(train, "i", "key", "y", period = 12))
limm <- data.frame(RMSE = sqrt(by_key(scm$squared)), RMSSE12 = sqrt(by_key(scm$rmsse)),
                   winkler80 = by_key(scm$interval_0.8),
                   pinball = by_key(rowMeans(scm[c("pinball_0.1", "pinball_0.5", "pinball_0.9")])))
mm <- c("RMSE", "RMSSE12", "winkler80")
cmpm <- data.frame(measure = paste0(c(mm, "quantile score / pinball"), " (mixture)"),
  max_rel_diff = c(vapply(mm, function(k) max(abs(limm[[k]] - accm[[k]]) /
                                                pmax(1, abs(accm[[k]]))), 0),
                   max(abs(accm$qscore / limm$pinball - 2) / 2)))
cmp <- rbind(cmp, cmpm)
print(cmp, row.names = FALSE)
## which quantiles are right? The exact mixture CDF at each one's quantiles
Fmix <- function(q) rowMeans(stats::pnorm(q, mu, s))
qres <- vapply(c(0.1, 0.5, 0.9), function(pr) c(
  limn = max(abs(Fmix(lmn_pred_quantile(predm, pr)) - pr)),
  distributional = max(abs(Fmix(unlist(stats::quantile(mix, pr))) - pr))), c(0, 0))
print(qres)
write.csv(cbind(source = "limn", lim), file.path(out_dir, "limn.csv"), row.names = FALSE)
write.csv(cbind(source = "fabletools", acc[c("key", measures, "qscore")]),
          file.path(out_dir, "fabletools.csv"), row.names = FALSE)
write.csv(cmp, file.path(out_dir, "comparison.csv"), row.names = FALSE)

md <- c(
  "# Accuracy measures against fabletools::accuracy()", "",
  sprintf("## limn %s (%s; fabletools %s)", version, format(Sys.Date()),
          utils::packageVersion("fabletools")), "",
  "Generated by `studies/scripts/accuracy_fabletools.R`. Do not edit by hand.", "",
  "Three monthly series, 36 months of history, 6-month normal forecasts that are biased and widen with the horizon; then the same forecasts as 20-draw mixtures of normals, whose quantiles limn finds by bisection and `distributional` by its own root-finding. Largest relative difference over the three series:", "",
  paste0("- ", cmp$measure, ": ", formatC(cmp$max_rel_diff, format = "e", digits = 1)), "",
  sprintf(paste0("The mixtures' quantile-based measures differ by about 1e-6, and that is `distributional`'s root-finding tolerance, not limn's: at the 10%%, 50%% and 90%% quantiles the exact mixture CDF is within %.1e of its target for limn's quantiles and within %.1e for `distributional`'s."),
          max(qres["limn", ]), max(qres["distributional", ])), "",
  "fabletools' quantile score is twice the pinball loss (Gneiting & Raftery's convention); limn reports the pinball loss itself, so the ratio is compared with 2.",
  "For a normal forecast the median is the mean, so the median-based MAE and MASE agree with fabletools' mean-based ones; for a skewed forecast they differ by design.")
writeLines(md, file.path("studies", "findings", "accuracy_fabletools.md"))
cat(md, sep = "\n")
