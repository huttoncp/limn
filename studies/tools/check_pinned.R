## check_pinned: before a long run, exercise every kind of model a study fits,
## on the pinned library alone (LIMN_STUDY_LIB), and in a worker too. Fails
## loudly if anything cannot load or fit. A null model from fable counts as a
## failure: fable reports a model's error as a warning and returns NULL.
##
## LIMN_STUDY_LIB=<lib> Rscript studies/tools/check_pinned.R
source(file.path("studies", "scripts", "pinned_library.R"))
suppressPackageStartupMessages({ pkgload::load_all(".", quiet = TRUE); library(tsibble); library(fable) })
stopifnot(startsWith(normalizePath(find.package("illume")), normalizePath(lmn_study_lib)),
          startsWith(normalizePath(find.package("fable")), normalizePath(lmn_study_lib)))
set.seed(1)
d <- data.frame(g = factor(rep(sprintf("s%d", 1:6), each = 16)), t = rep(1:16, 6))
d$y <- rpois(nrow(d), exp(1 + 0.3 * sin(d$t))); d$z <- d$y + rnorm(nrow(d))
for (fam in c("gaussian", "poisson", "nbinom")) {
  yy <- if (fam == "gaussian") "z" else "y"
  f <- suppressMessages(suppressWarnings(illume::ilm_model(stats::as.formula(paste(yy, "~ 1 + (1 | g)")),
         data = d, family = fam, ar = illume::ilm_ar1(~ t | g), verbose = FALSE)))
  fc <- ilm_forecast(f, data.frame(g = "s1", t = 17:18), uncertainty = "fixed", nsim = 100)
  stopifnot(all(is.finite(lmn_pred_scores(fc, c(3, 3), rules = c("crps", "log"))$crps)))
}
ts <- tsibble::as_tsibble(d, key = g, index = t)
m <- suppressWarnings(fabletools::model(ts, ets = fable::ETS(z), arima = fable::ARIMA(z)))
null <- vapply(c("ets", "arima"), function(k) any(fabletools::is_null_model(m[[k]])), TRUE)
if (any(null)) stop("fable returned null models for: ", paste(names(null)[null], collapse = ", "))
fcf <- fabletools::forecast(m, h = 2)
stopifnot(all(is.finite(mean(fcf$z))))
cl <- parallel::makePSOCKcluster(1); lmn_study_workers(cl)
w <- parallel::clusterCall(cl, function() { suppressPackageStartupMessages(library(fable)); find.package("fable") })[[1]]
parallel::stopCluster(cl)
stopifnot(startsWith(normalizePath(w), normalizePath(lmn_study_lib)))
cat("pinned library OK:", lmn_study_lib, "\n")
