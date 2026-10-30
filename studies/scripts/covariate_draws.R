## covariate_draws: when a covariate is not known in advance, do forecasts that
## mix over its draws cover, where forecasts conditional on its point forecast
## do not?
##
## Simulation. 20 series,
## each with a random intercept (SD 0.5), a latent AR(1) (rho 0.7, marginal SD
## 0.5) and a covariate x that is itself an AR(1) (phi 0.8, marginal SD 1),
## observed at times 1 to 24 and forecast at 25 to 30. Gaussian (noise SD 0.5,
## x's coefficient 0.8) and Poisson (x's coefficient 0.4). The model is the
## true one, y ~ x + (1 | g) with an AR(1) over time.
##
## Four ways of giving the future x, all with joint parameter draws:
##   known:  its realised values, an ex-post forecast; the reference, which
##           should cover at its nominal level;
##   point:  its point forecast from its own (known) process, phi^h x_24, as if
##           it were known: the usual ex-ante practice (fpp3 section 10.3);
##   draws:  C = 100 paths from its own process given x_24, through a .draw
##           column, which ilm_forecast() mixes over;
##   generated: C = 100 paths from ilm_covariate_paths(), from an illume model
##           of x fitted to its history (x ~ 1 with an AR(1) by series, the
##           right form), so its parameters are estimated, not known.
## "draws" knows the covariate's process, so it isolates what drawing it adds;
## "generated" is what a user would do. Coverage of central 50, 80 and 90% intervals by the randomised PIT,
## by horizon, and the median width of the 80% interval against "known".
##
## Run from the package root: Rscript studies/scripts/covariate_draws.R [R] [cores]
## Runs on a pinned library: set LIMN_STUDY_LIB (see studies/scripts/pinned_library.R).
## `... covariate_draws.R summarise` rewrites the findings from the kept replicates.
args <- commandArgs(TRUE)
summarise_only <- identical(args[1], "summarise")
R <- if (!summarise_only && length(args) >= 1) as.integer(args[1]) else 500L
cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
source(file.path("studies", "scripts", "pinned_library.R"))
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "covariate_draws")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
families <- c("gaussian", "poisson")

rep_cell <- function(rep, fam, G = 20L, Tn = 24L, H = 6L, C = 100L, S = 1000L, phi = 0.8) {
  set.seed(3e6 + 1e4 * match(fam, c("gaussian", "poisson")) + rep)
  n_t <- Tn + H
  ar1 <- function(r, sd_marg) as.numeric(stats::arima.sim(list(ar = r), n_t, sd = sd_marg * sqrt(1 - r^2)))
  full <- expand.grid(t = seq_len(n_t), g = factor(sprintf("s%02d", seq_len(G))))
  full$x <- unlist(lapply(seq_len(G), function(k) ar1(phi, 1)))
  full$a <- unlist(lapply(seq_len(G), function(k) ar1(0.7, 0.5)))
  b <- stats::rnorm(G, 0, 0.5)
  eta <- if (fam == "gaussian") 1 + 0.8 * full$x + b[full$g] + full$a
         else 0.5 + 0.4 * full$x + b[full$g] + full$a
  full$y <- if (fam == "gaussian") stats::rnorm(nrow(full), eta, 0.5) else stats::rpois(nrow(full), exp(eta))
  hist <- full[full$t <= Tn, ]; fut <- full[full$t > Tn, ]
  h <- fut$t - Tn
  x_last <- hist$x[hist$t == Tn][match(fut$g, hist$g[hist$t == Tn])]
  ## the covariate's future: point forecast, and C paths given its last value
  point <- fut; point$x <- phi^h * x_last
  paths <- do.call(rbind, lapply(seq_len(C), function(c) {
    z <- fut; innov <- matrix(stats::rnorm(G * H, 0, sqrt(1 - phi^2)), H)
    xs <- numeric(nrow(fut))
    for (k in seq_len(G)) {
      i <- which(fut$g == levels(fut$g)[k]); i <- i[order(fut$t[i])]
      prev <- x_last[i[1]]
      for (j in seq_along(i)) { prev <- phi * prev + innov[j, k]; xs[i[j]] <- prev }
    }
    z$x <- xs; cbind(.draw = c, z)
  }))
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressWarnings(suppressMessages(
    illume::ilm_model(y ~ x + (1 | g), data = hist, family = fam, ar = illume::ilm_ar1(~ t | g),
                      verbose = FALSE))), error = function(e) NULL)
  base <- data.frame(rep = rep, family = fam)
  if (is.null(f)) return(cbind(base, failed = "fit error"))
  fx <- tryCatch(suppressWarnings(suppressMessages(
    illume::ilm_model(x ~ 1, data = hist, family = "gaussian", ar = illume::ilm_ar1(~ t | g),
                      verbose = FALSE))), error = function(e) NULL)
  gen <- if (is.null(fx)) NULL else tryCatch(
    ilm_covariate_paths(fx, fut[c("g", "t", "x")], nsim = C, seed = rep), error = function(e) NULL)
  arms <- list(known = fut[c("g", "t", "x")], point = point[c("g", "t", "x")],
               draws = paths[c(".draw", "g", "t", "x")], generated = gen)
  out <- list()
  for (a in names(arms)) {
    if (is.null(arms[[a]])) {
      out[[a]] <- data.frame(arm = a, level = NA, h = NA, cover = NA, width80 = NA,
                             error = "covariate model or paths failed")
      next
    }
    out[[a]] <- tryCatch({
      fc <- suppressMessages(ilm_forecast(f, arms[[a]], nsim = S, seed = rep))
      pit <- ilm_pit(fc, fut$y, seed = rep)
      w80 <- lmn_pred_quantile(fc, 0.9) - lmn_pred_quantile(fc, 0.1)
      do.call(rbind, lapply(c(0.5, 0.8, 0.9), function(L)
        data.frame(arm = a, level = L, h = seq_len(H),
                   cover = vapply(seq_len(H), function(k) mean(abs(pit[h == k] - 0.5) <= L / 2), 0),
                   width80 = vapply(seq_len(H), function(k) stats::median(w80[h == k]), 0))))
    }, error = function(e) data.frame(arm = a, level = NA, h = NA, cover = NA, width80 = NA,
                                      error = conditionMessage(e)))
  }
  res <- do.call(rbind, lapply(out, function(z) { if (is.null(z$error)) z$error <- NA; z }))
  cbind(base, failed = NA, convergence = f$opt$convergence,
        boundary = length(f$boundary_terms) > 0,
        seconds = proc.time()[["elapsed"]] - t0, res, row.names = NULL)
}

rds <- file.path(out_dir, "replicates.rds")
if (!summarise_only) {
  jobs <- expand.grid(rep = seq_len(R), family = families, stringsAsFactors = FALSE)
  t_start <- Sys.time()
  cl <- parallel::makePSOCKcluster(cores)
  lmn_study_workers(cl)
  invisible(parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE),
                                  normalizePath(".")))
  reps <- parallel::clusterMap(cl, rep_cell, jobs$rep, jobs$family, .scheduling = "dynamic")
  parallel::stopCluster(cl)
  minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  all_cols <- unique(unlist(lapply(reps, names)))
  reps <- do.call(rbind, lapply(reps, function(x) { x[setdiff(all_cols, names(x))] <- NA; x[all_cols] }))
  attr(reps, "run") <- list(R = R, cores = cores, minutes = minutes,
                            commit = system("git rev-parse --short HEAD", intern = TRUE),
                            illume = as.character(utils::packageVersion("illume")),
                            pinned = lmn_study_header())
  saveRDS(reps, rds)
}
reps <- readRDS(rds)
run <- attr(reps, "run")

ok <- reps[is.na(reps$failed) & is.na(reps$error) & !is.na(reps$cover), ]
key <- c("family", "arm", "level", "h")
summ <- do.call(rbind, lapply(split(ok, ok[key], drop = TRUE), function(s) cbind(s[1, key],
  coverage = mean(s$cover), se = stats::sd(s$cover) / sqrt(nrow(s)), R = nrow(s))))
summ$arm <- factor(summ$arm, c("known", "point", "draws", "generated"))
summ <- summ[do.call(order, unname(summ[key])), ]; rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "coverage.csv"), row.names = FALSE)

## the 80% interval's median width over "known"'s, per replicate, then the median
w <- ok[ok$level == 0.8, c("rep", "family", "arm", "h", "width80")]
wk <- w[w$arm == "known", ]; names(wk)[names(wk) == "width80"] <- "known"
w <- merge(w[w$arm != "known", ], wk[c("rep", "family", "h", "known")], by = c("rep", "family", "h"))
w$ratio <- w$width80 / w$known
wr <- stats::aggregate(ratio ~ family + arm + h, data = w, FUN = stats::median)
failed <- stats::aggregate(cbind(fit = !is.na(failed)) ~ family, data = reps[!duplicated(reps[c("rep", "family")]), ],
                           FUN = sum)
errs <- reps[!is.na(reps$error), ]

tab <- function(L) {
  s <- summ[summ$level == L, ]
  wide <- stats::reshape(s[c("family", "h", "arm", "coverage", "se")], idvar = c("family", "h"),
                         timevar = "arm", direction = "wide")
  c(sprintf("**%d%% intervals.** Coverage (Monte Carlo SE):", round(100 * L)), "",
    "| family | horizon | known | point | draws | generated |", "|---|---|---|---|---|---|",
    sprintf("| %s | %d | %.3f (%.3f) | %.3f (%.3f) | %.3f (%.3f) | %.3f (%.3f) |", wide$family, wide$h,
            wide$coverage.known, wide$se.known, wide$coverage.point, wide$se.point,
            wide$coverage.draws, wide$se.draws, wide$coverage.generated, wide$se.generated), "")
}
wwide <- stats::reshape(wr, idvar = c("family", "h"), timevar = "arm", direction = "wide")
wwide <- wwide[order(wwide$family, wwide$h), ]
md <- c(
  "# Covariates given as draws: coverage against conditioning on a point forecast", "",
  sprintf("## limn %s (%s; illume %s; commit %s)", version, format(Sys.Date()), run$illume, run$commit), "",
  "Generated by `studies/scripts/covariate_draws.R`. Do not edit by hand.", "",
  if (!is.null(run$pinned)) c(run$pinned, "") else "Not run on a pinned library.",
  paste(sprintf("%d replicates per family, %d cores, %.0f minutes.", run$R, run$cores, run$minutes),
        "20 series with a random intercept and a latent AR(1), and a covariate x that is itself an",
        "AR(1) (phi 0.8, marginal SD 1), observed at times 1 to 24 and forecast 1 to 6 steps ahead with",
        "joint parameter draws. The future x is given four ways: its realised values (known, ex post);",
        "its point forecast from its own process (point); 100 paths from that process (draws, a",
        "`.draw` column); or 100 paths from `ilm_covariate_paths()` on an illume model of x fitted to",
        "its history (generated). Coverage is by the randomised PIT, its SE from replicate means."), "",
  tab(0.8), tab(0.9), tab(0.5),
  "**Width.** The median width of the 80% interval over \"known\"'s, per replicate and horizon, then the median over replicates:", "",
  "| family | horizon | point | draws | generated |", "|---|---|---|---|---|",
  sprintf("| %s | %d | %.3f | %.3f | %.3f |", wwide$family, wwide$h, wwide$ratio.point, wwide$ratio.draws,
          wwide$ratio.generated), "",
  sprintf("Fits that failed: %s. Forecasts that failed: %d.",
          paste(sprintf("%s %d", failed$family, failed$fit), collapse = ", "), nrow(errs)),
  if (any(errs$arm == "generated")) c("", paste(
    sprintf("In the generated arm, %d forecasts failed (%s).", sum(errs$arm == "generated"),
            paste(sprintf("%s %d", names(table(errs$family[errs$arm == "generated"])),
                          table(errs$family[errs$arm == "generated"])), collapse = ", ")),
    "The covariate is simulated without measurement noise, so its model's residual SD is estimated",
    "near zero (about 0.02 of the covariate's SD), above illume's hold threshold, and is not held; its",
    "joint draws run to about 1e37, and a few covariate paths then overflow the Poisson forecast. That",
    "gap is illume's to close (its hold study's AR-only panel arm). These failures are expected",
    "to clear once illume holds that SD, and the generated arm is to be rerun then.")))
writeLines(md, file.path("studies", "findings", "covariate_draws.md"))
cat(md, sep = "\n")
