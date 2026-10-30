## aggregate_coverage: do intervals for totals cover at their nominal level,
## and what do joint draws buy over summing each row's forecast on its own?
##
## Simulation. The panels are
## forecast_coverage's: 20 series with a random intercept (SD 0.5) and a
## covariate, AR(1) with rho 0.7 and marginal SD 0.8, observed at times 1 to 16
## and forecast 1 to 4 steps ahead. Families: gaussian (noise SD 0.5) and
## Poisson, whose totals are exact mixtures, and negative binomial (size 5),
## whose totals are simulated from each draw.
##
## Three kinds of total from each forecast: each series over the 4 steps (20
## totals), all series at each step (4), and the grand total (1). Two ways to
## form them:
##   - joint: ilm_aggregate() on the forecast, summing within each draw;
##   - row by row: each row's draws shuffled independently first, which keeps
##     every row's own forecast but removes the correlation between rows (a
##     series' path, the shared parameters and effects), as summing forecasts
##     made separately does.
## Coverage of the central 50, 80 and 90% intervals by the randomised PIT of
## the realised total; interval widths as the row-by-row width over the joint.
## The forecasts hold the variances at their estimates (uncertainty = "fixed"):
## joint draws of the variances fail on fits whose variance sits at its
## boundary until illume holds it, and the study is of how totals are formed,
## which the held variances leave intact (the paths and the shared parameters
## and effects are still drawn). To be rerun with joint draws after the fix.
## Monte Carlo SEs from replicate means.
##
## Run from the package root: Rscript studies/scripts/aggregate_coverage.R [R] [cores]
## `... aggregate_coverage.R summarise` rewrites the findings from the kept replicates.
args <- commandArgs(TRUE)
summarise_only <- identical(args[1], "summarise")
R <- if (!summarise_only && length(args) >= 1) as.integer(args[1]) else 500L
cores <- if (length(args) >= 2) as.integer(args[2]) else 6L
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "aggregate_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
FAMS <- c("gaussian", "poisson", "nbinom")

## each row's draws in an independent random order: the rows' own forecasts,
## without their correlation
row_by_row <- function(pred, seed) {
  set.seed(seed)
  n <- lmn_pred_n(pred); S <- length(pred$weights)
  perm <- t(vapply(seq_len(n), function(i) sample.int(S), integer(S)))    # n x S
  idx <- cbind(rep(seq_len(n), S), as.vector(perm))
  out <- pred
  out$lp <- lapply(pred$lp, function(M) matrix(M[idx], n, S))
  out$par <- lapply(pred$par, function(v) if (is.matrix(v)) matrix(v[idx], n, S) else matrix(v[perm], n, S))
  out
}

rep_cell <- function(rep, fam, fams, row_by_row, G = 20L, Tn = 16L, H = 4L, S = 1000L) {
  set.seed(3e6 + 1e5 * match(fam, fams) + rep)
  full <- expand.grid(t = seq_len(Tn + H), g = factor(sprintf("s%02d", seq_len(G))))
  full$x <- stats::rnorm(nrow(full))
  full$a <- unlist(lapply(seq_len(G), function(k)
    as.numeric(stats::arima.sim(list(ar = 0.7), Tn + H, sd = 0.8 * sqrt(1 - 0.49)))))
  b <- stats::rnorm(G, 0, 0.5)
  eta <- if (fam == "gaussian") 1 + 0.5 * full$x + b[full$g] + full$a
         else 1 + 0.3 * full$x + b[full$g] + full$a
  full$y <- switch(fam, gaussian = stats::rnorm(nrow(full), eta, 0.5),
                   poisson = stats::rpois(nrow(full), exp(eta)),
                   nbinom = stats::rnbinom(nrow(full), mu = exp(eta), size = 5))
  hist <- full[full$t <= Tn, ]; fut <- full[full$t > Tn, ]
  f <- tryCatch(suppressWarnings(suppressMessages(illume::ilm_model(
    y ~ x + (1 | g), data = hist, family = fam, ar = illume::ilm_ar1(~ t | g), verbose = FALSE))),
    error = function(e) NULL)
  base <- data.frame(rep = rep, family = fam)
  if (is.null(f)) return(cbind(base, failed = "fit error"))
  fc <- tryCatch(ilm_forecast(f, fut[c("g", "t", "x")], uncertainty = "fixed", nsim = S, seed = rep),
                 error = function(e) conditionMessage(e))
  if (is.character(fc)) return(cbind(base, failed = fc))
  res <- tryCatch(score_totals(fc, fut, rep, row_by_row), error = function(e) conditionMessage(e))
  if (is.character(res)) return(cbind(base, failed = res))
  cbind(base, failed = NA, res)
}

score_totals <- function(fc, fut, rep, row_by_row) {
  preds <- list(joint = fc, row_by_row = row_by_row(fc, seed = rep))
  kinds <- list(series = "g", step = "t", total = NULL)
  out <- list()
  for (m in names(preds)) for (k in names(kinds)) {
    a <- ilm_aggregate(preds[[m]], by = kinds[[k]])
    realised <- if (is.null(kinds[[k]])) sum(fut$y)
      else as.numeric(tapply(fut$y, fut[[kinds[[k]]]], sum)[as.character(a$rows[[kinds[[k]]]])])
    u <- ilm_pit(a, realised, seed = rep)
    w90 <- lmn_pred_quantile(a, 0.95) - lmn_pred_quantile(a, 0.05)
    for (L in c(0.5, 0.8, 0.9))
      out[[length(out) + 1L]] <- data.frame(method = m, total = k, level = L,
                                            cover = mean(abs(u - 0.5) <= L / 2), width90 = mean(w90),
                                            aggregate = a$meta$aggregate)
  }
  do.call(rbind, out)
}

rds <- file.path(out_dir, "replicates.rds")
if (!summarise_only) {
  jobs <- expand.grid(rep = seq_len(R), family = FAMS, stringsAsFactors = FALSE)
  t_start <- Sys.time()
  cl <- parallel::makePSOCKcluster(cores)
  invisible(parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE), normalizePath(".")))
  parallel::clusterExport(cl, "score_totals")
  reps <- parallel::clusterMap(cl, rep_cell, jobs$rep, jobs$family,
                               MoreArgs = list(fams = FAMS, row_by_row = row_by_row),
                               .scheduling = "dynamic")
  parallel::stopCluster(cl)
  minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  all_cols <- unique(unlist(lapply(reps, names)))
  reps <- do.call(rbind, lapply(reps, function(x) { x[setdiff(all_cols, names(x))] <- NA; x[all_cols] }))
  attr(reps, "run") <- list(R = R, cores = cores, minutes = minutes,
                            commit = system("git rev-parse --short HEAD", intern = TRUE))
  saveRDS(reps, rds)
}
reps <- readRDS(rds)
run <- attr(reps, "run"); R <- run$R

ok <- reps[is.na(reps$failed), ]
key <- c("family", "method", "total", "level")
summ <- do.call(rbind, lapply(split(ok, ok[key], drop = TRUE), function(s) cbind(s[1, key],
  coverage = mean(s$cover), se = stats::sd(s$cover) / sqrt(nrow(s)), R = nrow(s),
  width90 = mean(s$width90[is.finite(s$width90)]), width_nonfinite = sum(!is.finite(s$width90)))))
summ <- summ[do.call(order, unname(summ[key])), ]; rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "coverage.csv"), row.names = FALSE)
fr <- unique(reps[!is.na(reps$failed), c("rep", "family")])
fails <- if (nrow(fr)) aggregate(rep ~ family, FUN = length, data = fr) else
  data.frame(family = character(0), rep = integer(0))
how <- unique(ok[c("family", "aggregate")])

f3 <- function(v) sprintf("%.3f", v)
cellf <- function(fam, m, k, L) {
  s <- summ[summ$family == fam & summ$method == m & summ$total == k & summ$level == L, ]
  if (!nrow(s)) "" else paste0(f3(s$coverage), " (", f3(s$se), ")")
}
## within each replicate, the row-by-row width over the joint (both from the same
## forecast); the median over replicates, so a replicate with an unheld
## dispersion's enormous widths cannot dominate
ratio <- function(fam, k) {
  z <- ok[ok$family == fam & ok$total == k & ok$level == 0.9, ]
  j <- z[z$method == "joint", ]; r <- z[z$method == "row_by_row", ]
  rr <- r$width90[match(j$rep, r$rep)] / j$width90
  rr <- rr[is.finite(rr)]
  if (!length(rr)) "" else sprintf("%.2f", stats::median(rr))
}
labs <- c(series = "each series over 4 steps", step = "all series at a step", total = "the grand total")
tab <- function(L) c(
  "| family | total | joint | row by row | width, row by row / joint |", "|---|---|---|---|---|",
  unlist(lapply(FAMS, function(fam) vapply(names(labs), function(k)
    sprintf("| %s | %s | %s | %s | %s |", fam, labs[[k]], cellf(fam, "joint", k, L),
            cellf(fam, "row_by_row", k, L), ratio(fam, k)), ""))))
md <- c(
  "# Coverage of intervals for totals", "",
  sprintf("## limn %s (%s; illume %s)", version, format(Sys.Date()), utils::packageVersion("illume")), "",
  "Generated by `studies/scripts/aggregate_coverage.R`. Do not edit by hand.", "",
  sprintf(paste("20 series with a random intercept and a covariate, AR(1) with rho 0.7, observed at 16",
                "times and forecast 1 to 4 steps ahead, the variances held at their estimates",
                "(uncertainty = \"fixed\"; joint draws wait on illume's boundary fix). %d replicates per",
                "family, 1000 draws per forecast. Run on limn %s: %.0f minutes on %d cores."),
          R, run$commit, run$minutes, run$cores), "",
  "How the totals were formed, by family:", "",
  sprintf("- %s: %s", how$family, how$aggregate), "",
  paste0("Replicates whose fit, forecast or totals failed (left out): ",
         if (nrow(fails)) paste0(fails$family, " ", fails$rep, collapse = "; ") else "none", "."),
  if (nrow(fr)) c("", "Why, by kind of failure:", "", {
    fm <- reps[!is.na(reps$failed), c("rep", "family", "failed")]; fm <- unique(fm)
    fm$kind <- substr(sub("[0-9]+ of [0-9]+ ", "", fm$failed), 1, 90)
    k <- aggregate(rep ~ family + kind, data = fm, FUN = length)
    sprintf("- %s: %d, \"%s...\"", k$family, k$rep, k$kind) }), "",
  "Coverage, with its Monte Carlo SE in brackets; the last column is the median, over replicates, of the 90% intervals' width row by row over joint.", "",
  "### Central 90% intervals", "", tab(0.9), "",
  sprintf("Widths not finite, left out of the width ratios: %d replicate-totals.",
          sum(summ$width_nonfinite[summ$level == 0.9])), "",
  "### Central 80% intervals", "", tab(0.8), "",
  "### Central 50% intervals", "", tab(0.5))
writeLines(md, file.path("studies", "findings", "aggregate_coverage.md"))
cat(md, sep = "\n")
