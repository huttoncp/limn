## engine_coverage: do ilm_predict()'s intervals cover at their nominal level?
##
## Simulation (prototype P3, now through the engine, with three more
## families and the means' intervals added). Random-intercept models,
## y ~ x + (1 | g), with x ~ N(0, 1), 6 rows per group and an effect SD of 0.8,
## for 8 groups (where the sources of uncertainty should matter) and 40 (the
## control). Families: gaussian (residual SD 1), poisson, nbinom (size 2),
## binomial (0/1) and beta (precision 10).
##
## Each replicate fits once and predicts a new row in every seen group, and one
## row in each of 5 new groups (levels the fit has not seen, so groups =
## "fitted" predicts them as new units):
##   - a new observation, judged by its randomised PIT: the central L interval
##     covers when |PIT - 0.5| <= L / 2, for L = 0.5, 0.8 and 0.9;
##   - the mean, by whether its 90% and 95% intervals hold the true value: a
##     seen unit's mean and a new unit's (with their true effects), a typical
##     unit's (effect zero) and the population average (integrated exactly).
## Under uncertainty = "none", "fixed" and "joint". The settings share the
## PIT's random numbers, so their differences carry no randomisation noise. A
## plug-in ("none") has no interval for a mean, so means use "fixed" and
## "joint" only.
##
## Monte Carlo SEs come from replicate means, since a replicate's rows share a
## fit. Per-replicate results are kept as rds (not committed).
##
## Run from the package root: Rscript studies/scripts/engine_coverage.R [R] [cores]
## `Rscript studies/scripts/engine_coverage.R summarise` rewrites the findings
## from the kept replicates without rerunning them.
args <- commandArgs(TRUE)
summarise_only <- identical(args[1], "summarise")
R <- if (!summarise_only && length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 6L
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "engine_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
FAMS <- c("gaussian", "poisson", "nbinom", "binomial", "beta")

rep_cell <- function(rep, fam, G, fams, n = 6L, sb = 0.8, S = 1000L, n_new = 5L) {
  set.seed(1e6 + 1e5 * match(fam, fams) + 1e4 * (G == 40L) + rep)
  inv <- if (fam %in% c("binomial", "beta")) stats::plogis else if (fam == "gaussian") identity else exp
  gen <- function(mu) switch(fam,
    gaussian = stats::rnorm(length(mu), mu, 1),
    poisson = stats::rpois(length(mu), mu),
    nbinom = stats::rnbinom(length(mu), mu = mu, size = 2),
    binomial = stats::rbinom(length(mu), 1, mu),
    beta = pmin(pmax(stats::rbeta(length(mu), mu * 10, (1 - mu) * 10), 1e-10), 1 - 1e-10))
  pop <- function(eta) switch(fam, gaussian = eta, poisson = , nbinom = exp(eta + sb^2 / 2),
    vapply(eta, function(e) stats::integrate(function(z) stats::plogis(e + sb * z) *
                                               stats::dnorm(z), -Inf, Inf, rel.tol = 1e-10)$value, 0))
  b0 <- switch(fam, gaussian = 1, poisson = , nbinom = log(3), 0); b1 <- 0.5
  d <- data.frame(g = factor(sprintf("g%02d", rep(seq_len(G), each = n))), x = stats::rnorm(G * n))
  b <- stats::rnorm(G, 0, sb)
  d$y <- gen(inv(b0 + b1 * d$x + b[d$g]))
  ts <- data.frame(g = c(levels(d$g), sprintf("new%d", seq_len(n_new))),
                   x = stats::rnorm(G + n_new))
  eta0 <- b0 + b1 * ts$x
  mu_ts <- inv(eta0 + c(b, stats::rnorm(n_new, 0, sb)))
  ts$y <- gen(mu_ts)
  unit <- rep(c("seen", "new"), c(G, n_new))
  truth <- list(fitted = mu_ts, typical = inv(eta0), population = pop(eta0))

  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressWarnings(suppressMessages(
    illume::ilm_model(y ~ x + (1 | g), data = d, family = fam, verbose = FALSE))),
    error = function(e) NULL)
  if (is.null(f)) return(data.frame(rep = rep, family = fam, G = G, failed = TRUE))
  sd_hat <- sqrt(illume::ilm_varcorr(f)$re$g[1, 1])
  out <- list()
  add <- function(...) out[[length(out) + 1L]] <<- data.frame(rep = rep, family = fam, G = G, ...)
  for (u in c("none", "fixed", "joint")) {
    n0 <- length(out)                             # a failed setting keeps none of its rows
    ok <- tryCatch({
      r <- ilm_predict(f, ts, target = "response", groups = "fitted", uncertainty = u,
                       nsim = S, seed = rep)
      pit <- ilm_pit(r, ts$y, seed = rep)
      for (L in c(0.5, 0.8, 0.9)) for (k in c("seen", "new"))
        add(uncertainty = u, what = "response", target = k, level = L,
            cover = mean(abs(pit[unit == k] - 0.5) <= L / 2))
      if (u != "none") for (gr in c("fitted", "typical", "population")) {
        m <- as.data.frame(ilm_predict(f, ts, target = "mean", groups = gr, uncertainty = u,
                                       nsim = S, seed = rep), level = c(0.9, 0.95))
        for (L in c(0.9, 0.95)) {
          inside <- m[[paste0("lower_", L)]] <= truth[[gr]] & truth[[gr]] <= m[[paste0("upper_", L)]]
          if (gr == "fitted") {
            for (k in c("seen", "new"))
              add(uncertainty = u, what = "mean", target = k, level = L, cover = mean(inside[unit == k]))
          } else add(uncertainty = u, what = "mean", target = gr, level = L, cover = mean(inside))
        }
      }
      TRUE
    }, error = function(e) conditionMessage(e))
    if (!isTRUE(ok)) {
      out <- out[seq_len(n0)]
      add(uncertainty = u, what = "error", target = ok, level = NA, cover = NA)
    }
  }
  res <- do.call(rbind, out)
  res$failed <- FALSE; res$sd_hat <- sd_hat
  res$seconds <- proc.time()[["elapsed"]] - t0
  res
}

rds <- file.path(out_dir, "replicates.rds")
if (!summarise_only) {
jobs <- expand.grid(rep = seq_len(R), family = FAMS, G = c(8L, 40L), stringsAsFactors = FALSE)
t_start <- Sys.time()
cl <- parallel::makePSOCKcluster(cores)
invisible(parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE),
                                normalizePath(".")))
reps <- parallel::clusterMap(cl, rep_cell, jobs$rep, jobs$family, jobs$G,
                             MoreArgs = list(fams = FAMS), .scheduling = "dynamic")
parallel::stopCluster(cl)
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
all_cols <- unique(unlist(lapply(reps, names)))
reps <- do.call(rbind, lapply(reps, function(x) { x[setdiff(all_cols, names(x))] <- NA; x[all_cols] }))
attr(reps, "run") <- list(R = R, cores = cores, minutes = minutes,
                          commit = system("git rev-parse --short HEAD", intern = TRUE))
saveRDS(reps, rds)
}
reps <- readRDS(rds)
run <- attr(reps, "run"); R <- run$R; cores <- run$cores; minutes <- run$minutes

## per cell: coverage, the MC SE over replicate means, and how many replicates
fails <- aggregate(failed ~ family + G, data = unique(reps[c("rep", "family", "G", "failed")]), FUN = sum)
errs <- reps[!reps$failed & reps$what == "error", ]
ok <- reps[!reps$failed & reps$what != "error", ]
boundary <- aggregate(sd_hat ~ family + G, data = unique(ok[c("rep", "family", "G", "sd_hat")]),
                      FUN = function(v) mean(v < 0.01))
key <- c("family", "G", "uncertainty", "what", "target", "level")
cells <- split(ok, ok[key], drop = TRUE)
summ <- do.call(rbind, lapply(cells, function(s) cbind(s[1, key],
  coverage = mean(s$cover), se = stats::sd(s$cover) / sqrt(nrow(s)), R = nrow(s))))
summ <- summ[do.call(order, summ[key]), ]; rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "coverage.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

## findings
f3 <- function(v) sprintf("%.3f", v)
cellf <- function(fam, G, u, what, tg, L) {
  s <- summ[summ$family == fam & summ$G == G & summ$uncertainty == u & summ$what == what &
              summ$target == tg & summ$level == L, ]
  if (!nrow(s)) return("")
  paste0(f3(s$coverage), " (", f3(s$se), ")")
}
tab <- function(what, targets, L, us) {
  hdr <- c(paste0("| family | groups | target | ", paste(us, collapse = " | "), " |"),
           paste0("|---|---|---|", strrep("---|", length(us))))
  body <- character(0)
  for (fam in FAMS) for (G in c(8L, 40L)) for (tg in targets)
    body <- c(body, paste0("| ", fam, " | ", G, " | ", tg, " | ",
                           paste(vapply(us, function(u) cellf(fam, G, u, what, tg, L), ""),
                                 collapse = " | "), " |"))
  c(hdr, body)
}
low <- summ[summ$uncertainty == "joint" & summ$coverage < summ$level - 3 * summ$se, ]
md <- c(
  "# Coverage of ilm_predict()'s intervals", "",
  sprintf("## limn %s (%s; illume %s)", version, format(Sys.Date()), utils::packageVersion("illume")), "",
  "Generated by `studies/scripts/engine_coverage.R`. Do not edit by hand.", "",
  sprintf(paste("Random-intercept fits, `y ~ x + (1 | g)`, 6 rows per group, effect SD 0.8; %d replicates per",
                "family and group count; 1000 draws per prediction. Each replicate predicts a new row in",
                "every seen group and one row in each of 5 new groups. Run on limn %s: %.0f minutes on %d cores."),
          R, run$commit, minutes, cores), "",
  paste0("Fits that failed: ", paste0(fails$family, " ", fails$G, " groups: ", fails$failed, collapse = "; "),
         ". Replicates whose predictions failed under a setting (those replicates are left out of",
         " that setting's cells): ",
         if (!nrow(errs)) "none." else paste0(nrow(errs), ", as follows."), ""),
  if (nrow(errs)) {
    errs$kind <- ifelse(grepl("are not finite", errs$target),
                        "the variance's draws were not finite (an estimate at the boundary, not held)",
                        ifelse(grepl("missing values|could not be computed", errs$target),
                               "the population mean overflowed in some draws (a large drawn variance)",
                               substr(errs$target, 1, 80)))
    e <- aggregate(rep ~ family + G + uncertainty + kind, data = errs, FUN = length)
    e <- e[order(e$family, e$G), ]
    paste0("- ", e$family, ", ", e$G, " groups, uncertainty = \"", e$uncertainty, "\": ", e$rep,
           " of ", R, ": ", e$kind, ".")
  }, "",
  sprintf(paste("For reference: a z-interval covers %.3f at 90%% and %.3f at 95%% when the truth follows t",
                "with 7 df, and %.3f and %.3f with 39 df."),
          2 * stats::pt(stats::qnorm(0.95), 7) - 1, 2 * stats::pt(stats::qnorm(0.975), 7) - 1,
          2 * stats::pt(stats::qnorm(0.95), 39) - 1, 2 * stats::pt(stats::qnorm(0.975), 39) - 1), "",
  paste0("Fits with the effect's SD at the boundary (below 0.01): ",
         paste0(boundary$family, " ", boundary$G, ": ", sprintf("%.1f%%", 100 * boundary$sd_hat),
                collapse = "; "), "."), "",
  "Coverage, with its Monte Carlo SE in brackets.", "",
  "### A new observation, central 90% interval", "",
  tab("response", c("seen", "new"), 0.9, c("none", "fixed", "joint")), "",
  "### A new observation, central 50% interval", "",
  tab("response", c("seen", "new"), 0.5, c("none", "fixed", "joint")), "",
  "### The mean, 95% interval", "",
  tab("mean", c("seen", "new", "typical", "population"), 0.95, c("fixed", "joint")), "",
  "### The mean, 90% interval", "",
  tab("mean", c("seen", "new", "typical", "population"), 0.9, c("fixed", "joint")), "",
  "### Cells where joint uncertainty covers more than 3 SEs below nominal", "",
  if (nrow(low)) paste0("- ", low$family, ", ", low$G, " groups, ", low$what, " (", low$target,
                        "), ", low$level, ": ", f3(low$coverage), " (", f3(low$se), ")")
  else "None.")
writeLines(md, file.path("studies", "findings", "engine_coverage.md"))
cat(md, sep = "\n")
