## forecast_coverage: do forecast intervals cover at their nominal level, by
## horizon?
##
## Simulation (extending prototype P3's AR(1) design). 20 series, each with
## a random intercept (SD 0.5) and a covariate, observed at times 1 to T and
## forecast at T + 1 to T + 4 with the covariate known in advance. Cells:
##   - AR(1), true and fitted, rho 0.7 and marginal SD 0.8: gaussian (noise SD
##     0.5) and Poisson, T = 8, 16 and 32;
##   - AR(2) (phi 0.5 and 0.3, the same marginal SD), fitted as AR(1): gaussian
##     and Poisson, T = 16. A misspecified cell, reported whatever it shows;
##   - a random walk, true and fitted (step SD 0.3): Poisson, T = 8, 16 and 32.
## Coverage of the central 50, 80 and 90% intervals by the randomised PIT of
## the realised value, under uncertainty = "none", "fixed" and "joint". The
## settings share the PIT's random numbers. Monte Carlo SEs come from
## replicate means, since a replicate's series share a fit.
##
## Run from the package root: Rscript studies/scripts/forecast_coverage.R [R] [cores]
## `... forecast_coverage.R summarise` rewrites the findings from the kept replicates.
args <- commandArgs(TRUE)
summarise_only <- identical(args[1], "summarise")
R <- if (!summarise_only && length(args) >= 1) as.integer(args[1]) else 500L
cores <- if (length(args) >= 2) as.integer(args[2]) else 6L
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "forecast_coverage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- rbind(
  expand.grid(family = c("gaussian", "poisson"), truth = "ar1", fitted = "ar1", Tn = c(8L, 16L, 32L),
              stringsAsFactors = FALSE),
  expand.grid(family = c("gaussian", "poisson"), truth = "ar2", fitted = "ar1", Tn = 16L,
              stringsAsFactors = FALSE),
  expand.grid(family = "poisson", truth = "rw1", fitted = "rw1", Tn = c(8L, 16L, 32L),
              stringsAsFactors = FALSE))
cells$cell <- seq_len(nrow(cells))

rep_cell <- function(rep, cell, cells, G = 20L, H = 4L, S = 1000L) {
  ce <- cells[cell, ]; Tn <- ce$Tn; fam <- ce$family
  set.seed(2e6 + 1e4 * cell + rep)
  n_t <- Tn + H
  latent <- function() switch(ce$truth,
    ar1 = as.numeric(stats::arima.sim(list(ar = 0.7), n_t, sd = 0.8 * sqrt(1 - 0.7^2))),
    ar2 = { p1 <- 0.5; p2 <- 0.3; r1 <- p1 / (1 - p2); r2 <- p1 * r1 + p2
            as.numeric(stats::arima.sim(list(ar = c(p1, p2)), n_t,
                                        sd = 0.8 * sqrt(1 - p1 * r1 - p2 * r2))) },
    rw1 = c(0, cumsum(stats::rnorm(n_t - 1, 0, 0.3))))
  full <- expand.grid(t = seq_len(n_t), g = factor(sprintf("s%02d", seq_len(G))))
  full$x <- stats::rnorm(nrow(full))
  full$a <- unlist(lapply(seq_len(G), function(k) latent()))
  b <- stats::rnorm(G, 0, 0.5)
  eta <- if (fam == "gaussian") 1 + 0.5 * full$x + b[full$g] + full$a
         else 0.3 + 0.3 * full$x + b[full$g] + full$a
  full$y <- if (fam == "gaussian") stats::rnorm(nrow(full), eta, 0.5) else stats::rpois(nrow(full), exp(eta))
  hist <- full[full$t <= Tn, ]; fut <- full[full$t > Tn, ]
  ar <- if (ce$fitted == "ar1") illume::ilm_ar1(~ t | g) else illume::ilm_rw1(~ t | g)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressWarnings(suppressMessages(
    illume::ilm_model(y ~ x + (1 | g), data = hist, family = fam, ar = ar, verbose = FALSE))),
    error = function(e) NULL)
  base <- data.frame(rep = rep, cell = cell)
  if (is.null(f)) return(cbind(base, failed = "fit error"))
  conv <- f$opt$convergence
  boundary <- length(f$boundary_terms) > 0       # a variance or correlation at its edge
  out <- list()
  for (u in c("none", "fixed", "joint")) {
    res <- tryCatch({
      fc <- ilm_forecast(f, fut[c("g", "t", "x")], uncertainty = u, nsim = S, seed = rep)
      pit <- ilm_pit(fc, fut$y, seed = rep)
      h <- fut$t - Tn
      do.call(rbind, lapply(c(0.5, 0.8, 0.9), function(L)
        data.frame(uncertainty = u, level = L, h = seq_len(H),
                   cover = vapply(seq_len(H), function(k) mean(abs(pit[h == k] - 0.5) <= L / 2), 0))))
    }, error = function(e) data.frame(uncertainty = u, level = NA, h = NA, cover = NA,
                                      error = conditionMessage(e)))
    out[[u]] <- res
  }
  res <- do.call(rbind, lapply(out, function(z) { if (is.null(z$error)) z$error <- NA; z }))
  cbind(base, failed = NA, convergence = conv, boundary = boundary,
        seconds = proc.time()[["elapsed"]] - t0, res, row.names = NULL)
}

rds <- file.path(out_dir, "replicates.rds")
if (!summarise_only) {
  jobs <- expand.grid(rep = seq_len(R), cell = cells$cell)
  t_start <- Sys.time()
  cl <- parallel::makePSOCKcluster(cores)
  invisible(parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE),
                                  normalizePath(".")))
  reps <- parallel::clusterMap(cl, rep_cell, jobs$rep, jobs$cell, MoreArgs = list(cells = cells),
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

## per cell, setting, level and horizon: coverage and its MC SE over replicates
ok <- reps[is.na(reps$failed) & is.na(reps$error) & !is.na(reps$cover), ]
key <- c("cell", "uncertainty", "level", "h")
summ <- do.call(rbind, lapply(split(ok, ok[key], drop = TRUE), function(s) cbind(s[1, key],
  coverage = mean(s$cover), se = stats::sd(s$cover) / sqrt(nrow(s)), R = nrow(s))))
summ <- merge(cells, summ, by = "cell")
summ <- summ[do.call(order, summ[c("cell", "uncertainty", "level", "h")]), ]; rownames(summ) <- NULL
write.csv(summ, file.path(out_dir, "coverage.csv"), row.names = FALSE)
fit_info <- do.call(rbind, lapply(split(reps, reps$cell), function(s) {
  one <- s[!duplicated(s$rep), ]
  data.frame(cell = s$cell[1], fit_failed = sum(!is.na(one$failed)),
             not_converged = sum(one$convergence %in% 1:10, na.rm = TRUE),
             boundary = mean(one$boundary, na.rm = TRUE),
             pred_failed = sum(!is.na(s$error) & !duplicated(s[c("rep", "uncertainty")])))
}))
fit_info <- merge(cells, fit_info, by = "cell")

## the settings compared on the same replicates: those where every setting
## ran and the fit converged; and, apart, the replicates where joint could not
## run (a variance at its boundary, which illume does not yet hold)
all_ok <- aggregate(ok ~ rep + cell, FUN = all,
                    data = transform(reps[is.na(reps$failed), ],
                                     ok = is.na(error) & convergence %in% 0))
common <- merge(ok, all_ok[all_ok$ok, c("rep", "cell")], by = c("rep", "cell"))
jfail <- aggregate(f ~ rep + cell, FUN = any,
                   data = transform(reps[is.na(reps$failed), ], f = !is.na(error)))
bound <- merge(ok, jfail[jfail$f, c("rep", "cell")], by = c("rep", "cell"))
cov_by <- function(z) do.call(rbind, lapply(split(z, z[c("cell", "uncertainty", "h")], drop = TRUE),
  function(s) data.frame(cell = s$cell[1], uncertainty = s$uncertainty[1], h = s$h[1],
                         coverage = mean(s$cover), se = stats::sd(s$cover) / sqrt(nrow(s)),
                         R = nrow(s))))
common_s <- cov_by(common[common$level == 0.9 & common$h %in% c(1, 4), ])
bound_s <- cov_by(bound[bound$level == 0.9 & bound$h %in% c(1, 4), ])
write.csv(common_s, file.path(out_dir, "coverage_common.csv"), row.names = FALSE)
write.csv(bound_s, file.path(out_dir, "coverage_boundary.csv"), row.names = FALSE)

f3 <- function(v) sprintf("%.3f", v)
lab <- function(ce) sprintf("%s %s%s, T = %d", ce$family, toupper(ce$truth),
                            if (ce$truth != ce$fitted) paste0(" fitted as ", toupper(ce$fitted)) else "", ce$Tn)
cellf <- function(cl, u, L, h) {
  s <- summ[summ$cell == cl & summ$uncertainty == u & summ$level == L & summ$h == h, ]
  if (!nrow(s)) "" else paste0(f3(s$coverage), " (", f3(s$se), ")")
}
tab_h <- function(u, L) c(
  "| cell | h = 1 | h = 2 | h = 3 | h = 4 |", "|---|---|---|---|---|",
  vapply(cells$cell, function(cl) paste0("| ", lab(cells[cl, ]), " | ",
    paste(vapply(1:4, function(h) cellf(cl, u, L, h), ""), collapse = " | "), " |"), ""))
tab_u <- c("| cell | h | none | fixed | joint |", "|---|---|---|---|---|",
  unlist(lapply(cells$cell, function(cl) vapply(c(1, 4), function(h)
    paste0("| ", lab(cells[cl, ]), " | ", h, " | ",
           paste(vapply(c("none", "fixed", "joint"), function(u) cellf(cl, u, 0.9, h), ""),
                 collapse = " | "), " |"), ""))))
tab_s <- function(z, us) c(
  paste0("| cell | h | replicates | ", paste(us, collapse = " | "), " |"),
  paste0("|---|---|---|", strrep("---|", length(us))),
  unlist(lapply(cells$cell, function(cl) vapply(c(1, 4), function(h) {
    zz <- z[z$cell == cl & z$h == h, ]
    if (!nrow(zz)) return(NA_character_)
    paste0("| ", lab(cells[cl, ]), " | ", h, " | ", max(zz$R), " | ",
           paste(vapply(us, function(u) { q <- zz[zz$uncertainty == u, ]
             if (!nrow(q)) "" else paste0(f3(q$coverage), " (", f3(q$se), ")") }, ""),
             collapse = " | "), " |") }, "")))) 
low <- summ[which(summ$uncertainty == "joint" & summ$coverage < summ$level - 3 * summ$se), ]
md <- c(
  "# Coverage of forecast intervals, by horizon", "",
  sprintf("## limn %s (%s; illume %s)", version, format(Sys.Date()), utils::packageVersion("illume")), "",
  "Generated by `studies/scripts/forecast_coverage.R`. Do not edit by hand.", "",
  sprintf(paste("20 series with a random intercept and a covariate, observed at times 1 to T and forecast 1 to 4",
                "steps ahead with the covariate known. %d replicates per cell, 1000 draws per forecast. Run on",
                "limn %s: %.0f minutes on %d cores."), R, run$commit, run$minutes, run$cores), "",
  "Fits per cell: failed, not converged (nlminb), and the share with a variance at its boundary; and forecasts that failed:", "",
  sprintf("- %s: %d failed, %d not converged, %.1f%% at a boundary; %d forecasts failed",
          vapply(fit_info$cell, function(cl) lab(cells[cl, ]), ""), fit_info$fit_failed,
          fit_info$not_converged, 100 * fit_info$boundary, fit_info$pred_failed), "",
  "Coverage, with its Monte Carlo SE in brackets.", "",
  "### Central 90% intervals, joint uncertainty", "", tab_h("joint", 0.9), "",
  "### Central 90% intervals at the first and last horizon, by uncertainty", "",
  "Each setting over the replicates where it ran, so joint leaves out the boundary fits that the others include:", "",
  tab_u, "",
  "### The same, on the replicates where every setting ran and the fit converged", "",
  "The fair comparison of the settings:", "",
  stats::na.omit(tab_s(common_s, c("none", "fixed", "joint"))), "",
  "### Where joint could not run: a variance at its boundary", "",
  "illume does not yet hold a variance estimated at its boundary, so joint draws are not usable there; the other settings, on those replicates:", "",
  stats::na.omit(tab_s(bound_s, c("none", "fixed"))), "",
  "### Central 50% intervals, joint uncertainty", "", tab_h("joint", 0.5), "",
  "### Central 80% intervals, joint uncertainty", "", tab_h("joint", 0.8), "",
  "### Cells where joint uncertainty covers more than 3 SEs below nominal", "",
  if (nrow(low)) sprintf("- %s, h = %d, %s: %s (%s)", vapply(low$cell, function(cl) lab(cells[cl, ]), ""),
                         low$h, low$level, f3(low$coverage), f3(low$se)) else "None.")
writeLines(md, file.path("studies", "findings", "forecast_coverage.md"))
cat(md, sep = "\n")
