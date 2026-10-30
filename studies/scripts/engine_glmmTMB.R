## engine_glmmTMB: does ilm_predict() give an independent implementation's
## predictions, at the same parameters, and the exact population mean?
##
## Outside implementation: glmmTMB (not a limn dependency).
##   1. Shared parameters. glmmTMB is fitted with every parameter held (`map`)
##      at illume's estimates, so only the random effects are solved for, and
##      its predict() gives the fitted-unit mean (re.form = NULL) and the
##      typical-unit mean (re.form = NA) at exactly illume's parameters. Random
##      intercept models, gaussian, Poisson and logit binomial. glmmTMB and
##      illume share TMB's inner optimiser, so the unit's effect is also found
##      without TMB: Newton's method on its exact log posterior, in base R.
##   2. Parameter uncertainty. The standard deviation of the engine's joint
##      draws of the linear predictor, against glmmTMB's se.fit on its own fit:
##      for a typical unit (the fixed effects' delta method) and for a known
##      unit (TMB's generalised delta method over the fixed effects, the
##      variance and the unit's effect). Both are the same Gaussian
##      approximation, so they should agree within the draws' Monte Carlo
##      error (about 1% at 5000 draws).
##   3. Population means against exact forms: exp(eta + s^2/2) for the log
##      link, and integrate() over the effect's normal for the logit. illume's
##      binomial family is logit only, so the probit's closed form does not apply.
##   4. illume's own ilm_scenario(over = "sample"), which averages population
##      means over the observed rows: the same estimate, and interval limits
##      within Monte Carlo error, for random-intercept and random-slope fits.
##
## Run from the package root: Rscript studies/scripts/engine_glmmTMB.R
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(glmmTMB)
})
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "engine_glmmTMB")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260925)
G <- 20L; n <- 8L
d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(G)), each = n)), x = rnorm(G * n))
b <- rnorm(G, 0, 0.8)
eta <- 0.3 + 0.6 * d$x + b[d$g]
d$y_gaussian <- rnorm(nrow(d), eta, 0.9)
d$y_poisson <- rpois(nrow(d), exp(eta))
d$y_binomial <- rbinom(nrow(d), 1, plogis(eta))
nd <- data.frame(g = rep(c("g01", "g05", "g12"), each = 3), x = rep(c(-1, 0, 1.5), 3))
fams <- list(gaussian = gaussian(), poisson = poisson(), binomial = binomial())

## a unit's conditional mode given the parameters: Newton's method on
## sum log p(y | offset + b) - b^2 / (2 sd_g^2), with its exact derivatives
newton_mode <- function(y, off, fam, sd_g, sig) {
  b <- 0
  for (it in 1:100) {
    eta <- off + b
    if (fam == "gaussian") { s <- sum(y - eta) / sig^2; h <- length(y) / sig^2 }
    if (fam == "poisson") { m <- exp(eta); s <- sum(y - m); h <- sum(m) }
    if (fam == "binomial") { m <- stats::plogis(eta); s <- sum(y - m); h <- sum(m * (1 - m)) }
    step <- (s - b / sd_g^2) / (h + 1 / sd_g^2)
    b <- b + step
    if (abs(step) < 1e-14) break
  }
  b
}

rows <- list(); unc <- list(); popr <- list()
for (fam in names(fams)) {
  fml <- stats::as.formula(paste0("y_", fam, " ~ x + (1 | g)"))
  f <- suppressWarnings(suppressMessages(
    illume::ilm_model(fml, data = d, family = fam, verbose = FALSE)))
  vc <- illume::ilm_varcorr(f)
  beta <- unname(coef(f)); sd_g <- sqrt(vc$re$g[1, 1])
  start <- list(beta = beta, theta = log(sd_g))
  map <- list(beta = factor(c(NA, NA)), theta = factor(NA))
  if (fam == "gaussian") {
    start$betadisp <- log(vc$dispersion$value[1]); map$betadisp <- factor(NA)
  }
  m_shared <- suppressWarnings(glmmTMB(fml, data = d, family = fams[[fam]],
                                       start = start, map = map))
  e_fit <- ilm_predict(f, nd, groups = "fitted", uncertainty = "none")$estimate
  e_typ <- ilm_predict(f, nd, groups = "typical", uncertainty = "none")$estimate
  t_fit <- as.numeric(predict(m_shared, nd, re.form = NULL, type = "response"))
  t_typ <- as.numeric(predict(m_shared, nd, re.form = NA, type = "response"))
  sig <- if (fam == "gaussian") vc$dispersion$value[1] else NA
  b_hat <- vapply(unique(nd$g), function(k) {
    on <- d$g == k
    newton_mode(d[[paste0("y_", fam)]][on], beta[1] + beta[2] * d$x[on], fam, sd_g, sig)
  }, 0)
  n_fit <- fams[[fam]]$linkinv(beta[1] + beta[2] * nd$x + b_hat[nd$g])
  rows[[fam]] <- data.frame(family = fam,
    fitted_max_rel_diff = max(abs(e_fit - t_fit) / pmax(1, abs(t_fit))),
    fitted_newton_max_rel_diff = max(abs(e_fit - n_fit) / pmax(1, abs(n_fit))),
    typical_max_rel_diff = max(abs(e_typ - t_typ) / pmax(1, abs(t_typ))))

  ## 2. the spread of the linear predictor from the parameters' uncertainty
  m_own <- suppressWarnings(glmmTMB(fml, data = d, family = fams[[fam]]))
  for (unit in c("typical", "fitted")) {
    se_t <- as.numeric(predict(m_own, nd, re.form = if (unit == "fitted") NULL else NA,
                               se.fit = TRUE)$se.fit)
    lk <- ilm_predict(f, nd, target = "link", groups = unit, nsim = 5000L, seed = 3L)
    sd_e <- apply(lk$draws, 1L, stats::sd)
    unc[[paste(fam, unit)]] <- data.frame(family = fam, unit = unit,
                                          sd_ratio_min = min(sd_e / se_t),
                                          sd_ratio_max = max(sd_e / se_t))
  }

  ## 3. the population mean, against exact forms
  if (fam != "gaussian") {
    pop <- ilm_predict(f, nd, groups = "population", uncertainty = "none")$estimate
    e0 <- beta[1] + beta[2] * nd$x
    exact <- if (fam == "poisson") exp(e0 + sd_g^2 / 2) else
      vapply(e0, function(e) stats::integrate(function(z) stats::plogis(e + sd_g * z) *
                                                stats::dnorm(z), -Inf, Inf,
                                              rel.tol = 1e-12)$value, 0)
    popr[[fam]] <- data.frame(family = fam, exact = if (fam == "poisson")
      "exp(eta + s^2/2)" else "integrate()", max_rel_diff = max(abs(pop - exact) / exact))
  }
}
shared <- do.call(rbind, rows); unc <- do.call(rbind, unc); popr <- do.call(rbind, popr)

## 4. against ilm_scenario(over = "sample"): the population mean with x set for
## every observed row, averaged. Its interval limits come from draws, so they
## are compared over R independent pairs of seeds: the mean difference, in SDs
## of the draws, against its standard error across the pairs.
S <- 4000L; R <- 12L; xs <- c(-1, 0, 1)
d$y_slope <- rpois(nrow(d), exp(0.3 + 0.6 * d$x + b[d$g] + rnorm(G, 0, 0.4)[d$g] * d$x))
scen_pair <- function(fml, r, d, xs, S) {
  f <- suppressWarnings(suppressMessages(
    illume::ilm_model(stats::as.formula(fml), data = d, family = "poisson", verbose = FALSE)))
  s_i <- illume::ilm_scenario(f, x = xs, over = "sample", sims = S, level = 0.95,
                              seed = 100L + r)
  p <- ilm_predict(f, d, groups = "population", set = list(x = xs), nsim = S, seed = 200L + r)
  do.call(rbind, lapply(seq_along(xs), function(k) {
    on <- p$rows$x == xs[k]
    dr_k <- colMeans(p$draws[on, , drop = FALSE])
    q_e <- stats::quantile(dr_k, c(0.025, 0.975), names = FALSE)
    data.frame(model = sub("^y_\\w+ ~ ", "", fml), x = xs[k], pair = r,
               estimate_rel_diff = abs(mean(p$estimate[on]) - s_i$estimate[k]) / s_i$estimate[k],
               lower_diff_sd = (q_e[1] - s_i$lower[k]) / stats::sd(dr_k),
               upper_diff_sd = (q_e[2] - s_i$upper[k]) / stats::sd(dr_k))
  }))
}
jobs <- expand.grid(fml = c("y_poisson ~ x + (1 | g)", "y_slope ~ x + (1 + x | g)"),
                    r = seq_len(R), stringsAsFactors = FALSE)
cl <- parallel::makePSOCKcluster(6L)
invisible(parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE),
                                normalizePath(".")))
scen <- do.call(rbind, parallel::clusterMap(cl, scen_pair, jobs$fml, jobs$r,
                                            MoreArgs = list(d = d, xs = xs, S = S)))
parallel::stopCluster(cl)
scen_sum <- do.call(rbind, lapply(split(scen, list(scen$model, scen$x)), function(s)
  data.frame(model = s$model[1], x = s$x[1], estimate_rel_diff = max(s$estimate_rel_diff),
             lower_mean = mean(s$lower_diff_sd), lower_se = stats::sd(s$lower_diff_sd) / sqrt(R),
             upper_mean = mean(s$upper_diff_sd), upper_se = stats::sd(s$upper_diff_sd) / sqrt(R))))
print(shared, row.names = FALSE); print(unc, row.names = FALSE); print(popr, row.names = FALSE)
print(scen_sum, row.names = FALSE)
write.csv(shared, file.path(out_dir, "shared_parameters.csv"), row.names = FALSE)
write.csv(unc, file.path(out_dir, "uncertainty.csv"), row.names = FALSE)
write.csv(popr, file.path(out_dir, "population.csv"), row.names = FALSE)
write.csv(scen, file.path(out_dir, "scenario_pairs.csv"), row.names = FALSE)
write.csv(scen_sum, file.path(out_dir, "scenario.csv"), row.names = FALSE)

f1 <- function(v) formatC(v, format = "e", digits = 1)
md <- c(
  "# ilm_predict() against glmmTMB and exact integrals", "",
  sprintf("## limn %s (%s; glmmTMB %s; illume %s)", version, format(Sys.Date()),
          utils::packageVersion("glmmTMB"), utils::packageVersion("illume")), "",
  "Generated by `studies/scripts/engine_glmmTMB.R`. Do not edit by hand.", "",
  "20 groups of 8, random intercepts, predictions for 9 rows in 3 groups.", "",
  "**At shared parameters.** glmmTMB with every parameter held at illume's estimates, solving only for the random effects, and the unit's effect by Newton's method in base R. Largest relative difference:", "",
  paste0("- ", shared$family, ": fitted-unit mean ", f1(shared$fitted_max_rel_diff),
         " (glmmTMB), ", f1(shared$fitted_newton_max_rel_diff),
         " (Newton); typical-unit mean ", f1(shared$typical_max_rel_diff)), "",
  "**Parameter uncertainty.** The standard deviation of 5000 joint draws of the linear predictor, over glmmTMB's standard error on its own fit: the delta method for a typical unit, and TMB's generalised delta method (fixed effects, variance and the unit's effect) for a known unit. Both are the same Gaussian approximation; the Monte Carlo error is about 1%:", "",
  paste0("- ", unc$family, ", ", unc$unit, " unit: ratio ", sprintf("%.3f", unc$sd_ratio_min),
         " to ", sprintf("%.3f", unc$sd_ratio_max)), "",
  "**Population means** at the estimates, against exact forms (illume's binomial is logit only, so there is no probit case):", "",
  paste0("- ", popr$family, " (", popr$exact, "): ", f1(popr$max_rel_diff)), "",
  sprintf("**Against `ilm_scenario(over = \"sample\")`**, Poisson, x set to -1, 0 and 1 for all 160 rows. Largest relative difference in the estimate over %d pairs of seeds, and the 95%% interval limits' mean difference (engine minus `ilm_scenario()`, in SDs of the draws, %d draws each) with its standard error across the pairs:", R, S), "",
  sprintf("- `%s`, x = %g: estimate %s; lower limit %+.3f (%.3f); upper %+.3f (%.3f)",
          scen_sum$model, scen_sum$x, f1(scen_sum$estimate_rel_diff), scen_sum$lower_mean,
          scen_sum$lower_se, scen_sum$upper_mean, scen_sum$upper_se), "",
  sprintf(paste("The largest of these differences is %.1f of its standard errors. The three",
                "values of x share their draws, so each model's rows are not independent checks."),
          max(abs(c(scen_sum$lower_mean / scen_sum$lower_se,
                    scen_sum$upper_mean / scen_sum$upper_se)))))
writeLines(md, file.path("studies", "findings", "engine_glmmTMB.md"))
cat(md, sep = "\n")
