## aggregate_mint: limn's coherent totals against fable's reconciled forecasts,
## on Australian domestic tourism. Studies only; reported whichever way it
## goes.
##
## Data: tsibble::tourism, overnight trips (thousands) by quarter, 1998-2017,
## summed over purpose to 76 regions in 8 states and rounded to whole thousands,
## so that both sides forecast the same counts. Train 1998-2015 (72 quarters),
## test 2016-2017 (8 quarters), as fpp3 section 11.3 does.
##
## fable: ETS with additive errors on every series of the State / Region
## hierarchy (76 regions, 8 states, the total), whose forecasts are normal and
## so reconcile whole, reconciled bottom-up and by MinT with shrinkage (fpp3
## section 11.3).
## limn: negative binomial models for the 76 regions, each with an intercept
## and a random-walk level per region, forecast with joint draws; the state and
## national totals come from ilm_aggregate(), so they are coherent by
## construction. Two specifications (registered in advance in
## studies/designs/aggregate_mint_richer.md):
##   simple: a quarterly pattern common to all regions;
##   richer: that, plus a drift common to all (a trend in t) and a random
##           region-by-quarter effect, so each region's seasonal pattern may
##           depart from the common one.
##
## Scores, for each quarter's value at each level (national, state, region):
## the CRPS, and the coverage of central 80% and 95% intervals.
##
## Run from the package root, on a pinned library (studies/scripts/pinned_library.R):
##   LIMN_STUDY_LIB=<lib> Rscript studies/scripts/aggregate_mint.R
source(file.path("studies", "scripts", "pinned_library.R"))
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(tsibble); library(fabletools); library(fable); library(dplyr)
})
version <- as.character(utils::packageVersion("limn"))
out_dir <- file.path("studies", "runs", version, "aggregate_mint")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reg <- tsibble::tourism |> group_by(Region, State) |> summarise(Trips = round(sum(Trips)), .groups = "drop")
hier <- reg |> aggregate_key(State / Region, Trips = sum(Trips))
cut_q <- yearquarter("2015 Q4")

## ---- fable: ETS on every series, reconciled ----
t0 <- Sys.time()
fb <- hier |> filter(Quarter <= cut_q) |> model(base = ETS(Trips ~ error("A"))) |>
  reconcile(bu = bottom_up(base), mint = min_trace(base, method = "mint_shrink"))
fc_fb <- forecast(fb, h = 8)
fable_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lev <- function(state, region) ifelse(is_aggregated(state), "national",
                                      ifelse(is_aggregated(region), "state", "region"))
fbd <- data.frame(model = fc_fb$.model, level = lev(fc_fb$State, fc_fb$Region),
                  state = as.character(fc_fb$State), region = as.character(fc_fb$Region),
                  quarter = as.character(fc_fb$Quarter),
                  mu = mean(fc_fb$Trips), sd = sqrt(distributional::variance(fc_fb$Trips)))
test <- hier |> filter(Quarter > cut_q)
act <- data.frame(state = as.character(test$State), region = as.character(test$Region),
                  quarter = as.character(test$Quarter), y = test$Trips)
fbd <- merge(fbd, act, by = c("state", "region", "quarter"))
crps_norm <- function(y, mu, sd) { z <- (y - mu) / sd
  sd * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi)) }
fbd$crps <- crps_norm(fbd$y, fbd$mu, fbd$sd)
fbd$in80 <- abs(fbd$y - fbd$mu) <= stats::qnorm(0.9) * fbd$sd
fbd$in95 <- abs(fbd$y - fbd$mu) <= stats::qnorm(0.975) * fbd$sd

## ---- limn: two specifications for the regions, totals by ilm_aggregate() ----
rd <- as.data.frame(reg)
rd$quarter_label <- as.character(rd$Quarter)
qs <- sort(unique(rd$Quarter))
rd$t <- match(rd$Quarter, qs); rd$season <- factor(as.integer(format(as.Date(rd$Quarter), "%m")) %/% 3 + 1)
rd$Region <- factor(rd$Region)
rd$region_season <- interaction(rd$Region, rd$season, drop = TRUE)
train <- rd[rd$t <= 72, ]; testr <- rd[rd$t > 72, ]
specs <- list(simple = Trips ~ season + (1 | Region),
              richer = Trips ~ season + t + (1 | Region) + (1 | region_season))
rows_new <- testr[c("Region", "State", "t", "season", "region_season")]
limn_runs <- lapply(names(specs), function(k) {
  t0 <- Sys.time()
  fit <- suppressWarnings(illume::ilm_model(specs[[k]], data = train, family = "nbinom",
                                            ar = illume::ilm_rw1(~ t | Region), verbose = FALSE))
  fc <- tryCatch(suppressWarnings(ilm_forecast(fit, rows_new, nsim = 1000, seed = 1)),
                 error = function(e) e)
  used <- "joint"
  if (inherits(fc, "error")) {                     # an unheld boundary variance: hold the variances
    used <- paste0("fixed (joint draws failed: ", conditionMessage(fc), ")")
    fc <- suppressWarnings(ilm_forecast(fit, rows_new, uncertainty = "fixed", nsim = 1000, seed = 1))
  }
  ## illume's checks on the fit (?ilm_model, "Always read the checks"): the fit
  ## runs under suppressWarnings(), so they are reported here instead
  ch <- fit$checks
  bad <- ch[ch$status %in% c("WARN", "FAIL"), , drop = FALSE]
  checks <- sprintf("%d FAIL, %d WARN%s", sum(ch$status == "FAIL"), sum(ch$status == "WARN"),
                    if (nrow(bad)) paste0(" (", paste0(bad$check, ": ", bad$status, collapse = "; "), ")") else "")
  list(spec = k, fit = fit, fc = fc, used = used, checks = checks, loglik = as.numeric(logLik(fit)),
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
})
names(limn_runs) <- names(specs)
score_limn <- function(p, y, level, model) {
  s <- lmn_pred_scores(p, y, rules = "crps")
  q <- function(pr) lmn_pred_quantile(p, pr)
  data.frame(model = model, level = level, crps = s$crps,
             in80 = y >= q(0.1) & y <= q(0.9), in95 = y >= q(0.025) & y <= q(0.975))
}
nat_rows <- list(); lm_rows <- list()
for (k in names(limn_runs)) {
  fc <- limn_runs[[k]]$fc; m <- paste0("limn_", k)
  by_t <- ilm_aggregate(fc, by = "t")
  y_nat <- as.numeric(tapply(testr$Trips, testr$t, sum)[as.character(by_t$rows$t)])
  by_st <- ilm_aggregate(fc, by = c("State", "t"))
  y_st <- vapply(seq_len(nrow(by_st$rows)), function(i)
    sum(testr$Trips[testr$State == by_st$rows$State[i] & testr$t == by_st$rows$t[i]]), 0)
  lm_rows[[k]] <- rbind(score_limn(by_t, y_nat, "national", m), score_limn(by_st, y_st, "state", m),
                        score_limn(fc, testr$Trips, "region", m))
  dn <- as.data.frame(by_t, level = 0.8)
  nat_rows[[k]] <- data.frame(mean = dn$mean, sd = apply(by_t$draws, 1L, stats::sd))
}
last_year <- mean(tail(as.numeric(tapply(train$Trips, train$t, sum)), 4))

## ---- summary ----
all_rows <- rbind(do.call(rbind, lm_rows), fbd[c("model", "level", "crps", "in80", "in95")])
summ <- do.call(rbind, lapply(split(all_rows, all_rows[c("model", "level")], drop = TRUE), function(s)
  data.frame(model = s$model[1], level = s$level[1], n = nrow(s), crps = mean(s$crps),
             cover80 = mean(s$in80), cover95 = mean(s$in95))))
summ$model <- factor(summ$model, c("limn_simple", "limn_richer", "mint", "bu", "base"))
summ$level <- factor(summ$level, c("national", "state", "region"))
summ <- summ[order(summ$level, summ$model), ]; rownames(summ) <- NULL
print(summ)
write.csv(summ, file.path(out_dir, "scores.csv"), row.names = FALSE)

## the national totals quarter by quarter: centre and spread against the outcome
nat_fb <- fbd[fbd$level == "national" & fbd$model == "mint", ]; nat_fb <- nat_fb[order(nat_fb$quarter), ]
nat <- data.frame(quarter = sort(unique(testr$quarter_label)),
                  actual = as.numeric(tapply(testr$Trips, testr$t, sum)),
                  simple_mean = nat_rows$simple$mean, simple_sd = nat_rows$simple$sd,
                  richer_mean = nat_rows$richer$mean, richer_sd = nat_rows$richer$sd,
                  mint_mean = nat_fb$mu, mint_sd = nat_fb$sd)
write.csv(nat, file.path(out_dir, "national.csv"), row.names = FALSE)
vc_r <- illume::ilm_varcorr(limn_runs$richer$fit)
drift <- coef(limn_runs$richer$fit)[["t"]]

names_m <- c(limn_simple = "limn, simple: one NB model, a common season, random-walk levels",
             limn_richer = "limn, richer: plus a common drift and region-specific seasons",
             mint = "fable: ETS, MinT (shrink) reconciled", bu = "fable: ETS, bottom-up",
             base = "fable: ETS on each series, unreconciled")
md <- c(
  "# Totals for Australian domestic tourism: limn against fable's reconciliation", "",
  sprintf("## limn %s (%s; illume %s; fable %s; fabletools %s; commit %s)", version, format(Sys.Date()),
          utils::packageVersion("illume"), utils::packageVersion("fable"),
          utils::packageVersion("fabletools"), system("git rev-parse --short HEAD", intern = TRUE)), "",
  "Generated by `studies/scripts/aggregate_mint.R`. Do not edit by hand. The richer limn specification was registered in advance in `studies/designs/aggregate_mint_richer.md`.", "",
  lmn_study_header(), "",
  paste("Scored with the discrete CRPS by running sums of the mass function (limn's faster scoring,",
        "the same result to rounding: 2.6e-15 relative on these 608 regional forecasts); the",
        "pinned library, data, models and seeds are those of the first run."), "",
  "illume's checks on each limn fit (the fits run under suppressWarnings(), so the checks are reported here):", "",
  sprintf("- %s: %s; log-likelihood %.3f.", names(limn_runs),
          vapply(limn_runs, `[[`, "", "checks"), vapply(limn_runs, `[[`, 0, "loglik")), "",
  paste("Trips summed over purpose to 76 regions in 8 states and rounded to whole thousands;",
        "trained on 1998-2015 and scored on the 8 quarters of 2016-2017, as fpp3 section 11.3 does.",
        sprintf("limn's forecasts used %s uncertainty (simple) and %s (richer); they took %.1f and %.1f minutes, fable's %.1f.",
                limn_runs$simple$used, limn_runs$richer$used, limn_runs$simple$minutes,
                limn_runs$richer$minutes, fable_minutes)), "",
  sprintf(paste("The richer fit's drift is %.4f on the log scale a quarter (%.1f%% a year), and its",
                "region-by-quarter SD %.3f, beside a region SD of %.3f."),
          drift, 100 * (exp(4 * drift) - 1),
          sqrt(c(vc_r$re$region_season)), sqrt(c(vc_r$re$Region))), "",
  "Mean CRPS (thousands of trips, lower is better) and the coverage of central 80% and 95% intervals, by level; n is the number of quarter-by-series values scored:", "",
  "| level | model | n | CRPS | 80% covered | 95% covered |", "|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %.1f | %.2f | %.2f |", summ$level, names_m[as.character(summ$model)],
          summ$n, summ$crps, summ$cover80, summ$cover95), "",
  "One split of one dataset: a benchmark, not a test. The national level has 8 values, so its coverage is coarse.", "",
  sprintf("The national totals quarter by quarter (thousands of trips; the last training year averaged %.0f):", last_year), "",
  "| quarter | actual | simple mean | simple SD | richer mean | richer SD | MinT mean | MinT SD |",
  "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f |", nat$quarter, nat$actual, nat$simple_mean,
          nat$simple_sd, nat$richer_mean, nat$richer_sd, nat$mint_mean, nat$mint_sd))
writeLines(md, file.path("studies", "findings", "aggregate_mint.md"))
cat(md, sep = "\n")
