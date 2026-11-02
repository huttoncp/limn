# Validation studies

The evidence behind `limn`'s validity claims, kept the way `illume` keeps its
own. These are **not** part of the R package: `.Rbuildignore` leaves this
directory out of `R CMD build`. It is versioned so that results stay with the
code that produced them.

```
studies/
  scripts/     the study scripts
  runs/        <version>/<study>/  raw per-replicate .rds and summary .csv
  findings/    <study>.md  the permanent record, newest version first
```

**Every claim is checked against something `limn` did not write.** That is an
independent implementation where one exists, such as `rsample`, `forecast` or
`scoringRules`. Where none exists, it is a simulation with a known answer. The
packages used here are not dependencies of `limn`, so a study runs only where
they are installed.

**What is committed.**

- Committed: scripts, per-cell summary `.csv` files, and findings.
- Not committed: raw per-replicate `.rds`, which the scripts' fixed seeds
  reproduce exactly. It stays on local disk, as in `illume`.

**Commit IDs in these files.** Some findings, designs and scripts here cite
`limn` commits by their IDs: the commit a study was run at, as in the version
line at the top of a findings file, and commits for fixes, changes, earlier
runs and code kept for reference. Every `limn` commit ID in these files
refers to `limn`'s development history before its first public commit. That
history is not published, so those IDs cannot be looked up in this
repository. They are kept as they were recorded.

Results that go against `limn` are reported as plainly as those that do not.

**Where the outside packages come from.** Most comparators live in the
ordinary R library. `scoringRules` 1.1.3 and `hypergeo` 1.2.14 are the
exception: `hypergeo` (with `elliptic`, `contfrac` and `deSolve`) is what
`scoringRules`' negative binomial CRPS needs. Both sit in a scratch library
beside the repository, installed with the author's agreement on 2026-09-24.
- The tests and `R CMD check` are run with that library first on the path:
  `.libPaths(c(scratch, .libPaths()))`, then `rcmdcheck(libpath = .libPaths())`.
- Without it, `R CMD check` reports `scoringRules` as a missing suggested
  package, and the comparisons in `tests/testthat/test-scores.R` are skipped
  rather than run.
- Run with it, on 2026-09-24:
  - limn's negative binomial CRPS matched the exact energy form, summed over
    the probability mass function, to 4.5e-14 in 45 cases;
  - `scoringRules`' closed form agreed wherever it was finite. It returned -Inf
    in 10 of the 45 cases, size 2.5 with a mean of 2 or more, so it is
    compared only where finite.
  - NEWS records both.

## Scripts

| script | what it answers | outside implementation |
|---|---|---|
| `folds_rsample.R` | Are `ilm_folds()`'s rolling origins and grouped folds the standard ones? | `rsample::rolling_origin()`, `rsample::group_vfold_cv()` |
| `compare_dm.R` | Does `ilm_compare(by = "origin")` give the Diebold-Mariano test? | `forecast::dm.test()` |
| `compare_dm_pooled.R` | Does the pooled Diebold-Mariano test of `ilm_compare(by = "origin")` hold its size at 5%, against the single-horizon tests, and what is its power? | none: a simulation with known truth |
| `accuracy_fabletools.R` | Do MASE, RMSSE, MAE, RMSE, the Winkler score and the quantile score agree, for normal forecasts and mixtures? | `fabletools::accuracy()` |
| `engine_glmmTMB.R` | Does `ilm_predict()` give an independent implementation's predictions at the same parameters, with the same spread, and exact population means? | `glmmTMB` with its parameters held at `illume`'s; Newton's method in base R; closed forms; `illume::ilm_scenario()` |
| `forecast_exact.R` | Are forecasts, and values between observed times, the exact conditional distribution for AR(1), CAR(1) and a random walk, and how do joint draws compare with the delta method? | exact multivariate-normal conditioning; `glmmTMB` `ar1()` and `ou()` with the new rows at weight zero |
| `forecast_coverage.R` | Do forecast intervals cover at their nominal level by horizon, for AR(1), a misspecified AR(2) and a random walk, gaussian and Poisson, with 8 to 32 times per series? | none: a simulation with known truth, 11 cells of 500 replicates |
| `aggregate_coverage.R` | Do totals from `ilm_aggregate()` cover at their nominal level, over a series' horizon, across series and in all, and what do joint draws buy over summing rows separately? | none: a simulation with known truth, 3 families of 500 replicates |
| `aggregate_mint.R` | How do limn's coherent totals compare with fable's reconciled forecasts on real hierarchical data? A benchmark, reported whichever way it goes. | `fable` ETS with MinT and bottom-up reconciliation on `tsibble::tourism` |
| `engine_coverage.R` | Do the engine's intervals cover at their nominal level, for new observations and for means, for seen and new units? | none: a simulation with known truth, 5 families, 8 and 40 groups |
