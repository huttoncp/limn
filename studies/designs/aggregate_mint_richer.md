# Design: a richer limn specification for the tourism benchmark

Written before the run, so that the result cannot shape the specification.
Settled on 25 September 2026: the benchmark of `aggregate_mint` adds one
richer limn specification, chosen in advance, and reports both limn
specifications whichever way they go.

## Why

In the first benchmark (`studies/findings/aggregate_mint.md`) limn's simple
model scored 1.6 to 2 times MinT's CRPS at every level, and its national totals
were biased low: 2016-17 grew, while a random-walk level with no trend and one
seasonal pattern for every region lags a rise. The richer specification adds
exactly those two things, and nothing tuned to the test quarters.

## The two limn specifications

Both are negative binomial `illume::ilm_model()` fits to the 76 regions, with a
random-walk level per region (`ar = ilm_rw1(~ t | Region)`), forecast with
joint draws (1,000, seed 1). Where joint draws fail, the forecast holds the
variances (`uncertainty = "fixed"`) and the findings say so.

- **Simple** (the first benchmark's): `Trips ~ season + (1 | Region)`.
- **Richer**: `Trips ~ season + t + (1 | Region) + (1 | region_season)`.
  - `t` is the quarter index, so the model has a drift common to all regions:
    the random-walk levels move around a trend rather than around a constant.
  - `region_season` is the region-by-quarter factor, so each region's seasonal
    pattern may depart from the common one, shrunk towards it by its variance.

## Everything else as the first benchmark

- Trips summed over purpose to 76 regions in 8 states and rounded to whole
  thousands; trained on 1998-2015 (72 quarters), scored on 2016-2017 (8
  quarters), as fpp3 section 11.3 does.
- fable: ETS on every series of the State / Region hierarchy, reconciled
  bottom-up and by MinT with shrinkage.
- Totals at state and national level from `ilm_aggregate()`.
- Scores: CRPS, and the coverage of central 80% and 95% intervals, by level.
- Both limn specifications are rerun on the same pinned library, so that they
  compare on one build.
- The findings report both limn specifications beside fable's, and the national
  totals quarter by quarter, whichever way they go.

## Changes

- 2026-09-28, after the first run and before any decision rested on it:
  - The findings now report illume's checks on each limn fit: the counts,
    and every WARN or FAIL by name, with the log-likelihood. The fits run
    under `suppressWarnings()`, and the first run's findings did not show
    their checks. Both fits carry WARNs, for the latent budget and the
    gradient.
  - Scoring uses limn's faster discrete CRPS, which gives the same result
    to rounding.
  - Nothing else changes: the pinned library (illume 0.0.8.9000), the
    data, the models and the seeds are all as before. The rerun is
    checked to reproduce every score and log-likelihood of the first run.
  - The rerun reproduced the first run exactly: every CRPS, coverage and
    count, the national table, and both log-likelihoods (the scoring
    change is below the CSV's 15 significant digits).
  - The fits' own warnings, beyond their checks, were not recorded in this
    run. The same fits run alone, on the
    study's data build, give "NA/NaN function evaluation" while fitting
    the richer model, and no warning while fitting the simple one.
