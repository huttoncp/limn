# limn 0.0.0.9000

* limn is licensed under the GNU General Public License, version 2 or (at
  your option) any later version: GPL (>= 2).
* A prediction knows its time and series columns without a correlation over
  time (`ilm_predict()`, *Time and series columns*):
  - rows from `ilm_future()` say which they are, so a model with time only in
    its formula, a trend or a random slope, gets the trend check, and
    `ilm_scores(history = )` works for it;
  - `time` and `group` in `ilm_predict()` and `ilm_forecast()` now name columns
    of `newdata`, for rows made otherwise or changed by `merge()` or
    `subset()`. They used to take the times and groups themselves, for a
    correlation over time built from vectors; such a model now takes them as
    columns, which follow the rows through `set` and scenarios;
  - a correlation over time built by name keeps its own columns, and naming a
    different one is an error.
* `ilm_forecast()` takes `time` and `group` as arguments of its own. Passed
  on through `...`, `group =` was taken as `groups =`.
* The trend check also covers a random slope in time, which gives each unit a
  trend of its own, and rows past the data that did not come from
  `ilm_future()`. Grouping by time, as in `(1 | week)`, is not a trend.
* The help for cross-validation and model comparison cites its sources by
  section, equation and page, and the comparison plots say how to read them.
  The few-origins warning of `ilm_compare()` now cites Coroneo and Iacone
  (2020) by section (4.2 and 5) instead of by page.
* The help for prediction and forecasting lists every element of a
  prediction (`ilm_predict()`, *The prediction object*), cites its sources by
  section, explains how to read a fan chart, and gives `as.data.frame()` and
  `ilm_plot_fan()` examples. Claims the sources did not support were removed:
  that ETS intervals leave out parameter uncertainty, and that glmmTMB and
  lme4 give the same point forecast.
* The package help page (`?limn`) gives the workflow, explains how checks
  and remedies work, and has a worked example.
* Every check's source was read again at its page:
  - the few-groups check for means cites Kenward and Roger (1997) at
    pp. 983-984, not p. 983 alone;
  - the few-groups check for new units says that Maas and Hox (2005) cite
    Browne (1998) for restricted maximum likelihood being at least as good;
  - the typical-unit coverage the few-groups check quotes now gives binary
    outcomes' 0.92 beside the 0.84 to 0.86 of the other families;
  - the trend check labels its point about intervals as limn's own reasoning;
  - the Diebold-Mariano check pages Harvey, Leybourne and Newbold's Table 1.
* Checks whose remedy is a change to the model now name
  `illume::ilm_apply_remedy()` beside `illume::ilm_remedies()`.
* The help for scores and calibration defines every scoring rule as an R
  expression or formula. It cites Gneiting and Raftery (2007) by section,
  equation or table, and page for the log, CRPS, Brier, pinball, interval
  and squared-error scores, and fpp3 section 5.8 for the scaled errors; the
  absolute error is twice the pinball loss at 1/2. The
  PIT, its histogram and the calibration plot say how to read them, and
  `ilm_pit()` and `ilm_plot_calibration()` gain examples. The interval score
  is no longer called "Winkler's": Gneiting and Raftery trace it to Dunsmore,
  Winkler, and Winkler and Murphy. The calibration verdicts say which of
  their thresholds are limn's own choices.
* A test reproduces Gneiting and Raftery's closed-form gaussian CRPS, their
  interval score and their quantile score.
* A vignette, "Getting started with limn" (`vignette("limn")`), works through
  one example from fit to comparison: predictions for known and new units, a
  forecast, weekly totals, scores, and a calibration check whose verdict is
  acted on, the model refitted and checked again. The models are then
  compared on held-back weeks and over rolling origins. knitr and rmarkdown
  join Suggests, to build it.
* A second vignette, "Predictions: which unit, and how uncertain"
  (`vignette("predictions")`), works through `ilm_predict()`'s choices: which
  unit (a fitted, new, typical or average one), the mean or a new
  observation, averages over a population with `ilm_average()`, and how much
  of the fit's uncertainty to draw. It ends with the few-groups check and its
  remedy.
* The few-groups check for a new unit's observation no longer suggests
  refitting with `reml = TRUE` when the model was already fitted by restricted
  maximum likelihood. It still warns, and gives instead what limn's prototype
  study found with REML and 8 groups: 90% intervals covered 0.892 (Poisson)
  and 0.903 (gaussian).
* The help for `ilm_predict()`'s uncertainty choices now says that a plug-in
  prediction (`uncertainty = "none"`) still draws a new unit's effects, and a
  series' cells at times it was not fitted at. The few-groups coverage it
  gives for `"joint"` is now said to be for a new unit's observation, against
  the other methods the prototype tried, with its simulation error.
* A third vignette, "Forecasting series" (`vignette("forecasting")`),
  forecasts several series at once: the rows a forecast needs, the trend
  check and its remedy (a random walk in place of an autoregressive term),
  a covariate whose future is drawn from its own model or taken from another
  tool, fan charts, and totals across series and over the horizon.
* **`ilm_compare(method = "dm")` no longer drops to lag 0 when its long-run
  variance is negative.**
  - It follows Diebold and Mariano's (1995) own rule: the variance is taken as
    0 and the null hypothesis rejected, so the statistic is infinite and the
    p-value 0. It warns, naming `method = "wpe"`, the default, which needs no
    lag.
  - Results change only in that case.
* **Two generics let another package's fits use limn.**
  - `ilm_predict_draws()` is the step of `ilm_predict()` that draws the
    model's parameters, and `ilm_cv_refit()` the step of `ilm_cv()` that
    refits the model to a fold's training rows.
  - A package that fits illume's models another way, by sampling a posterior
    say, registers methods for its own class of fit. Predictions, forecasts,
    scores, calibration, comparison and cross-validation then work on its
    fits unchanged.
  - The methods for illume's fits are what limn did before, so every result
    is the same.
  - Draws in the wrong shape stop with a message saying what is wrong. Each
    `ilm_cv(ncores = )` worker also loads the package whose methods the fit
    needs.
* **A prediction's quantiles are found faster**, with the same answers.
  - The mixture quantile behind calibration, the pinball and interval scores
    and every interval now closes its bracket by a capped Illinois secant,
    not by bisection alone.
  - It never takes more than two passes of the mixture's CDF beyond the
    bisection's count, and its answer agrees with the bisection's to 1e-12.
  - On 2,000 rows by 1,000 draws, calibration took 0.62 of its earlier time
    on a gaussian model and 0.24 on a weibull survival model. The pinball and
    interval scores took 0.72.
  - Coverage, interval widths and calibration verdicts are unchanged.
* **The shared display rules write very large and very small numbers in
  scientific notation**, as in illumex 0.0.8.9003.
  - A number is written in full up to 15 digits before the decimal point and
    down to 6 decimal places. Beyond either, it uses R's form: 1.23e+300,
    1.23e-30.
  - The shared file also gains a rule for a group of values written to common
    decimals.
  - None of limn's usual prints changes. Only an estimate beyond those limits,
    such as a standard error of 1.23e-08, now reads that way, rather than as a
    long run of zeros.
* **Numbers are written by the family's shared display rules**, in
  `R/shared-helpers.R`, the same file byte for byte in each of the family's
  packages. Each number declares what it is:
  - an estimate (a score, a difference, a standard error, a width, a test
    statistic) to three significant figures with its trailing zeros kept:
    0.290, not 0.29;
  - a count (rows, draws, folds, clusters) whole, grouped from 1,000:
    "1,200 rows";
  - a share (a coverage, an observed or predicted rate) as a percentage to
    three significant figures, in prints and verdicts alike: "70.4% of
    outcomes", not "70%". A coverage's standard error is an estimate in
    percentage points: "(SE 1.77 percentage points)". An interval's level
    ("95%") is a setting and is written as given;
  - a p-value by the family's own rule.

  Rounding is half away from zero on the value's decimal form, not on its
  binary one. `ilm_cv()`'s print now shows its mean scores to three
  significant figures, where it printed seven.
* **Fix: the CRPS of a continuous prediction is now exact, or integrated to a
  stated error, for every continuous family.**
  - **What was wrong.** It was integrated with 16 equal panels over the range
    of all the draws. The panels were sized by the widest component, while
    the integrand's detail is as narrow as the narrowest one.
    - Gaussian mixtures whose components differ widely in spread were off by
      up to 0.77% relative: sds from 0.01 to 10 across draws.
    - So were tiny-sd components with spread means, such as a forecast whose
      latent path varies far more than its noise: 0.1%.
    - On the positive families (lognormal, weibull, loglogistic), mixtures
      whose scale parameters spread up to 60-fold across draws were off by
      14% to 100% relative at small and middling outcomes.
    - A loglogistic's upper tail beyond its 1 - 1e-10 quantile was dropped:
      6.6e-4 of the CRPS at s = 1.67 and 17% at s = 1.92.
    - On the shapes typical predictions have, the error was at most about
      1e-11. Beta was within 1e-7.
  - **Now.**
    - A gaussian mixture's CRPS is the exact closed form (Grimit et al.
      2006), as `scoringRules::crps_mixnorm` gives it. The old quadrature is
      kept only for rows it has been shown to integrate to rounding.
    - Lognormal, weibull and loglogistic are integrated on log y, and beta on
      logit y, with panels sized by the narrowest component.
    - The upper tail comes from the survival function, never 1 - F. A
      loglogistic is integrated out to a stated bound on what is left.
  - **Edge cases.**
    - A loglogistic component with scale 2 or more makes the CRPS infinite,
      and it is returned as Inf, with a message.
    - A row whose integration fails is NA, named in a message, never a number.
    - A continuous family without a CRPS method is an error, not a fallback.
  - **Cost**, at 2,000 rows by 1,000 draws on one core: 1.19 times before
    for a gaussian prediction, 1.68 times for a 480-row gaussian forecast,
    and 1.52 times for a weibull AFT prediction.
  - Validated in `studies/findings/crps_continuous.md`, which reports all
    four runs. The first three found a bug in the study script, R's `qbeta()`
    failing at tiny shape parameters, and a reference that silently missed
    narrow mass, each fixed before the next run.
* **Fix: quantiles no longer trust a family's own quantile function blindly.**
  - The bisection behind intervals, the pinball and interval scores,
    calibration widths and forecast bands took its bracket from the draws'
    own quantiles.
  - Where those are wrong, as R's `qbeta()` is at tiny shape parameters (off
    by up to 0.5 in probability), it returned a wrong quantile with no
    warning.
  - The bracket is now checked against the mixture's CDF and widened until it
    holds. A row that cannot be bracketed is NA, with a message.
* `ilm_calibration()` on a prediction is now public, as a method of
  illume's generic, as `ilm_scores()` is. Until now it was internal.
  - It gives coverage and width by level, bias on the outcome's scale, PIT
    uniformity, and for a binary outcome the calibration intercept and slope.
  - It ends with a verdict and a remedy, such as `family = "nbinom"` when
    Poisson intervals are too narrow.
  - `ilm_plot_calibration()` draws it.
  - Two verdicts are reworded:
    - A biased mean is now reported first. A shortfall in coverage beside it
      is reported as one the bias can cause, and the remedies that widen
      intervals wait until the mean is fixed.
    - A calibration slope whose interval lies below 0 is reported as
      reversed: the predictions order the outcomes the wrong way. Before, it
      was only called overconfident.
* Predictions, totals, comparisons and cross-validation results keep the
  diagnostics they raise (`$checks`, or the `checks` attribute), and a
  prediction keeps the draws asked for and its seed with the generator's kind.
  A forecast names its time and series columns and keeps each row's origin:
  the series' last observed time. The warnings are raised as before.
* `ilm_scores()` scores a prediction by proper rules, and the scaled errors
  (MASE, RMSSE) take their scale either way:
  - `history = train`, read through the columns the prediction now records
    (its response, time and series);
  - `naive = ilm_naive_scale(...)`, now exported, for a prediction with no
    time column, or one scale shared across models.

  The two give identical scores. When `history =` can't work, the error
  gives the `naive =` call to use. The help of both says when to use each.
* The discrete CRPS on wide supports is much faster.
  - A row wider than 2,000 values (was 10,000) is summed exactly within 50 of
    the outcome, and elsewhere by Euler-Maclaurin with Simpson's rule, to an
    error estimate of at most 1e-6 of the score.
  - Each row records its error estimate, a flag where it's exceeded, and the
    rigorous bracket width.
  - Checked against the exact sum on 105 rows, including uneven mixtures,
    outcomes in the far tails, steep near-Poisson rows and widths either side
    of 2,000: the largest relative error is 1.1e-7, and every error lies
    within its estimate.
  - On large counts (PBS scripts) it is 3 to 27 times faster than the block
    bracket it replaces.
* A prediction whose distribution cannot be computed in some draws now stops
  with limn's "draws not computed" message and the count of draws affected,
  in two cases that slipped through before:
  - a zero-inflated model whose zero part is not a number, where the family's
    quantile function stopped with an unexplained error;
  - an ordered outcome, whose quantile is NA rather than NaN in such a draw,
    and which passed unflagged.
* The discrete CRPS is faster again: 2.5 times on small counts and 1.9 times on
  tourism's large counts (medians of 5 alternated runs), with the same result to
  rounding. The
  support now runs between the mixture's own 1e-6 and 1 - 1e-6 quantiles,
  where the terms left out are at most 1e-12 each. Before, it was the
  widest of every draw's 1e-12 quantiles. A row whose outcome falls outside
  that support keeps the old bounds, so it is scored exactly as before.
  Rows summed exactly agree with the old code to 1.2e-11 relative, and
  rows too wide to sum exactly now more often need no bracket.
* `ilm_average()` takes weighted averages of predicted means, draw by draw,
  within groups of rows: standardisation over a reference population, or
  post-stratification with population cell sizes as weights (as in MRP).
  Weights are a number per row or a column name, normalised within each
  group, and scenarios are never averaged together. It takes a prediction of
  the mean. `ilm_aggregate(fun = "mean")` remains the verb for unweighted
  averages of new observations.
* `ilm_plot_fan()` passes `facet.args = list(free = TRUE)` and `ylim` on to
  tinyplot, so each series' panel can have its own y scale, which suits series
  of very different sizes. Before, the scale was fixed, and a `ylim` given
  through `...` was an error. Every scale holds the history as well as the
  bands.
* The CRPS of a count forecast is about 8 times faster for large counts, with
  the same result to rounding (2.6e-15 relative on 608 tourism forecasts;
  7,232 s against 917 s, measured on a loaded machine). The mixture's CDF
  comes from running sums of its mass function, and rows are grouped by
  support width, so a narrow row is no longer padded to the widest.
* `ilm_cv()` scores a model fold by fold on data it was not fitted to: it
  refits the model to each fold's training rows (`illume::ilm_refit()`,
  E6), predicts the test rows and scores them. Grouped folds predict their
  held-out units as new ones; other folds predict each unit as itself.
  Rolling folds report origin and horizon, and whether the evaluation is
  ex-ante or ex-post, from the covariates that vary within a series and are
  not named in `ilm_folds(known = )`. Random or blocked folds on a model
  with a correlation over time warn that they score interpolation. Folds
  can be shared out to worker processes (`ncores`), with the same results
  whatever the number. The results go straight into `ilm_compare()` and
  `ilm_plot_horizon()`. It needs illume 0.0.8.9001 or later for
  `ilm_refit()`, which DESCRIPTION now requires.
* `ilm_compare()` gives its test as data beside the `method` description: a
  `test` code (`"independent"`, `"clustered"`, `"wpe"` or `"dm"`), the WPE
  test's `bandwidth` m, and the DM test's `lag` as used.
* `ilm_compare(by = "origin")` uses the fixed-smoothing Diebold-Mariano
  test of Coroneo and Iacone (2020) by default (`method = "wpe"`): the
  periodogram at the first floor(T^(1/3)) Fourier frequencies of the origin
  means, referred to t on 2m degrees of freedom. The usual test stays as
  `method = "dm"`, and warns with fewer than 100 origins over several
  horizons. Comparisons over several horizons now give a row pooled over them
  and a row for each horizon.
* Covariate paths without building `.draw` blocks by hand:
  - `ilm_covariate_paths()` simulates a covariate's future paths from an
    illume model of it. Each path is one parameter draw, the series' latent
    path continued, and a new observation at each time, so the paths carry the
    covariate model's own uncertainty. It returns the rows stacked into
    `.draw` blocks, ready for `ilm_forecast()`. Several covariates are drawn
    from their own models, independently.
  - `ilm_as_draws()` turns paths simulated elsewhere (fable's `generate()`,
    an ensemble) into the same blocks, matching rows on `by` columns.
* Covariates can be given as draws: a `.draw` column in
  `newdata` marks one block of rows per draw of the covariates, every block
  listing the same rows in the same order. Parameter draw s uses covariate
  draw ((s - 1) mod C) + 1, so the prediction mixes over both, and `nsim` is
  rounded up to a multiple of C, with a message.
  - Blocks that differ in their scenario, units, series or times are refused.
  - The prediction has one row per row of a block; drawn covariates show as
    `NA` where they vary, `meta` records C and the drawn columns, and a drawn
    row's estimate is its predictive median.
  - With every parameter held, a gaussian row's predictive distribution is the
    mixture over the covariate draws to 1e-8.
* Scoring and quantiles of predictive mixtures are faster: the family is
  called once per block of draws rather than once per draw, which matters for
  small evaluations such as the quantile search. A predictive quantile of 160
  rows from 1,000 draws is about 10 times faster. Results are unchanged, to
  1e-12 in tests across Poisson, negative binomial, zero-inflated and
  gaussian-with-dispersion fits.
* A new observation whose family overflows in some draws now stops, in the
  registry's words (`nan_family`), instead of giving NaN quantiles and scores.
  It happens when a parameter sits at the edge of its range and the fit did not
  hold it there: a residual SD near zero, or an AR correlation near 1. The
  remedy named is `uncertainty = "fixed"` or `"latent"`.
* The discrete CRPS stays an exact sum up to 10^4 support values per row.
  - Beyond that, it sums each of the integrand's two monotone pieces in blocks
    of d values and takes the midpoint of each block's bracket. The error is at
    most d, reported per row in the scores' `crps_bound` attribute.
  - On a Poisson mixture with mean 10^6 (about 1.3 x 10^5 support values,
    d = 13), the error was 2e-4, or 1e-7 relative, against a stated bound of
    1.27.
  - The support runs between each draw's 1e-12 and 1 - 1e-12 quantiles. The
    tails beyond it add at most 1e-12 times the outcome's expected distance
    past its end.
* `ilm_aggregate()` totals or averages a prediction over groups of its rows:
  a series' forecasts over a horizon, every series at each time, or all rows.
  - Totals are taken draw by draw. A series' path and the shared parameters
    and effects carry into them, and totals at different levels add up.
  - A total of new observations from a Poisson or gaussian model is an exact
    mixture, as in the decision "exact where the family allows".
  - Other families' totals are summed from four outcomes simulated per draw
    and held as samples. They are scored by the fair CRPS, quantiles and PITs,
    with no log score. The printout says which method was used.
  - An average is the total over the number of rows and is scored through it,
    so an average of counts keeps its whole-number steps.
  - Scenarios are never added together.
* `uncertainty = "latent"` in `ilm_predict()` and `ilm_forecast()` holds
  every parameter at its estimate, the dispersion included, and draws only the
  random effects and the cells of a correlation over time, given the
  parameters (`illume::ilm_draws(given = "parameters")`). For a gaussian model
  with a correlation over time, that is the Kalman filter's forecast, which
  treats the parameters as known, as ETS and ARIMA intervals do.
* Every check limn makes now takes its words from one registry.
  - Each check says what was found, what it means and why the remedy helps,
    the remedy and how to apply it, and a reference where one exists.
  - Warnings, errors, the calibration verdicts and the help pages (the Checks
    section of `?ilm_predict`) all compose their text from the same entries.
    The entries' fields are `illume`'s check fields plus a reference.
  - Warnings and errors carry the filled-in fields, and a class naming the
    check (`lmn_diagnostic`, `lmn_<check>`), so code can catch and read them.
  - A printed prediction now notes rows that name a level the fit has not
    seen.
* Forecasting: `ilm_forecast()`, and `ilm_predict()` at new times. A
  correlation over time (`ilm_ar1()`, `ilm_car1()` or `ilm_rw1()`) is carried
  to new times from the fit's own parameters, with no refit.
  - A seen series continues from its fitted cells: forward past its last cell,
    as a bridge between two cells, and backward before its first.
  - A new series is a fresh path. AR(1) and CAR(1) paths start from their
    stationary distribution, and a random walk starts from zero at its first
    time.
  - A series' new times are drawn as one path per draw, and rows naming the
    same series and time share a value. So totals over a horizon, and
    contrasts between scenarios, keep the path's correlation.
  - `ilm_forecast()` is `ilm_predict()` with `target = "response"` and
    `groups = "fitted"`.
  - For a mean or the linear predictor at a seen series' new time, the
    estimate is the point forecast at the estimates: the fitted cells carried
    forward with no new shocks, as glmmTMB, lme4 and a Kalman filter give.
  - A forecast past the data from a model whose fixed part trends in time
    warns. Such a trend is assumed to carry on unchanged, so the intervals
    leave out any change in it (fpp3 section 10.4). The warning names a random
    walk, `ilm_rw1()`, as the stochastic trend that lets them widen. Periodic
    terms such as `ilm_fourier()` are safe to forecast and don't count as
    trends. The time column comes from the correlation over time, or from
    `ilm_future()`, which now records its time and series columns.
  - `as.data.frame()` puts the expected value beside the estimate as `mean`.
    This is what fable reports as its point forecast, and it is what totals
    and squared error need. For a Poisson AR(1) series with latent SD 0.77,
    it was 31% to 42% above the estimate at horizons 1 to 8.
  - Each conditional step equals conditioning the process's own covariance to
    1e-12, and the paths' draws have the joint conditional distribution of the
    new times given the cells.
* `ilm_predict()`, the prediction engine. It gives the predictive distribution
  of the mean, a new observation or the linear predictor for each row of new
  data.
  - Parameter uncertainty is drawn jointly by `illume::ilm_draws()`, or with the
    variances held, or left out.
  - `groups` names how the random effects enter, and has no default for a
    mixed model:
    - `"fitted"`: a unit's own effects;
    - `"new"`: an unseen unit's, drawn afresh in every draw and shared by that
      unit's rows;
    - `"typical"`: effects at zero;
    - `"population"`: the mean averaged over units, by `illume`'s own
      quadrature.
  - A level not seen in the fit is predicted as a new unit and marked.
  - A variance at the edge of its range, whose draws are not finite, stops a
    prediction that needs them. The message names `uncertainty = "fixed"` and
    `illume`'s `boundary = "avoid"`. So does a draw that cannot be computed.
  - With fewer than about 20 levels of a random term, a prediction warns that
    its intervals cover below their nominal level, for means too.
    - For means and linear predictors, whatever the unit, the warning says they
      cover as z-intervals do.
    - For a new unit's observation, it names `reml = TRUE`, which narrows the
      gap.
  - `newdata` may be a named list of scenarios, and `set` gives every row
    several values of a variable.
  - Response predictions are the same objects the scores, PITs, calibration and
    fan chart already use. Forecasting past a series' observed times comes
    with the forecasting engine.
* The first development version: the package skeleton, and the parts that need
  nothing from `illume` beyond what it exports today.
* `ilm_folds()` builds cross-validation folds as row indices: random, grouped
  (whole groups held out, for predictions of new units), blocked (contiguous
  time, with an optional gap either side), and rolling origins (expanding or
  sliding windows, a horizon, a gap, and the covariates known in advance).
* `ilm_future()` builds the rows a forecast needs: the next `h` time points
  after each series' last one.
  - Columns constant within a series are copied forward.
  - Columns that vary over time are left missing and named, not guessed.
  - Time steps come from the data. For dates they are calendar units, so a
    month-end stays a month-end.
* `ilm_plot_fan()` draws predictive distributions over time as a fan chart.
  - The median is a line, with nested central-interval bands.
  - Each series gets its own panel, with the observed history drawn first.
  - It returns the data behind the plot.
* Diagnostic plots, each returning the data behind it:
  - `ilm_plot_horizon()`: each model's mean score by forecast horizon, with
    standard errors clustered by origin, as in fpp3's accuracy-by-horizon
    plots.
  - `ilm_plot_compare()`: each comparison's difference with an interval from
    the same standard error as its test.
  - `ilm_plot_calibration()`: coverage against nominal level or, for a binary
    outcome, the reliability diagram with Wilson intervals, with the
    calibration intercept and slope beneath it.
* `ilm_pit()` gives probability integral transforms. For discrete outcomes they
  are randomised, so they are uniform exactly when the predictive distribution
  is right.
* `ilm_compare()` compares two or more models' scores on the same rows, with
  standard errors for independent rows, clustered by a column (the HC1
  cluster-robust estimator), or over rolling origins by the Diebold-Mariano
  test with the Harvey-Leybourne-Newbold correction.
* `ilm_plot_pit()` draws the PIT histogram, with a band showing where each bar
  falls if the predictions are right. It uses tinyplot, as illumex does.
* Calibration of a predictive distribution (internal until `illume` makes
  `ilm_calibration()` a generic):
  - PIT uniformity, and coverage and width by level. Coverage is read from the
    PITs, so it is exactly nominal for counts too when the model is right.
  - Bias, on the outcome's own scale.
  - For binary outcomes, the calibration intercept and slope (Van Calster et
    al. 2019).
  - A verdict with a remedy that exists, in `illume` where most are: for
    example `family = "nbinom"` when Poisson intervals are too narrow.
* Two scaled errors join the scores:
  - `mase`: the absolute error over the series' in-sample naive error.
  - `rmsse`: the squared error over the mean squared naive error.

  The naive error can be seasonal. Following fpp3, these make accuracy
  comparable across series of different scales.
* Predictive distributions are scored from their mixture over parameter draws,
  with no response simulated. The scores are log, CRPS, Brier, pinball,
  interval, squared and absolute, and all are losses. Scoring is
  internal until `illume` makes `ilm_scores()` a generic, which `limn` will
  then extend.

## Validation

* A richer limn specification, registered in advance, closes most of the gap
  to fable's MinT on Australian tourism (`studies/findings/aggregate_mint.md`).
  It adds a common drift and region-specific seasons (a random
  region-by-quarter effect) to the simple model.
  - Mean CRPS against MinT's: 27.9 against 27.8 over the regions, 188.5 against
    208.3 over the states (better), and 1,551 against 1,494 for the national
    total. The simple model had 45.0, 352.9 and 3,040.
  - It covered better than MinT at every level: at 80%, 0.81, 0.73 and 0.25,
    against 0.71, 0.45 and 0.25.
  - The simple model's national totals were biased low. The richer model's sit
    with MinT's, both a little below the outcomes.
  - fable's unreconciled ETS was still best at the state and national levels.
    One split of one dataset: a benchmark, not a test.
* The fixed-smoothing Diebold-Mariano tests hold their size with few origins
  (`studies/findings/compare_dm_fixed.md`, R = 1000 per cell, T from 10 to
  200). Over 40 null cells the 5% rejection rate of the pooled test was
  0.027 to 0.119 for the usual test, 0.027 to 0.068 for WCE and 0.027 to 0.059
  for WPE. At a single horizon: up to 0.126, 0.081 and 0.068. The statistics
  equal `sandwich::lrvar()` (WCE) and `stats::spec.pgram()` (WPE), and the
  critical values reproduce Coroneo and Iacone's published ones at T = 20.
* Covariate paths drawn by `ilm_covariate_paths()` from an illume model of the
  covariate cover as well as paths from its true process
  (`studies/findings/covariate_draws.md`, "generated" arm; R = 500 per family,
  on a pinned library). At 80% they covered 0.807 to 0.821 (gaussian) and
  0.798 to 0.806 (Poisson). They run slightly wide at short gaussian
  horizons, since they also carry the covariate model's own parameter
  uncertainty. With every parameter held, the paths are the covariate's Kalman
  filter forecast.
  - 25 of 500 Poisson replicates failed in that arm. The covariate is
    noise-free, so its model's residual SD (about 0.02 of the covariate's SD)
    is not held by illume. Its joint draws reach 1e37, and a few paths
    overflow the Poisson forecast. The gap is illume's to close.
* Forecasts that mix over covariate draws cover at their nominal level, where
  forecasts conditional on the covariate's point forecast do not
  (`studies/findings/covariate_draws.md`). The design is 20 series whose
  covariate is itself an AR(1), forecast 1 to 6 steps ahead with joint
  parameter draws, 500 replicates per family.
  - At 80%, drawn covariates (100 paths) covered 0.799 to 0.806 (gaussian) and
    0.799 to 0.805 (Poisson), as the realised covariates did (0.806 to 0.811,
    0.796 to 0.807).
  - Conditioning on the point forecast covered 0.708 falling to 0.634
    (gaussian) and 0.785 to 0.767 (Poisson), worst at longer horizons, where
    the covariate is least known.
  - Drawn covariates' intervals were 22% to 44% wider (gaussian) and up to 11%
    (Poisson), which is what the covariate's uncertainty adds.
  - 2 of 1,000 replicates failed in every arm: a random-intercept SD at 0.02
    that the fit did not hold, and a Hessian that was not positive definite.
* The pooled Diebold-Mariano test of `ilm_compare(by = "origin")`
  (`studies/findings/compare_dm_pooled.md`). The design is an
  AR(1) series forecast 1 to H steps ahead from T origins, with 1000
  replicates per cell.
  - With one horizon it holds its size: 0.038 to 0.063 at 5%.
  - With 2 or 4 horizons and 50 origins or fewer it over-rejects, at up to
    0.12. By 100 or more origins it is at 0.037 to 0.069.
  - The single-horizon test at the longest horizon, which is
    `forecast::dm.test()`, over-rejects by as much, so the excess is the
    Diebold-Mariano test's own with few origins, not the pooling.
  - The pooled test's power was much higher: 0.53 against about 0.20 for
    either single horizon, at 200 origins with one forecaster's noise 25%
    larger.
* Totals from `ilm_aggregate()` cover at their nominal level
  (`studies/findings/aggregate_coverage.md`).
  - The design is 20 AR(1) series forecast 4 steps ahead, gaussian, Poisson
    and negative binomial, with 500 replicates per family. The variances are
    held, since joint draws wait on `illume`'s boundary holds.
  - At 90%, joint totals covered 0.888 to 0.892 for each series over the
    horizon, 0.891 to 0.904 across series at a step, and 0.873 to 0.905 for
    the grand total.
  - Summing each row's forecast on its own under-covers where rows are
    correlated: 0.74 to 0.79 per series and 0.70 to 0.74 for the grand total,
    with intervals 22% to 34% narrower. Across series at one step it makes
    little difference (0.87 to 0.90).
  - When outcomes cannot be simulated (a dispersion at its limit that
    `illume` does not hold), a total now stops and names
    `family = "poisson"`.
* Against fable on Australian tourism (`studies/findings/aggregate_mint.md`,
  a benchmark, not a test), limn's simple global model does worse, reported as
  found. The model is one negative binomial model with random-walk levels and
  a quarterly pattern common to all regions.
  - Its CRPS was 1.6 to 2 times MinT's at every level.
  - Its national totals were biased low. They sat below the last training
    year while 2016-17 grew, so they covered 0 of 8 quarters at 80%. Its
    intervals were wider than MinT's, so the misses are bias, not
    overconfidence.
  - The model has no trend, and a smoothed level lags the rise.

* With every parameter held, forecasts are the Kalman filter's
  (`studies/findings/forecast_kalman.md`). Gaussian AR(1), CAR(1) at
  irregular times, and random-walk panels, forecast 1 to 6 steps ahead with
  `uncertainty = "latent"` and compared with `stats::KalmanForecast()` at
  illume's estimates.
  - Point forecasts agree to 9e-16.
  - With 40,000 draws, the draws' means are within 1.5 Monte Carlo SEs, their
    SDs within 0.994 to 1.006 of the filter's, and the correlation between
    horizons 1 and 6 within 0.008 of the filter's.
  - Against the filter, "none" is too narrow at short horizons (0.79 to 0.92
    of its SD at horizon 1), since it leaves out the last cell's uncertainty;
    "fixed" and "joint" are up to 5% wider, adding the parameters'.
* Forecast coverage by horizon, by simulation
  (`studies/findings/forecast_coverage.md`). The design is 20 series with
  random intercepts, observed at 8, 16 or 32 times and forecast 1 to 4 steps
  ahead, with 500 replicates per cell. The cells are AR(1) for gaussian and
  Poisson, an AR(2) fitted as AR(1), and a Poisson random walk. On the
  replicates where every setting ran, 90% intervals covered:
  - **Plug-in:** 0.850 to 0.896, under-covering everywhere.
  - **With the latent's variances held:** 0.884 to 0.921.
  - **With joint draws:** 0.889 to 0.917, except for gaussian series of 8
    times, where joint over-covers at 0.943 to 0.953. Eight times barely
    separate the latent from the noise, and the Gaussian approximation to
    their variances puts too much weight on large ones.
  - **With 32 times,** both settings are at 0.895 to 0.905.
  - **The misspecified AR(2) cell** covers 0.893 to 0.904 under joint draws.
  - **Boundary fits** stop joint draws until `illume` holds them: 36% to 40%
    of fits with 8 times, 13% to 19% with 16, and 0% to 4% with 32. On those
    replicates the held-variance intervals cover 0.85 to 0.91.

* Forecasts have the exact conditional distribution
  (`studies/findings/forecast_exact.md`). The test data are gaussian panels of
  8 series, with AR(1), CAR(1) or a random walk.
  - **Latent parameters held.** Forecasts at horizons 1 to 4, and values
    between observed times, match exact multivariate-normal conditioning at
    `illume`'s estimates within Monte Carlo error. Means are within 2.6 MC SEs,
    and SDs 0.992 to 1.014 times the exact ones, from 40,000 draws.
  - **The reference itself.** glmmTMB, with its variances mapped to the same
    estimates and the new rows at weight zero, reproduces it to 1.7e-15 for
    AR(1) and 3.9e-10 for CAR(1).
  - **Joint draws.** At the last observed time, limn and glmmTMB's delta
    method agree, with SD ratios of 0.995 to 1.007. Ahead of it they part: at
    one step, limn's spread is up to 4% wider and its centre up to 0.1 SD
    away, fading by four steps. limn draws each new time from its exact
    conditional given each draw's parameters, while the delta method
    linearises in them.

* `ilm_predict()` gives an independent implementation's predictions at the
  same parameters (`studies/findings/engine_glmmTMB.md`). The comparison is
  glmmTMB, with every parameter held at `illume`'s estimates, for gaussian,
  Poisson and binomial random-intercept fits.
  - Known units' means agree to 5e-12 or better, 2e-16 for a gaussian, and
    typical units' exactly. A unit's effect solved by Newton's method in base R
    agrees to the same, so the remainder is the fit's own tolerance.
  - The joint draws' spread matches glmmTMB's standard errors within Monte
    Carlo error: 0.980 to 1.012 times, for typical units (the delta method) and
    known units (TMB's generalised delta method), from 5000 draws.
  - Population means equal their exact forms to 1.2e-15: exp(eta + s^2/2) for a
    log link, and `integrate()` for a logit.
  - Averaged over the sample, they reproduce `illume::ilm_scenario(over =
    "sample")`, for random intercepts and random slopes. The estimates are the
    same, and the interval limits agree within Monte Carlo error over 12 pairs
    of seeds.
* The coverage of `ilm_predict()`'s intervals, by simulation
  (`studies/findings/engine_coverage.md`). The fits have random intercepts, 6
  rows per group and 8 or 40 groups, across five families, with 1000
  replicates per cell.
  - **New observations.** With joint uncertainty, central 90% intervals
    covered:
    - 0.891 to 0.900 with 40 groups, for seen and new units in every family;
    - 0.887 to 0.904 for seen units with 8 groups;
    - 0.872 to 0.896 for new units with 8 groups. Poisson is lowest, which is
      B2's few-groups shortfall.

    Plug-in intervals covered 0.845 to 0.898, and holding the variances fell in
    between, so the parameters' uncertainty is needed.
  - **Means with 40 groups.** 95% intervals covered 0.933 to 0.950, except a
    binomial seen unit's mean at 0.914. With six 0/1 rows per unit, the
    Gaussian approximation to a unit's effect is poor.
  - **Means with 8 groups.** Intervals under-cover as z-intervals do, which is
    also how `illume`'s `confint()` works. A typical unit's 90% interval
    covered 0.838 to 0.861 in four families; a z-interval covers 0.856 when the
    truth follows t with 7 df. A binomial new unit's mean over-covers (0.973
    at 90%), because its variance is barely determined.
  - **Boundary fits left out.** The joint cells exclude fits whose variance sat
    at the boundary: 27% of binomial fits with 8 groups, and 4% or fewer
    elsewhere. `illume` is fixing the hold that makes their draws usable, and
    the study will be rerun then.

* Scores match `scoringRules`, measured as the largest relative difference over
  200 outcomes:
  - the log score to machine precision (1.4e-15 or better);
  - the CRPS to 1.2e-12 for a gaussian and 2.8e-15 for a Poisson predictive
    distribution;
  - for 100-draw mixtures of normals, from tight to wide, against its closed
    forms: the log score and the CRPS to 4.9e-12 or better.

* The negative binomial CRPS equals the exact energy form, summed over the
  probability mass function, to 4.5e-14 in 45 cases (3 sizes × 3 means × 5
  outcomes). That is a different computation from limn's sum over the CDF.
  - `scoringRules`' own negative binomial CRPS agrees wherever it is finite.
  - That closed form (1.1.3, through `hypergeo` 1.2.14) returns −Inf in 10 of
    the 45 cases: size 2.5 with a mean of 2 or more. So it is compared only at
    sizes where it is finite.
  - The negative binomial log score matches `scoringRules` to machine
    precision.
* A clustered standard error equals `sandwich::vcovCL()` to 1e-12, and one over
  independent rows equals the paired t-test.
* `ilm_compare(by = "origin")` equals `forecast::dm.test()` over 600
  comparisons with correlated errors: to 3.6e-15 in the statistic and 3.3e-16
  in the p-value. The comparisons covered 20 to 200 origins, horizons 1 to 4,
  and absolute and squared losses (`studies/findings/compare_dm.md`).
* Rolling-origin folds are identical, row for row, to
  `rsample::rolling_origin()` in all 48 configurations tried: expanding and
  sliding windows, and several horizons and steps.
* The point and interval measures agree with `fabletools::accuracy()` on the
  same forecasts (`studies/findings/accuracy_fabletools.md`).
  - For normal forecasts, all agree exactly: MASE and RMSSE (seasonal and
    not), MAE, RMSE, the 80% Winkler score, and the quantile score, which is
    twice limn's pinball loss.
  - For 20-draw normal mixtures, the mean-based measures agree to 8e-16. The
    quantile-based ones differ by 1.5e-6, which is `distributional`'s
    root-finding tolerance: at the quantiles, the exact mixture CDF is within
    3.9e-12 of its target for limn and within 4.1e-6 for `distributional`.
* Grouped folds have `rsample::group_vfold_cv()`'s properties in all 180 random
  datasets: no group split, every group held out once, and group counts within
  one of each other (`studies/findings/folds_rsample.md`).
