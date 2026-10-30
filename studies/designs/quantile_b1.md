# Design: a faster quantile root-finder, never slower than bisection (quantile_b1)

Written before the method is built and before any run. It was settled on
1 October 2026, with its cap; a looser tolerance for calibration's widths was
not taken, so the stopping rule stays at 1e-12.

**Why.** After the CRPS fix (main b28b739), calibration at the profiling size
(gaussian, 2,000 rows x 1,000 draws, 1 core) takes 263 s. 99.7% of that is
`lmn_pred_quantile()`: eight calls of about 34 s, each a bisection of 41 or
42 mixture-CDF passes. The pinball and
interval scores (155 s) are all quantiles too.

A prototype of the Illinois method, run on 1 October, found:
- agreement with the bisection to 8.1e-13 everywhere;
- 8 to 27 passes against 41 on smooth mixtures;
- up to 92 passes against 43 on flat or staircase CDFs, which is what the cap
  removes.

## The method

For a continuous family, `lmn_pred_quantile()` keeps:
- today's bracket: the draws' own quantiles, checked and widened against the
  mixture's CDF (7c52c68);
- today's stopping rule: hi - lo <= 1e-12 max(1, |lo|, |hi|);
- today's answer: the bracket's midpoint, with hi always at F >= p and lo at
  F < p, so the answer is inf{x: F(x) >= p} to that tolerance.

Only the choice of the next point changes:

1. **Illinois.** The secant through the bracket's ends, with the value at the
   end that stayed put halved, so both ends converge. Where the secant point
   is not finite or not strictly inside the bracket, the midpoint is used.
2. **The end step.** A secant point within one tolerance of an end is moved one
   tolerance inside it, so the bracket closes on that side instead of
   creeping. This was the prototype's fix.
3. **Stalls.** A bisection step follows three steps that have not halved the
   bracket.
4. **The cap.** Bisection would take B = ceiling(log2(w0 / tol)) passes for a
   row, where w0 is its initial bracket width and tol its stopping width.
   - After each pass, if the passes so far plus ceiling(log2(w / tol)) for the
     current width w reach B, the row bisects from then on.
   - Bisection halves the bracket every pass, so the row finishes within B
     passes: never more than today.
   - The passes are taken across rows together, so a call never takes more
     passes than the bisection would.

The discrete families keep their integer bisection (`lmn_mix_qint()`), which
is not changed.

## Validation (pre-registered)

The script is `studies/scripts/quantile_b1.R`. It runs on the pinned library
mrp-442639f, on 1 core. The old function is kept as `lmn_pred_quantile_bisect()`
in the script only, copied verbatim from main 3841634, as the reference.

1. **Agreement with the bisection.**
   - The mixtures:
     - the six mixture shapes of crps_continuous (baseline, bimodal, scale
       spread, tail, near-degenerate, spiky);
     - check 4's added shapes there (weibull shape 0.3, beta(0.3, 0.3), the
       heavy loglogistic tails, sdlog 3, sdlog-15 spread);
     - 200 random mixtures per family, generated as in crps_continuous;
     - every continuous family (gaussian, lognormal, weibull, loglogistic,
       beta);
     - at S = 100, 1,000 and 4,000 draws.
   - At p = 1e-6, 0.001, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975,
     0.999 and 1 - 1e-6.
   - Seeds: 2950 for the shapes, 2951 for the random mixtures.
   - Rule: on every row, |q_B1 - q_bisect| <= 2e-12 max(1, |q_bisect|), twice
     the stopping tolerance, since each lies within one tolerance of the same
     point.
   - An NA in either is a failed case unless both are NA.
2. **The cap holds.** On every call, the new method's mixture-CDF passes are at
   most the bisection's, counted by tracing `lmn_mix()`.
   - Rule: never more, on any call.
3. **The bracket property.** On every row, F(q - 2 tol) < p <= F(q + 2 tol).
   - Rule: on every row where F is continuous there.
   - Rows where F is flat at p, where the quantile is the left end of the flat
     stretch, are reported separately.
4. **Downstream results are unchanged.** On the profiling data (the gaussian
   at 2,000 x 1,000, and the weibull AFT at 2,000 x 1,000):
   - `ilm_calibration()`'s coverage, widths and verdict text;
   - `ilm_scores()`'s pinball and interval scores.
   - Rule: the widths and scores within 1e-10 relative of today's, and the
     printed verdict text identical.
5. **limn's test suite** passes, with NOT_CRAN=true.
6. **Cost**, reported, at the profiling sizes against today, one run each,
   since all are over 5 minutes:
   - `ilm_calibration()` and pinball plus interval on the gaussian;
   - `ilm_calibration()` on the weibull AFT.

   The cap makes "slower" impossible in passes. Time is reported, with any
   arm above 1.05 times today named.

## The decision rule (approved with the design, before any build)

B1 merges only if all of these pass:
- checks 1, 2, 4 and 5;
- check 3 on every row except the flat-at-p rows, which are reported apart.

If any of them fails, B1 does not merge.
- The failing rows are reported.
- The method is not tuned and rerun under this registration. A revised method
  needs an addendum, committed before its run.

Check 6 (the cost) is reported, not a gate: the cap already rules out more
passes than today.

**The pilot.** If its times are off the estimate below by more than 1.5x,
they are reported before the full run.

## Addendum (2 October, before any build)

**The cap as registered could not be faster than bisection.** Its budget,
B = ceiling(log2(w0 / tol)), is exactly the number of bisection passes still
needed at the start. So "passes so far plus ceiling(log2(w / tol)) reaches B"
holds on the first pass and on every pass after: the row always bisects. A
strict per-row "never more passes than bisection" leaves no slack for a step
that might not pay.

**The new cap.** Each row's budget is bisection's count plus 2 passes:
B = ceiling(log2(w0 / tol)) + 2.
- A row takes a secant step only if, should that step gain nothing,
  bisection could still finish within the budget: passes so far + 1 +
  ceiling(log2(w / tol)) <= B.
- Otherwise the row bisects from then on.
- So no row takes more than 2 passes beyond the bisection's count.
- The ends' CDF values come from the bracket check, which evaluates them
  already. So the secant spends none of the slack on them.

**Check 2's rule** becomes: on every call, the new method's passes are at most
the bisection's plus 2.

Checks 1 and 3 to 6, and the decision rule, stay as registered.

## Pilot first, and the estimate

The largest pieces are piloted before the run: the gaussian calibration at
2,000 x 1,000 (about 263 s today, so about 80 to 110 s expected), and one call
at S = 4,000.

The full run is estimated at about 1 core-hour:
- check 1's quantile calls, about 25 minutes, mostly at S = 4,000;
- checks 4 and 6, about 25 minutes;
- the rest, small.

It stops at twice that, 2 core-hours.

## What is reported

- Each check's worst case, its rule and pass or fail.
- The pass counts per call: median and worst, against the bisection's.
- The costs.

NEWS: calibration and the quantile-based scores are faster, with the factor,
and their results are unchanged to 1e-10.
