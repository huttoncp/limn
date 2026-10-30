# Design: the capped quantile root-finder, again (quantile_b1_retry)

Written before any run of the retry. On 5 October 2026 it was settled to try the method again,
under a new registration on new seeds, with the three causes of the first run's failure fixed, and the
same decision rule.

The first registration is `studies/designs/quantile_b1.md`, with its decision rule and addendum (the cap
at the bisection's count plus 2). Its findings are `studies/findings/quantile_b1.md`. Everything
there stands unless this design changes it below.

## What the first run found, and what changed

The first run met check 1 on every row: agreement with the bisection to 9.9e-13. It failed three checks,
each for a reason found in the run:

1. **Check 2, the cap: 41 calls took the bisection's passes plus 3.**
   - The cause was in the method. The guard asked whether bisection could still finish within the
     budget, counting passes from the bracket's current stopping width, 1e-12 max(1, |lo|, |hi|). A
     bracket closing towards zero has a shrinking stopping width, so that count could fall one pass short.
   - **The fix** (a14dca6): the guard counts from the smallest stopping width reachable inside the
     bracket, 1e-12 max(1, min |x| over [lo, hi]), with 1e-3 added to the pass count for the midpoints'
     rounding.
   - The budget is unchanged: the starting width's count plus 2. The starting width is the largest the
     bracket allows, so that count is a lower bound on the bisection's own.
   - So no row can take more than the bisection's passes plus 2. Each pass a row spends is either within
     a guard that held, or a bisection step that the guard's count already covers.
   - A test (`tests/testthat/test-pred.R`) holds two three-draw mixtures that overran under the old
     guard, by 45 against 42 and 43 against 40. It requires the bisection's passes plus 2 and the
     bisection's answer.
   - In a probe on 3,146 calls (on a probe seed, not this study's), no call passed the cap. The total
     passes were 0.778 of the bisection's, against 0.734 for the old guard.
2. **Check 3, the bracket property: 4 rows could not be read.**
   - The cause was in the harness. The check evaluates F at q - 2 tol. For a loglogistic quantile near
     0 that falls below 0, where the harness's loglogistic CDF took log(q) and gave NaN.
   - **The fix** (30482dd): the harness's loglogistic CDF is 0 below 0, as plnorm() and pweibull() are.
     A self-check stops the run unless every harness CDF is a number on the whole line, 0 below its
     support and 1 above.
   - illume's own loglogistic family has the same flaw, NaN with a warning below 0. That is reported to
     illume and is not part of this study.
3. **Check 4, downstream: the gaussian pinball and interval scores differed by 1.3e-9 relative, against
   1e-10.**
   - The cause was the rule's form. The quantiles agreed to 1e-13, but a pinball loss near 0 makes any
     difference large relative to itself.
   - **The fix** (b984d36): scores are compared as |difference| / max(1, |score|), the form of the
     quantiles' own rule. The calibration widths and coverage keep the relative rule. A self-check stops
     the run unless the rule passes a 1.3e-13 difference near 0 and fails 1e-9 at every size.

Nothing else changes: the method, the checks, the sizes and the decision rule. The script is
`studies/scripts/quantile_b1_retry.R`, a copy of the first run's script with the two harness fixes, the
new seeds and the stop. The first run's script stays as its record.

Not part of this registration: a budget that grows as the bracket closes, which would keep more of the
slack. It was measured in the same probe (0.596 of the bisection's passes) and is not in the method.

## The seeds

Drawn on 6 October 2026 at 10:41 with Python's `secrets.randbelow(2^31 - 1)`, before any fit:

| use | first run | this run |
|---|---|---|
| the shapes, checks 1 to 3, and the pilot's S = 4,000 call | 2950 | 248791001 |
| the 200 random mixtures per family | 2951 | 382963453 |
| the profiling data for checks 4 and 6 (gaussian, then the weibull AFT from the same stream) | 1 | 1289645533 |

## The sizes

As registered for the first run:
- the six shapes and check 4's added shapes, at S = 100, 1,000 and 4,000;
- 200 random mixtures per family;
- 5 continuous families and 13 probabilities;
- the profiling data at 2,000 rows by 1,000 draws, gaussian and weibull AFT.

The library is mrp-442639f (illume 442639f), as before, on 1 core. The run uses a frozen copy of the script
as a detached process, one long run at a time.

## The decision rule (unchanged, word for word from the first registration and its addendum)

B1 merges only if all of these pass:
- checks 1, 2, 4 and 5;
- check 3 on every row except the flat-at-p rows, which are reported apart.

If any of them fails, B1 does not merge.
- The failing rows are reported.
- The method is not tuned and rerun under this registration. A revised method needs an addendum, committed before its run.

Check 6 (the cost) is reported, not a gate: the cap already rules out more passes than today.

**Check 2's rule** becomes: on every call, the new method's passes are at most the bisection's plus 2.

## The pilot, the estimate and the stop

The largest pieces are piloted first, on this run's seeds:
- the gaussian calibration at 2,000 x 1,000, both ways;
- one call at S = 4,000.

The estimate and the stop, at twice the estimate, are registered from the pilot in a commit before the run.
The first run took 26.2 minutes.

**The pilot** (6 October 2026, frozen 426cf48, 1 core):
- The gaussian calibration at 2,000 x 1,000 took 225.4 s today and 187.8 s with B1, 1.20 times faster.
  The widths agree to 1.5e-13, and the verdict text is identical.
- The first run's pilot of the same piece was 2.07 times faster. The fixed guard bisects more often from
  wide brackets, as the probe found.
- One call at S = 4,000 (13 probabilities, both ways) took 6.3 s, with a median of 43 passes each way.
  All agree.

**The estimate: 35 minutes.**
- The first run took 26.2 minutes: 11.7 for the six timed check-6 arms, and the rest for checks 1 to 3
  and the fits.
- Today's calibration took 1.17 times as long in this pilot as in the first run's check 6, so both parts
  are scaled by that.
- B1's arms are scaled for the smaller gain: about 16 minutes for checks 4 and 6, and 17 for checks 1
  to 3.

**The stop: 70 minutes**, checked between families in checks 1 to 3. A stopped run records what it
reached, and the decision rule cannot pass on it.

## Addendum (6 October 2026, before any run): the budget grows as the bracket closes

**Why.** The registered pilot above measured the fixed guard at 1.20 times faster than today on the
gaussian calibration, against 2.07 in the first run's pilot. The guard is right but conservative: it
counts from the smallest stopping width the bracket can reach, while the budget stayed at the starting
width's count. So rows whose bracket starts wide, or near zero, bisected from the first pass.

Speed is what B1 is for. In the probe of 3,146 calls (a probe seed, not this study's), a budget that
grows as the bracket closes took 0.596 of the bisection's passes, against 0.778 for the registered
budget, and no call passed the cap.

**The change** (72ef45b):
- **The budget**, computed at every pass for each row:
  B(lo, hi) = ceiling(log2(w0 / (1e-12 max(1, |lo|, |hi|))) - 1e-3) + 2,
  where w0 is the starting bracket's width.
- **The guard**: a row takes a secant step only if used + 1 + N(lo, hi) <= B(lo, hi). Here N is the
  passes bisection needs from the current bracket down to the smallest stopping width reachable inside
  it:
  N(lo, hi) = ceiling(log2((hi - lo) / (1e-12 max(1, min |x| over [lo, hi]))) + 1e-3).
  It is unchanged from a14dca6.
- **The guard is asked again at every pass.** A row that bisects one pass may take a secant step on a
  later one.

**The proof that no row passes the bisection's count plus 2.** Write q for the quantile,
inf{x: F(x) >= p}.

1. *Every bracket holds q.*
   - The bracket check gives F(lo) < p <= F(hi).
   - Every pass keeps it: a point with F >= p becomes hi, and any other point becomes lo.
   - So q is in (lo, hi] for every bracket the root keeps.
   - The bisection's brackets, from the same start, keep the same property.
2. *The budget is a lower bound on the bisection's count, plus 2.*
   - Let the bisection stop after R passes, at [a, b]. It stops because b - a <= 1e-12 max(1, |a|, |b|).
   - Since q is in (a, b], max(|a|, |b|) <= |q| + (b - a).
   - For any bracket [lo, hi] of the root, |q| <= M = max(|lo|, |hi|).
   - So (b - a)(1 - 1e-12) <= 1e-12 max(1, M).
   - The bisection halves the width each pass: b - a = w0 / 2^R, up to the midpoints' rounding.
     That rounding is at most a relative 1e-4 even in the last passes, where the width is near
     1e-12 |x|.
   - So R >= log2(w0 / (1e-12 max(1, M))) - 1.5e-4, and as R is a whole number,
     R >= ceiling(log2(w0 / (1e-12 max(1, M))) - 1e-3).
   - So B(lo, hi) <= R + 2 for every bracket of the root.
3. *N counts enough.*
   - Bisection from [lo, hi] stops by the time the width reaches the stopping width of a bracket inside
     [lo, hi].
   - That stopping width is never below 1e-12 max(1, min |x| over [lo, hi]).
   - So it takes at most N(lo, hi) passes, the 1e-3 again covering rounding.
   - N cannot grow as the bracket shrinks inside itself.
4. *The cap.*
   - Suppose a row's guard last held at some pass, with u passes used and bracket [lo, hi]. Every later
     pass is a bisection from a bracket inside [lo, hi], so it finishes within N(lo, hi) passes.
   - The row's total is then at most u + 1 + N(lo, hi) <= B(lo, hi) <= R + 2.
   - If the guard never held, the row is the bisection itself, pass for pass: R.
   - A call's passes are the most any of its rows takes, so a call takes at most the bisection's passes
     plus 2.

Where F is not a number at a point, the property in step 1 cannot be read there. The bracket check then
stops the row with an NA quantile, as today.

**The test** (`tests/testthat/test-pred.R`, 72ef45b):
- the two mixtures that overran the first build;
- a bracket from -3 to 1,000 closing towards a quantile near 0.001, where the stopping width shrinks
  1,000-fold.

Each must take at most the bisection's passes plus 2 and give the bisection's answer. It passes, and with
the growing budget but the first build's guard it fails on two of the cases.

**Unchanged:** the checks, the seeds, the sizes, the library, the frozen-copy and detached run, and the
decision rule word for word. Check 6, the cost, stays reported and is not a gate.

**The pilot is run again** from a frozen copy at this addendum's commit, and the estimate and stop are
registered again from it before the run. The 35- and 70-minute figures above are superseded.

**The re-pilot** (6 October 2026, frozen 67e7832, 1 core):
- The gaussian calibration at 2,000 x 1,000 took 208.9 s today and 129.6 s with B1, 1.61 times faster.
  The widths agree to 2.1e-13, and the verdict text is identical.
- One call at S = 4,000 (13 probabilities, both ways) took 5.7 s, with medians of 43 passes today and
  14 for B1. All agree.

**The estimate: 30 minutes.** The first run took 26.2 minutes, and B1's gain here is close to its gain
then (1.61 against the first run's check-6 ratio of 1/0.63). Today's calibration took 1.08 times as long
as the first run's, which gives 28 minutes, rounded up.

**The stop: 60 minutes**, checked between families in checks 1 to 3, replacing the 70 above.
