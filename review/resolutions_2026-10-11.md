# What was done about the audit

Companion to `audit_2026-10-11.md`. One entry per finding, saying what changed
and where. Three findings turned into new analysis rather than new wording, and
one of them changed a conclusion, so those are marked.

---

## A. Substantive

**A1. Exposure-response curve has zero degrees of freedom and is extrapolated.**
Methods 2.9 now states the zero degrees of freedom outright, says that the Hill
coefficient of 1.22 is a consequence of the fit rather than a measurement of
cooperativity, and lists the exposures at which the curve is used outside its
7.5 to 25 nmol g⁻¹ calibration span. Figure 6A shades the calibrated range so
every later exposure can be read against it, and its title now says "2 points,
2 parameters". Added as Limitations item 5.

**A2. The primate figure is a different measurement than the one it is used as.**
The source was retrieved and the target identified: the oligonucleotide was
directed against *Malat1*, a nuclear long non-coding RNA used as a reporter of
intrathecal distribution, and the readout is residual RNA. Methods 2.9 and
Results 3.8 both name the transcript, state that MSH3 equivalence is assumed,
and state that residual RNA is standing in for residual protein. The Table V
column heading changed from "Knockdown achieved" to "Target knockdown measured
in primate", the suppression column now reads "if MSH3 behaves the same way",
and the Figure 6B legend changed to match. Added as Limitations item 6.

**A3. The first-order MSH3 result is close to circular as tested.**
`01_calibrate.py` now computes a non-circular statistic, the worst error against
the two measured points, and selects the best order on that rather than on the
RMSE against the regression line. Results 3.2 leads with the two measurements
(worst error 4.1 points for k = 1 against 20.1 and 34.4 for k = 2 and 3),
reports the RMSE table second and explains in plain terms why its margin is
partly an artefact of its own construction. The conclusion was narrowed from
"first-order" to "not a cooperative assembly", which is what two points support.
Methods 2.6 rewritten to match.

**A4. The sharpness of the switch is assumed. (new analysis)**
`10_profile.py` now holds w_sw at 19 values from 3 to 12 repeats and refits all
four parameters at each, writing `results/wsw_profile.csv` and a fifth panel on
Figure 11. The anchors turn out to carry more information than expected: an
exact fit exists only for widths between about 4 and 7.5 repeats, and the fit
misses by 5.2% at 3 repeats and 13.5% at 12. But the switch position moves with
the assumed width, from L_sw = 74.6 to 101.6 repeats, so Results 3.1 now reports
the two jointly and says the qualitative claim is robust while the specific
figure of 80.2 is conditional on a width of 6. Added as Limitations item 7.

**A5. Three assumed parameters are among the strongest drivers.**
Results 3.13 now names them, contraction rate at −0.68, contraction step at
−0.45 and switch width at −0.42, says they sit above every elimination constant
and that they are the main reason the interval on the benefit is as wide as it
is. The reassuring sentence was narrowed to the elimination submodel, where it
is true.

**A6. "Every assumed value is varied" was not true.**
Methods 2.1 now distinguishes the constants swept globally and simultaneously
from those swept one at a time, and says why. `L_max` was added to the global
Latin hypercube, swept 600 to 1,400 repeats, which answers C6 at the same time.
Table III was split so each group of assumed constants states its own treatment
rather than all carrying "all swept".

**A7. A falsifiable prediction rests on an undisclosed assumed number.**
The assumed parental age at onset of 45 years is now stated in Results 3.5 where
the regression is first used, added to Table III, and swept from 35 to 55 years
in `06_trial.py`, writing `results/parent_aao_sensitivity.csv`. The sweep shows
the rejection of the faster model survives the whole range, while the 18 to 31%
onset prediction widens, and section 4.3 now carries that caveat.

---

## B. Missing explanations

**B1. The mouse and primate arms treat the same physics differently.**
Methods 2.9 names the asymmetry, explains why the mouse calculation is kept in
whole-brain form, and notes that the resulting error in EC₅₀ is conservative for
the paper's conclusion.

**B2. The solver check is noise-limited. (new analysis)**
`03_coupling.py` now compares two exact runs differing only in seed, giving a
Monte Carlo noise floor of 0.013 to 0.024, which is the same size as the
tau-leaping differences it was being used to bound. Methods 2.5 states the
resolution limit: no bias larger than about two to three percentage points of
survival, and nothing smaller is detectable by this test.

**B3. The shape of the step-size distribution is unexplained.**
Methods 2.2 now justifies the geometric choice as the memoryless law on the
positive integers, the weakest assumption available given one measured mean, and
points out that it is deliberately the thin-tailed option, so the narrow
terminal tail in Results 3.3 follows from a conservative choice.

**B4. Population size varies between scripts.**
Methods 2.11 now accounts for every size used and the reason for each.

**B5. Calibration fits a mean-field integral to stochastic means. (new analysis,
changed a conclusion)**
This was the most consequential finding. `hdmodel.mean_first_passage` was added
and `03_coupling.py` now compares the two estimators directly on the calibration
legs: the stochastic process is 21 to 50% slower than the integral, about two
thirds of it from the jump process and the rest from contraction. The claim in
Methods 2.3 that the integral was "checked against the stochastic simulation"
was false and has been replaced with the measured discrepancy.

A new script, `17_estimator.py`, refits all four parameters against stochastic
mean first-passage times. The refit converges, raises λ(40) from 0.070 to 0.092
per year and R from 17.7 to 24.3, and leaves the switch position at 80.0 repeats
and β unidentifiable, so the structural conclusions of section 3.1 are not
artefacts of the estimator. Mean loss at CAP 600 rises from 39.8% to 62.4%,
which closes more than half of the gap that sections 3.4 and 3.5 could not
explain, and the therapeutic benefit moves only from 13.7 to 14.7 points.
Results 3.5 now carries this as a third candidate explanation, and it is the
only one of the three that partly works. The mean-field fit is kept as primary
because it reproduces the anchors exactly as published and the source does not
say which estimator its numbers represent. Added as Limitations item 8.

**B6. The subtype rate difference is a rate-constant multiplier.**
Methods 2.7 now states the substitution, names the alternative (a shift in L_sw
rather than in A), says the two cannot be distinguished from available data, and
says the multiplier is the more conservative of them.

---

## C. Minor

- **C1.** `hdmodel.py` docstring now points at `03_coupling.py`.
- **C2.** Methods 2.1 now says forty minutes, matching section 6 and the README.
- **C3.** `ASSUMED["L_switch"]` and `ASSUMED["R_cluster"]` were removed. Both
  `clustered_fraction` and `lambda_L` now raise if the fitted value is not
  supplied, so the silent substitution is gone.
- **C4.** Results 3.3 no longer counts the first row as a pass. The heading says
  three of five, the prose says the 99.9% figure is outside the published 95 to
  98%, and the missing protected fraction is Limitations item 9.
- **C5.** The branch is not unreachable, which the audit got wrong: it is reached
  at complete MSH3 ablation below the elimination threshold, where λ and the
  hazard are both zero. No current caller does that, so the branch is kept and
  now carries a comment saying when it fires.
- **C6.** Answered by the `L_max` sweep described under A6.

---

## What is not fixed, and why

Three findings describe limits of the published data rather than defects that
could be engineered away, and they are fixed by being stated accurately rather
than by disappearing:

- the exposure-response curve still has two points behind it, because two is how
  many doses were published;
- the MSH3 order still rests on two measurements from a regression with
  R² = 0.65;
- the parental age at onset is still assumed, because the source reports no
  cohort mean.

In each case the paper now says what the evidence is, what it can carry and
what it cannot, and the sensitivity of the conclusion to it is computed rather
than asserted.

---

## One correction to the audit itself

A1 says the exposure-response curve is used at steady-state peaks up to
157 nmol g⁻¹ in the dosing analysis. That is right for the schedule the paper
headlines, the 75% peak-knockdown regimen. Across the whole grid in
`results/dosing_timeaverage.csv`, which includes a 90% target, the maximum is
388 nmol g⁻¹, fifteen times the highest calibration point rather than six. Both
figures are correct about different subsets of the grid, so A1 was not an
error; the paper now states the larger one, because the extrapolation
limitation should be stated against the worst case in the analysis rather than
against the headline case.

C5 was wrong, as noted above: the `tot <= 0` branch is reachable at complete
MSH3 ablation, so it was documented rather than deleted.
