"""17 -- Does the choice of transit-time estimator change any conclusion?

Calibration (section 2.3) minimises residuals on the mean-field integral

    T(L1 -> L2) = int dL / (s(L) * lambda(L)),

which follows the mean trajectory. The published anchors are mean
first-passage times of a jump process. For a nonlinear lambda these are
different quantities, and the integral additionally ignores contractions, so
the two estimators do not have to agree. `results/transit_time_check.csv`
shows they do not: the stochastic process is 15 to 41% slower than the
integral before contractions are counted and 22 to 53% slower after.

Which estimator matches the published anchors depends on how those anchors
were defined in the source, which is not recoverable from the published text.
This script therefore refits the same four parameters directly against
stochastic mean first-passage times and re-runs the headline predictions under
both parameterisations, so the reader can see the size of the dependence
instead of being asked to trust one choice.

Outputs
-------
results/estimator_calibration.csv
results/estimator_predictions.csv
"""
import csv
import json
import os

import numpy as np
from scipy.optimize import least_squares

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

TARGETS = [hd.LIT["t_40_60_years"], hd.LIT["t_60_80_years"],
           hd.LIT["t_80_150_years"], hd.LIT["rate_rise_70_90_min"]]
LEGS = ((40.0, 60.0), (60.0, 80.0), (80.0, 150.0))
N_FIT = 3000
CAGS = [42, 46, 50]
N_CELLS = 20000
DT = 1 / 24
RESIDUAL_ACHIEVABLE = 0.71
START_AGE = 30.0


def stochastic_predictions(vec, n=N_FIT, seed=99):
    p = hd._params_from_vector(vec)
    out = []
    for (a, b), sd in zip(LEGS, (seed, seed + 1, seed + 2)):
        m, _, frac = hd.mean_first_passage(a, b, p, n_cells=n, seed=sd)
        out.append(m if frac > 0.99 else 1e6)
    out.append(float(hd.lam(90.0, p) / hd.lam(70.0, p)))
    return np.array(out)


# ------------------------------------------------- the two parameterisations
par_mf = hd.calibrate()
x_mf = np.array([np.log(par_mf["lambda40"]), par_mf["beta"],
                 np.log(par_mf["R"]), par_mf["L_sw"]])

best = None
for jitter in (0.0, 0.25, -0.20):
    start = x_mf.copy()
    start[0] += jitter
    r = least_squares(lambda v: np.log(stochastic_predictions(v) / TARGETS),
                      start, bounds=([-12, 1e-4, 0.0, 50.0],
                                     [5.0, 1.0, np.log(1e5), 140.0]),
                      diff_step=0.02, xtol=1e-8, ftol=1e-8)
    if best is None or r.cost < best.cost:
        best = r
par_st = hd._params_from_vector(best.x)

cal_rows = []
for name, par, vec in (("mean-field integral", par_mf, x_mf),
                       ("stochastic first passage", par_st, best.x)):
    pred = stochastic_predictions(vec, n=8000, seed=501)
    cal_rows.append({
        "calibration_estimator": name,
        "lambda40_per_year": round(par["lambda40"], 5),
        "beta_per_repeat": round(par["beta"], 6),
        "R": round(par["R"], 2),
        "L_sw_repeats": round(par["L_sw"], 2),
        "stochastic_T_40_60_y": round(float(pred[0]), 1),
        "stochastic_T_60_80_y": round(float(pred[1]), 1),
        "stochastic_T_80_150_y": round(float(pred[2]), 2),
        "lambda90_over_lambda70": round(float(pred[3]), 2)})

print("calibration under the two estimators (published anchors: "
      f"{TARGETS[0]:.0f}, {TARGETS[1]:.0f}, {TARGETS[2]:.0f} y, "
      f"{TARGETS[3]:.0f}-fold)")
for r in cal_rows:
    print("  ", r)


# ------------------------------------------------------- headline predictions
def loss_at_cap(cag, times, surviving, cap=600.0):
    age = cap / (cag - hd.LIT["cap_offset"])
    return 1.0 - float(np.interp(age, times, surviving))


pred_rows = []
for name, par in (("mean-field integral", par_mf),
                  ("stochastic first passage", par_st)):
    caps, a25s, a50s = [], [], []
    for cag in CAGS:
        _, tt, ss = hd.simulate_population_tau(float(cag), par,
                                               n_cells=N_CELLS, dt=DT,
                                               seed=300 + cag)
        caps.append(100 * loss_at_cap(cag, tt, ss))
        a25s.append(float(np.interp(0.75, ss[::-1], tt[::-1])))
        a50s.append(float(np.interp(0.50, ss[::-1], tt[::-1])))

    # benefit of the measured achievable knockdown, 42 repeats, started at 30
    _, t0, s0 = hd.simulate_population_tau(42.0, par, n_cells=N_CELLS, dt=DT,
                                           seed=777)
    _, t1, s1 = hd.simulate_population_tau(42.0, par, n_cells=N_CELLS, dt=DT,
                                           t_treat=START_AGE,
                                           msh3_treat=RESIDUAL_ACHIEVABLE,
                                           seed=777)
    untreated = 100 * (1 - float(np.interp(70.0, t0, s0)))
    treated = 100 * (1 - float(np.interp(70.0, t1, s1)))

    # loss at the age the external onset regression predicts
    onset_loss = []
    for cag in CAGS:
        _, tt, ss = hd.simulate_population_tau(float(cag), par,
                                               n_cells=N_CELLS, dt=DT,
                                               seed=300 + cag)
        age = hd.ranen_predicted_onset(cag)
        onset_loss.append(100 * (1 - float(np.interp(age, tt, ss))))

    pred_rows.append({
        "calibration_estimator": name,
        "mean_loss_at_CAP600_pct": round(float(np.mean(caps)), 1),
        "age_25pct_loss_CAG42": round(a25s[0], 1),
        "age_50pct_loss_CAG42": round(a50s[0], 1),
        "age_50pct_loss_CAG50": round(a50s[2], 1),
        "loss_at_predicted_onset_CAG42_pct": round(onset_loss[0], 1),
        "loss_at_predicted_onset_CAG50_pct": round(onset_loss[2], 1),
        "untreated_loss_at_70_pct": round(untreated, 1),
        "treated_loss_at_70_pct": round(treated, 1),
        "benefit_points": round(untreated - treated, 1)})

print("\nheadline predictions under the two parameterisations")
for r in pred_rows:
    print("  ", r)


def _write(name, rws):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rws[0].keys()))
        w.writeheader()
        w.writerows(rws)


_write("estimator_calibration.csv", cal_rows)
_write("estimator_predictions.csv", pred_rows)
with open(os.path.join(RES, "estimator_fit.json"), "w") as f:
    json.dump({"mean_field": {k: (None if v == np.inf else v)
                              for k, v in par_mf.items()},
               "stochastic": {k: (None if v == np.inf else v)
                              for k, v in par_st.items()}}, f, indent=2)
print("\nwrote results/estimator_calibration.csv and "
      "results/estimator_predictions.csv")
