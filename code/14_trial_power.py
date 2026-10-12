"""14 -- Could a trial of this size and length actually see the effect?

A predicted benefit is only useful if it is detectable. This script turns the
lifetime model into a conventional sample-size calculation for a randomised,
placebo-controlled trial whose endpoint is the change in striatal projection
neuron survival over the trial period, the quantity that striatal volume on MRI
is standing in for.

Per-subject endpoint: loss(age0 + T) - loss(age0), in percentage points of the
SPN population. Between-subject variation comes from two sources that are kept
separate: the spread of inherited CAG lengths among enrolled participants, and
measurement error on the imaging endpoint. The latter is not a property of this
model, so it is swept.

n per arm = 2 (z_(1-alpha/2) + z_(1-beta))^2 * sigma^2 / delta^2

Outputs
-------
results/trial_power.csv
figures/fig15_power.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

N_CELLS = 12000
DT = 1 / 24
AGES = np.linspace(0, 110, 2201)
CAGS = list(range(40, 51))
# Enrolment distribution: triangular, centred on 44, which is roughly where
# adult-onset trial cohorts sit. Assumed, and swept in the reporting below.
WEIGHTS = np.array([1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1], dtype=float)
WEIGHTS /= WEIGHTS.sum()

ALPHA, POWER = 0.05, 0.80
Z = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)

MEASUREMENT_SD = [0.0, 1.0, 2.0, 4.0]   # percentage points, per measurement
START_AGES = [30, 40, 50]
DURATIONS = [2, 3, 5]
ARMS = {"29% knockdown (measured achievable)": 0.71,
        "55% knockdown (best caudate figure)": 0.45,
        "75% knockdown (aspirational)": 0.25}

par = hd.calibrate()
cache = {}


def loss_curve(cag, residual, t_treat):
    key = (cag, residual, t_treat)
    if key not in cache:
        _, tt, ss = hd.simulate_population_tau(
            float(cag), par, n_cells=N_CELLS, dt=DT,
            t_treat=float(t_treat), msh3_treat=float(residual),
            seed=4000 + cag + int(1000 * residual)
            + (0 if not np.isfinite(t_treat) else int(t_treat)))
        cache[key] = 100 * (1 - np.interp(AGES, tt, ss))
    return cache[key]


rows = []
for start in START_AGES:
    placebo = {c: loss_curve(c, 1.0, np.inf) for c in CAGS}
    for arm, residual in ARMS.items():
        treated = {c: loss_curve(c, residual, start) for c in CAGS}
        for T in DURATIONS:
            d_pl = np.array([np.interp(start + T, AGES, placebo[c])
                             - np.interp(start, AGES, placebo[c]) for c in CAGS])
            d_tr = np.array([np.interp(start + T, AGES, treated[c])
                             - np.interp(start, AGES, treated[c]) for c in CAGS])
            mean_pl = float(WEIGHTS @ d_pl)
            mean_tr = float(WEIGHTS @ d_tr)
            delta = mean_pl - mean_tr
            # variance across the enrolled CAG distribution, pooled over arms
            var_cag = float(WEIGHTS @ (d_pl - mean_pl) ** 2
                            + WEIGHTS @ (d_tr - mean_tr) ** 2) / 2.0
            for m_sd in MEASUREMENT_SD:
                # two measurements per subject (baseline and end of trial)
                sigma = np.sqrt(var_cag + 2 * m_sd ** 2)
                n = (2 * Z ** 2 * sigma ** 2 / delta ** 2) if delta > 1e-9 else np.inf
                rows.append({
                    "treatment_start_age": start,
                    "arm": arm,
                    "trial_duration_years": T,
                    "placebo_loss_over_trial_points": round(mean_pl, 2),
                    "treated_loss_over_trial_points": round(mean_tr, 2),
                    "effect_delta_points": round(delta, 2),
                    "between_subject_sd_from_CAG": round(np.sqrt(var_cag), 2),
                    "measurement_sd_points": m_sd,
                    "total_sd_points": round(float(sigma), 2),
                    "n_per_arm_80pct_power": (None if not np.isfinite(n)
                                              else int(np.ceil(n)))})

# ---------------------------------------------------------------------------
# The variance above counts only one source of between-subject difference: the
# spread of inherited CAG lengths across the cohort. Real participants also
# differ in the modifier genotypes of section 3.7, which move age at 25% loss
# by more years than the CAG spread does, and that variation is absent from
# the simulation. It enters the sample-size formula exactly as an extra
# between-subject standard deviation, so its effect can be computed without
# re-simulating anything.
EXTRA_BETWEEN_SD = [0.0, 0.5, 1.0, 2.0]   # percentage points over the trial
mod_rows = []
for r in rows:
    if (r["treatment_start_age"] != 40 or r["trial_duration_years"] != 3
            or r["measurement_sd_points"] != 2.0):
        continue
    for extra in EXTRA_BETWEEN_SD:
        sigma = float(np.sqrt(r["total_sd_points"] ** 2 + extra ** 2))
        delta = r["effect_delta_points"]
        n = 2 * Z ** 2 * sigma ** 2 / delta ** 2
        mod_rows.append({
            "arm": r["arm"],
            "trial_duration_years": r["trial_duration_years"],
            "treatment_start_age": r["treatment_start_age"],
            "measurement_sd_points": r["measurement_sd_points"],
            "extra_between_subject_sd_points": extra,
            "total_sd_points": round(sigma, 2),
            "n_per_arm_80pct_power": int(np.ceil(n))})

print("\neffect of between-subject variation the model does not simulate "
      "(modifier genotype), 3-year trial enrolling at 40, imaging SD 2")
for r in mod_rows:
    print(f"  {r['arm']:<40} extra SD {r['extra_between_subject_sd_points']:>4} "
          f"-> n = {r['n_per_arm_80pct_power']}")

print("sample size per arm, 80% power, alpha 0.05")
for r in rows:
    if (r["arm"].startswith("29%") and r["measurement_sd_points"] == 2.0):
        print(f"  start {r['treatment_start_age']}, {r['trial_duration_years']} y: "
              f"effect {r['effect_delta_points']:>5} pts, SD {r['total_sd_points']:>5}"
              f"  -> n = {r['n_per_arm_80pct_power']} per arm")

for name, rws in (("trial_power.csv", rows),
                  ("modifier_variance_power.csv", mod_rows)):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rws[0].keys()))
        w.writeheader()
        w.writerows(rws)

fig, ax = plt.subplots(1, 3, figsize=(10.4, 2.95))

for arm, c in zip(ARMS, ("#c0504d", "#8064a2", "#1f4e79")):
    sub = [r for r in rows if r["arm"] == arm and r["treatment_start_age"] == 40
           and r["measurement_sd_points"] == 2.0]
    ax[0].plot([r["trial_duration_years"] for r in sub],
               [r["effect_delta_points"] for r in sub], "o-", color=c, lw=1.6,
               ms=4, label=arm.split(" (")[0])
ax[0].set_xlabel("trial duration (years)")
ax[0].set_ylabel("effect size (percentage points\nof striatum preserved)")
ax[0].set_title("A  effect size, treated from 40")
ax[0].legend(frameon=False, fontsize=6.3)

for arm, c in zip(ARMS, ("#c0504d", "#8064a2", "#1f4e79")):
    sub = [r for r in rows if r["arm"] == arm and r["treatment_start_age"] == 40
           and r["measurement_sd_points"] == 2.0
           and r["n_per_arm_80pct_power"] is not None]
    ax[1].semilogy([r["trial_duration_years"] for r in sub],
                   [r["n_per_arm_80pct_power"] for r in sub], "o-", color=c,
                   lw=1.6, ms=4, label=arm.split(" (")[0])
ax[1].axhline(200, color="0.5", lw=0.8, ls=":")
ax[1].annotate("200 per arm", xy=(2.05, 215), fontsize=6.3)
ax[1].set_xlabel("trial duration (years)")
ax[1].set_ylabel("n per arm for 80% power")
ax[1].set_title("B  sample size")
ax[1].legend(frameon=False, fontsize=6.3)

sub3 = [r for r in rows if r["arm"].startswith("29%")
        and r["trial_duration_years"] == 3
        and r["n_per_arm_80pct_power"] is not None]
for start, c in zip(START_AGES, ("#1f4e79", "#8064a2", "#c0504d")):
    s = [r for r in sub3 if r["treatment_start_age"] == start]
    ax[2].semilogy([r["measurement_sd_points"] for r in s],
                   [r["n_per_arm_80pct_power"] for r in s], "o-", color=c,
                   lw=1.6, ms=4, label=f"enrol at {start}")
ax[2].set_xlabel("imaging measurement SD (points)")
ax[2].set_ylabel("n per arm for 80% power")
ax[2].set_title("C  3-year trial, achievable dose")
ax[2].legend(frameon=False, fontsize=6.3)
fig.savefig(os.path.join(FIG, "fig15_power.png"))
plt.close(fig)
print("\nwrote figure 15 and result table")
