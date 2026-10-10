"""12 -- Genetic modifiers as perturbations of the same rate law.

The strongest modifiers of age at onset in Huntington's disease sit in DNA
mismatch-repair genes and in the CAG tract's own sequence. Each maps onto a
parameter of the fitted rate law:

  MSH3 expression        -> the first-order factor theta
  FAN1 activity          -> the protective divisor f
  CAA interruption loss  -> a downward shift in the clustering switch L_sw,
                            because an uninterrupted tract reaches the
                            clustered regime at a shorter total length

Running each perturbation through the population model turns a genotype into a
predicted shift in age at 25% striatal neuron loss, in years. These are
predictions, not fits: no modifier effect size was used in calibration.

Outputs
-------
results/modifier_effects.csv
figures/fig13_modifiers.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

N = 12000
DT = 1 / 24
AGES = np.linspace(0, 110, 1101)
CAG = 42.0

par = hd.calibrate()


def age_at_loss(frac, msh3=1.0, fan1=1.0, dL_sw=0.0, seed=0):
    p = dict(par)
    p["L_sw"] = par["L_sw"] + dL_sw
    _, tt, ss = hd.simulate_population_tau(CAG, p, n_cells=N, dt=DT,
                                           msh3=msh3, fan1=fan1, seed=seed)
    if ss.min() > 1 - frac:
        return np.nan
    return float(np.interp(1 - frac, ss[::-1], tt[::-1]))


base = age_at_loss(0.25, seed=1)
print(f"baseline age at 25% SPN loss (42 CAG): {base:.1f} y")

rows = []
scans = {
    "MSH3 expression (theta)": [("msh3", v) for v in
                                (0.6, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5)],
    "FAN1 activity (f)": [("fan1", v) for v in
                          (0.5, 0.75, 0.9, 1.0, 1.25, 1.5, 2.0)],
    "clustering switch shift (repeats)": [("dL_sw", v) for v in
                                          (-15, -10, -5, 0, 5, 10, 15)],
}

curves = {}
for label, entries in scans.items():
    xs, ys = [], []
    for i, (kind, v) in enumerate(entries):
        kw = {"msh3": 1.0, "fan1": 1.0, "dL_sw": 0.0}
        kw[kind] = v
        a = age_at_loss(0.25, seed=100 * (hash(label) % 97) + i, **kw)
        shift = (a - base) if np.isfinite(a) else np.nan
        xs.append(v)
        ys.append(shift)
        rows.append({"parameter": label, "value": v,
                     "age_at_25pct_loss": None if not np.isfinite(a) else round(a, 1),
                     "shift_vs_baseline_years": None if not np.isfinite(shift)
                     else round(shift, 1)})
        print(f"  {label:<34} {v:<6} -> {shift:+.1f} y")
    curves[label] = (xs, ys)

with open(os.path.join(RES, "modifier_effects.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

fig, ax = plt.subplots(1, 3, figsize=(10.2, 2.9))
for a, (label, (xs, ys)) in zip(ax, curves.items()):
    a.plot(xs, ys, "o-", color="#1f4e79", lw=1.7, ms=4)
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel(label, fontsize=8)
    a.set_ylabel("shift in age at 25% SPN loss (y)", fontsize=8)
ax[0].set_title("A  MSH3 dose", fontsize=9)
ax[1].set_title("B  FAN1 protection", fontsize=9)
ax[2].set_title("C  loss of CAA interruption", fontsize=9)
fig.savefig(os.path.join(FIG, "fig13_modifiers.png"))
plt.close(fig)
print("\nwrote figure 13 and result table")
