"""11 -- Direct and indirect pathway neurons: a held-out test of the model.

Indirect-pathway SPNs are reported to expand about 15% faster than
direct-pathway SPNs, and to make up a smaller share of surviving SPNs in
Huntington's disease donors than in controls (47 +/- 6% in controls). Neither
number was used in calibration, so the subtype composition over time is a
prediction.

Two questions:
  1. Does a 15% rate difference, on its own, produce a measurable shift in the
     surviving iSPN fraction, and how large is it?
  2. What rate difference would be needed to produce a given depletion?

Outputs
-------
results/subtype_trajectory.csv
results/subtype_rate_scan.csv
figures/fig12_subtypes.png
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

ISPN_BASELINE = 0.47     # reported control fraction
ISPN_RATE_RATIO = 1.15   # reported iSPN / dSPN somatic expansion rate
N = 15000
DT = 1 / 24
AGES = np.linspace(0, 100, 1001)

par = hd.calibrate()


def survival(cag, rate_mult, seed):
    p = dict(par)
    p["A"] = par["A"] * rate_mult
    _, tt, ss = hd.simulate_population_tau(float(cag), p, n_cells=N, dt=DT,
                                           seed=seed)
    return np.interp(AGES, tt, ss)


def ispn_fraction(cag, ratio, seed=0):
    s_d = survival(cag, 1.0, seed)
    s_i = survival(cag, ratio, seed + 1)
    tot = ISPN_BASELINE * s_i + (1 - ISPN_BASELINE) * s_d
    frac = np.where(tot > 1e-9, ISPN_BASELINE * s_i / np.maximum(tot, 1e-9),
                    np.nan)
    overall = ISPN_BASELINE * s_i + (1 - ISPN_BASELINE) * s_d
    return frac, overall, s_d, s_i


traj_rows = []
curves = {}
for cag in (42, 46, 50):
    frac, overall, s_d, s_i = ispn_fraction(cag, ISPN_RATE_RATIO, seed=600 + cag)
    curves[cag] = (frac, overall)
    for target_loss in (0.10, 0.25, 0.50, 0.75):
        loss = 1 - overall
        if loss.max() < target_loss:
            continue
        k = int(np.argmax(loss >= target_loss))
        traj_rows.append({"inherited_CAG": cag,
                          "overall_SPN_loss_pct": int(100 * target_loss),
                          "age_years": round(float(AGES[k]), 1),
                          "iSPN_fraction_of_survivors_pct": round(
                              100 * float(frac[k]), 1),
                          "control_fraction_pct": round(100 * ISPN_BASELINE, 1),
                          "depletion_points": round(
                              100 * (ISPN_BASELINE - float(frac[k])), 2)})

print("subtype trajectory (15% faster iSPN expansion)")
for r in traj_rows:
    print("  ", r)

# How large a rate difference would a given depletion require?
scan_rows = []
for ratio in (1.00, 1.05, 1.10, 1.15, 1.25, 1.50, 2.00):
    frac, overall, _, _ = ispn_fraction(42, ratio, seed=900 + int(100 * ratio))
    loss = 1 - overall
    k = int(np.argmax(loss >= 0.50)) if loss.max() >= 0.50 else len(AGES) - 1
    scan_rows.append({"iSPN_dSPN_rate_ratio": ratio,
                      "age_at_50pct_loss": round(float(AGES[k]), 1),
                      "iSPN_fraction_at_50pct_loss_pct": round(
                          100 * float(frac[k]), 1),
                      "depletion_points": round(
                          100 * (ISPN_BASELINE - float(frac[k])), 2)})

print("\nrate-ratio scan (42 CAG, measured at 50% overall SPN loss)")
for r in scan_rows:
    print("  ", r)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("subtype_trajectory.csv", traj_rows, list(traj_rows[0].keys()))
_write("subtype_rate_scan.csv", scan_rows, list(scan_rows[0].keys()))

fig, ax = plt.subplots(1, 3, figsize=(10.2, 2.9))

for cag, c in zip((42, 46, 50), ("#1f4e79", "#8064a2", "#c0504d")):
    frac, overall = curves[cag]
    ax[0].plot(AGES, 100 * frac, color=c, lw=1.7, label=f"{cag} CAG")
ax[0].axhline(100 * ISPN_BASELINE, color="k", lw=1.0, ls="--")
ax[0].axhspan(41, 53, color="0.9", lw=0, zorder=0)
ax[0].annotate("control 47 +/- 6%", xy=(3, 100 * ISPN_BASELINE + 0.6),
               fontsize=6.5)
ax[0].set_xlabel("age (years)")
ax[0].set_ylabel("iSPNs as % of surviving SPNs")
ax[0].set_ylim(30, 50)
ax[0].set_title("A  predicted subtype composition")
ax[0].legend(frameon=False, fontsize=6.5, loc="lower left")

for cag, c in zip((42, 46, 50), ("#1f4e79", "#8064a2", "#c0504d")):
    frac, overall = curves[cag]
    ax[1].plot(100 * (1 - overall), 100 * frac, color=c, lw=1.7,
               label=f"{cag} CAG")
ax[1].axhline(100 * ISPN_BASELINE, color="k", lw=1.0, ls="--")
ax[1].set_xlabel("overall SPN loss (%)")
ax[1].set_ylabel("iSPNs as % of survivors")
ax[1].set_ylim(30, 50)
ax[1].set_title("B  depletion tracks total loss")
ax[1].legend(frameon=False, fontsize=6.5, loc="lower left")

ratios = [r["iSPN_dSPN_rate_ratio"] for r in scan_rows]
dep = [r["depletion_points"] for r in scan_rows]
ax[2].plot(ratios, dep, "o-", color="#1f4e79", lw=1.6, ms=4)
ax[2].axvline(ISPN_RATE_RATIO, color="0.4", ls=":", lw=1.0)
ax[2].annotate("reported 1.15", xy=(ISPN_RATE_RATIO + 0.02, 0.5), fontsize=6.5)
ax[2].set_xlabel("iSPN / dSPN expansion-rate ratio")
ax[2].set_ylabel("iSPN depletion at 50% loss\n(percentage points)")
ax[2].set_title("C  how big a difference is needed")
fig.savefig(os.path.join(FIG, "fig12_subtypes.png"))
plt.close(fig)
print("\nwrote figure 12 and result tables")
