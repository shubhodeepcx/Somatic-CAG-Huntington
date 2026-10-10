"""04 -- Reconciling two independent published constraints.

The fitted rate law reproduces the published per-cell dwell times exactly, but
it under-predicts end-stage striatal loss: donors with CAP > 600 are reported to
have lost 80-99% of their SPNs, whereas the model loses ~43% at CAP = 600.
Script 03 showed that accelerating post-threshold elimination (even 40-fold)
cannot close the gap, because reaching 150 CAG, not dying afterwards, is
rate-limiting.

Here the two remaining structural possibilities are tested quantitatively:
  (i)  the expansion rate is uniformly faster than the published dwell times
       imply, by a factor g;
  (ii) the length at which an SPN becomes committed to elimination is lower
       than 150 CAG.

A 2-D feasibility map is computed, and two bracketing parameterisations are
carried forward to the therapeutic simulations:
  M1  "dwell-time faithful"  g = 1,      L_commit = 150
  M2  "loss faithful"        g = fitted, L_commit = 150

Outputs
-------
results/reconciliation_map.csv
results/reconciled_models.csv
figures/fig5_reconciliation.png
"""
import csv
import json
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

par0 = hd.calibrate()
CAGS = [40, 42, 44, 46, 48, 50]
DT = 1 / 24
N_SCAN = 6000
N_FINAL = 20000


def scaled(par, g):
    p = dict(par)
    p["A"] = par["A"] * g
    return p


def loss_at_cap(cag, times, surviving, cap=600.0):
    age = cap / (cag - hd.LIT["cap_offset"])
    return 1.0 - float(np.interp(age, times, surviving))


def mean_loss_at_cap600(g, L_commit, n=N_SCAN, seed0=1000):
    saved = hd.LIT["L_phaseC"]
    hd.LIT["L_phaseC"] = L_commit
    try:
        vals = []
        for cag in CAGS:
            _, tt, ss = hd.simulate_population_tau(
                float(cag), scaled(par0, g), n_cells=n, dt=DT,
                hazard_kw={"L_death_50": hd.ASSUMED["L_death_50"] *
                           (L_commit / 150.0)},
                seed=seed0 + cag)
            vals.append(loss_at_cap(cag, tt, ss))
        return float(np.mean(vals))
    finally:
        hd.LIT["L_phaseC"] = saved


# ------------------------------------------------------- feasibility map
g_grid = [1.0, 1.25, 1.5, 2.0, 3.0]
Lc_grid = [150.0, 130.0, 110.0, 90.0]
rows = []
grid_vals = np.zeros((len(Lc_grid), len(g_grid)))
for i, Lc in enumerate(Lc_grid):
    for j, g in enumerate(g_grid):
        v = mean_loss_at_cap600(g, Lc)
        grid_vals[i, j] = 100 * v
        rows.append({"rate_multiplier_g": g, "L_commit": Lc,
                     "mean_SPN_loss_at_CAP600_pct": round(100 * v, 1),
                     "implied_t_40_60_years": round(hd.LIT["t_40_60_years"] / g, 1)})
        print(f"  g={g:<5} L_commit={Lc:<6} loss@CAP600={100*v:5.1f}%")

# ------------------------------------------- fit g at L_commit = 150 (model M2)
target = 0.90
lo, hi = 1.0, 6.0
for _ in range(14):
    mid = 0.5 * (lo + hi)
    if mean_loss_at_cap600(mid, 150.0) < target:
        lo = mid
    else:
        hi = mid
g_fit = 0.5 * (lo + hi)
print(f"\nfitted rate multiplier g = {g_fit:.2f} "
      f"(implied 40->60 CAG transit {hd.LIT['t_40_60_years']/g_fit:.0f} y)")

models = {"M1_dwell_faithful": 1.0, "M2_loss_faithful": g_fit}
mrows = []
curves = {}
for name, g in models.items():
    p = scaled(par0, g)
    for cag in (42, 46):
        arr, tt, ss = hd.simulate_population_tau(float(cag), p, n_cells=N_FINAL,
                                                 dt=DT, seed=7000 + cag)
        curves[(name, cag)] = (tt, ss)
    vals, a50s = [], []
    for cag in CAGS:
        arr, tt, ss = hd.simulate_population_tau(float(cag), p, n_cells=N_FINAL,
                                                 dt=DT, seed=8000 + cag)
        vals.append(100 * loss_at_cap(cag, tt, ss))
        a50s.append(float(np.interp(0.5, ss[::-1], tt[::-1])))
        if cag == 42:
            died = arr["died"]
            life = np.where(died, arr["t_death"], 110.0)
            flb = 100 * float(np.mean(arr["years_below_150"] / life))
    mrows.append({"model": name, "rate_multiplier_g": round(g, 3),
                  "implied_t_40_60_y": round(hd.LIT["t_40_60_years"] / g, 1),
                  "implied_t_60_80_y": round(hd.LIT["t_60_80_years"] / g, 1),
                  "implied_t_80_150_y": round(hd.LIT["t_80_150_years"] / g, 2),
                  "mean_loss_at_CAP600_pct": round(float(np.mean(vals)), 1),
                  "age_50pct_loss_CAG42": round(a50s[CAGS.index(42)], 1),
                  "age_50pct_loss_CAG46": round(a50s[CAGS.index(46)], 1),
                  "frac_life_below_150_CAG42_pct": round(flb, 1),
                  "ranen_onset_CAG42": round(hd.ranen_predicted_onset(42), 1)})

print("\nbracketing models")
for r in mrows:
    print("  ", r)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("reconciliation_map.csv", rows, list(rows[0].keys()))
_write("reconciled_models.csv", mrows, list(mrows[0].keys()))
with open(os.path.join(RES, "g_fit.json"), "w") as f:
    json.dump({"g_fit": g_fit, "target_loss_at_CAP600": target}, f, indent=2)

# -------------------------------------------------------------------- Fig. 5
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.1),
                       gridspec_kw={"width_ratios": [1.15, 1]})

im = ax[0].imshow(grid_vals, origin="lower", aspect="auto", cmap="magma",
                  vmin=20, vmax=100,
                  extent=[-0.5, len(g_grid) - 0.5, -0.5, len(Lc_grid) - 0.5])
ax[0].set_xticks(range(len(g_grid)))
ax[0].set_xticklabels([f"{g:g}" for g in g_grid])
ax[0].set_yticks(range(len(Lc_grid)))
ax[0].set_yticklabels([f"{int(l)}" for l in Lc_grid])
ax[0].set_xlabel("expansion-rate multiplier $g$")
ax[0].set_ylabel("commitment length $L_{commit}$ (CAG)")
for i in range(len(Lc_grid)):
    for j in range(len(g_grid)):
        v = grid_vals[i, j]
        ax[0].text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                   color="white" if v < 70 else "black")
cb = fig.colorbar(im, ax=ax[0], fraction=0.046)
cb.set_label("SPN loss at CAP 600 (%)", fontsize=7)
cb.ax.tick_params(labelsize=7)
ax[0].set_title("A  what would reconcile the two constraints", fontsize=8.5)

for (name, cag), (tt, ss) in curves.items():
    style = "-" if name.startswith("M2") else "--"
    col = "#1f4e79" if cag == 42 else "#c0504d"
    ax[1].plot(tt, 100 * ss, style, color=col, lw=1.6,
               label=f"{cag} CAG, {'M2' if name.startswith('M2') else 'M1'}")
for cag, col in ((42, "#1f4e79"), (46, "#c0504d")):
    age600 = 600.0 / (cag - hd.LIT["cap_offset"])
    ax[1].axvspan(age600 - 0.4, age600 + 0.4, color=col, alpha=0.25, lw=0)
ax[1].axhspan(1, 20, color="0.9", lw=0, zorder=0)
ax[1].set_xlabel("age (years)")
ax[1].set_ylabel("surviving SPNs (%)")
ax[1].set_xlim(0, 100)
ax[1].set_title("B  bracketing models\n(bands: age at CAP 600)", fontsize=8.5)
ax[1].legend(frameon=False, fontsize=6.5, loc="lower left")
fig.savefig(os.path.join(FIG, "fig5_reconciliation.png"))
plt.close(fig)
print("\nwrote figure 5 and result tables")
