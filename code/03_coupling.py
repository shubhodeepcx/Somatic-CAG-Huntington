"""03 -- Does a purely cell-autonomous DNA-repeat clock explain the observed
steepness of striatal degeneration?

Step 1: verify the fast tau-leaping population solver against the exact
        Gillespie solver (k_nca = 0).
Step 2: compare the cell-autonomous model with the reported relationship
        between CAP score and SPN loss (CAP > 600 -> 80-99% of SPNs lost).
Step 3: if it under-predicts, fit a single non-cell-autonomous coupling
        constant k_nca and re-validate.

Outputs
-------
results/solver_agreement.csv
results/knca_fit.csv
results/population_validation_coupled.csv
figures/fig4_coupling.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

par = hd.calibrate()
CAGS = [40, 42, 44, 46, 48, 50]
N = 20000
TARGET_LOSS_AT_CAP600 = 0.90   # midpoint of the reported 80-99% range


def loss_at_cap(cag, times, surviving, cap=600.0):
    age = cap / (cag - hd.LIT["cap_offset"])
    return 1.0 - float(np.interp(age, times, surviving))


# ----------------------------------------------- step 1: solver agreement
agree = []
for cag in (42, 46):
    exact = hd.simulate_population(float(cag), par, n_cells=4000, seed=5)
    ages = np.linspace(0, 110, 1101)
    s_exact = hd.survival_curve(exact, ages)
    for dt in (1 / 12, 1 / 24, 1 / 48):
        arr, tt, ss = hd.simulate_population_tau(float(cag), par, n_cells=4000,
                                                 dt=dt, seed=5)
        s_tau = np.interp(ages, tt, ss)
        agree.append({"inherited_CAG": cag, "dt_years": round(dt, 5),
                      "max_abs_survival_difference": round(
                          float(np.max(np.abs(s_tau - s_exact))), 4),
                      "age_50pct_loss_exact": round(hd.age_at_loss(exact, 0.5, ages), 2),
                      "age_50pct_loss_tau": round(
                          float(np.interp(0.5, ss[::-1], tt[::-1])), 2)})

print("solver agreement (exact Gillespie vs tau-leaping)")
for a in agree:
    print("  ", a)

DT = 1 / 24

# ------------------------------- step 2: cell-autonomous model vs CAP anchor
ca_rows = []
curves0 = {}
for cag in CAGS:
    arr, tt, ss = hd.simulate_population_tau(float(cag), par, n_cells=N, dt=DT,
                                             k_nca=0.0, seed=200 + cag)
    curves0[cag] = (tt, ss)
    ca_rows.append({"inherited_CAG": cag, "k_nca": 0.0,
                    "SPN_loss_at_CAP600_pct": round(100 * loss_at_cap(cag, tt, ss), 1)})
mean_loss_ca = np.mean([r["SPN_loss_at_CAP600_pct"] for r in ca_rows])
print(f"\ncell-autonomous model: mean SPN loss at CAP 600 = {mean_loss_ca:.1f}% "
      f"(reported 80-99%)")

# ----------------------------------------- step 3: fit k_nca to the anchor
def mean_loss(k_nca, n=8000):
    vals = []
    for cag in CAGS:
        arr, tt, ss = hd.simulate_population_tau(float(cag), par, n_cells=n,
                                                 dt=DT, k_nca=k_nca,
                                                 seed=900 + cag)
        vals.append(loss_at_cap(cag, tt, ss))
    return float(np.mean(vals))


grid = [0.0, 2.0, 5.0, 10.0, 20.0, 40.0]
scan = [{"k_nca": k, "mean_SPN_loss_at_CAP600_pct": round(100 * mean_loss(k), 1)}
        for k in grid]
print("\nk_nca scan")
for s in scan:
    print("  ", s)

try:
    k_fit = brentq(lambda k: mean_loss(k) - TARGET_LOSS_AT_CAP600, 0.0, 60.0,
                   xtol=0.2, rtol=1e-3, maxiter=40)
except Exception as exc:  # pragma: no cover
    print("brentq failed:", exc)
    k_fit = max(s["k_nca"] for s in scan)
print(f"\nfitted k_nca = {k_fit:.2f}")

# --------------------------------------- re-validate the coupled model
curves1 = {}
coupled_rows = []
for cag in CAGS:
    arr, tt, ss = hd.simulate_population_tau(float(cag), par, n_cells=N, dt=DT,
                                             k_nca=k_fit, seed=400 + cag)
    curves1[cag] = (tt, ss)
    a50 = float(np.interp(0.5, ss[::-1], tt[::-1]))
    a25 = float(np.interp(0.75, ss[::-1], tt[::-1]))
    died = arr["died"]
    life = np.where(died, arr["t_death"], 110.0)
    coupled_rows.append({
        "inherited_CAG": cag,
        "age_25pct_loss": round(a25, 1),
        "age_50pct_loss": round(a50, 1),
        "SPN_loss_at_CAP600_pct": round(100 * loss_at_cap(cag, tt, ss), 1),
        "frac_life_below_150_pct": round(
            100 * float(np.mean(arr["years_below_150"] / life)), 1),
        "ranen_predicted_onset_age": round(hd.ranen_predicted_onset(cag), 1)})

print("\ncoupled model")
for r in coupled_rows:
    print("  ", r)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("solver_agreement.csv", agree, list(agree[0].keys()))
_write("knca_fit.csv", scan + [{"k_nca": round(k_fit, 3),
                                "mean_SPN_loss_at_CAP600_pct": "fitted"}],
       ["k_nca", "mean_SPN_loss_at_CAP600_pct"])
_write("population_validation_coupled.csv", coupled_rows,
       list(coupled_rows[0].keys()))
with open(os.path.join(RES, "k_nca_value.txt"), "w") as f:
    f.write(f"{k_fit:.4f}\n")

# -------------------------------------------------------------------- Fig. 4
fig, ax = plt.subplots(1, 3, figsize=(10.2, 2.95))

for cag, c in zip((42, 46, 50), ("#1f4e79", "#4f81bd", "#c0504d")):
    tt, ss = curves0[cag]
    ax[0].plot(tt, 100 * ss, color=c, lw=1.3, ls="--")
    tt, ss = curves1[cag]
    ax[0].plot(tt, 100 * ss, color=c, lw=1.8, label=f"{cag} CAG")
    age600 = 600.0 / (cag - hd.LIT["cap_offset"])
    ax[0].scatter([age600], [10], marker="v", color=c, s=26, zorder=5)
ax[0].plot([], [], color="0.4", ls="--", lw=1.2, label="cell-autonomous only")
ax[0].plot([], [], color="0.4", lw=1.8, label="with coupling")
ax[0].axhspan(1, 20, color="0.9", lw=0, zorder=0)
ax[0].set_xlabel("age (years)")
ax[0].set_ylabel("surviving SPNs (%)")
ax[0].set_xlim(0, 100)
ax[0].set_title("A  SPN survival\n(triangles: age at CAP 600)", fontsize=8.5)
ax[0].legend(frameon=False, fontsize=6.2, loc="lower left")

ks = [s["k_nca"] for s in scan]
ls_ = [s["mean_SPN_loss_at_CAP600_pct"] for s in scan]
ax[1].plot(ks, ls_, "o-", color="#1f4e79", lw=1.6, ms=4)
ax[1].axhspan(80, 99, color="#c0504d", alpha=0.18, lw=0, label="reported 80-99%")
ax[1].axvline(k_fit, color="k", ls=":", lw=1.0)
ax[1].annotate(f"$k_{{nca}}$ = {k_fit:.1f}", xy=(k_fit, 40), fontsize=7,
               xytext=(k_fit + 3, 35))
ax[1].set_xlabel("non-cell-autonomous coupling $k_{nca}$")
ax[1].set_ylabel("mean SPN loss at CAP 600 (%)")
ax[1].set_title("B  coupling required")
ax[1].legend(frameon=False, fontsize=6.5, loc="lower right")

ax[2].plot(CAGS, [r["age_50pct_loss"] for r in coupled_rows], "o-",
           color="#1f4e79", lw=1.6, ms=4, label="model: 50% SPN loss")
ax[2].plot(CAGS, [r["age_25pct_loss"] for r in coupled_rows], "^-",
           color="#4f81bd", lw=1.3, ms=4, label="model: 25% SPN loss")
ax[2].plot(CAGS, [r["ranen_predicted_onset_age"] for r in coupled_rows], "s--",
           color="#c0504d", lw=1.4, ms=4, label="Ranen regression (onset)")
ax[2].set_xlabel("inherited CAG length")
ax[2].set_ylabel("age (years)")
ax[2].set_title("C  external check")
ax[2].legend(frameon=False, fontsize=6.5)
fig.savefig(os.path.join(FIG, "fig4_coupling.png"))
plt.close(fig)
print("\nwrote figure 4 and result tables")
