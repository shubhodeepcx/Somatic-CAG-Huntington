"""02 -- Stochastic SPN population simulation and validation against held-out
observations that were not used in calibration.

Outputs
-------
results/population_validation.csv
results/survival_by_cag.csv
results/onset_vs_cag.csv
figures/fig3_population.png
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

N_CELLS = 20000
AGES = np.linspace(0, 110, 1101)

par = hd.calibrate()

# ------------------------------------------------- reference cohort: 42 CAG
ref = hd.simulate_population(42.0, par, n_cells=N_CELLS, seed=11)

frac_expanded = float(np.mean(ref["L_max"] > 42.0))
died = ref["died"]
life = np.where(died, ref["t_death"], 110.0)
frac_life_below = float(np.mean(ref["years_below_150"] / life))
# time from crossing 150 CAG to elimination, among cells that died
crossed = ~np.isnan(ref["t_C"]) & died
t_C_to_death = ref["t_death"][crossed] - ref["t_C"][crossed]
# de-repression proxy: maximal length reached before elimination
max_len = ref["L_max"][died]

val = [
    {"observation": "SPNs expanding beyond inherited length",
     "published": "95-98%", "model": f"{100*frac_expanded:.1f}%"},
    {"observation": "fraction of SPN life spent at 42-150 CAG",
     "published": "96.4% (SD 2.0)", "model": f"{100*frac_life_below:.1f}%"},
    {"observation": "time from 150 CAG to elimination (median)",
     "published": "months, not decades",
     "model": f"{np.median(t_C_to_death)*12:.0f} months"},
    {"observation": "median maximal CAG length reached before elimination",
     "published": "100 to >500 observed (max 842)",
     "model": f"{np.median(max_len):.0f}"},
    {"observation": "97.5th percentile of maximal CAG length",
     "published": "max observed 842",
     "model": f"{np.percentile(max_len, 97.5):.0f}"},
]

# ---------------------------------------------- survival / loss across CAG
cags = [40, 42, 44, 46, 48, 50]
surv = {}
onset_rows = []
for cag in cags:
    arr = hd.simulate_population(float(cag), par, n_cells=N_CELLS,
                                 seed=100 + cag)
    s = hd.survival_curve(arr, AGES)
    surv[cag] = s
    a50 = hd.age_at_loss(arr, 0.50, AGES)
    a25 = hd.age_at_loss(arr, 0.25, AGES)
    a80 = hd.age_at_loss(arr, 0.80, AGES)
    onset_rows.append({
        "inherited_CAG": cag,
        "age_25pct_SPN_loss": round(a25, 1),
        "age_50pct_SPN_loss": round(a50, 1),
        "age_80pct_SPN_loss": round(a80, 1),
        "CAP_at_50pct_loss": round(hd.cap_score(a50, cag), 0) if np.isfinite(a50) else None,
        "CAP_at_80pct_loss": round(hd.cap_score(a80, cag), 0) if np.isfinite(a80) else None,
        "ranen_predicted_onset_age": round(hd.ranen_predicted_onset(cag), 1),
    })

# Handsaker: donors with CAP > 600 had lost an estimated 80-99% of SPNs.
# Check the model's SPN loss at the age where CAP = 600, for each CAG length.
cap_rows = []
for cag in cags:
    age_cap600 = 600.0 / (cag - hd.LIT["cap_offset"])
    arr_s = surv[cag]
    loss = 1.0 - float(np.interp(age_cap600, AGES, arr_s))
    cap_rows.append({"inherited_CAG": cag, "age_at_CAP_600": round(age_cap600, 1),
                     "model_SPN_loss_at_CAP_600_pct": round(100 * loss, 1)})

val.append({"observation": "SPN loss at CAP = 600 (mean over CAG 40-50)",
            "published": "80-99%",
            "model": f"{np.mean([r['model_SPN_loss_at_CAP_600_pct'] for r in cap_rows]):.1f}%"})

# ------------------------------------------------------------------ outputs
def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("population_validation.csv", val, ["observation", "published", "model"])
_write("onset_vs_cag.csv", onset_rows, list(onset_rows[0].keys()))
_write("cap600_check.csv", cap_rows, list(cap_rows[0].keys()))

with open(os.path.join(RES, "survival_by_cag.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["age_years"] + [f"surviving_fraction_CAG{c}" for c in cags])
    for i, a in enumerate(AGES):
        w.writerow([round(a, 2)] + [round(surv[c][i], 5) for c in cags])

print("held-out validation (none of these were fitted)")
for v in val:
    print(f"  {v['observation']:<50} published={v['published']:<22} model={v['model']}")
print("\nonset proxies")
for r in onset_rows:
    print("  ", r)

# -------------------------------------------------------------------- Fig. 3
fig, ax = plt.subplots(1, 3, figsize=(10.2, 2.95))

cmap = plt.get_cmap("viridis")
for i, c in enumerate(cags):
    ax[0].plot(AGES, 100 * surv[c], color=cmap(i / (len(cags) - 1)), lw=1.6,
               label=f"{c} CAG")
ax[0].set_xlabel("age (years)")
ax[0].set_ylabel("surviving SPNs (%)")
ax[0].set_xlim(0, 100)
ax[0].set_title("A  simulated SPN survival")
ax[0].legend(frameon=False, fontsize=6.5, ncol=2)

ages_m = [r["age_50pct_SPN_loss"] for r in onset_rows]
ranen = [r["ranen_predicted_onset_age"] for r in onset_rows]
ax[1].plot(cags, ages_m, "o-", color="#1f4e79", lw=1.6, ms=4,
           label="model: age at 50% SPN loss")
ax[1].plot(cags, ranen, "s--", color="#c0504d", lw=1.4, ms=4,
           label="Ranen regression (clinical onset)")
ax[1].set_xlabel("inherited CAG length")
ax[1].set_ylabel("age (years)")
ax[1].set_title("B  external check on onset")
ax[1].legend(frameon=False, fontsize=6.5)

ax[2].hist(max_len, bins=np.arange(100, 1050, 25), color="#8064a2",
           edgecolor="white", linewidth=0.3)
ax[2].axvline(hd.LIT["max_observed_CAG"], color="k", ls="--", lw=1.0)
ax[2].annotate("largest length\nreported (842)", xy=(hd.LIT["max_observed_CAG"], 0),
               xytext=(600, ax[2].get_ylim()[1] * 0.6), fontsize=6.5,
               arrowprops=dict(arrowstyle="->", lw=0.7))
ax[2].set_xlabel("maximal CAG length before elimination")
ax[2].set_ylabel("SPNs")
ax[2].set_title("C  predicted terminal lengths (42 CAG)")
fig.savefig(os.path.join(FIG, "fig3_population.png"))
plt.close(fig)
print("\nwrote figure 3 and result tables")
