"""13 -- From an idealised step change to a real dosing schedule.

The in silico trial in script 06 assumed that treatment sets striatal MSH3 to a
new level instantly and holds it there for life. A real intrathecal regimen is
a series of boluses that decay between doses, so the achieved knockdown
oscillates between a peak and a trough.

Because the expansion rate is first-order in MSH3 protein (script 01), the
cumulative effect over decades depends only on the time-averaged MSH3 level,
not on its waveform. That makes the correction exact rather than approximate:
compute the time average, feed it back into the trial as an effective residual
MSH3, and the lifetime prediction follows.

Neither the CNS half-life of a stabilised oligonucleotide nor the eventual
dosing interval is settled, so both are swept rather than assumed.

Outputs
-------
results/dosing_timeaverage.csv
results/dosing_benefit.csv
figures/fig14_dosing.png
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

EC50, HILL = 10.12, 1.215     # fitted in script 05 (nmol/g, dimensionless)
HALF_LIVES_M = [1, 2, 3, 6, 9, 12]        # months, tissue elimination
INTERVALS_M = [3, 4, 6, 12]               # months between doses
PEAK_EXPOSURES = [10.12, 25.0, 61.8]      # nmol/g: 50, 75 and 90% peak knockdown


def msh3_from_exposure(x):
    return 100.0 / (1.0 + (np.asarray(x, float) / EC50) ** HILL)


def steady_state_average(peak_dose, half_life_m, interval_m, n_pts=4000):
    """Time-averaged MSH3 (% of control) over one steady-state dosing cycle.

    Accumulation is geometric: with elimination constant k and interval tau,
    the steady-state peak is D / (1 - exp(-k*tau)).
    """
    k = np.log(2.0) / half_life_m
    tau = interval_m
    peak = peak_dose / (1.0 - np.exp(-k * tau))
    t = np.linspace(0.0, tau, n_pts)
    conc = peak * np.exp(-k * t)
    msh3 = msh3_from_exposure(conc)
    return (float(np.trapezoid(msh3, t) / tau), float(msh3.min()),
            float(msh3.max()), float(peak))


rows = []
for dose in PEAK_EXPOSURES:
    single_peak_msh3 = float(msh3_from_exposure(dose))
    for hl in HALF_LIVES_M:
        for iv in INTERVALS_M:
            avg, best, worst, peak_conc = steady_state_average(dose, hl, iv)
            rows.append({
                "per_dose_exposure_nmol_per_g": dose,
                "single_dose_peak_knockdown_pct": round(100 - single_peak_msh3, 1),
                "tissue_half_life_months": hl,
                "dosing_interval_months": iv,
                "steady_state_peak_conc_nmol_per_g": round(peak_conc, 2),
                "msh3_at_peak_pct": round(best, 1),
                "msh3_at_trough_pct": round(worst, 1),
                "time_averaged_msh3_pct": round(avg, 1),
                "time_averaged_knockdown_pct": round(100 - avg, 1),
                "expansion_suppression_pct": round(
                    float(hd.prevention_from_msh3(avg)), 1)})

print("time-averaged knockdown, per-dose exposure 25 nmol/g")
for r in rows:
    if r["per_dose_exposure_nmol_per_g"] == 25.0 and r["dosing_interval_months"] in (4, 6):
        print(f"  t1/2 {r['tissue_half_life_months']:>2} mo, every "
              f"{r['dosing_interval_months']:>2} mo -> peak "
              f"{100 - r['msh3_at_peak_pct']:5.1f}%, trough "
              f"{100 - r['msh3_at_trough_pct']:5.1f}%, average "
              f"{r['time_averaged_knockdown_pct']:5.1f}%")

# Feed representative time-averaged levels back into the lifetime model.
par = hd.calibrate()
N, DT = 9000, 1 / 24
AGES = np.linspace(0, 110, 1101)


def loss_at_70(residual, t_treat=30.0, seed=0):
    _, tt, ss = hd.simulate_population_tau(
        42.0, par, n_cells=N, dt=DT,
        t_treat=t_treat, msh3_treat=float(residual), seed=seed)
    return 100 * (1 - float(np.interp(70.0, tt, ss)))


untreated = loss_at_70(1.0, t_treat=np.inf, seed=11)
bench_rows = []
cases = [("idealised step, 75% lowering", 0.25),
         ("6-monthly, 6-month half-life", None),
         ("6-monthly, 3-month half-life", None),
         ("12-monthly, 3-month half-life", None),
         ("4-monthly, 9-month half-life", None)]
lookup = {(6, 6): None, (6, 3): None, (12, 3): None, (4, 9): None}
for r in rows:
    if r["per_dose_exposure_nmol_per_g"] == 25.0:
        key = (r["dosing_interval_months"], r["tissue_half_life_months"])
        if key in lookup:
            lookup[key] = r["time_averaged_msh3_pct"] / 100.0

residuals = [0.25, lookup[(6, 6)], lookup[(6, 3)], lookup[(12, 3)],
             lookup[(4, 9)]]
# Seeds come from the case index, not from hash(label), because Python
# randomises string hashing per process and the numbers have to reproduce.
for k, ((label, _), res) in enumerate(zip(cases, residuals)):
    l70 = loss_at_70(res, seed=7100 + 17 * k)
    bench_rows.append({"regimen": label,
                       "effective_residual_MSH3": round(float(res), 3),
                       "effective_knockdown_pct": round(100 * (1 - res), 1),
                       "SPN_loss_at_70_pct": round(l70, 1),
                       "untreated_loss_at_70_pct": round(untreated, 1),
                       "benefit_points": round(untreated - l70, 1)})

print("\nlifetime effect of real schedules (42 CAG, from age 30)")
for r in bench_rows:
    print("  ", r)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("dosing_timeaverage.csv", rows, list(rows[0].keys()))
_write("dosing_benefit.csv", bench_rows, list(bench_rows[0].keys()))

fig, ax = plt.subplots(1, 3, figsize=(10.4, 2.9))

# A: one steady-state cycle for three half-lives, 6-monthly dosing
t = np.linspace(0, 12, 600)
for hl, c in zip((3, 6, 12), ("#c0504d", "#1f4e79", "#4f6228")):
    k = np.log(2) / hl
    tau = 6.0
    peak = 25.0 / (1 - np.exp(-k * tau))
    conc = peak * np.exp(-k * (t % tau))
    ax[0].plot(t, 100 - msh3_from_exposure(conc), color=c, lw=1.6,
               label=f"half-life {hl} mo")
ax[0].set_xlabel("months within a 6-monthly schedule")
ax[0].set_ylabel("MSH3 knockdown (%)")
ax[0].set_ylim(0, 100)
ax[0].set_title("A  what the patient actually gets")
ax[0].legend(frameon=False, fontsize=6.5, loc="lower right")

# B: time-averaged knockdown map
sub = [r for r in rows if r["per_dose_exposure_nmol_per_g"] == 25.0]
grid = np.zeros((len(INTERVALS_M), len(HALF_LIVES_M)))
for r in sub:
    i = INTERVALS_M.index(r["dosing_interval_months"])
    j = HALF_LIVES_M.index(r["tissue_half_life_months"])
    grid[i, j] = r["time_averaged_knockdown_pct"]
im = ax[1].imshow(grid, origin="lower", aspect="auto", cmap="viridis",
                  vmin=0, vmax=100)
ax[1].set_xticks(range(len(HALF_LIVES_M)))
ax[1].set_xticklabels(HALF_LIVES_M)
ax[1].set_yticks(range(len(INTERVALS_M)))
ax[1].set_yticklabels(INTERVALS_M)
ax[1].set_xlabel("tissue half-life (months)")
ax[1].set_ylabel("dosing interval (months)")
for i in range(len(INTERVALS_M)):
    for j in range(len(HALF_LIVES_M)):
        ax[1].text(j, i, f"{grid[i, j]:.0f}", ha="center", va="center",
                   fontsize=6.5, color="white" if grid[i, j] < 60 else "black")
cb = fig.colorbar(im, ax=ax[1], fraction=0.046)
cb.set_label("time-averaged knockdown (%)", fontsize=7)
cb.ax.tick_params(labelsize=7)
ax[1].set_title("B  75% peak knockdown per dose", fontsize=9)

labels = [r["regimen"] for r in bench_rows]
vals = [r["benefit_points"] for r in bench_rows]
cols = ["#8064a2"] + ["#1f4e79"] * (len(vals) - 1)
ax[2].barh(range(len(vals)), vals, color=cols, height=0.52)
# Regimen names go above their bars; as tick labels they would run out of the
# axes and over the colour bar of panel B.
ax[2].set_yticks(range(len(labels)))
ax[2].set_yticklabels([])
ax[2].tick_params(axis="y", length=0)
for k, lab in enumerate(labels):
    ax[2].annotate(lab, xy=(0.6, k + 0.32), fontsize=6.1, ha="left",
                   va="bottom")
ax[2].set_xlim(0, max(vals) * 1.12)
ax[2].set_ylim(-0.6, len(vals) - 0.3)
ax[2].set_xlabel("striatum preserved at age 70\n(percentage points)")
ax[2].set_title("C  lifetime cost of intermittency", fontsize=9)
fig.subplots_adjust(wspace=0.34)
fig.savefig(os.path.join(FIG, "fig14_dosing.png"))
plt.close(fig)
print("\nwrote figure 14 and result tables")
