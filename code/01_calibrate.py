"""01 -- Calibrate the biochemical rate law and test it against held-out anchors.

Outputs
-------
results/fitted_parameters.json
results/calibration_check.csv
results/model_comparison_msh3.csv
results/msh3_measured_points.csv
results/identifiability_profiles.csv
figures/fig1_rate_law.png
figures/fig2_msh3_order.png
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
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

RT = 0.593  # kcal/mol at 310 K

# ---------------------------------------------------------------- calibration
par = hd.calibrate(verbose=True)

rows = []
for (L1, L2, obs) in [(40, 60, hd.LIT["t_40_60_years"]),
                      (60, 80, hd.LIT["t_60_80_years"]),
                      (80, 150, hd.LIT["t_80_150_years"])]:
    pred = hd.mean_transit_time(L1, L2, par)
    rows.append({"quantity": f"mean time {L1}->{L2} CAG (y)",
                 "published": obs, "model": round(pred, 2), "role": "calibration"})

fold = float(hd.lam(90.0, par) / hd.lam(70.0, par))
rows.append({"quantity": "lambda(90)/lambda(70)",
             "published": f">={hd.LIT['rate_rise_70_90_min']:.0f}",
             "model": round(fold, 2), "role": "calibration"})

# --- held-out validation anchors -------------------------------------------
lam_max_low = float(np.max(hd.lam(np.arange(36, 56, dtype=float), par)))
rows.append({"quantity": "max lambda, 36-55 CAG (events/cell/y)",
             "published": "<1", "model": round(lam_max_low, 4),
             "role": "held-out"})

t_150_350 = hd.mean_transit_time(150, 350, par)
rows.append({"quantity": "mean time 150->350 CAG (months)",
             "published": "months, not decades", "model": round(t_150_350 * 12, 1),
             "role": "held-out"})

# ------------------------------------------------- MSH3 order model comparison
msh3_grid = np.linspace(0, 100, 21)
obs_pd = hd.prevention_from_msh3(msh3_grid)
cmp_rows = []
for order in (1.0, 1.5, 2.0, 3.0):
    pred = hd.model_prevention_from_msh3(msh3_grid, order=order)
    cmp_rows.append({"binding_order_k": order,
                     "rmse_vs_published_regression": round(
                         float(np.sqrt(np.mean((pred - obs_pd) ** 2))), 2)})
best_order = min(cmp_rows, key=lambda r: r["rmse_vs_published_regression"])

meas = [(hd.LIT["msh3_protein_10nmol_pct"], hd.LIT["expansion_prevented_10nmol"]),
        (hd.LIT["msh3_protein_3nmol_pct"], hd.LIT["expansion_prevented_3nmol"])]
meas_rows = [{"msh3_protein_pct_of_control": m,
              "measured_pct_expansion_prevented": y,
              "model_k1": round(float(hd.model_prevention_from_msh3(m, 1.0)), 1),
              "model_k2": round(float(hd.model_prevention_from_msh3(m, 2.0)), 1),
              "model_k3": round(float(hd.model_prevention_from_msh3(m, 3.0)), 1)}
             for m, y in meas]

# ------------------------------------------------------- identifiability scan
prof_rows = []
for tag, kwargs in [("baseline", {}),
                    ("t(80->150) = 2 y", {"targets": {"t3": 2.0}}),
                    ("t(80->150) = 8 y", {"targets": {"t3": 8.0}}),
                    ("fold-rise = 10", {"fold_target": 10.0}),
                    ("fold-rise = 20", {"fold_target": 20.0}),
                    ("switch width 3 repeats", {"assumed": {"w_switch": 3.0}}),
                    ("switch width 12 repeats", {"assumed": {"w_switch": 12.0}}),
                    ("phi ceiling L50 = 250", {"assumed": {"L50_ceiling": 250.0}}),
                    ("phi ceiling L50 = 800", {"assumed": {"L50_ceiling": 800.0}})]:
    p = hd.calibrate(**kwargs)
    prof_rows.append({"variant": tag,
                      "lambda40_per_year": round(p["lambda40"], 5),
                      "beta_per_repeat": round(p["beta"], 4),
                      "alpha_kcal_mol_repeat": round(p["beta"] * RT, 4),
                      "R_cluster": round(p["R"], 2),
                      "L_switch": round(p["L_sw"], 1),
                      "fit_cost": round(p["cost"], 6)})

# ------------------------------------------------------------------- outputs
with open(os.path.join(RES, "fitted_parameters.json"), "w") as f:
    json.dump({"fitted": {k: (None if v == np.inf else v) for k, v in par.items()},
               "derived": {"alpha_kcal_per_mol_per_repeat": round(par["beta"] * RT, 4),
                           "lambda_at_40": round(float(hd.lam(40.0, par)), 5),
                           "lambda_at_150": round(float(hd.lam(150.0, par)), 3)},
               "interpretation": {
                   "A": "events per year per MutSbeta site at unit slip-out occupancy",
                   "beta": "d(dG_slip)/dL divided by RT, per repeat unit",
                   "R": "fold rise in expansion probability for a clustered vs "
                        "isolated slip-out",
                   "L_sw": "tract length at which half of slip-outs are clustered"},
               "assumed_structural_constants": hd.ASSUMED,
               "best_msh3_binding_order": best_order}, f, indent=2)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("calibration_check.csv", rows, ["quantity", "published", "model", "role"])
_write("model_comparison_msh3.csv", cmp_rows,
       ["binding_order_k", "rmse_vs_published_regression"])
_write("msh3_measured_points.csv", meas_rows, list(meas_rows[0].keys()))
_write("identifiability_profiles.csv", prof_rows, list(prof_rows[0].keys()))

print("\ncalibration / validation")
for r in rows:
    print(f"  {r['quantity']:<42} published={str(r['published']):>18}  "
          f"model={str(r['model']):<9} [{r['role']}]")
print("\nMSH3 binding order (RMSE vs published regression):")
for r in cmp_rows:
    print("   k =", r["binding_order_k"], "->", r["rmse_vs_published_regression"])
print("\nalpha =", round(par["beta"] * RT, 4), "kcal/mol per repeat")
print("\nidentifiability:")
for r in prof_rows:
    print("  ", r)

# -------------------------------------------------------------------- Fig. 1
Lg = np.linspace(36, 420, 1500)
fig, ax = plt.subplots(1, 3, figsize=(10.0, 2.95))

C = hd.clustered_fraction(Lg, par["L_sw"], par["w_sw"])
phi = hd.phi_slipout(Lg, par["beta"], par["L50"])
nst = hd.n_sites(Lg)
ax[0].plot(Lg, C, color="#c0504d", lw=1.8, label="$C(L)$ clustered slip-outs")
ax[0].plot(Lg, phi / phi.max(), color="#1f4e79", lw=1.6,
           label=r"$\phi(L)/\phi_{\max}$ slip-out occupancy")
ax[0].plot(Lg, nst / nst.max(), color="#4f6228", lw=1.4, ls="--",
           label=r"$n(L)/n_{\max}$ MutS$\beta$ sites")
ax[0].axvline(hd.LIT["L_phaseC"], color="0.5", lw=0.8, ls=":")
ax[0].set_xlabel("CAG repeat length $L$")
ax[0].set_ylabel("normalised")
ax[0].set_title("A  substrate layer")
ax[0].legend(frameon=False, fontsize=6.5, loc="center right")

ax[1].semilogy(Lg, hd.lam(Lg, par), color="#1f4e79", lw=1.9, label="untreated")
for m, c, lab in [(0.5, "#4f81bd", "50% MSH3 lowering"),
                  (0.25, "#c0504d", "75% MSH3 lowering")]:
    ax[1].semilogy(Lg, hd.lam(Lg, par, msh3=m), color=c, lw=1.3, ls="--", label=lab)
ax[1].axvline(hd.LIT["L_phaseC"], color="0.5", lw=0.8, ls=":")
ax[1].axhline(1.0, color="0.75", lw=0.8, ls="--")
ax[1].set_xlabel("CAG repeat length $L$")
ax[1].set_ylabel(r"$\lambda$ (events cell$^{-1}$ y$^{-1}$)")
ax[1].set_title("B  fitted mutation rate")
ax[1].legend(frameon=False, fontsize=6.5, loc="lower right")

seg = [(40, 60, hd.LIT["t_40_60_years"]), (60, 80, hd.LIT["t_60_80_years"]),
       (80, 150, hd.LIT["t_80_150_years"])]
pred = [hd.mean_transit_time(a, b, par) for a, b, _ in seg]
x = np.arange(3)
ax[2].bar(x - 0.18, [s[2] for s in seg], width=0.34, color="#8064a2",
          label="published estimate")
ax[2].bar(x + 0.18, pred, width=0.34, color="#1f4e79", label="model")
ax[2].set_xticks(x)
ax[2].set_xticklabels([f"{a}→{b}" for a, b, _ in seg])
ax[2].set_ylabel("mean transit time (years)")
ax[2].set_title("C  calibration")
ax[2].legend(frameon=False, fontsize=7)
fig.savefig(os.path.join(FIG, "fig1_rate_law.png"))
plt.close(fig)

# -------------------------------------------------------------------- Fig. 2
m = np.linspace(0, 110, 400)
fig, ax = plt.subplots(figsize=(4.1, 3.2))
base = hd.prevention_from_msh3(m)
ax.fill_between(m, np.clip(base - 12, 0, 100), np.clip(base + 12, 0, 100),
                color="0.86", lw=0, label="regression $\\pm$ 12% band")
ax.plot(m, base, color="k", lw=2.0, label="published regression (mouse striatum)")
for order, c in zip((1.0, 2.0, 3.0), ("#1f4e79", "#4f81bd", "#c0504d")):
    ax.plot(m, hd.model_prevention_from_msh3(m, order), lw=1.5, ls="--", color=c,
            label=f"model, $k$ = {order:.0f}")
ax.scatter([p[0] for p in meas], [p[1] for p in meas], s=38, color="#e36c0a",
           zorder=5, edgecolor="k", linewidth=0.4,
           label="measured (3 and 10 nmol ICV)")
ax.set_xlabel("striatal MSH3 protein (% of control)")
ax.set_ylabel("% somatic expansion prevented")
ax.set_xlim(0, 112)
ax.set_ylim(-2, 102)
ax.legend(frameon=False, fontsize=6.5, loc="lower left")
ax.set_title("MutS$\\beta$ stoichiometry inferred from the\nMSH3 dose-response",
             fontsize=9)
fig.savefig(os.path.join(FIG, "fig2_msh3_order.png"))
plt.close(fig)
print("\nwrote figures 1-2 and result tables")
