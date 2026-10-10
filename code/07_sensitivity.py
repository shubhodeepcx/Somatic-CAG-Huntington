"""07 -- Global sensitivity analysis of the therapeutic prediction.

Every quantity that was assumed rather than measured, plus the published
calibration anchors themselves (which carry real uncertainty), is varied
simultaneously over a Latin hypercube. For each sample the rate law is
re-calibrated to the (perturbed) anchors and the whole pipeline is re-run, so
the analysis propagates calibration uncertainty, not just forward uncertainty.

Outcome measures
----------------
y1  SPN loss at age 70 with no treatment (%)
y2  SPN loss at age 70 with the *measured achievable* striatal MSH3 knockdown
    (29%) started at age 30 (%)
y3  absolute benefit, y1 - y2 (percentage points of striatum preserved)

Sensitivity index: partial rank correlation coefficient (PRCC) of each input
against each outcome, which is the standard choice for monotonic but nonlinear
model responses.

Outputs
-------
results/sensitivity_samples.csv
results/sensitivity_prcc.csv
figures/fig10_sensitivity.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

N_SAMPLES = 160
N_CELLS = 3000
DT = 1 / 24
CAG = 42.0
RESIDUAL_ACHIEVABLE = 0.71   # measured putamen SPN residual MSH3
START_AGE = 30.0

# input name -> (low, high, kind)
RANGES = {
    # published anchors, varied over their stated or plausible uncertainty
    "t_40_60_years":      (40.0, 60.0, "anchor"),
    "t_60_80_years":      (9.0, 15.0, "anchor"),
    "t_80_150_years":     (2.0, 8.0, "anchor"),
    "fold_rise_70_90":    (6.0, 20.0, "anchor"),
    "L_phaseC":           (120.0, 180.0, "anchor"),
    "step_repeats":       (2.0, 8.0, "anchor"),
    # structural constants of the rate law
    "w_switch":           (3.0, 12.0, "assumed"),
    "L50_ceiling":        (250.0, 800.0, "assumed"),
    "mutsb_footprint":    (9.0, 25.0, "assumed"),
    # elimination submodel
    "L_death_50":         (200.0, 450.0, "assumed"),
    "death_width":        (20.0, 80.0, "assumed"),
    "h_death_max":        (2.0, 12.0, "assumed"),
    # contraction submodel
    "contraction_ratio":  (0.0, 0.30, "assumed"),
    "contraction_step":   (1.0, 4.0, "assumed"),
    # pharmacodynamic coupling: order of the MSH3 dependence of lambda
    "msh3_order":         (0.8, 1.3, "anchor"),
}

rng = np.random.default_rng(20261010)
names = list(RANGES)
k = len(names)

# Latin hypercube
X = np.empty((N_SAMPLES, k))
for j, nm in enumerate(names):
    lo, hi = RANGES[nm][0], RANGES[nm][1]
    perm = rng.permutation(N_SAMPLES)
    X[:, j] = lo + (hi - lo) * (perm + rng.random(N_SAMPLES)) / N_SAMPLES

Y = np.full((N_SAMPLES, 3), np.nan)

lit_keys = ["t_40_60_years", "t_60_80_years", "t_80_150_years", "L_phaseC",
            "step_repeats", "mutsb_footprint_repeats"]
lit_backup = {kk: hd.LIT[kk] for kk in lit_keys}
assumed_backup = dict(hd.ASSUMED)

for i in range(N_SAMPLES):
    s = dict(zip(names, X[i]))
    try:
        hd.LIT["t_40_60_years"] = s["t_40_60_years"]
        hd.LIT["t_60_80_years"] = s["t_60_80_years"]
        hd.LIT["t_80_150_years"] = s["t_80_150_years"]
        hd.LIT["L_phaseC"] = s["L_phaseC"]
        hd.LIT["step_repeats"] = s["step_repeats"]
        hd.LIT["mutsb_footprint_repeats"] = s["mutsb_footprint"]
        hd.ASSUMED["w_switch"] = s["w_switch"]
        hd.ASSUMED["L50_ceiling"] = s["L50_ceiling"]
        hd.ASSUMED["contraction_ratio"] = s["contraction_ratio"]
        hd.ASSUMED["contraction_step"] = s["contraction_step"]

        p = hd.calibrate(fold_target=s["fold_rise_70_90"], quick=True)
        hazard = {"L_death_50": s["L_death_50"], "width": s["death_width"],
                  "h_max": s["h_death_max"]}

        # untreated
        _, tt0, ss0 = hd.simulate_population_tau(
            CAG, p, n_cells=N_CELLS, dt=DT, hazard_kw=hazard, seed=1000 + i)
        y1 = 100 * (1 - float(np.interp(70.0, tt0, ss0)))

        # treated: residual MSH3 enters lambda with the sampled order
        m_eff = RESIDUAL_ACHIEVABLE ** s["msh3_order"]
        _, tt1, ss1 = hd.simulate_population_tau(
            CAG, p, n_cells=N_CELLS, dt=DT, hazard_kw=hazard,
            t_treat=START_AGE, msh3_treat=m_eff, seed=1000 + i)
        y2 = 100 * (1 - float(np.interp(70.0, tt1, ss1)))
        Y[i] = (y1, y2, y1 - y2)
    except Exception as exc:
        print(f"sample {i} failed: {exc}")
    finally:
        hd.LIT.update(lit_backup)
        hd.ASSUMED.clear()
        hd.ASSUMED.update(assumed_backup)
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{N_SAMPLES} samples done")

ok = np.all(np.isfinite(Y), axis=1)
print(f"\n{ok.sum()}/{N_SAMPLES} samples usable")
Xo, Yo = X[ok], Y[ok]

print("\noutcome distributions (median [2.5th, 97.5th percentile])")
labels = ["untreated SPN loss at 70 (%)",
          "treated SPN loss at 70 (%)",
          "absolute benefit (percentage points)"]
summary = []
for j, lab in enumerate(labels):
    q = np.percentile(Yo[:, j], [2.5, 50, 97.5])
    print(f"  {lab:<42} {q[1]:6.1f}  [{q[0]:.1f}, {q[2]:.1f}]")
    summary.append({"outcome": lab, "median": round(q[1], 2),
                    "p2.5": round(q[0], 2), "p97.5": round(q[2], 2)})


def prcc(Xm, y):
    """Partial rank correlation coefficient of every column of Xm against y."""
    Rx = np.apply_along_axis(rankdata, 0, Xm).astype(float)
    ry = rankdata(y).astype(float)
    n, m = Rx.shape
    out = np.zeros(m)
    for j in range(m):
        others = np.delete(np.arange(m), j)
        Z = np.column_stack([np.ones(n), Rx[:, others]])
        # residualise both the input of interest and the outcome on the rest
        bx, *_ = np.linalg.lstsq(Z, Rx[:, j], rcond=None)
        by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
        rx = Rx[:, j] - Z @ bx
        res_y = ry - Z @ by
        denom = np.sqrt(np.sum(rx ** 2) * np.sum(res_y ** 2))
        out[j] = float(np.sum(rx * res_y) / denom) if denom > 0 else 0.0
    return out


prcc_rows = []
P = np.zeros((3, k))
for j in range(3):
    P[j] = prcc(Xo, Yo[:, j])
for idx, nm in enumerate(names):
    rho_s = spearmanr(Xo[:, idx], Yo[:, 2]).statistic
    prcc_rows.append({"input": nm, "kind": RANGES[nm][2],
                      "low": RANGES[nm][0], "high": RANGES[nm][1],
                      "PRCC_untreated_loss": round(P[0, idx], 3),
                      "PRCC_treated_loss": round(P[1, idx], 3),
                      "PRCC_benefit": round(P[2, idx], 3),
                      "spearman_benefit": round(float(rho_s), 3)})

prcc_rows.sort(key=lambda r: -abs(r["PRCC_benefit"]))
print("\nPRCC against absolute benefit (ranked)")
for r in prcc_rows:
    print(f"  {r['input']:<20} {r['PRCC_benefit']:+.3f}   [{r['kind']}]")


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("sensitivity_prcc.csv", prcc_rows, list(prcc_rows[0].keys()))
_write("sensitivity_summary.csv", summary, list(summary[0].keys()))
with open(os.path.join(RES, "sensitivity_samples.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(names + ["untreated_loss_70", "treated_loss_70", "benefit"])
    for xi, yi in zip(Xo, Yo):
        w.writerow([round(v, 5) for v in xi] + [round(v, 3) for v in yi])

# ------------------------------------------------------------------- Fig. 10
fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.1))

ax[0].hist(Yo[:, 2], bins=24, color="#4f81bd", edgecolor="white", linewidth=0.3)
ax[0].axvline(np.median(Yo[:, 2]), color="k", ls="--", lw=1.1)
ax[0].set_xlabel("striatum preserved at age 70\n(percentage points)")
ax[0].set_ylabel("samples")
ax[0].set_title("A  benefit of achievable\nMSH3 lowering", fontsize=8.5)

vals = [r["PRCC_benefit"] for r in prcc_rows][::-1]
labs = [r["input"] for r in prcc_rows][::-1]
kinds = [r["kind"] for r in prcc_rows][::-1]
cols = ["#1f4e79" if kk == "anchor" else "#c0504d" for kk in kinds]
ax[1].barh(range(len(vals)), vals, color=cols, height=0.65)
ax[1].set_yticks(range(len(labs)))
ax[1].set_yticklabels(labs, fontsize=5.8)
ax[1].axvline(0, color="k", lw=0.7)
ax[1].set_xlabel("PRCC vs benefit")
ax[1].set_title("B  sensitivity\n(blue: published anchor, red: assumed)",
                fontsize=8.5)

top = prcc_rows[0]["input"]
jtop = names.index(top)
ax[2].scatter(Xo[:, jtop], Yo[:, 2], s=11, color="#8064a2", alpha=0.75)
ax[2].set_xlabel(top)
ax[2].set_ylabel("benefit (percentage points)")
ax[2].set_title(f"C  strongest driver\nPRCC = {prcc_rows[0]['PRCC_benefit']:+.2f}",
                fontsize=8.5)
fig.savefig(os.path.join(FIG, "fig10_sensitivity.png"))
plt.close(fig)
print("\nwrote figure 10 and result tables")
