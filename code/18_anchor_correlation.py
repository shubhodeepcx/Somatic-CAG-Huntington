"""18 -- Does treating the calibration anchors as independent misstate the
spread of the therapeutic prediction?

`07_sensitivity.py` samples the three published dwell-time anchors
(40->60, 60->80, 80->150 CAG) independently over their plausible ranges. All
three come from the same single-cell study on the same post-mortem cohort, so a
systematic error in that study's clock would move all three together. Sampling
them independently lets such errors cancel across the ensemble, which could make
the reported spread of the benefit too narrow.

This script repeats the global analysis with the three dwell-time anchors
coupled through a Gaussian copula, at several strengths of coupling, and compares
the resulting spread of the benefit with the independent case.

Design
------
A Gaussian copula is used rather than a shared multiplicative factor so that
every anchor keeps *exactly* the marginal distribution it has in
`07_sensitivity.py` (uniform over the same range). Only the dependence between
them changes, so the comparison isolates the effect of correlation and nothing
else. Pearson rho on the latent normals maps to a Spearman rank correlation of
(6/pi)*arcsin(rho/2), which is reported alongside.

rho = 0 reproduces the independent design (up to sampling: this script uses
copula draws rather than a Latin hypercube, so the rho = 0 arm is also the
control for that change of sampler).

Outputs
-------
results/anchor_correlation.csv     spread of each outcome at each rho
figures/fig16_anchor_correlation.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, spearmanr

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
RESIDUAL_ACHIEVABLE = 0.71
START_AGE = 30.0

# Identical to 07_sensitivity.py, so the two analyses are comparable.
RANGES = {
    "t_40_60_years":      (40.0, 60.0, "anchor"),
    "t_60_80_years":      (9.0, 15.0, "anchor"),
    "t_80_150_years":     (2.0, 8.0, "anchor"),
    "fold_rise_70_90":    (6.0, 20.0, "anchor"),
    "L_phaseC":           (120.0, 180.0, "anchor"),
    "step_repeats":       (2.0, 8.0, "anchor"),
    "w_switch":           (3.0, 12.0, "assumed"),
    "L50_ceiling":        (250.0, 800.0, "assumed"),
    "mutsb_footprint":    (9.0, 25.0, "assumed"),
    "L_death_50":         (200.0, 450.0, "assumed"),
    "death_width":        (20.0, 80.0, "assumed"),
    "h_death_max":        (2.0, 12.0, "assumed"),
    "contraction_ratio":  (0.0, 0.30, "assumed"),
    "contraction_step":   (1.0, 4.0, "assumed"),
    "L_max":              (600.0, 1400.0, "assumed"),
    "msh3_order":         (0.8, 1.3, "anchor"),
}

# The anchors that share a source and are therefore coupled here.
COUPLED = ["t_40_60_years", "t_60_80_years", "t_80_150_years"]

RHOS = [0.0, 0.5, 0.9]

names = list(RANGES)
k = len(names)
coupled_idx = [names.index(nm) for nm in COUPLED]


def sample_inputs(rho, seed):
    """Copula sample: uniform marginals, latent correlation rho among COUPLED."""
    rng = np.random.default_rng(seed)
    # latent covariance: identity, with an equicorrelated block on the coupled
    # anchors (equicorrelation with rho >= 0 is positive definite for any size)
    S = np.eye(k)
    for a in coupled_idx:
        for b in coupled_idx:
            if a != b:
                S[a, b] = rho
    Lc = np.linalg.cholesky(S)
    Z = rng.standard_normal((N_SAMPLES, k)) @ Lc.T
    U = norm.cdf(Z)
    X = np.empty_like(U)
    for j, nm in enumerate(names):
        lo, hi = RANGES[nm][0], RANGES[nm][1]
        X[:, j] = lo + (hi - lo) * U[:, j]
    return X


def run_ensemble(X):
    """Return the (n, 3) outcome matrix for a matrix of sampled inputs."""
    lit_keys = ["t_40_60_years", "t_60_80_years", "t_80_150_years", "L_phaseC",
                "step_repeats", "mutsb_footprint_repeats"]
    lit_backup = {kk: hd.LIT[kk] for kk in lit_keys}
    assumed_backup = dict(hd.ASSUMED)
    Y = np.full((len(X), 3), np.nan)

    for i in range(len(X)):
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
            hd.ASSUMED["L_max"] = s["L_max"]

            p = hd.calibrate(fold_target=s["fold_rise_70_90"], quick=True)
            hazard = {"L_death_50": s["L_death_50"], "width": s["death_width"],
                      "h_max": s["h_death_max"]}

            _, tt0, ss0 = hd.simulate_population_tau(
                CAG, p, n_cells=N_CELLS, dt=DT, hazard_kw=hazard, seed=1000 + i)
            y1 = 100 * (1 - float(np.interp(70.0, tt0, ss0)))

            m_eff = RESIDUAL_ACHIEVABLE ** s["msh3_order"]
            _, tt1, ss1 = hd.simulate_population_tau(
                CAG, p, n_cells=N_CELLS, dt=DT, hazard_kw=hazard,
                t_treat=START_AGE, msh3_treat=m_eff, seed=1000 + i)
            y2 = 100 * (1 - float(np.interp(70.0, tt1, ss1)))
            Y[i] = (y1, y2, y1 - y2)
        except Exception as exc:
            print(f"    sample {i} failed: {exc}")
        finally:
            hd.LIT.update(lit_backup)
            hd.ASSUMED.clear()
            hd.ASSUMED.update(assumed_backup)
        if (i + 1) % 40 == 0:
            print(f"    {i+1}/{len(X)}")
    return Y


LABELS = ["untreated SPN loss at 70 (%)",
          "treated SPN loss at 70 (%)",
          "absolute benefit (percentage points)"]

rows = []
benefit_by_rho = {}

for rho in RHOS:
    print(f"\nrho = {rho}")
    X = sample_inputs(rho, seed=20261101)
    # realised rank correlation among the coupled anchors, as a check
    rs = [spearmanr(X[:, a], X[:, b]).statistic
          for ii, a in enumerate(coupled_idx) for b in coupled_idx[ii + 1:]]
    rank_target = (6.0 / np.pi) * np.arcsin(rho / 2.0)
    print(f"  realised Spearman among coupled anchors: "
          f"{np.mean(rs):+.3f} (target {rank_target:+.3f})")

    Y = run_ensemble(X)
    ok = np.all(np.isfinite(Y), axis=1)
    Yo = Y[ok]
    print(f"  {ok.sum()}/{N_SAMPLES} usable")
    benefit_by_rho[rho] = Yo[:, 2]

    for j, lab in enumerate(LABELS):
        q = np.percentile(Yo[:, j], [2.5, 50, 97.5])
        rows.append({
            "latent_rho": rho,
            "realised_spearman_coupled": round(float(np.mean(rs)), 3),
            "outcome": lab,
            "n_usable": int(ok.sum()),
            "median": round(float(q[1]), 2),
            "p2.5": round(float(q[0]), 2),
            "p97.5": round(float(q[2]), 2),
            "central95_width": round(float(q[2] - q[0]), 2),
            "min": round(float(Yo[:, j].min()), 2),
            "benefit_positive_in_all_samples": (
                bool(np.all(Yo[:, 2] > 0)) if j == 2 else ""),
        })
        print(f"    {lab:<42} {q[1]:6.1f}  [{q[0]:.1f}, {q[2]:.1f}]  "
              f"width {q[2]-q[0]:.1f}")

with open(os.path.join(RES, "anchor_correlation.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# --------------------------------------------------------------- Figure 16
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1))

cols = ["#1f4e79", "#4f81bd", "#c0504d"]
for c, rho in zip(cols, RHOS):
    b = benefit_by_rho[rho]
    ax[0].hist(b, bins=20, histtype="step", lw=1.6, color=c,
               label=f"$\\rho$ = {rho}")
ax[0].set_xlabel("striatum preserved at age 70\n(percentage points)")
ax[0].set_ylabel("samples")
ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title("A  benefit, by anchor coupling", fontsize=8.5)

widths = [np.percentile(benefit_by_rho[r], 97.5)
          - np.percentile(benefit_by_rho[r], 2.5) for r in RHOS]
mins = [benefit_by_rho[r].min() for r in RHOS]
ax[1].plot(RHOS, widths, "o-", color="#1f4e79", lw=1.6, label="central 95% width")
ax[1].plot(RHOS, mins, "s--", color="#c0504d", lw=1.4, label="minimum benefit")
ax[1].axhline(0, color="k", lw=0.7)
ax[1].set_xlabel("latent correlation among dwell-time anchors")
ax[1].set_ylabel("percentage points")
ax[1].legend(fontsize=7, frameon=False)
ax[1].set_title("B  spread vs coupling", fontsize=8.5)

fig.savefig(os.path.join(FIG, "fig16_anchor_correlation.png"))
plt.close(fig)
print("\nwrote results/anchor_correlation.csv and figure 16")
