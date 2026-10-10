"""10 -- Feasible intervals for the four fitted parameters (profile analysis).

The calibration has four free parameters and four published anchors, so the fit
is exact and no residual variance is left to build a conventional confidence
interval from. Instead we ask a question that can be answered: for each
parameter, over what range can it be held fixed while the remaining three are
refitted and all four anchors are still matched to within a stated tolerance?

Tolerance is the maximum relative deviation allowed on any single anchor. The
anchors themselves are model-derived point estimates reported without
uncertainty, so three tolerances are reported: 10%, 25% and 50%.

Outputs
-------
results/profile_intervals.csv
results/profile_curves.csv
figures/fig11_profile.png
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

ANCHORS = [hd.LIT["t_40_60_years"], hd.LIT["t_60_80_years"],
           hd.LIT["t_80_150_years"], hd.LIT["rate_rise_70_90_min"]]
NAMES = ["lambda40", "beta", "R", "L_sw"]
LABELS = {"lambda40": r"$\lambda(40)$  (events y$^{-1}$)",
          "beta": r"$\beta$  (repeat$^{-1}$)",
          "R": "$R$  (clustered / isolated)",
          "L_sw": r"$L_{sw}$  (repeats)"}
BOUNDS_LO = np.array([-12.0, 1e-4, 0.0, 50.0])
BOUNDS_HI = np.array([5.0, 1.0, np.log(1e5), 140.0])


def predictions(vec):
    p = hd._params_from_vector(vec)
    return np.array([hd.mean_transit_time(40, 60, p),
                     hd.mean_transit_time(60, 80, p),
                     hd.mean_transit_time(80, 150, p),
                     float(hd.lam(90.0, p) / hd.lam(70.0, p))])


def worst_relative_deviation(vec):
    pred = predictions(vec)
    return float(np.max(np.abs(pred - ANCHORS) / np.array(ANCHORS)))


def refit_with_fixed(idx, value, x_start):
    """Hold component `idx` at `value`; refit the other three."""
    free = [j for j in range(4) if j != idx]

    def unpack(z):
        v = np.empty(4)
        v[idx] = value
        v[free] = z
        return v

    def resid(z):
        pred = predictions(unpack(z))
        return np.log(pred / np.array(ANCHORS))

    best = None
    for jitter in (0.0, 0.4, -0.4):
        z0 = np.array(x_start)[free] * (1.0 + jitter)
        z0 = np.clip(z0, BOUNDS_LO[free] + 1e-9, BOUNDS_HI[free] - 1e-9)
        try:
            r = least_squares(resid, z0,
                              bounds=(BOUNDS_LO[free], BOUNDS_HI[free]),
                              xtol=1e-14, ftol=1e-14, gtol=1e-14)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return unpack(best.x) if best is not None else None


par = hd.calibrate()
x_hat = np.array([np.log(par["lambda40"]), par["beta"], np.log(par["R"]),
                  par["L_sw"]])
print("optimum worst deviation:", f"{worst_relative_deviation(x_hat):.2e}")

# Grids are chosen on the natural scale of each parameter.
GRIDS = {
    "lambda40": np.exp(np.linspace(np.log(0.01), np.log(0.5), 25)),
    "beta": np.concatenate([[1e-4], np.linspace(0.005, 0.12, 24)]),
    "R": np.exp(np.linspace(np.log(2.0), np.log(3000.0), 25)),
    "L_sw": np.linspace(55.0, 135.0, 25),
}
TO_VEC = {"lambda40": np.log, "beta": lambda v: v, "R": np.log,
          "L_sw": lambda v: v}

curve_rows, interval_rows = [], []
curves = {}
for idx, name in enumerate(NAMES):
    devs, vals = [], []
    for v in GRIDS[name]:
        vec = refit_with_fixed(idx, float(TO_VEC[name](v)), x_hat)
        d = worst_relative_deviation(vec) if vec is not None else np.inf
        vals.append(float(v))
        devs.append(float(d))
        curve_rows.append({"parameter": name, "value": round(float(v), 6),
                           "worst_relative_deviation": round(float(d), 5)})
    devs = np.array(devs)
    vals = np.array(vals)
    curves[name] = (vals, devs)
    row = {"parameter": name, "best_fit": round(float(
        {"lambda40": par["lambda40"], "beta": par["beta"], "R": par["R"],
         "L_sw": par["L_sw"]}[name]), 5)}
    for tol in (0.10, 0.25, 0.50):
        ok = vals[devs <= tol]
        if ok.size:
            lo, hi = float(ok.min()), float(ok.max())
            # flag intervals that run into the edge of the grid
            edge = (ok.min() <= vals.min() + 1e-12) or (ok.max() >= vals.max() - 1e-12)
            row[f"interval_{int(tol*100)}pct"] = (
                f"{lo:.4g} to {hi:.4g}" + (" (grid-limited)" if edge else ""))
        else:
            row[f"interval_{int(tol*100)}pct"] = "empty"
    interval_rows.append(row)
    print(f"  {name:<10} {row}")


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("profile_intervals.csv", interval_rows, list(interval_rows[0].keys()))
_write("profile_curves.csv", curve_rows, list(curve_rows[0].keys()))

fig, ax = plt.subplots(1, 4, figsize=(11.0, 2.7))
for a, name in zip(ax, NAMES):
    vals, devs = curves[name]
    if name in ("lambda40", "R"):
        a.semilogx(vals, 100 * devs, color="#1f4e79", lw=1.7)
    else:
        a.plot(vals, 100 * devs, color="#1f4e79", lw=1.7)
    for tol, c in ((10, "#4f81bd"), (25, "#8064a2"), (50, "#c0504d")):
        a.axhline(tol, color=c, lw=0.8, ls=":")
    best = {"lambda40": par["lambda40"], "beta": par["beta"], "R": par["R"],
            "L_sw": par["L_sw"]}[name]
    a.axvline(best, color="k", lw=0.9, ls="--")
    a.set_ylim(0, 80)
    a.set_xlabel(LABELS[name], fontsize=8)
    a.set_ylabel("worst anchor error (%)" if name == "lambda40" else "")
ax[0].set_title("A", loc="left")
ax[1].set_title("B", loc="left")
ax[2].set_title("C", loc="left")
ax[3].set_title("D", loc="left")
fig.suptitle("Feasible ranges of the fitted parameters "
             "(dotted lines: 10, 25 and 50% tolerance)", fontsize=9, y=1.04)
fig.savefig(os.path.join(FIG, "fig11_profile.png"))
plt.close(fig)
print("\nwrote figure 11 and result tables")
