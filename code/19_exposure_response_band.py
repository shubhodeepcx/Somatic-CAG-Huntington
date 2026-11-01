"""19 -- How much do two mouse doses actually constrain the exposure-response
curve, and what does that uncertainty do to the delivery conclusions?

`05_delivery.py` fits an Emax curve with two free parameters, EC50 and a Hill
coefficient, to exactly two published mouse ICV doses. Two parameters on two
points is an exact interpolation with zero degrees of freedom, so the fitted
Hill coefficient of 1.215 is not a measurement of cooperativity: it is whatever
value makes the curve pass through both points. Everything in the delivery and
dosing analysis passes through this curve.

Limitation 5 of the paper states that. This script asks the quantitative
question that statement leaves open: *how badly does the curve miss the two
measured points as the Hill coefficient is moved away from 1.215?* That is the
same profile question `10_profile.py` asks of the rate-law parameters, applied
to the pharmacodynamic curve, and it is answerable from the same two points.

Method
------
For each of a grid of fixed Hill coefficients, EC50 is refitted by least squares
through the two measured points and the worst relative deviation across them is
recorded. An admissible window is read off at a 10% tolerance, matching the
tolerance used for the rate law. The resulting EC50 range is then propagated
into the three numbers the paper quotes downstream:

  * the striatal exposure required for 75% MSH3 knockdown (the headline
    25 nmol/g figure, and the ~7 mg of drug it implies);
  * the exposure required for 50% and 90% knockdown.

The systemic fold shortfall is not recomputed here: the required systemic dose
is strictly proportional to the required striatal exposure (see
`systemic_dose_nmol` in `05_delivery.py`), so the relative band on the exposure
carries across to the shortfall unchanged. Reporting it that way avoids
duplicating the receptor-transit and endosomal-escape constants in a second
place, which is how they went stale once before (audit pass 2, finding 3).

Outputs
-------
results/hill_profile.csv          worst miss and refitted EC50 vs Hill
results/exposure_response_band.csv  downstream quantities across the window
figures/fig17_exposure_response_band.png
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

# Same constants as 05_delivery.py.
MOUSE_BRAIN_G = 0.40
HUMAN_STRIATUM_G = 20.0
ASSUMED_MW_KDA = 14.0

dose_mouse_nmol = np.array([3.0, 10.0])
dose_per_g = dose_mouse_nmol / MOUSE_BRAIN_G
msh3_obs = np.array([hd.LIT["msh3_protein_3nmol_pct"],
                     hd.LIT["msh3_protein_10nmol_pct"]])

TOLERANCE = 0.10   # same tolerance as the rate-law profiles in 10_profile.py


def msh3_from_exposure(x, ec50, hill):
    x = np.asarray(x, float)
    return 100.0 / (1.0 + (x / ec50) ** hill)


def exposure_for_msh3(target_pct, ec50, hill):
    return ec50 * ((100.0 / target_pct) - 1.0) ** (1.0 / hill)


def refit_ec50(hill):
    """Best EC50 at a fixed Hill, and the worst relative miss on the two points."""
    r = least_squares(
        lambda p: msh3_from_exposure(dose_per_g, np.exp(p[0]), hill) - msh3_obs,
        [np.log(10.0)], xtol=1e-14, ftol=1e-14)
    ec50 = float(np.exp(r.x[0]))
    pred = msh3_from_exposure(dose_per_g, ec50, hill)
    worst = float(np.max(np.abs(pred - msh3_obs) / msh3_obs))
    return ec50, worst


# ------------------------------------------------------- free (2-point) fit
free = least_squares(
    lambda p: msh3_from_exposure(dose_per_g, np.exp(p[0]), np.exp(p[1])) - msh3_obs,
    [np.log(10.0), np.log(1.0)], xtol=1e-14, ftol=1e-14)
EC50_FREE, HILL_FREE = float(np.exp(free.x[0])), float(np.exp(free.x[1]))
print(f"free fit (2 params, 2 points): EC50 = {EC50_FREE:.2f} nmol/g, "
      f"Hill = {HILL_FREE:.3f}  [zero degrees of freedom]")

# ------------------------------------------------------- profile over Hill
hills = np.round(np.arange(0.5, 3.01, 0.05), 2)
prof = []
for h in hills:
    ec50, worst = refit_ec50(float(h))
    prof.append({"hill": float(h), "refitted_EC50_nmol_per_g": round(ec50, 3),
                 "worst_relative_deviation": round(worst, 5),
                 "within_10pct_tolerance": bool(worst <= TOLERANCE)})

with open(os.path.join(RES, "hill_profile.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(prof[0].keys()))
    w.writeheader()
    w.writerows(prof)

adm = [r for r in prof if r["within_10pct_tolerance"]]
h_lo, h_hi = adm[0]["hill"], adm[-1]["hill"]
ec_at = {r["hill"]: r["refitted_EC50_nmol_per_g"] for r in prof}
print(f"\nadmissible Hill window at {TOLERANCE:.0%} tolerance: "
      f"{h_lo} to {h_hi}")
print(f"  EC50 across that window: {ec_at[h_lo]:.2f} to {ec_at[h_hi]:.2f} nmol/g")

# ------------------------------------------- propagate to downstream numbers
corners = [("lower Hill bound", h_lo, ec_at[h_lo]),
           ("free 2-point fit", round(HILL_FREE, 3), EC50_FREE),
           ("upper Hill bound", h_hi, ec_at[h_hi])]

rows = []
for label, h, ec in corners:
    for target_pct, kd in [(50.0, 50), (25.0, 75), (10.0, 90)]:
        x_need = exposure_for_msh3(target_pct, ec, h)
        striatal_nmol = x_need * HUMAN_STRIATUM_G
        striatal_mg = striatal_nmol * ASSUMED_MW_KDA / 1000.0
        rows.append({
            "scenario": label,
            "hill": h,
            "EC50_nmol_per_g": round(ec, 3),
            "MSH3_knockdown_pct": kd,
            "required_striatal_exposure_nmol_per_g": round(float(x_need), 2),
            "required_striatal_dose_mg": round(float(striatal_mg), 2),
            "relative_to_point_estimate": "",
        })
        print(f"  {label:<18} Hill {h:<5} {kd}% knockdown -> "
              f"{x_need:7.2f} nmol/g, {striatal_mg:6.2f} mg striatal")

# express each corner relative to the point estimate at the same knockdown
# level: the required systemic dose is proportional to the required striatal
# exposure, so this ratio is also the factor on the reported fold shortfall
point = {r["MSH3_knockdown_pct"]: r["required_striatal_exposure_nmol_per_g"]
         for r in rows if r["scenario"] == "free 2-point fit"}
for r in rows:
    r["relative_to_point_estimate"] = round(
        r["required_striatal_exposure_nmol_per_g"] / point[r["MSH3_knockdown_pct"]], 3)

with open(os.path.join(RES, "exposure_response_band.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("\nrelative band on required exposure (= band on the systemic fold shortfall,"
      "\nsince required systemic dose is proportional to required exposure):")
for q in [50, 75, 90]:
    rel = [r["relative_to_point_estimate"] for r in rows
           if r["MSH3_knockdown_pct"] == q]
    print(f"  {q}% knockdown: x{min(rel):.2f} to x{max(rel):.2f}")
print("  applied to the reported 18-36-fold systemic shortfall at 75% knockdown: "
      f"{18 * min(r['relative_to_point_estimate'] for r in rows if r['MSH3_knockdown_pct'] == 75):.0f}"
      f" to "
      f"{36 * max(r['relative_to_point_estimate'] for r in rows if r['MSH3_knockdown_pct'] == 75):.0f}"
      "-fold")

# headline band for the 75% knockdown case
hl = [r for r in rows if r["MSH3_knockdown_pct"] == 75]
x_lo = min(r["required_striatal_exposure_nmol_per_g"] for r in hl)
x_hi = max(r["required_striatal_exposure_nmol_per_g"] for r in hl)
mg_lo = min(r["required_striatal_dose_mg"] for r in hl)
mg_hi = max(r["required_striatal_dose_mg"] for r in hl)
print(f"\nheadline band, 75% knockdown:")
print(f"  required striatal exposure {x_lo:.1f} to {x_hi:.1f} nmol/g "
      f"(point estimate {[r for r in hl if r['scenario']=='free 2-point fit'][0]['required_striatal_exposure_nmol_per_g']:.1f})")
print(f"  required striatal dose     {mg_lo:.1f} to {mg_hi:.1f} mg "
      f"(point estimate {[r for r in hl if r['scenario']=='free 2-point fit'][0]['required_striatal_dose_mg']:.1f} mg)")

# --------------------------------------------------------------- Figure 17
fig, ax = plt.subplots(1, 3, figsize=(11.2, 3.1))

hs = [r["hill"] for r in prof]
ws = [100 * r["worst_relative_deviation"] for r in prof]
ax[0].plot(hs, ws, color="#1f4e79", lw=1.7)
ax[0].axhline(100 * TOLERANCE, color="#c0504d", ls="--", lw=1.1,
              label=f"{TOLERANCE:.0%} tolerance")
ax[0].axvspan(h_lo, h_hi, color="#4f81bd", alpha=0.17)
ax[0].axvline(HILL_FREE, color="k", ls=":", lw=1.1, label="free fit")
ax[0].set_xlabel("Hill coefficient (fixed)")
ax[0].set_ylabel("worst miss on the two\nmeasured points (%)")
ax[0].set_ylim(-1, 40)
ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title(f"A  the two doses admit\nHill {h_lo}-{h_hi}", fontsize=8.5)

xx = np.logspace(-0.6, 2.9, 400)
for label, h, ec in corners:
    style = "-" if label == "free 2-point fit" else "--"
    ax[1].semilogx(xx, msh3_from_exposure(xx, ec, h), ls=style, lw=1.6,
                   label=f"Hill {h}")
ax[1].axvspan(dose_per_g.min(), dose_per_g.max(), color="#4f81bd", alpha=0.17)
ax[1].scatter(dose_per_g, msh3_obs, s=34, color="#c0504d", zorder=5,
              edgecolor="k", lw=0.4)
ax[1].set_xlabel("tissue exposure (nmol g$^{-1}$)")
ax[1].set_ylabel("MSH3 protein (% of control)")
ax[1].legend(fontsize=7, frameon=False, loc="lower left")
ax[1].set_title("B  curves through the\nsame two points", fontsize=8.5)

kds = [50, 75, 90]
lo = [min(r["required_striatal_exposure_nmol_per_g"]
          for r in rows if r["MSH3_knockdown_pct"] == q) for q in kds]
hi = [max(r["required_striatal_exposure_nmol_per_g"]
          for r in rows if r["MSH3_knockdown_pct"] == q) for q in kds]
pt = [[r["required_striatal_exposure_nmol_per_g"] for r in rows
       if r["MSH3_knockdown_pct"] == q and r["scenario"] == "free 2-point fit"][0]
      for q in kds]
xpos = np.arange(len(kds))
ax[2].vlines(xpos, lo, hi, color="#1f4e79", lw=7, alpha=0.4)
ax[2].plot(xpos, pt, "o", color="#1f4e79", ms=6)
ax[2].set_yscale("log")
ax[2].set_xticks(xpos)
ax[2].set_xticklabels([f"{q}%" for q in kds])
ax[2].set_xlabel("MSH3 knockdown target")
ax[2].set_ylabel("required striatal exposure\n(nmol g$^{-1}$)")
ax[2].set_title("C  downstream band\n(point estimate, and range)", fontsize=8.5)

fig.tight_layout(w_pad=2.2)
fig.savefig(os.path.join(FIG, "fig17_exposure_response_band.png"))
plt.close(fig)
print("\nwrote results/hill_profile.csv, results/exposure_response_band.csv, figure 17")
