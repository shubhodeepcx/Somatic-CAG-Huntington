"""05 -- Delivery layer: what striatal MSH3 lowering is actually achievable?

Two routes are modelled and compared on a common axis (nmol of oligonucleotide
per gram of striatal tissue):

  Route A  CSF-administered oligonucleotide (intrathecal / intracerebroventricular)
  Route B  systemically administered receptor-mediated-transcytosis (RMT)
           nanocarrier

Anchors
-------
*  Dose -> MSH3 protein: an Emax curve fitted to the two published mouse ICV
   di-siRNA doses (3 nmol -> 59% of control; 10 nmol -> 25% of control).
*  Route A regional penalty: measured intrathecal ASO activity in macaque deep
   brain (putamen SPN residual target RNA 0.71; caudate bulk 0.45-0.57), versus
   the much deeper knockdown achieved in cortex.
*  Route B efficiency: parenchymal accumulation of TfR-targeted nanoparticles
   after intravenous dosing (0.04-0.23 %ID/g, valency- and affinity-dependent).

Outputs
-------
results/emax_fit.csv
results/route_comparison.csv
results/rmt_design_space.csv
figures/fig6_delivery.png
figures/fig7_rmt_design.png
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

# ---------------------------------------------------------------------------
# Anatomical / physiological constants used for scaling, with sources or an
# explicit "assumed" flag. Assumed values are swept in script 07.
# ---------------------------------------------------------------------------
SCALE = {
    "mouse_brain_g": 0.40,        # assumed; standard adult C57BL/6 brain mass
    "human_brain_g": 1300.0,      # assumed; adult human brain mass
    "human_striatum_g": 20.0,     # assumed; combined caudate + putamen, one side ~10 g
    "human_plasma_volume_L": 3.0, # assumed
}

# ------------------------------------------------- Emax fit: dose -> MSH3
# Published mouse ICV di-siRNA doses, expressed per gram of brain
dose_mouse_nmol = np.array([3.0, 10.0])
dose_per_g = dose_mouse_nmol / SCALE["mouse_brain_g"]          # nmol/g brain
msh3_obs = np.array([hd.LIT["msh3_protein_3nmol_pct"],
                     hd.LIT["msh3_protein_10nmol_pct"]])        # % of control


def msh3_from_exposure(x, ec50, hill):
    """Striatal MSH3 protein (% of control) at tissue exposure x (nmol/g)."""
    x = np.asarray(x, float)
    return 100.0 / (1.0 + (x / ec50) ** hill)


fit = least_squares(
    lambda p: msh3_from_exposure(dose_per_g, np.exp(p[0]), np.exp(p[1])) - msh3_obs,
    [np.log(10.0), np.log(1.0)], xtol=1e-14, ftol=1e-14)
EC50, HILL = float(np.exp(fit.x[0])), float(np.exp(fit.x[1]))
print(f"Emax fit: EC50 = {EC50:.2f} nmol/g brain, Hill = {HILL:.2f}")

emax_rows = [{"dose_nmol_ICV_mouse": d, "exposure_nmol_per_g": round(x, 2),
              "measured_MSH3_pct": o,
              "fitted_MSH3_pct": round(float(msh3_from_exposure(x, EC50, HILL)), 1)}
             for d, x, o in zip(dose_mouse_nmol, dose_per_g, msh3_obs)]

# exposure required for a target MSH3 level
def exposure_for_msh3(target_pct):
    return EC50 * ((100.0 / target_pct) - 1.0) ** (1.0 / HILL)


# ------------------------------------------------- Route A: CSF delivery
# The macaque study gives regional *activity*; converting it to an equivalent
# tissue exposure through the Emax curve gives the striatal exposure that an
# intrathecal regimen actually achieves, without needing a CSF transport model.
nhp_residual = {"putamen (SPNs)": hd.LIT["nhp_putamen_spn_residual"],
                "putamen (bulk)": hd.LIT["nhp_putamen_residual_bulk"],
                "caudate (snRNA-seq)": hd.LIT["nhp_caudate_residual_snrna"],
                "caudate (qPCR)": hd.LIT["nhp_caudate_residual_qpcr"]}

routeA = []
for region, resid in nhp_residual.items():
    achieved_msh3 = 100.0 * resid
    x_eq = exposure_for_msh3(achieved_msh3)
    routeA.append({"route": "A: intrathecal oligonucleotide",
                   "region": region,
                   "achieved_knockdown_pct": round(100 * (1 - resid), 1),
                   "implied_tissue_exposure_nmol_per_g": round(float(x_eq), 2),
                   "predicted_expansion_suppression_pct": round(
                       float(hd.prevention_from_msh3(achieved_msh3)), 1)})

# Human penetration penalty: the observed >10-fold cortex-to-putamen gradient
# implies a tissue penetration length. With an exponential depth profile
# c(d) = c0 * exp(-d / Lpen), a gradient ratio r over a depth difference dd
# gives Lpen = dd / ln(r).
DEPTH_CORTEX_MM = 1.0      # assumed: cortical ribbon is adjacent to CSF
DEPTH_PUTAMEN_MM = 10.0    # assumed: 8-12 mm from the nearest CSF surface
GRADIENT = 10.0            # "more than an order of magnitude" lower in deep brain
Lpen_mm = (DEPTH_PUTAMEN_MM - DEPTH_CORTEX_MM) / np.log(GRADIENT)
print(f"\nimplied tissue penetration length Lpen = {Lpen_mm:.2f} mm")

# For the human brain the same penetration length must cover a greater depth.
# Depths are assumed to scale with the cube root of brain mass (macaque ~ 90 g,
# human ~ 1300 g), so the human putamen sits proportionally further from CSF.
MACAQUE_BRAIN_G = 90.0   # assumed
depth_human_mm = DEPTH_PUTAMEN_MM * (SCALE["human_brain_g"] / MACAQUE_BRAIN_G) ** (1 / 3)
extra_attenuation = float(np.exp(-(depth_human_mm - DEPTH_PUTAMEN_MM) / Lpen_mm))
print(f"macaque->human putamen depth {DEPTH_PUTAMEN_MM:.1f} -> "
      f"{depth_human_mm:.1f} mm; extra attenuation x{extra_attenuation:.3g}")

# ------------------------------------------------- Route B: RMT nanocarrier
# Transcytosis efficiency as a function of anti-TfR affinity and valency:
#   eta(KD, v) = f_v * [KD_b / (KD + KD_b)] * [1 / (1 + K_r / KD)]
# term 1 = surface capture (favours high affinity),
# term 2 = abluminal release (favours low affinity) -> an optimum in between.
AU = hd.LIT["aunp_parenchyma_pctID_g"]
f_biv = AU["bivalent_highaff"] / AU["monovalent"]       # valency factor at matched KD
KD_pts = np.array([21.0, 149.0])                        # nM, bivalent variants
eta_pts = np.array([AU["bivalent_highaff"], AU["bivalent_lowaff"]])


def eta(KD, f_v, KD_b, K_r):
    KD = np.asarray(KD, float)
    return f_v * (KD_b / (KD + KD_b)) * (1.0 / (1.0 + K_r / KD))


fitB = least_squares(
    lambda p: eta(KD_pts, 1.0, np.exp(p[0]), np.exp(p[1])) - eta_pts / f_biv,
    [np.log(300.0), np.log(100.0)], xtol=1e-14, ftol=1e-14)
KD_B, K_R = float(np.exp(fitB.x[0])), float(np.exp(fitB.x[1]))
KD_opt = float(np.sqrt(KD_B * K_R))   # maximum of the product form
eta_opt_mono = float(eta(KD_opt, 1.0, KD_B, K_R))
eta_mono_obs = AU["monovalent"]
print(f"\nRMT design: KD_b = {KD_B:.1f} nM, K_r = {K_R:.1f} nM, "
      f"optimal KD = {KD_opt:.1f} nM")
print(f"monovalent parenchymal %ID/g at optimum = {eta_opt_mono:.3f} "
      f"(measured at KD 22 nM: {eta_mono_obs:.3f}); "
      f"headroom x{eta_opt_mono/eta_mono_obs:.2f}")

design_rows = []
for KD in [1, 3, 10, 21, 30, 50, 100, 149, 300, 1000, 3000]:
    for vlabel, fv in (("monovalent", 1.0), ("bivalent", f_biv)):
        design_rows.append({"KD_nM": KD, "valency": vlabel,
                            "parenchymal_pctID_per_g": round(
                                float(eta(KD, fv, KD_B, K_R)), 4)})

# Systemic dose required to match a target striatal exposure
def systemic_dose_nmol(target_exposure_nmol_per_g, pct_id_per_g,
                       endosomal_escape=0.02):
    """Injected dose (nmol) needed to reach `target_exposure_nmol_per_g` of
    *cytosolic* oligonucleotide in striatal tissue, given parenchymal
    accumulation `pct_id_per_g` (%ID/g) and an endosomal escape fraction.
    """
    frac_per_g = pct_id_per_g / 100.0
    return target_exposure_nmol_per_g / (frac_per_g * endosomal_escape)


targets = [50.0, 25.0, 10.0]   # target MSH3 protein, % of control
routeB = []
for tgt in targets:
    x_need = exposure_for_msh3(tgt)
    for vlabel, fv, KD in (("monovalent, KD 22 nM", 1.0, 22.0),
                           ("monovalent, KD at optimum", 1.0, KD_opt),
                           ("bivalent, KD 21 nM", f_biv, 21.0)):
        e = float(eta(KD, fv, KD_B, K_R))
        D = systemic_dose_nmol(x_need, e)
        routeB.append({"route": "B: systemic RMT nanocarrier",
                       "design": vlabel,
                       "target_MSH3_pct": tgt,
                       "required_striatal_exposure_nmol_per_g": round(x_need, 2),
                       "parenchymal_pctID_per_g": round(e, 4),
                       "required_systemic_dose_nmol": round(D, 1),
                       "fold_vs_direct_CSF_dose": round(
                           D / (x_need * SCALE["human_striatum_g"]), 1)})

for tgt in targets:
    x_need = exposure_for_msh3(tgt)
    routeA.append({"route": "A: intrathecal oligonucleotide",
                   "region": f"target MSH3 {tgt:.0f}% -> required exposure",
                   "achieved_knockdown_pct": round(100 - tgt, 1),
                   "implied_tissue_exposure_nmol_per_g": round(float(x_need), 2),
                   "predicted_expansion_suppression_pct": round(
                       float(hd.prevention_from_msh3(tgt)), 1)})


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("emax_fit.csv", emax_rows, list(emax_rows[0].keys()))
_write("route_comparison.csv", routeA + routeB,
       ["route", "region", "design", "target_MSH3_pct",
        "achieved_knockdown_pct", "implied_tissue_exposure_nmol_per_g",
        "required_striatal_exposure_nmol_per_g", "parenchymal_pctID_per_g",
        "required_systemic_dose_nmol", "fold_vs_direct_CSF_dose",
        "predicted_expansion_suppression_pct"])
_write("rmt_design_space.csv", design_rows, list(design_rows[0].keys()))
with open(os.path.join(RES, "delivery_constants.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["constant", "value", "role"])
    for k, v in SCALE.items():
        w.writerow([k, v, "assumed scaling constant"])
    w.writerow(["EC50_nmol_per_g", round(EC50, 3), "fitted to mouse ICV dose-response"])
    w.writerow(["Hill", round(HILL, 3), "fitted to mouse ICV dose-response"])
    w.writerow(["Lpen_mm", round(Lpen_mm, 3), "inferred from NHP regional gradient"])
    w.writerow(["KD_b_nM", round(KD_B, 2), "fitted to AuNP affinity series"])
    w.writerow(["K_r_nM", round(K_R, 2), "fitted to AuNP affinity series"])
    w.writerow(["KD_optimal_nM", round(KD_opt, 2), "derived"])
    w.writerow(["valency_factor_bivalent", round(f_biv, 3), "from AuNP data"])

print("\nroute A (what intrathecal delivery actually achieves)")
for r in routeA[:4]:
    print("  ", r)
print("\nroute B (systemic RMT dose requirement)")
for r in routeB:
    print("  ", r)

# -------------------------------------------------------------------- Fig. 6
fig, ax = plt.subplots(1, 3, figsize=(10.4, 2.95))

xx = np.logspace(-1, 3, 300)
# Two parameters fitted to two points: the curve is an exact interpolation
# inside this span and an extrapolation outside it. Shade the span so that
# every exposure quoted elsewhere can be read against it.
ax[0].axvspan(float(dose_per_g.min()), float(dose_per_g.max()),
              color="#c6d9f1", alpha=0.55, lw=0, zorder=0)
ax[0].annotate("calibrated\nrange", xy=(float(np.sqrt(dose_per_g.min() * dose_per_g.max())), 92),
               ha="center", va="top", fontsize=6.0, color="#1f4e79")
ax[0].semilogx(xx, msh3_from_exposure(xx, EC50, HILL), color="#1f4e79", lw=1.9)
ax[0].scatter(dose_per_g, msh3_obs, s=42, color="#e36c0a", zorder=5,
              edgecolor="k", lw=0.4, label="mouse ICV (3, 10 nmol)")
for tgt, c in zip(targets, ("#8064a2", "#4f81bd", "#c0504d")):
    ax[0].axhline(tgt, color=c, lw=0.8, ls=":")
    ax[0].annotate(f"{100-tgt:.0f}% KD", xy=(0.12, tgt + 2), fontsize=6.2, color=c)
ax[0].set_xlabel("tissue exposure (nmol g$^{-1}$)")
ax[0].set_ylabel("MSH3 protein (% of control)")
ax[0].set_title("A  exposure-response (2 points, 2 parameters)")
ax[0].legend(frameon=False, fontsize=6.5)

regions = list(nhp_residual.keys())
kd_vals = [100 * (1 - nhp_residual[r]) for r in regions]
supp = [float(hd.prevention_from_msh3(100 * nhp_residual[r])) for r in regions]
y = np.arange(len(regions))
ax[1].barh(y - 0.19, kd_vals, height=0.36, color="#4f81bd",
           label="target knockdown measured in primate")
ax[1].barh(y + 0.19, supp, height=0.36, color="#c0504d",
           label="predicted expansion suppression")
ax[1].set_yticks(y)
ax[1].set_yticklabels(regions, fontsize=6.8)
ax[1].set_xlabel("%")
ax[1].set_xlim(0, 100)
ax[1].set_title("B  intrathecal route, deep brain")
ax[1].legend(frameon=False, fontsize=6.3, loc="lower right")

labels = [r["design"] for r in routeB if r["target_MSH3_pct"] == 25.0]
doses = [r["required_systemic_dose_nmol"] for r in routeB if r["target_MSH3_pct"] == 25.0]
direct = [r["required_striatal_exposure_nmol_per_g"] * SCALE["human_striatum_g"]
          for r in routeB if r["target_MSH3_pct"] == 25.0][0]
ax[2].barh(range(len(labels)), doses, color="#8064a2", height=0.46)
ax[2].axvline(direct, color="k", ls="--", lw=1.1)
ax[2].set_xscale("log")
# Long design names are written above each bar rather than as tick labels,
# which would otherwise run left out of the axes and over panel B.
ax[2].set_yticks(range(len(labels)))
ax[2].set_yticklabels([])
ax[2].tick_params(axis="y", length=0)
for k, lab in enumerate(labels):
    ax[2].annotate(lab.replace("KD", "$K_D$"), xy=(direct * 1.25, k + 0.29),
                   fontsize=6.3, ha="left", va="bottom")
ax[2].annotate("dose delivered directly to the striatum",
               xy=(direct * 1.25, -0.60), fontsize=6.2, ha="left", va="center")
ax[2].set_ylim(-1.0, len(labels) - 0.35)
ax[2].set_xlabel("systemic dose required (nmol)")
ax[2].set_title("C  75% MSH3 lowering, systemic route")
fig.subplots_adjust(wspace=0.30)
fig.savefig(os.path.join(FIG, "fig6_delivery.png"))
plt.close(fig)

# -------------------------------------------------------------------- Fig. 7
kds = np.logspace(-0.5, 3.7, 400)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
ax[0].semilogx(kds, eta(kds, 1.0, KD_B, K_R), color="#1f4e79", lw=1.9,
               label="monovalent")
ax[0].semilogx(kds, eta(kds, f_biv, KD_B, K_R), color="#c0504d", lw=1.6, ls="--",
               label="bivalent")
ax[0].scatter([21, 149], [AU["bivalent_highaff"], AU["bivalent_lowaff"]], s=36,
              color="#c0504d", edgecolor="k", lw=0.4, zorder=5)
ax[0].scatter([22], [AU["monovalent"]], s=36, color="#1f4e79", edgecolor="k",
              lw=0.4, zorder=5)
ax[0].axvline(KD_opt, color="0.4", ls=":", lw=1.0)
ax[0].annotate(f"optimum\n$K_D \\approx$ {KD_opt:.0f} nM", xy=(KD_opt, 0.02),
               fontsize=6.5, ha="center")
ax[0].set_xlabel("anti-TfR $K_D$ (nM)")
ax[0].set_ylabel("parenchymal accumulation (%ID g$^{-1}$)")
ax[0].set_title("A  affinity-valency optimum", fontsize=9)
ax[0].legend(frameon=False, fontsize=7)

eff = np.logspace(np.log10(AU["bivalent_highaff"]), 2, 200)
x25 = exposure_for_msh3(25.0)
for esc, c, lab in ((0.005, "#c0504d", "2% of 2% escape"),
                    (0.02, "#1f4e79", "2% endosomal escape"),
                    (0.20, "#4f6228", "20% endosomal escape")):
    ax[1].loglog(eff, [systemic_dose_nmol(x25, e, esc) for e in eff], lw=1.6,
                 color=c, label=lab)
ax[1].axvline(eta_opt_mono, color="0.4", ls=":", lw=1.0)
ax[1].axhline(x25 * SCALE["human_striatum_g"], color="k", ls="--", lw=1.0)
ax[1].annotate("direct striatal dose", xy=(eff[0], x25 * SCALE["human_striatum_g"] * 1.4),
               fontsize=6.3)
ax[1].annotate("best measured\nmonovalent design", xy=(eta_opt_mono, 1e8),
               fontsize=6.3, ha="center")
ax[1].set_xlabel("parenchymal accumulation (%ID g$^{-1}$)")
ax[1].set_ylabel("systemic dose for 75% MSH3 lowering (nmol)")
ax[1].set_title("B  dose penalty of the systemic route", fontsize=9)
ax[1].legend(frameon=False, fontsize=6.5)
fig.savefig(os.path.join(FIG, "fig7_rmt_design.png"))
plt.close(fig)
print("\nwrote figures 6-7 and result tables")
