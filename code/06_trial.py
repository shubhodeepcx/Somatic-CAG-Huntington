"""06 -- In silico trial: how much striatal neuron loss does MSH3 lowering
actually prevent, and when does treatment have to start?

Design
------
Virtual cohorts of SPNs carrying an inherited CAG length, treated at a given
age with a step change in striatal MSH3 protein that is held for life. The
MSH3 level is mapped onto the expansion rate first-order (script 01), so a
treatment that leaves a fraction m of MSH3 protein multiplies lambda(L) by m.

Scenarios are anchored to what the delivery layer (script 05) says is
achievable: 71% residual MSH3 in putamen SPNs and 45-57% in caudate under a
supraclinical intrathecal regimen, versus aspirational 25% and 10% residual.

Both bracketing kinetic parameterisations from script 04 are run.

Outputs
-------
results/trial_grid.csv
results/trial_scenarios.csv
results/loss_vs_cap.csv
results/parent_aao_sensitivity.csv
figures/fig8_trial.png
figures/fig9_start_age.png
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
with open(os.path.join(RES, "g_fit.json")) as f:
    G_FIT = json.load(f)["g_fit"]

DT = 1 / 24
N = 9000
MODELS = {"M1": 1.0, "M2": G_FIT}

# residual MSH3 protein levels (fraction of control)
RESIDUALS = [1.00, 0.75, 0.71, 0.57, 0.45, 0.25, 0.10, 0.05]
START_AGES = [20, 30, 40, 50]
LABELS = {0.71: "putamen SPNs, intrathecal (measured)",
          0.57: "caudate, intrathecal (snRNA-seq)",
          0.45: "caudate, intrathecal (qPCR)",
          0.25: "aspirational 75% lowering",
          0.10: "aspirational 90% lowering",
          0.05: "near-complete lowering"}


def _fmt_age(a, t_max):
    return f">{t_max:.0f}" if not np.isfinite(a) else round(a, 1)


def _fmt_delay(a, a_base, t_max):
    """Years gained. If the treated arm is censored, report a lower bound."""
    if not np.isfinite(a_base):
        return "n/a"
    if not np.isfinite(a):
        return f">{t_max - a_base:.1f}"
    return round(a - a_base, 1)


def _delay_numeric(a, a_base, t_max):
    """Numeric version for plotting: censored values become the lower bound."""
    if not np.isfinite(a_base):
        return np.nan
    return (t_max - a_base) if not np.isfinite(a) else (a - a_base)


def run(cag, g, residual, t_treat, n=N, seed=0):
    p = dict(par0)
    p["A"] = par0["A"] * g
    arr, tt, ss = hd.simulate_population_tau(
        float(cag), p, n_cells=n, dt=DT, t_treat=float(t_treat),
        msh3_treat=float(residual), seed=seed)
    # Right-censoring matters: under strong treatment the cohort never reaches
    # 25% or 50% loss within the simulated lifespan, and interpolation would
    # silently return the final time point. Censored endpoints are returned as
    # NaN and reported as "> t_max".
    def age_at(level):
        if ss.min() > level:
            return np.nan
        return float(np.interp(level, ss[::-1], tt[::-1]))

    return {"age_25pct_loss": age_at(0.75), "age_50pct_loss": age_at(0.50),
            "loss_at_age_70": 1.0 - float(np.interp(70.0, tt, ss)),
            "loss_at_age_80": 1.0 - float(np.interp(80.0, tt, ss)),
            "t_max": float(tt[-1]), "times": tt, "surviving": ss}


# ------------------------------------------------------------- full grid (M1)
grid_rows = []
curves = {}
base = {}
for mname, g in MODELS.items():
    b = run(42, g, 1.0, np.inf, seed=1)
    base[mname] = b
    curves[(mname, "untreated")] = (b["times"], b["surviving"])

seed = 100
for mname, g in MODELS.items():
    ages = START_AGES if mname == "M1" else [20, 40]
    resids = RESIDUALS if mname == "M1" else [1.0, 0.71, 0.25, 0.05]
    for t0 in ages:
        for r in resids:
            seed += 1
            out = run(42, g, r, t0, seed=seed)
            grid_rows.append({
                "model": mname, "inherited_CAG": 42,
                "treatment_start_age": t0,
                "residual_MSH3_fraction": r,
                "MSH3_knockdown_pct": round(100 * (1 - r), 1),
                "expansion_suppression_pct": round(
                    float(hd.prevention_from_msh3(100 * r)), 1),
                "age_25pct_SPN_loss": _fmt_age(out["age_25pct_loss"], out["t_max"]),
                "age_50pct_SPN_loss": _fmt_age(out["age_50pct_loss"], out["t_max"]),
                "delay_25pct_years": _fmt_delay(out["age_25pct_loss"],
                                                base[mname]["age_25pct_loss"],
                                                out["t_max"]),
                "delay_50pct_years": _fmt_delay(out["age_50pct_loss"],
                                                base[mname]["age_50pct_loss"],
                                                out["t_max"]),
                "SPN_loss_at_age_70_pct": round(100 * out["loss_at_age_70"], 1),
                "SPN_loss_at_70_untreated_pct": round(
                    100 * base[mname]["loss_at_age_70"], 1),
                "SPN_loss_at_age_80_pct": round(100 * out["loss_at_age_80"], 1),
                "SPN_loss_at_80_untreated_pct": round(
                    100 * base[mname]["loss_at_age_80"], 1),
                "delay_25pct_numeric": round(_delay_numeric(
                    out["age_25pct_loss"], base[mname]["age_25pct_loss"],
                    out["t_max"]), 2),
                "delay_25pct_censored": not np.isfinite(out["age_25pct_loss"])})
            if mname == "M1" and t0 == 30:
                curves[("M1", r)] = (out["times"], out["surviving"])

# --------------------------------------------- CAG sensitivity of the headline
scen_rows = []
for cag in (40, 42, 45, 48):
    b = run(cag, 1.0, 1.0, np.inf, seed=500 + cag)
    for r in (0.71, 0.45, 0.25, 0.05):
        out = run(cag, 1.0, r, 30, seed=600 + cag + int(100 * r))
        scen_rows.append({
            "model": "M1", "inherited_CAG": cag, "treatment_start_age": 30,
            "scenario": LABELS.get(r, f"residual {r}"),
            "residual_MSH3_fraction": r,
            "delay_25pct_years": _fmt_delay(out["age_25pct_loss"],
                                            b["age_25pct_loss"], out["t_max"]),
            "delay_25pct_numeric": round(_delay_numeric(
                out["age_25pct_loss"], b["age_25pct_loss"], out["t_max"]), 2),
            "delay_25pct_censored": not np.isfinite(out["age_25pct_loss"]),
            "SPN_loss_at_70_treated_pct": round(100 * out["loss_at_age_70"], 1),
            "SPN_loss_at_70_untreated_pct": round(100 * b["loss_at_age_70"], 1),
            "age_25pct_untreated": _fmt_age(b["age_25pct_loss"], b["t_max"]),
            "age_25pct_treated": _fmt_age(out["age_25pct_loss"], out["t_max"])})

# ---------------------------------- loss-vs-CAP curve for the untreated model
cap_rows = []
for mname, g in MODELS.items():
    for cag in (42, 46, 50):
        p = dict(par0); p["A"] = par0["A"] * g
        _, tt, ss = hd.simulate_population_tau(float(cag), p, n_cells=N, dt=DT,
                                               seed=3000 + cag)
        for cap in (300, 400, 500, 600, 700, 800, 900, 1000):
            age = cap / (cag - hd.LIT["cap_offset"])
            if age > 100:
                continue
            cap_rows.append({"model": mname, "inherited_CAG": cag, "CAP": cap,
                             "age_years": round(age, 1),
                             "SPN_loss_pct": round(
                                 100 * (1 - float(np.interp(age, tt, ss))), 1)})
        # SPN loss at the externally predicted clinical onset age
        onset = hd.ranen_predicted_onset(cag)
        cap_rows.append({"model": mname, "inherited_CAG": cag,
                         "CAP": "at Ranen onset", "age_years": round(onset, 1),
                         "SPN_loss_pct": round(
                             100 * (1 - float(np.interp(onset, tt, ss))), 1)})


# ------------------- sensitivity to the assumed parental age at onset
# The external onset regression takes a parental age at onset, which the
# source does not report as a cohort mean; 45 years is assumed here. Both
# conclusions that use the regression (the rejection of M2 and the predicted
# neuron loss at clinical onset) therefore need it swept.
# The survival curves do not depend on the assumed parental onset age, so
# simulate each model and length once and reuse them across the sweep.
aao_curves = {}
for mname, g in MODELS.items():
    for cag in (42, 46, 50):
        p = dict(par0); p["A"] = par0["A"] * g
        aao_curves[(mname, cag)] = hd.simulate_population_tau(
            float(cag), p, n_cells=N, dt=DT, seed=3000 + cag)[1:]

aao_rows = []
for paao in (35.0, 40.0, 45.0, 50.0, 55.0):
    row = {"assumed_parental_onset_age_y": paao}
    for mname, g in MODELS.items():
        tag = "dwell_faithful" if mname.startswith("M1") else "loss_faithful"
        for cag in (42, 46, 50):
            tt, ss = aao_curves[(mname, cag)]
            onset = hd.ranen_predicted_onset(cag, parent_aao=paao)
            row[f"onset_CAG{cag}_y"] = round(onset, 1)
            row[f"loss_at_onset_CAG{cag}_{tag}_pct"] = round(
                100 * (1 - float(np.interp(onset, tt, ss))), 1)
    aao_rows.append(row)

print("\nsensitivity to the assumed parental age at onset")
for r in aao_rows:
    print("  ", r)


def _write(name, rws, fields):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rws)


_write("trial_grid.csv", grid_rows, list(grid_rows[0].keys()))
_write("trial_scenarios.csv", scen_rows, list(scen_rows[0].keys()))
_write("loss_vs_cap.csv", cap_rows, list(cap_rows[0].keys()))
_write("parent_aao_sensitivity.csv", aao_rows, list(aao_rows[0].keys()))

print("headline grid (model M1, CAG 42)")
for r in grid_rows:
    if r["model"] == "M1" and r["treatment_start_age"] in (20, 30, 40):
        print(f"  start {r['treatment_start_age']}  KD {r['MSH3_knockdown_pct']:>5}%"
              f"  loss@70 {r['SPN_loss_at_age_70_pct']:>5}%"
              f"  (untreated {r['SPN_loss_at_70_untreated_pct']}%)"
              f"  delay(25% loss) {str(r['delay_25pct_years']):>7} y")
print("\nSPN loss at externally predicted clinical onset")
for r in cap_rows:
    if r["CAP"] == "at Ranen onset":
        print("  ", r)

# -------------------------------------------------------------------- Fig. 8
fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.0))

order = [1.00, 0.71, 0.45, 0.25, 0.05]
cmap = plt.get_cmap("viridis")
for i, r in enumerate(order):
    key = ("M1", "untreated") if r == 1.0 else ("M1", r)
    if key not in curves:
        continue
    tt, ss = curves[key]
    lab = "untreated" if r == 1.0 else f"{100*(1-r):.0f}% MSH3 lowering"
    ax[0].plot(tt, 100 * ss, color=cmap(i / (len(order) - 1)), lw=1.7, label=lab)
ax[0].axvline(30, color="0.5", ls=":", lw=0.9)
ax[0].annotate("treatment\nstarts", xy=(32, 72), fontsize=6.3)
ax[0].set_xlabel("age (years)")
ax[0].set_ylabel("surviving SPNs (%)")
ax[0].set_xlim(0, 100)
ax[0].set_title("A  42 CAG, treated from age 30")
ax[0].legend(frameon=False, fontsize=6.3, loc="lower left")

m1 = [r for r in grid_rows if r["model"] == "M1"]
for t0, c in zip(START_AGES, ("#1f4e79", "#4f81bd", "#8064a2", "#c0504d")):
    sub = sorted([r for r in m1 if r["treatment_start_age"] == t0],
                 key=lambda r: r["MSH3_knockdown_pct"])
    ax[1].plot([r["MSH3_knockdown_pct"] for r in sub],
               [r["delay_25pct_numeric"] for r in sub], "o-", color=c, lw=1.5,
               ms=3.4, label=f"from age {t0}")
    cens = [r for r in sub if r["delay_25pct_censored"]]
    if cens:
        ax[1].scatter([r["MSH3_knockdown_pct"] for r in cens],
                      [r["delay_25pct_numeric"] for r in cens], marker="^",
                      s=30, color=c, zorder=5)
ax[1].axvspan(25, 55, color="0.88", lw=0, zorder=0)
ax[1].annotate("measured intrathecal\nrange (striatum)", xy=(40, 1.0),
               fontsize=6.0, ha="center")
ax[1].set_xlabel("striatal MSH3 knockdown (%)")
ax[1].set_ylabel("delay to 25% SPN loss (years)")
ax[1].set_title("B  dose and timing\n(triangles: censored, lower bound)", fontsize=8.5)
ax[1].legend(frameon=False, fontsize=6.3)

sub = [r for r in scen_rows if r["residual_MSH3_fraction"] in (0.71, 0.25, 0.05)]
cags = sorted({r["inherited_CAG"] for r in sub})
w = 0.26
for k, (res, c, lab) in enumerate(((0.71, "#c0504d", "29% KD (measured)"),
                                   (0.25, "#4f81bd", "75% KD"),
                                   (0.05, "#1f4e79", "95% KD"))):
    vals = [next(r["delay_25pct_numeric"] for r in sub
                 if r["inherited_CAG"] == c0 and r["residual_MSH3_fraction"] == res)
            for c0 in cags]
    cens = [next(r["delay_25pct_censored"] for r in sub
                 if r["inherited_CAG"] == c0 and r["residual_MSH3_fraction"] == res)
            for c0 in cags]
    ax[2].bar(np.arange(len(cags)) + (k - 1) * w, vals, width=w, color=c,
              label=lab, hatch=["//" if cc else "" for cc in cens],
              edgecolor="white", linewidth=0.4)
ax[2].set_xticks(range(len(cags)))
ax[2].set_xticklabels(cags)
ax[2].set_xlabel("inherited CAG length")
ax[2].set_ylabel("delay to 25% SPN loss (years)")
ax[2].set_title("C  treated from age 30\n(hatched: censored lower bound)",
                fontsize=8.5)
ax[2].legend(frameon=False, fontsize=6.3)
fig.savefig(os.path.join(FIG, "fig8_trial.png"))
plt.close(fig)

# -------------------------------------------------------------------- Fig. 9
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
for res, c in zip((0.71, 0.45, 0.25, 0.05),
                  ("#c0504d", "#8064a2", "#4f81bd", "#1f4e79")):
    sub = sorted([r for r in m1 if r["residual_MSH3_fraction"] == res],
                 key=lambda r: r["treatment_start_age"])
    if not sub:
        continue
    ax[0].plot([r["treatment_start_age"] for r in sub],
               [r["delay_25pct_numeric"] for r in sub], "o-", color=c, lw=1.6,
               ms=4, label=f"{100*(1-res):.0f}% MSH3 lowering")
ax[0].set_xlabel("age at start of treatment (years)")
ax[0].set_ylabel("delay to 25% SPN loss (years)")
ax[0].set_title("A  value of starting early (42 CAG)", fontsize=9)
ax[0].legend(frameon=False, fontsize=6.5)

for mname, c in (("M1", "#1f4e79"), ("M2", "#c0504d")):
    for cag, ls in ((42, "-"), (46, "--"), (50, ":")):
        pts = [(r["CAP"], r["SPN_loss_pct"]) for r in cap_rows
               if r["model"] == mname and r["inherited_CAG"] == cag
               and isinstance(r["CAP"], int)]
        if pts:
            ax[1].plot([p[0] for p in pts], [p[1] for p in pts], ls, color=c,
                       lw=1.4, label=f"{mname}, {cag} CAG")
ax[1].axhspan(80, 99, color="0.85", lw=0, zorder=0)
ax[1].axvline(600, color="k", ls=":", lw=0.9)
ax[1].annotate("reported: >80% loss\nfor CAP > 600", xy=(620, 55), fontsize=6.2)
ax[1].set_xlabel("CAP score")
ax[1].set_ylabel("SPN loss (%)")
ax[1].set_title("B  model vs the CAP-loss constraint", fontsize=9)
ax[1].legend(frameon=False, fontsize=6.0, ncol=2)
fig.savefig(os.path.join(FIG, "fig9_start_age.png"))
plt.close(fig)
print("\nwrote figures 8-9 and result tables")
