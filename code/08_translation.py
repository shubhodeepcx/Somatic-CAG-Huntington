"""08 -- Translating the delivery requirement into a dose a patient could receive.

Converts the systemic doses computed in script 05 into mg/kg, compares them
with the dose range that systemically administered oligonucleotide conjugates
actually occupy, and reports the fold improvement in parenchymal delivery that
would be needed to close the gap. Molecular weight and the feasible dose
ceiling are both swept, because neither is a property of this model.

Outputs
-------
results/translation.csv
"""
import csv
import os

import numpy as np

import hdmodel as hd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")


def _fitted_delivery_constants():
    """Read EC50 and the Hill coefficient from the table script 05 wrote.

    They are duplicated nowhere: hard-coding them here would let the two
    scripts drift apart silently if the Emax fit were ever re-done.
    """
    path = os.path.join(RES, "delivery_constants.csv")
    with open(path) as f:
        vals = {r["constant"]: float(r["value"]) for r in csv.DictReader(f)}
    return vals["EC50_nmol_per_g"], vals["Hill"]


EC50, HILL = _fitted_delivery_constants()
STRIATUM_G = 20.0
BODY_KG = 70.0
ESCAPE = 0.02              # assumed endosomal escape fraction
MW_KDA = [7.0, 14.0, 28.0]  # assumed: single-strand, di-valent, conjugated
FEASIBLE_MGKG = [1.0, 3.0, 10.0]   # assumed ceilings for a systemic oligonucleotide

ETA = {"bivalent, KD 21 nM": 0.040,
       "monovalent, KD 22 nM (best measured)": 0.230,
       "monovalent, KD at modelled optimum": 0.461}


def exposure_for_msh3(target_pct):
    return EC50 * ((100.0 / target_pct) - 1.0) ** (1.0 / HILL)


rows = []
for target in (50.0, 25.0):
    x_need = exposure_for_msh3(target)
    for design, eta in ETA.items():
        dose_nmol = x_need / ((eta / 100.0) * ESCAPE)
        for mw in MW_KDA:
            mg = dose_nmol * mw * 1e-3          # nmol * kDa = ug; /1000 -> mg
            mgkg = mg / BODY_KG
            for ceiling in FEASIBLE_MGKG:
                rows.append({
                    "target_MSH3_pct_of_control": target,
                    "MSH3_knockdown_pct": round(100 - target, 1),
                    "expansion_suppression_pct": round(
                        float(hd.prevention_from_msh3(target)), 1),
                    "design": design,
                    "parenchymal_pctID_per_g": eta,
                    "assumed_MW_kDa": mw,
                    "required_systemic_dose_mg": round(mg, 1),
                    "required_dose_mg_per_kg": round(mgkg, 1),
                    "feasible_ceiling_mg_per_kg": ceiling,
                    "fold_improvement_needed": round(mgkg / ceiling, 1)})

with open(os.path.join(RES, "translation.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("systemic RMT nanocarrier, 75% striatal MSH3 lowering, 14 kDa payload")
for r in rows:
    if (r["target_MSH3_pct_of_control"] == 25.0 and r["assumed_MW_kDa"] == 14.0
            and r["feasible_ceiling_mg_per_kg"] == 3.0):
        print(f"  {r['design']:<40} {r['required_dose_mg_per_kg']:>9} mg/kg"
              f"   -> needs x{r['fold_improvement_needed']} better delivery")

print("\nrange of required improvement across all assumptions "
      f"({min(r['fold_improvement_needed'] for r in rows):.0f}x to "
      f"{max(r['fold_improvement_needed'] for r in rows):.0f}x)")

# direct comparison: dose actually delivered to the striatum
for target in (50.0, 25.0):
    x_need = exposure_for_msh3(target)
    print(f"\ntarget MSH3 {target:.0f}% of control -> striatal exposure "
          f"{x_need:.1f} nmol/g; mass in a {STRIATUM_G:.0f} g human striatum = "
          f"{x_need * STRIATUM_G:.0f} nmol "
          f"({x_need * STRIATUM_G * 14e-3:.1f} mg at 14 kDa)")
