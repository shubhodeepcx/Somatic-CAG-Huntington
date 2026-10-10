# Somatic CAG-repeat expansion as a druggable rate process in Huntington's disease

A self-contained computational research project: a multiscale model of somatic
CAG-repeat expansion in striatal projection neurons, calibrated to published
summary statistics, used to ask how much striatal neuron loss MSH3 lowering can
prevent, whether current delivery technology can achieve it, and whether a trial
of plausible size could see the effect.

**Read the paper:** [`paper.md`](paper.md), or the built versions in
[`dist/`](dist) as Word and PDF.

## Layout

| Path | Contents |
|---|---|
| `paper.md` | The full write-up: abstract, introduction, methods, results, discussion, limitations, references |
| `code/hdmodel.py` | Core model: rate law, calibration, exact Gillespie and vectorised tau-leaping solvers |
| `code/01_calibrate.py` | Calibration, identifiability scan, MutSβ stoichiometry inference, Figures 1 and 2 |
| `code/02_population.py` | Population simulation, held-out validation, Figure 3 |
| `code/03_coupling.py` | Solver verification, non-cell-autonomous coupling test, Figure 4 |
| `code/04_reconcile.py` | Reconciliation of conflicting published constraints, Figure 5 |
| `code/05_delivery.py` | Exposure-response, CSF route, nanocarrier design space, Figures 6 and 7 |
| `code/06_trial.py` | In silico trial, Figures 8 and 9 |
| `code/07_sensitivity.py` | Latin hypercube sensitivity analysis with recalibration, Figure 10 |
| `code/08_translation.py` | Dose translation to mg/kg |
| `code/09_build_documents.py` | Builds the Word and PDF versions of the paper |
| `code/10_profile.py` | Parameter profiles and feasible ranges, Figure 11 |
| `code/11_subtypes.py` | Direct and indirect pathway prediction, Figure 12 |
| `code/12_modifiers.py` | Genetic modifiers as parameter perturbations, Figure 13 |
| `code/13_dosing.py` | Intermittent dosing and time-averaged knockdown, Figure 14 |
| `code/14_trial_power.py` | Sample size and power for a real trial, Figure 15 |
| `results/` | Every table in the paper, as CSV or JSON, written by the script that produced it |
| `figures/` | All fifteen figures |
| `dist/` | The paper as `.docx` and `.pdf` |

## Running it

```bash
python3 -m pip install numpy scipy matplotlib
cd code
for f in 01_*.py 02_*.py 03_*.py 04_*.py 05_*.py 06_*.py 07_*.py 08_*.py \
         10_*.py 11_*.py 12_*.py 13_*.py 14_*.py; do python3 "$f"; done
python3 09_build_documents.py   # needs pandoc and LibreOffice
```

Total runtime is about forty minutes on a single core, most of it in the
sensitivity analysis, which recalibrates the model 160 times. All random seeds
are fixed, so the figures and tables regenerate exactly.

Tested with Python 3, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, pandoc 3.1.3
and LibreOffice 24.2.7.2.

## Provenance

No experimental data are redistributed here. Every calibration target is a
published summary statistic, listed with its source in Table 1 of the paper.
Parameters fitted in this work and parameters simply assumed are labelled as
such, and all of the latter are varied in the global sensitivity analysis.
Four published numbers were fitted; everything else, including five held-out
validations, three modifier predictions, a cell-subtype prediction and the trial
power calculation, follows without further fitting.
