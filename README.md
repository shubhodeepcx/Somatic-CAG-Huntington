# Somatic CAG-repeat expansion as a druggable rate process in Huntington's disease

A self-contained computational research project: a multiscale model of somatic
CAG-repeat expansion in striatal projection neurons, calibrated to published
summary statistics, used to ask how much striatal neuron loss MSH3 lowering can
prevent, whether current delivery technology can achieve it, and whether a trial
of plausible size could see the effect.

**Read the paper:** [`paper.md`](paper.md), or the built versions in
[`dist/`](dist). Four builds are produced from the same source: a plain Word
and PDF pair for reading, and an IEEE two-column Word and PDF pair in the
conference format, alongside the LaTeX source that produced it.

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
| `code/15_build_ieee.py` | Builds the IEEE two-column LaTeX source and PDF |
| `code/16_build_ieee_docx.py` | Builds the IEEE two-column Word version |
| `code/17_estimator.py` | Mean-field against stochastic calibration, section 3.5 |
| `review/` | An internal audit of the model, the code and the paper, and what was done about each finding |
| `results/` | Every table in the paper, as CSV or JSON, written by the script that produced it |
| `figures/` | All fifteen figures |
| `dist/` | The paper as `.docx` and `.pdf`, plus the IEEE versions and `paper_ieee.tex` |

## Running it

```bash
python3 -m pip install numpy scipy matplotlib
cd code
for f in 01_*.py 02_*.py 03_*.py 04_*.py 05_*.py 06_*.py 07_*.py 08_*.py \
         10_*.py 11_*.py 12_*.py 13_*.py 14_*.py 17_*.py; do python3 "$f"; done
python3 09_build_documents.py   # plain Word and PDF; needs pandoc and LibreOffice
python3 15_build_ieee.py        # IEEE PDF; needs pdflatex and the IEEEtran class
python3 16_build_ieee_docx.py   # IEEE Word; needs python-docx
```

Total runtime is about forty minutes on a single core, most of it in the
sensitivity analysis, which recalibrates the model 160 times. All random seeds
are fixed, so the figures and tables regenerate exactly.

Tested with Python 3, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, pandoc 3.1.3,
LibreOffice 24.2.7.2, python-docx 1.2.0 and TeX Live 2023 with `texlive-publishers`
for the IEEEtran class.

`review/` carries an internal audit of this work and a log of what was done
about each finding. It is kept in the repository rather than discarded because
two of the findings changed results, and the second of them is worth knowing
about before reading the paper: calibration minimises a mean-field integral
whose stochastic counterpart runs 21 to 50% slower, and section 3.5 reports
what refitting against the right quantity does to each conclusion.

The IEEE build is not a reformat of the Word version. It renumbers the figures
into order of first citation, as IEEE requires, and rewrites every in-text
reference to match, so a figure called 11 in `paper.md` may print as Fig. 2.

## Provenance

No experimental data are redistributed here. Every calibration target is a
published summary statistic, listed with its source in Table 1 of the paper.
Parameters fitted in this work and parameters simply assumed are labelled as
such. Every assumed value is varied, but not all in the same way: the constants
of the expansion and elimination model are swept simultaneously in the global
sensitivity analysis, which recalibrates the rate law inside each sample, while
the delivery, translation and onset-regression constants are swept one at a
time. Table 1 of the paper states which treatment each one receives.
Four published numbers were fitted; everything else, including five held-out
validations, three modifier predictions, a cell-subtype prediction and the trial
power calculation, follows without further fitting.
