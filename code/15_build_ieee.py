"""15 -- Build the IEEE two-column version of the paper.

Converts paper.md into an IEEEtran LaTeX source and compiles it with pdflatex.
The conversion is deliberately written out here rather than handed to a generic
markdown converter, because the IEEE layout imposes things a generic converter
gets wrong: figures must be numbered in order of first citation, tables need
captions above them, wide multi-panel figures have to span both columns, and
unicode mathematics has to become real LaTeX mathematics rather than a font
substitution.

What this script does, in order:

  1. splits paper.md into its title block, abstract and body
  2. strips the manual section numbering, since IEEEtran numbers sections itself
  3. turns each image plus its following bold caption paragraph into a float,
     and rewrites every in-text "Figure N" into a cross-reference, so the
     printed numbers follow citation order as IEEE requires
  4. converts pipe tables into tabularx floats with captions
  5. converts unicode mathematics into LaTeX mathematics
  6. reformats the reference list into IEEE citation style
  7. compiles three times so cross-references and the reference list settle

Outputs
-------
build/ieee/paper_ieee.tex
dist/Huntington_somatic_expansion_model_IEEE.pdf
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
BUILD = os.path.join(ROOT, "build", "ieee")
DIST = os.path.join(ROOT, "dist")
SRC = os.path.join(ROOT, "paper.md")
TEXNAME = "paper_ieee"
PDFOUT = os.path.join(DIST, "Huntington_somatic_expansion_model_IEEE.pdf")

# ---------------------------------------------------------------- abstract ----
# The source abstract is a structured, six-paragraph summary. IEEE wants a
# single paragraph of roughly 150 to 250 words, so this is the condensed
# version. Every figure quoted here appears in the body with its source.
ABSTRACT = (
    "Huntington's disease has long been modelled as a polyglutamine dose "
    "problem, but single-cell sequencing of post-mortem striatum shows that "
    "striatal projection neurons expand their CAG tract somatically for "
    "decades, tolerate it until roughly 150 repeats, and are then eliminated "
    "within months. That makes the expansion rate, not the inherited allele, "
    "the quantity a drug should target, and MSH3 the most direct way to reach "
    "it. This work builds a three-layer stochastic model of that clock, fits "
    "four parameters to four published kinetic anchors, and keeps five further "
    "observations back for validation. The fit requires a switch-like rather "
    "than a graded length dependence, centred on 80.2 repeats, and a "
    "first-order rather than cooperative dependence on MSH3, so suppression "
    "tracks target engagement proportionally. The calibrated model reproduces "
    "four of five held-out observations and predicts, without being fitted to "
    "them, indirect-pathway depletion, three genetic-modifier effect sizes and "
    "the loss at which clinical onset should occur. Twenty-nine percent "
    "striatal MSH3 lowering, the best figure measured in primates, begun at "
    "age 30, preserves 11.4 percentage points of striatum at age 70 across "
    "full parameter uncertainty; starting at 50 removes the benefit. The "
    "limiting factor is target engagement rather than pharmacology: the "
    "striatum needs 7 mg of drug and a systemic transferrin-receptor "
    "nanocarrier would require 7.6 g. A three-year trial enrolling at 40 would "
    "need about 133 participants per arm."
)

INDEX_TERMS = (
    "Huntington's disease, somatic repeat expansion, mismatch repair, MSH3, "
    "stochastic simulation, Gillespie algorithm, parameter identifiability, "
    "antisense oligonucleotide delivery, blood-brain barrier, clinical trial "
    "design"
)

# ------------------------------------------------------------ table captions --
# One per pipe table, in order of appearance in paper.md.
TABLE_CAPTIONS = [
    "Published kinetic anchors used for calibration.",
    "Mapping of known genetic modifiers onto parameters of the fitted rate law.",
    "Provenance of every parameter in the model. Entries marked \\emph{fitted "
    "here} were estimated in this work; entries marked \\emph{assumed} carry no "
    "supporting measurement and are all varied in the sensitivity analysis.",
    "Held-out observations, none of which was used in calibration, against the "
    "calibrated population model.",
    "Target engagement achieved by intrathecal delivery in primate striatum, "
    "and the expansion suppression it implies.",
    "Systemic dose required for 75\\% striatal MSH3 lowering, by nanocarrier "
    "design.",
    "Predicted outcome of lifelong treatment from age 30 in a 42-repeat "
    "carrier. Censored entries are lower bounds set by the simulated lifespan.",
    "Lifetime effect of representative dosing regimens, against the idealised "
    "step change, for a 42-repeat carrier treated from age 30.",
    "Participants per arm for 80\\% power at $\\alpha = 0.05$, enrolling at age "
    "40 with 2-point imaging measurement noise.",
]

# --------------------------------------------------------------- references ---
# Reformatted into IEEE style: author initials first, title in quotation marks,
# publication in italics. No volume or page range is invented; where the source
# page did not expose one, the DOI stands in its place.
REFERENCES = [
    ("handsaker", "R. E. Handsaker, S. Kashin, N. M. Reed, \\emph{et al.}, "
     "``Long somatic DNA-repeat expansion drives neurodegeneration in "
     "Huntington's disease,'' \\emph{Cell}, 2025, PMID 39824182."),
    ("breakpoint", "``The breaking point where repeat expansion triggers "
     "neuronal collapse in Huntington's disease,'' \\emph{Cell Genomics}, "
     "vol.~5, no.~4, 2025. Commentary on~\\cite{handsaker}."),
    ("tandem", "J. Li, H. Wang, and W. Yang, ``Tandem MutS$\\beta$ binding to "
     "long extruded DNA trinucleotide repeats underpins pathogenic "
     "expansions,'' \\emph{bioRxiv} 2023.12.12.571350, 2023, "
     "doi:10.1101/2023.12.12.571350. Preprint, not peer reviewed."),
    ("disirna", "``Dose-dependent reduction of somatic expansions but not Htt "
     "aggregates by di-valent siRNA-mediated silencing of MSH3 in HdhQ111 "
     "mice,'' \\emph{Scientific Reports}, 2024, doi:10.1038/s41598-024-52667-3."),
    ("senoussi", "I. Senoussi, V. Mengoli, A. Cerana, \\emph{et al.}, "
     "``Mechanism of trinucleotide repeat expansion by MutS$\\beta$ and "
     "MutL$\\gamma$ and contraction by FAN1,'' \\emph{Nature Communications}, "
     "vol.~16, p.~9445, 2025, doi:10.1038/s41467-025-64485-w."),
    ("panigrahi", "G. B. Panigrahi \\emph{et al.}, ``Isolated short CTG/CAG DNA "
     "slip-outs are repaired efficiently by hMutS$\\beta$, but clustered "
     "slip-outs are poorly repaired,'' \\emph{Proceedings of the National "
     "Academy of Sciences USA}, 2010, doi:10.1073/pnas.0909087107."),
    ("frei", "J. A. Frei, J. E. Gentile, Y. Lian, M. A. Mortberg, J. Capitanio, "
     "P. Jafar-nejad, S. M. Vallabh, H. T. Zhao, and E. V. Minikel, ``Cell type "
     "distribution of intrathecal antisense oligonucleotide activity in deep "
     "brain regions of non-human primates,'' \\emph{Nucleic Acid Therapeutics}, "
     "2025, doi:10.1177/21593337251371594."),
    ("johnsen", "K. B. Johnsen, M. Bak, P. J. Kempen, F. Melander, A. Burkhart, "
     "M. S. Thomsen, M. S. Nielsen, T. Moos, and T. L. Andresen, ``Antibody "
     "affinity and valency impact brain uptake of transferrin "
     "receptor-targeted gold nanoparticles,'' \\emph{Theranostics}, vol.~8, "
     "no.~12, pp.~3416--3436, 2018, doi:10.7150/thno.25228."),
    ("oreilly", "D. O'Reilly \\emph{et al.}, ``Di-valent siRNA-mediated "
     "silencing of MSH3 blocks somatic repeat expansion in mouse models of "
     "Huntington's disease,'' \\emph{Molecular Therapy}, 2023."),
    ("orth", "M. Orth and C. Schwenke, ``Age-at-onset in Huntington disease,'' "
     "\\emph{PLOS Currents Huntington Disease}, ed.~1, Jul. 2011, "
     "doi:10.1371/currents.RRN1258."),
    ("aso", "``Antisense oligonucleotide-mediated MSH3 suppression reduces "
     "somatic CAG repeat expansion in Huntington's disease iPSC-derived "
     "striatal neurons,'' \\emph{Science Translational Medicine}, 2025, "
     "doi:10.1126/scitranslmed.adn4600."),
    ("enroll", "Enroll-HD, ``Special characteristics of HD data,'' citing "
     "J.~H. Warner \\emph{et al.}, ``Standardizing the CAP score in "
     "Huntington's disease.'' [Online]. Available: "
     "\\url{https://enroll-hd.org/for-researchers/analyzing-data/}"),
]
CITEKEYS = [k for k, _ in REFERENCES]

# ------------------------------------------------------------------- unicode --
# Compound tokens first, so that related characters end up inside one set of
# maths delimiters rather than several.
COMPOUND = [
    ("z\u2081\u208b\u03b1/\u2082 + z\u2081\u208b\u03b2",
     "$z_{1-\\alpha/2} + z_{1-\\beta}$"),
    ("\u03bb(L)", "$\\lambda(L)$"),
    ("\u03bb(40)", "$\\lambda(40)$"),
    ("\u03bb(90)/\u03bb(70)", "$\\lambda(90)/\\lambda(70)$"),
    ("\u03bb(150)", "$\\lambda(150)$"),
    ("\u03c6(L)", "$\\varphi(L)$"),
    ("EC\u2085\u2080", "$\\mathrm{EC}_{50}$"),
    ("L\u2085\u2080", "$L_{50}$"),
    ("age\u2080", "$\\mathrm{age}_{0}$"),
    ("nmol g\u207b\u00b9", "nmol\\,g$^{-1}$"),
    ("y\u207b\u00b9", "y$^{-1}$"),
    ("repeat\u207b\u00b9", "repeat$^{-1}$"),
    ("kcal mol\u207b\u00b9", "kcal\\,mol$^{-1}$"),
    ("mol\u207b\u00b9", "mol$^{-1}$"),
    ("\u0394G_slip", "$\\Delta G_{\\mathrm{slip}}$"),
    ("\u03c3\u00b2/\u03b4\u00b2", "$\\sigma^{2}/\\delta^{2}$"),
    ("10\u207b\u00b3\u00b3", "$10^{-33}$"),
]
SINGLE = {
    "\u03b1": "$\\alpha$", "\u03b2": "$\\beta$", "\u03b3": "$\\gamma$",
    "\u03b4": "$\\delta$", "\u03b7": "$\\eta$", "\u03b8": "$\\theta$",
    "\u03bb": "$\\lambda$", "\u03c3": "$\\sigma$", "\u03c6": "$\\varphi$",
    "\u0394": "$\\Delta$",
    "\u2212": "$-$", "\u00d7": "$\\times$", "\u00b1": "$\\pm$",
    "\u221a": "$\\sqrt{\\ }$", "\u222b": "$\\int$",
    "\u00b7": "$\\cdot$",
    "\u2080": "$_{0}$", "\u2081": "$_{1}$", "\u2082": "$_{2}$",
    "\u2085": "$_{5}$", "\u208b": "$_{-}$",
    "\u00b9": "$^{1}$", "\u00b2": "$^{2}$", "\u00b3": "$^{3}$",
    "\u207b": "$^{-}$",
}

LATEX_SPECIALS = {
    "\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$",
    "&": r"\&", "#": r"\#", "_": r"\_", "%": r"\%", "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def escape(text):
    return "".join(LATEX_SPECIALS.get(c, c) for c in text)


# Variables written with an underscore in the source become real subscripts.
# These run on already-escaped text, so the underscore appears as "\\_" and the
# caret as "\\textasciicircum{}". Longest first, so K_D,b is not eaten by K_D.
SUBSCRIPTS = [
    ("K\\_D,b", "$K_{D,b}$"),
    ("K\\_D", "$K_{D}$"),
    ("K\\_r", "$K_{r}$"),
    ("L\\_sw", "$L_{\\mathrm{sw}}$"),
    ("w\\_sw", "$w_{\\mathrm{sw}}$"),
    ("L\\_pen", "$L_{\\mathrm{pen}}$"),
    ("L\\_death", "$L_{\\mathrm{death}}$"),
    ("L\\_commit", "$L_{\\mathrm{commit}}$"),
    ("L\\_C", "$L_{C}$"),
    ("k\\_nca", "$k_{\\mathrm{nca}}$"),
    ("f\\_lost", "$f_{\\mathrm{lost}}$"),
    ("f\\_v", "$f_{v}$"),
    ("h\\_max", "$h_{\\mathrm{max}}$"),
    ("m\\textasciicircum{}k", "$m^{k}$"),
    ("2\\textasciicircum{}k", "$2^{k}$"),
]


def to_math(text):
    for src, dst in SUBSCRIPTS:
        text = text.replace(src, dst)
    for src, dst in COMPOUND:
        text = text.replace(src, dst)
    text = "".join(SINGLE.get(c, c) for c in text)
    # merge adjacent maths groups: $a$$b$ -> $ab$
    text = re.sub(r"\$\$", "", text)
    return text


def inline(text, allow_cite=True):
    """Markdown inline markup to LaTeX, for one paragraph of prose."""
    spans = []

    def stash(m):
        spans.append(m.group(1))
        return "\x00%d\x00" % (len(spans) - 1)

    text = re.sub(r"`([^`]+)`", stash, text)
    text = escape(text)
    text = to_math(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"(?<![\w*])\*([^*\n]+?)\*(?![\w*])", r"\\textit{\1}", text)
    if allow_cite:
        text = re.sub(r"\[(\d+(?:\s*,\s*\d+)*)\]", _cite, text)
    text = re.sub(r"\bsections?\s+(\d+)\.(\d+)\b",
                  lambda m: "Section~\\ref{sec:%s.%s}" % (m.group(1), m.group(2)),
                  text, flags=re.I)
    # bare "Section 6", after the N.M form so it cannot eat one
    text = re.sub(r"\bsections?\s+(\d+)(?![.~\d])",
                  lambda m: "Section~\\ref{sec:%s}" % m.group(1),
                  text, flags=re.I)
    # "Table 1" in the source means the parameter-provenance table, which the
    # IEEE numbering reaches third
    text = re.sub(r"\bTable\s+1\b",
                  lambda m: "Table~\\ref{tab:3}", text)
    text = re.sub(r"\bFigures?~?\s*(\d+)([A-C])?\b", _figref, text)
    # straight quotation marks to LaTeX opening and closing pairs
    text = re.sub(r'"([^"]*)"', r"``\1''", text)
    for i, s in enumerate(spans):
        text = text.replace("\x00%d\x00" % i,
                            "\\texttt{%s}" % escape(s).replace("-", "-\\/"))
    return text


def _cite(m):
    nums = [int(x) for x in m.group(1).split(",")]
    return "\\cite{%s}" % ",".join(CITEKEYS[n - 1] for n in nums)


def _figref(m):
    panel = m.group(2) or ""
    return "Fig.~\\ref{fig:%s}%s" % (m.group(1), panel)


# ------------------------------------------------------------------- parsing --
def split_source(md):
    lines = md.split("\n")
    # abstract runs from '## Abstract' to the next '## '
    try:
        a = next(i for i, l in enumerate(lines) if l.strip() == "## Abstract")
    except StopIteration:
        sys.exit("no abstract heading found")
    b = next(i for i in range(a + 1, len(lines)) if lines[i].startswith("## "))
    body = [l for l in lines[b:]]
    return body


def convert_table(rows, caption, number):
    header = [c.strip() for c in rows[0].strip().strip("|").split("|")]
    data = [[c.strip() for c in r.strip().strip("|").split("|")]
            for r in rows[2:]]
    n = len(header)
    # Weight the column widths by the longest cell in each column, so a column
    # of short numbers does not get the same space as a column of prose.
    widths = []
    for j in range(n):
        longest = max([len(header[j])] + [len(r[j]) for r in data if j < len(r)])
        widths.append(max(longest, 6))
    total = float(sum(widths))
    spec = "".join(">{\\raggedright\\arraybackslash}p{%.3f\\colw}%s"
                   % (w / total, "" if j == n - 1 else "")
                   for j, w in enumerate(widths))
    wide = total > 70 or n >= 4
    env = "table*" if wide else "table"
    colw = ("\\dimexpr\\textwidth-%d\\tabcolsep\\relax" % (2 * n)) if wide \
        else ("\\dimexpr\\columnwidth-%d\\tabcolsep\\relax" % (2 * n))
    out = ["\\begin{%s}[!t]" % env,
           "\\caption{%s}" % caption,
           "\\label{tab:%d}" % number,
           "\\centering",
           "\\footnotesize",
           "\\setlength{\\colw}{%s}" % colw,
           "\\begin{tabular}{%s}" % spec,
           "\\toprule"]
    out.append(" & ".join("\\textbf{%s}" % inline(h) for h in header)
               + " \\\\")
    out.append("\\midrule")
    for r in data:
        cells = [inline(c) for c in r] + [""] * (n - len(r))
        out.append(" & ".join(cells[:n]) + " \\\\")
    out += ["\\bottomrule", "\\end{tabular}", "\\end{%s}" % env, ""]
    return out


def build_body(body):
    out, i = [], 0
    table_n = 0
    while i < len(body):
        line = body[i]
        s = line.strip()

        if s in ("---", ""):
            i += 1
            continue

        # headings ------------------------------------------------------------
        m = re.match(r"^###\s+(\d+)\.(\d+)\s+(.*)$", s)
        if m:
            out.append("\\subsection{%s}" % inline(m.group(3)))
            out.append("\\label{sec:%s.%s}" % (m.group(1), m.group(2)))
            i += 1
            continue
        m = re.match(r"^##\s+(\d+)\.\s+(.*)$", s)
        if m:
            out.append("\\section{%s}" % inline(m.group(2)))
            out.append("\\label{sec:%s}" % m.group(1))
            i += 1
            continue
        m = re.match(r"^##\s+(.*)$", s)
        if m:
            out.append("\\section{%s}" % inline(m.group(1)))
            i += 1
            continue

        # figure: image line, then a blank, then the bold caption paragraph ---
        m = re.match(r"^!\[[^\]]*\]\(([^)]+\.png)\)$", s)
        if m:
            path = m.group(1)
            num = re.search(r"fig(\d+)", os.path.basename(path)).group(1)
            j = i + 1
            while j < len(body) and not body[j].strip():
                j += 1
            cap = ""
            if j < len(body) and body[j].strip().startswith("**Figure"):
                cap = body[j].strip()
                cap = re.sub(r"^\*\*Figure\s+\d+\.\s*", "**", cap)
                i = j
            out += ["\\begin{figure*}[!t]",
                    "\\centering",
                    "\\includegraphics[width=\\textwidth]{%s}" % path,
                    "\\caption{%s}" % inline(cap),
                    "\\label{fig:%s}" % num,
                    "\\end{figure*}", ""]
            i += 1
            continue

        # display equation ----------------------------------------------------
        if s.startswith("> "):
            eq = s[2:].strip().strip("*").strip()
            out += ["\\begin{equation}",
                    "\\begin{split}",
                    _equation(eq),
                    "\\end{split}",
                    "\\end{equation}", ""]
            i += 1
            continue

        # fenced block: the file listing in the reproducibility section -------
        if s.startswith("```"):
            j = i + 1
            block = []
            while j < len(body) and not body[j].strip().startswith("```"):
                block.append(body[j])
                j += 1
            out += _listing(block)
            i = j + 1
            continue

        # pipe table ----------------------------------------------------------
        if s.startswith("|"):
            j = i
            rows = []
            while j < len(body) and body[j].strip().startswith("|"):
                rows.append(body[j])
                j += 1
            cap = (TABLE_CAPTIONS[table_n] if table_n < len(TABLE_CAPTIONS)
                   else "Results.")
            table_n += 1
            out += convert_table(rows, cap, table_n)
            i = j
            continue

        # numbered list -------------------------------------------------------
        if re.match(r"^\d+\.\s", s):
            items = []
            while i < len(body) and re.match(r"^\d+\.\s", body[i].strip()):
                items.append(re.sub(r"^\d+\.\s+", "", body[i].strip()))
                i += 1
            out.append("\\begin{enumerate}")
            out += ["\\item %s" % inline(t) for t in items]
            out += ["\\end{enumerate}", ""]
            continue

        # bulleted list -------------------------------------------------------
        if s.startswith("- "):
            items = []
            while i < len(body) and body[i].strip().startswith("- "):
                items.append(body[i].strip()[2:])
                i += 1
            out.append("\\begin{itemize}")
            out += ["\\item %s" % inline(t) for t in items]
            out += ["\\end{itemize}", ""]
            continue

        # reference list: handled separately, so stop here --------------------
        out += [inline(s), ""]
        i += 1
    return out


def _equation(eq):
    """The two display equations, written out as real LaTeX."""
    if eq.startswith("\u03bb(L)"):
        return (r"\lambda(L) = A \cdot n(L) \cdot \varphi(L) \cdot "
                r"\bigl[\,(1 - C(L)) + R\,C(L)\,\bigr] \cdot \theta / f")
    if eq.startswith("\u03b7"):
        return (r"\eta(K_D, v) = f_v \cdot "
                r"\frac{K_{D,b}}{K_D + K_{D,b}} \cdot "
                r"\frac{1}{1 + K_r/K_D}")
    return escape(eq)


def _listing(block):
    rows = []
    for ln in block:
        if not ln.strip():
            continue
        parts = re.split(r"\s{2,}", ln.strip(), maxsplit=1)
        rows.append((parts[0], parts[1] if len(parts) > 1 else ""))
    out = ["\\begin{table*}[!t]",
           "\\caption{Contents of the code repository. Scripts run in "
           "numerical order; each writes its own result tables.}",
           "\\label{tab:code}",
           "\\centering", "\\footnotesize",
           "\\begin{tabular}{@{}ll@{}}", "\\toprule",
           "\\textbf{Path} & \\textbf{Contents} \\\\", "\\midrule"]
    for a, b in rows:
        out.append("\\texttt{%s} & %s \\\\"
                   % (escape(a).replace("-", "-\\/"),
                      inline(b.replace("->", "gives"))))
    out += ["\\bottomrule", "\\end{tabular}", "\\end{table*}", ""]
    return out


PREAMBLE = r"""\documentclass[conference]{IEEEtran}
\IEEEoverridecommandlockouts
\usepackage{cite}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{array}
\usepackage{url}
\usepackage[hidelinks]{hyperref}
\usepackage{balance}
\newlength{\colw}
\graphicspath{{%s/}}
\def\BibTeX{{\rm B\kern-.05em{\sc i\kern-.025em b}\kern-.08em T\kern-.1667em\lower.7ex\hbox{E}\kern-.125emX}}
\begin{document}

\title{Somatic CAG-Repeat Expansion as a Druggable Rate Process in
Huntington's Disease\\[0.7\baselineskip]
\parbox{0.88\textwidth}{\centering\large\itshape A Multiscale In Silico Model
Linking MutS$\beta$ Biochemistry, Striatal Neuron Loss, and the Delivery Gap
for MSH3-Lowering Oligonucleotides\par}\\[0.25\baselineskip]}

\author{\IEEEauthorblockN{Shubho}
\IEEEauthorblockA{\textit{Independent Research} \\
shubhodeepcx@gmail.com}}

\maketitle

\begin{abstract}
%s
\end{abstract}

\begin{IEEEkeywords}
%s
\end{IEEEkeywords}
"""

POSTAMBLE = r"""
\balance
\begin{thebibliography}{99}
%s
\end{thebibliography}

\end{document}
"""


def main():
    os.makedirs(BUILD, exist_ok=True)
    os.makedirs(DIST, exist_ok=True)
    md = open(SRC, encoding="utf-8").read()
    body = split_source(md)

    # drop the markdown reference list; it is replaced by thebibliography
    try:
        cut = next(i for i, l in enumerate(body)
                   if re.match(r"^##\s+\d+\.\s+References", l.strip()))
        body = body[:cut]
    except StopIteration:
        print("warning: no references heading found, keeping whole body")

    tex = [PREAMBLE % (ROOT.replace("\\", "/"), ABSTRACT, INDEX_TERMS)]
    tex += build_body(body)
    bib = "\n".join("\\bibitem{%s} %s" % (k, v) for k, v in REFERENCES)
    tex.append(POSTAMBLE % bib)

    path = os.path.join(BUILD, TEXNAME + ".tex")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(tex) + "\n")
    print("wrote %s (%d lines)" % (path, len(tex)))

    leftover = [c for c in open(path, encoding="utf-8").read() if ord(c) > 127]
    if leftover:
        print("warning: %d non-ASCII characters remain: %s"
              % (len(leftover), sorted(set(leftover))))

    for run in (1, 2, 3):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                            "-halt-on-error", "-file-line-error",
                            TEXNAME + ".tex"],
                           cwd=BUILD, capture_output=True, text=True)
        if p.returncode != 0:
            print("pdflatex failed on pass %d" % run)
            for ln in p.stdout.split("\n"):
                if re.search(r"^.*:\d+:|^!|Undefined|Overfull.*(1[0-9]|[2-9])"
                             r"[0-9]\.[0-9]+pt", ln):
                    print("  " + ln)
            sys.exit(1)
        print("pdflatex pass %d ok" % run)

    shutil.copy(os.path.join(BUILD, TEXNAME + ".pdf"), PDFOUT)
    shutil.copy(path, os.path.join(DIST, TEXNAME + ".tex"))
    print("wrote %s" % PDFOUT)

    log = open(os.path.join(BUILD, TEXNAME + ".log"), encoding="utf-8",
               errors="replace").read()
    for kind in ("Undefined control sequence", "Undefined citation",
                 "There were undefined references",
                 "Reference .* undefined", "Citation .* undefined"):
        for m in re.finditer(kind, log):
            print("LOG: " + m.group(0))
    over = re.findall(r"Overfull \\hbox \((\d+\.\d+)pt", log)
    bad = [float(x) for x in over if float(x) > 10]
    print("overfull boxes: %d total, %d worse than 10pt" % (len(over), len(bad)))


if __name__ == "__main__":
    main()
