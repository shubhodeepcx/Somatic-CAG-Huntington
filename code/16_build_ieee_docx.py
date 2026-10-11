"""16 -- Build the IEEE two-column Word version of the paper.

The LaTeX build in script 15 is the authoritative IEEE version. This one
produces the same layout as a .docx, for people who need to edit or comment in
Word. It is written against python-docx directly rather than going through a
generic converter, because the IEEE layout needs things a converter will not
do: a single-column title block above a two-column body, wide figures that
break out of the columns and back again, tables numbered in Roman numerals
with the caption above, and cross-references that follow the printed numbering
rather than the numbering in the source file.

Numbering is worked out in a pre-pass over paper.md, so "Figure 11" in the
source becomes whatever number that figure actually receives once figures are
ordered by first citation, and the same for tables and sections.

Outputs
-------
dist/Huntington_somatic_expansion_model_IEEE.docx
"""
import os
import re

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

import importlib
ieee = importlib.import_module("15_build_ieee") if False else None

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "paper.md")
OUT = os.path.join(ROOT, "dist",
                   "Huntington_somatic_expansion_model_IEEE.docx")

FONT = "Times New Roman"
BODY_PT = 10
SMALL_PT = 8
ROMAN = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV"]

SUBSCRIPTS = ["K_D,b", "K_D", "K_r", "L_sw", "w_sw", "L_pen", "L_death",
              "L_commit", "L_C", "k_nca", "f_lost", "f_v", "h_max"]
SUPERSCRIPTS = ["m^k", "2^k"]

# Pulled from the LaTeX build so the two versions say the same thing.
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
TABLE_CAPTIONS = [
    "Published kinetic anchors used for calibration.",
    "Mapping of known genetic modifiers onto parameters of the fitted rate law.",
    "Provenance of every parameter in the model. Entries marked fitted here "
    "were estimated in this work; entries marked assumed carry no supporting "
    "measurement and are all varied in the sensitivity analysis.",
    "Held-out observations, none of which was used in calibration, against the "
    "calibrated population model.",
    "Target engagement achieved by intrathecal delivery in primate striatum, "
    "and the expansion suppression it implies.",
    "Systemic dose required for 75% striatal MSH3 lowering, by nanocarrier "
    "design.",
    "Predicted outcome of lifelong treatment from age 30 in a 42-repeat "
    "carrier. Censored entries are lower bounds set by the simulated lifespan.",
    "Lifetime effect of representative dosing regimens, against the idealised "
    "step change, for a 42-repeat carrier treated from age 30.",
    "Participants per arm for 80% power at alpha = 0.05, enrolling at age 40 "
    "with 2-point imaging measurement noise.",
]
REFERENCES = [
    "R. E. Handsaker, S. Kashin, N. M. Reed, et al., “Long somatic "
    "DNA-repeat expansion drives neurodegeneration in Huntington's disease,"
    "” Cell, 2025, PMID 39824182.",
    "“The breaking point where repeat expansion triggers neuronal "
    "collapse in Huntington's disease,” Cell Genomics, vol. 5, no. 4, "
    "2025. Commentary on [1].",
    "J. Li, H. Wang, and W. Yang, “Tandem MutSβ binding to long "
    "extruded DNA trinucleotide repeats underpins pathogenic expansions,” "
    "bioRxiv 2023.12.12.571350, 2023, doi:10.1101/2023.12.12.571350. "
    "Preprint, not peer reviewed.",
    "“Dose-dependent reduction of somatic expansions but not Htt "
    "aggregates by di-valent siRNA-mediated silencing of MSH3 in HdhQ111 "
    "mice,” Scientific Reports, 2024, doi:10.1038/s41598-024-52667-3.",
    "I. Senoussi, V. Mengoli, A. Cerana, et al., “Mechanism of "
    "trinucleotide repeat expansion by MutSβ and MutLγ and "
    "contraction by FAN1,” Nature Communications, vol. 16, p. 9445, "
    "2025, doi:10.1038/s41467-025-64485-w.",
    "G. B. Panigrahi et al., “Isolated short CTG/CAG DNA slip-outs are "
    "repaired efficiently by hMutSβ, but clustered slip-outs are poorly "
    "repaired,” Proceedings of the National Academy of Sciences USA, "
    "2010, doi:10.1073/pnas.0909087107.",
    "J. A. Frei, J. E. Gentile, Y. Lian, M. A. Mortberg, J. Capitanio, "
    "P. Jafar-nejad, S. M. Vallabh, H. T. Zhao, and E. V. Minikel, “Cell "
    "type distribution of intrathecal antisense oligonucleotide activity in "
    "deep brain regions of non-human primates,” Nucleic Acid "
    "Therapeutics, 2025, doi:10.1177/21593337251371594.",
    "K. B. Johnsen, M. Bak, P. J. Kempen, F. Melander, A. Burkhart, "
    "M. S. Thomsen, M. S. Nielsen, T. Moos, and T. L. Andresen, "
    "“Antibody affinity and valency impact brain uptake of transferrin "
    "receptor-targeted gold nanoparticles,” Theranostics, vol. 8, "
    "no. 12, pp. 3416–3436, 2018, doi:10.7150/thno.25228.",
    "D. O'Reilly et al., “Di-valent siRNA-mediated silencing of MSH3 "
    "blocks somatic repeat expansion in mouse models of Huntington's "
    "disease,” Molecular Therapy, 2023.",
    "M. Orth and C. Schwenke, “Age-at-onset in Huntington disease,” "
    "PLOS Currents Huntington Disease, ed. 1, Jul. 2011, "
    "doi:10.1371/currents.RRN1258.",
    "“Antisense oligonucleotide-mediated MSH3 suppression reduces somatic "
    "CAG repeat expansion in Huntington's disease iPSC-derived striatal "
    "neurons,” Science Translational Medicine, 2025, "
    "doi:10.1126/scitranslmed.adn4600.",
    "Enroll-HD, “Special characteristics of HD data,” citing "
    "J. H. Warner et al., “Standardizing the CAP score in Huntington's "
    "disease.” [Online]. Available: https://enroll-hd.org/for-researchers/"
    "analyzing-data/",
]


# ------------------------------------------------------------- docx helpers --
def set_columns(section, num, space=Inches(0.25)):
    sectPr = section._sectPr
    cols = sectPr.xpath("./w:cols")
    cols = cols[0] if cols else OxmlElement("w:cols")
    if cols.getparent() is None:
        sectPr.append(cols)
    cols.set(qn("w:num"), str(num))
    cols.set(qn("w:space"), str(int(space.twips)))
    cols.set(qn("w:equalWidth"), "1")


def continuous(doc, columns):
    s = doc.add_section(WD_SECTION.CONTINUOUS)
    page_margins(s)
    set_columns(s, columns)
    return s


def page_margins(s):
    s.top_margin = Inches(0.75)
    s.bottom_margin = Inches(1.0)
    s.left_margin = Inches(0.625)
    s.right_margin = Inches(0.625)


def para(doc, align=WD_ALIGN_PARAGRAPH.JUSTIFY, size=BODY_PT, space_after=0,
         indent=0.0, space_before=0):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = align
    pf.space_after = Pt(space_after)
    pf.space_before = Pt(space_before)
    pf.first_line_indent = Inches(indent)
    pf.line_spacing = 1.0
    p._size = size
    return p


def run(p, text, bold=False, italic=False, size=None, sub=False, sup=False,
        caps=False, mono=False):
    r = p.add_run(text)
    r.font.name = "Consolas" if mono else FONT
    r.font.size = Pt(size if size else getattr(p, "_size", BODY_PT))
    r.bold = bold
    r.italic = italic
    r.font.subscript = sub
    r.font.superscript = sup
    r.font.all_caps = caps
    rpr = r._element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts")
        rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs"):
        rf.set(qn(a), "Consolas" if mono else FONT)
    return r


# -------------------------------------------------------------- inline text --
def tokenise(text):
    """Split prose into (text, style) pieces."""
    out = []
    pattern = re.compile(r"(\*\*.+?\*\*|(?<![\w*])\*[^*\n]+?\*(?![\w*])|`[^`]+`)")
    for piece in pattern.split(text):
        if not piece:
            continue
        if piece.startswith("**") and piece.endswith("**"):
            out += _subscripts(piece[2:-2], {"bold": True})
        elif piece.startswith("`") and piece.endswith("`"):
            out.append((piece[1:-1], {"mono": True}))
        elif piece.startswith("*") and piece.endswith("*") and len(piece) > 2:
            out += _subscripts(piece[1:-1], {"italic": True})
        else:
            out += _subscripts(piece, {})
    return out


def _subscripts(text, base):
    """Turn L_sw and m^k into real subscript and superscript runs."""
    tokens = sorted(SUBSCRIPTS + SUPERSCRIPTS, key=len, reverse=True)
    pat = re.compile("(" + "|".join(re.escape(t) for t in tokens) + ")")
    out = []
    for piece in pat.split(text):
        if not piece:
            continue
        if piece in SUBSCRIPTS:
            stem, tail = piece.split("_", 1)
            out.append((stem, dict(base, italic=True)))
            out.append((tail, dict(base, sub=True, italic=True)))
        elif piece in SUPERSCRIPTS:
            stem, tail = piece.split("^", 1)
            out.append((stem, dict(base, italic=True)))
            out.append((tail, dict(base, sup=True, italic=True)))
        else:
            out.append((piece, base))
    return out


def emit(p, text, size=None):
    for txt, st in tokenise(text):
        run(p, txt, size=size, **st)


# ----------------------------------------------------------------- numbering --
def prepass(lines):
    """Work out the printed numbering before anything is written."""
    figs, tables, secs = {}, 0, {}
    sec_i, sub_i = 0, 0
    for ln in lines:
        s = ln.strip()
        m = re.match(r"^!\[[^\]]*\]\(([^)]+\.png)\)$", s)
        if m:
            src = re.search(r"fig(\d+)", os.path.basename(m.group(1))).group(1)
            if src not in figs:
                figs[src] = len(figs) + 1
            continue
        m = re.match(r"^##\s+(\d+)\.\s", s)
        if m:
            sec_i += 1
            sub_i = 0
            secs[m.group(1)] = ROMAN[sec_i]
            continue
        m = re.match(r"^###\s+(\d+)\.(\d+)\s", s)
        if m:
            sub_i += 1
            secs["%s.%s" % (m.group(1), m.group(2))] = \
                "%s-%s" % (ROMAN[sec_i], chr(ord("A") + sub_i - 1))
            continue
        if s.startswith("|---"):
            tables += 1
    return figs, secs


def rewrite_refs(text, figs, secs):
    text = re.sub(r"\bsections?\s+(\d+)\.(\d+)\b",
                  lambda m: "Section %s" % secs.get(
                      "%s.%s" % (m.group(1), m.group(2)), m.group(0)),
                  text, flags=re.I)
    text = re.sub(r"\bsections?\s+(\d+)(?![.\d])",
                  lambda m: "Section %s" % secs.get(m.group(1), m.group(0)),
                  text, flags=re.I)
    text = re.sub(r"\bTable\s+1\b", "Table III", text)
    text = re.sub(r"\bFigures?\s*(\d+)([A-C])?\b",
                  lambda m: "Fig. %s%s" % (figs.get(m.group(1), m.group(1)),
                                           m.group(2) or ""), text)
    return text


# --------------------------------------------------------------------- build --
def main():
    md = open(SRC, encoding="utf-8").read().split("\n")
    start = next(i for i, l in enumerate(md) if l.strip() == "## Abstract")
    body_start = next(i for i in range(start + 1, len(md))
                      if md[i].startswith("## "))
    cut = next(i for i, l in enumerate(md)
               if re.match(r"^##\s+\d+\.\s+References", l.strip()))
    body = md[body_start:cut]
    figs, secs = prepass(body)

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(BODY_PT)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)

    s0 = doc.sections[0]
    page_margins(s0)
    set_columns(s0, 1)

    p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 24, space_after=2)
    run(p, "Somatic CAG-Repeat Expansion as a Druggable Rate Process in "
           "Huntington's Disease", size=24)
    p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 13, space_after=10)
    run(p, "A Multiscale In Silico Model Linking MutSβ Biochemistry, "
           "Striatal Neuron Loss, and the Delivery Gap for MSH3-Lowering "
           "Oligonucleotides", size=13, italic=True)
    for txt, it in (("Shubho", False), ("Independent Research", True),
                    ("shubhodeepcx@gmail.com", False)):
        p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 11, space_after=1)
        run(p, txt, size=11, italic=it)
    para(doc, WD_ALIGN_PARAGRAPH.CENTER, 10, space_after=6)

    continuous(doc, 2)

    p = para(doc, space_after=6, indent=0.2)
    run(p, "Abstract—", size=9, bold=True, italic=True)
    run(p, ABSTRACT, size=9, bold=True)
    p = para(doc, space_after=10, indent=0.2)
    run(p, "Index Terms—", size=9, bold=True, italic=True)
    run(p, INDEX_TERMS, size=9, bold=True)

    write_body(doc, body, figs, secs)

    p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 10, space_before=8, space_after=4)
    run(p, "References", size=10, caps=True)
    for i, ref in enumerate(REFERENCES, 1):
        p = para(doc, WD_ALIGN_PARAGRAPH.JUSTIFY, SMALL_PT, space_after=2)
        p.paragraph_format.left_indent = Inches(0.22)
        p.paragraph_format.first_line_indent = Inches(-0.22)
        run(p, "[%d] " % i, size=SMALL_PT)
        run(p, ref, size=SMALL_PT)

    doc.core_properties.title = ("Somatic CAG-Repeat Expansion as a Druggable "
                                "Rate Process in Huntington's Disease")
    doc.core_properties.author = "Shubho"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    doc.save(OUT)
    print("wrote %s (%d figures, %d sections)"
          % (OUT, len(figs), len(secs)))


def write_body(doc, body, figs, secs):
    i, table_n = 0, 0
    while i < len(body):
        s = body[i].strip()
        if s in ("---", ""):
            i += 1
            continue

        m = re.match(r"^###\s+(\d+)\.(\d+)\s+(.*)$", s)
        if m:
            p = para(doc, WD_ALIGN_PARAGRAPH.LEFT, 10, space_before=6,
                     space_after=2)
            run(p, "%s.  %s" % (secs["%s.%s" % (m.group(1), m.group(2))]
                                .split("-")[1], m.group(3)),
                size=10, italic=True)
            i += 1
            continue
        m = re.match(r"^##\s+(\d+)\.\s+(.*)$", s)
        if m:
            p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 10, space_before=8,
                     space_after=3)
            run(p, "%s.  %s" % (secs[m.group(1)], m.group(2)), size=10,
                caps=True)
            i += 1
            continue

        m = re.match(r"^!\[[^\]]*\]\(([^)]+\.png)\)$", s)
        if m:
            rel = m.group(1)
            src = re.search(r"fig(\d+)", os.path.basename(rel)).group(1)
            j = i + 1
            while j < len(body) and not body[j].strip():
                j += 1
            cap = ""
            if j < len(body) and body[j].strip().startswith("**Figure"):
                cap = re.sub(r"^\*\*Figure\s+\d+\.\s*", "**",
                             body[j].strip())
                i = j
            continuous(doc, 1)
            p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, space_after=2)
            p.add_run().add_picture(os.path.join(ROOT, rel),
                                    width=Inches(7.1))
            p = para(doc, WD_ALIGN_PARAGRAPH.JUSTIFY, SMALL_PT, space_after=8)
            run(p, "Fig. %d.  " % figs[src], size=SMALL_PT)
            emit(p, rewrite_refs(cap, figs, secs), size=SMALL_PT)
            continuous(doc, 2)
            i += 1
            continue

        if s.startswith("> "):
            p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, 10, space_before=4,
                     space_after=4)
            run(p, s[2:].strip().strip("*").strip(), size=10, italic=True)
            i += 1
            continue

        if s.startswith("```"):
            j = i + 1
            while j < len(body) and not body[j].strip().startswith("```"):
                p = para(doc, WD_ALIGN_PARAGRAPH.LEFT, 7, space_after=0)
                run(p, body[j].rstrip(), size=7, mono=True)
                j += 1
            i = j + 1
            continue

        if s.startswith("|"):
            j = i
            rows = []
            while j < len(body) and body[j].strip().startswith("|"):
                rows.append(body[j])
                j += 1
            table_n += 1
            add_table(doc, rows, table_n, figs, secs)
            i = j
            continue

        if re.match(r"^\d+\.\s", s) or s.startswith("- "):
            n = 0
            while i < len(body) and (re.match(r"^\d+\.\s", body[i].strip())
                                     or body[i].strip().startswith("- ")):
                t = body[i].strip()
                numbered = bool(re.match(r"^\d+\.\s", t))
                n += 1
                t = re.sub(r"^(\d+\.|-)\s+", "", t)
                p = para(doc, WD_ALIGN_PARAGRAPH.JUSTIFY, space_after=2)
                p.paragraph_format.left_indent = Inches(0.2)
                p.paragraph_format.first_line_indent = Inches(-0.12)
                run(p, ("%d) " % n) if numbered else "•  ")
                emit(p, rewrite_refs(t, figs, secs))
                i += 1
            continue

        p = para(doc, space_after=3, indent=0.2)
        emit(p, rewrite_refs(s, figs, secs))
        i += 1


def add_table(doc, rows, number, figs, secs):
    header = [c.strip() for c in rows[0].strip().strip("|").split("|")]
    data = [[c.strip() for c in r.strip().strip("|").split("|")]
            for r in rows[2:]]
    wide = len(header) >= 4 or max(
        len(" ".join(r)) for r in [header] + data) > 70
    if wide:
        continuous(doc, 1)
    p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, SMALL_PT, space_before=6,
             space_after=1)
    run(p, "TABLE %s" % ROMAN[number], size=SMALL_PT)
    p = para(doc, WD_ALIGN_PARAGRAPH.CENTER, SMALL_PT, space_after=3)
    run(p, TABLE_CAPTIONS[number - 1] if number <= len(TABLE_CAPTIONS)
        else "Results.", size=SMALL_PT, caps=True)

    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for k, h in enumerate(header):
        cell = t.rows[0].cells[k]
        cell.text = ""
        pp = cell.paragraphs[0]
        pp.paragraph_format.space_after = Pt(0)
        pp._size = SMALL_PT
        run(pp, re.sub(r"\*\*", "", h), bold=True, size=SMALL_PT)
    for r in data:
        cells = t.add_row().cells
        for k in range(len(header)):
            pp = cells[k].paragraphs[0]
            pp.paragraph_format.space_after = Pt(0)
            pp._size = SMALL_PT
            emit(pp, rewrite_refs(r[k] if k < len(r) else "", figs, secs),
                 size=SMALL_PT)
    para(doc, space_after=6, size=4)
    if wide:
        continuous(doc, 2)


if __name__ == "__main__":
    main()
