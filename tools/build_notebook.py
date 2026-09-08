"""Generate the CLIPPR Colab notebook.

Written as a generator rather than hand-edited JSON so cell order, titles and the form
metadata stay consistent, and so the whole notebook can be regenerated after an API change
instead of patched.
"""
import json
from pathlib import Path

OUT = Path("notebooks/CLIPPR_designer.ipynb")
REPO = "SyedZainAliShah/clippr"
COLAB = f"https://colab.research.google.com/github/{REPO}/blob/main/notebooks/CLIPPR_designer.ipynb"
DIAGRAM = Path("notebooks/pipeline.svg")
#: Served from the repository rather than inlined -- see pipeline_svg().
DIAGRAM_URL = f"https://raw.githubusercontent.com/{REPO}/main/{DIAGRAM.as_posix()}"


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text, title=None, form=True):
    meta = {"cellView": "form"} if (title and form) else {}
    body = (f'#@title {title} {{display-mode: "form"}}\n' if title else "") + text.strip("\n")
    return {"cell_type": "code", "execution_count": None, "metadata": meta,
            "outputs": [], "source": body.splitlines(True)}


cells = []


#: Single ink colour for the pipeline diagram, chosen to read on both Colab themes.
#: Contrast is 3.3:1 on white and 4.9:1 on Colab's #202124 -- above the 3:1 the WCAG
#: non-text threshold asks for, on both. `currentColor` cannot be used: Colab strips
#: inline SVG from markdown, so the diagram has to be an external image, and an image has
#: no inherited colour to take.
INK = "#7d918a"


def pipeline_svg():
    """The pipeline diagram, as a standalone SVG file.

    Written to `notebooks/pipeline.svg` and referenced by URL rather than inlined: Colab's
    markdown sanitiser removes inline <svg> entirely, so an inlined diagram renders as
    nothing at all.
    """
    stages = [
        ("ppr", "RNA -> protein"),
        ("arelf", "where to cut"),
        ("overhangs", "will it join?"),
        ("codons", "make it real"),
        ("assembly", "fragments"),
        ("qc", "worth ordering?"),
    ]
    x0, box_w, gap, y = 96, 118, 16, 34
    parts = []
    for i, (name, sub) in enumerate(stages):
        x = x0 + i * (box_w + gap)
        parts.append(
            f'<rect x="{x}" y="{y}" width="{box_w}" height="52" rx="4" fill="none" '
            f'stroke="{INK}" stroke-opacity=".5"/>'
            f'<text x="{x + box_w / 2}" y="{y + 21}" text-anchor="middle" '
            f'font-family="ui-monospace,monospace" font-size="13" font-weight="600" '
            f'fill="{INK}">{name}</text>'
            f'<text x="{x + box_w / 2}" y="{y + 38}" text-anchor="middle" '
            f'font-family="ui-sans-serif,system-ui" font-size="10.5" '
            f'fill="{INK}" fill-opacity=".72">{sub}</text>')
        if i < len(stages) - 1:
            ax = x + box_w + 3
            parts.append(f'<path d="M{ax} {y + 26} l9 0 m-3 -3 l3 3 l-3 3" fill="none" '
                         f'stroke="{INK}" stroke-opacity=".6" stroke-width="1.3"/>')

    end_x = x0 + len(stages) * (box_w + gap) - gap
    # Right margin sized to the widest right-hand label ("4-7 fragments" at 12.5px),
    # which overflowed a 96px margin and was clipped.
    right_margin = 120
    label = ('font-family="ui-sans-serif,system-ui" font-size="11" '
             f'fill="{INK}" fill-opacity=".85"')
    return (
        f'<?xml version="1.0" encoding="UTF-8"?>'
        f'<svg viewBox="0 0 {end_x + right_margin} 122" '
        f'width="{end_x + right_margin}" height="122" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'role="img" aria-label="CLIPPR pipeline: target RNA through six stages to '
        f'orderable DNA">'
        f'<text x="0" y="{y + 21}" {label}>target</text>'
        f'<text x="0" y="{y + 36}" font-family="ui-monospace,monospace" font-size="12.5" '
        f'fill="{INK}">AAAAUGUGG</text>'
        f'<path d="M78 {y + 26} l11 0 m-4 -3.5 l4 3.5 l-4 3.5" fill="none" '
        f'stroke="{INK}" stroke-opacity=".6" stroke-width="1.3"/>'
        + "".join(parts) +
        f'<text x="{end_x + 14}" y="{y + 21}" {label}>order</text>'
        f'<text x="{end_x + 14}" y="{y + 36}" font-family="ui-sans-serif,system-ui" '
        f'font-size="12.5" fill="{INK}">4-7 fragments</text>'
        f'<text x="{x0}" y="112" font-family="ui-sans-serif,system-ui" font-size="10.5" '
        f'fill="{INK}" fill-opacity=".7">'
        f'overhangs are chosen before the sequence is optimised, then locked into it'
        f'</text>'
        f'</svg>')


# ---------------------------------------------------------------- header
cells.append(md(f"""
<div align="center">

# CLIPPR

### Design a PPR protein that binds any RNA sequence you choose

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB})
[![License: MIT](https://img.shields.io/badge/License-MIT-1a7f5a.svg)](https://github.com/{REPO}/blob/main/LICENSE)
[![Tests](https://img.shields.io/badge/tests-283%20passing-1a7f5a.svg)](https://github.com/{REPO})

**iGEM Marburg 2026**

</div>

<div align="center">

<img src="{DIAGRAM_URL}" alt="CLIPPR pipeline: a target RNA passes through ppr, arelf, overhangs, codons, assembly and qc to become orderable DNA fragments" width="100%" style="max-width:980px">

</div>

---

Pentatricopeptide repeat proteins are built from tandem ~31-residue repeats, and **each
repeat reads exactly one RNA base** through two specificity residues:

| 5th + last residue | reads |     | 5th + last residue | reads |
|:---:|:---:|---|:---:|:---:|
| `T` `N` | **A** |  | `T` `D` | **G** |
| `N` `N` | **C** |  | `N` `D` | **U** |

So the protein is a deterministic function of your target — no catalogue to search. Give it
nine bases and you get a nine-repeat protein, a synthesisable coding sequence, a Golden Gate
assembly plan, and the fragments to order.

**Run the cells top to bottom.** Only the *Design parameters* cell normally needs editing.

> ##### Before you read any number
> **Predicted fidelity** comes from published ligation-count matrices (Pryor *et al.* 2020) —
> it is not a measured assembly efficiency in your hands. **QC** is a sequence-complexity
> check, not calibrated against vendor outcomes. **Cost** is a list price, not a quote.
> Nothing here has been validated at the bench.
"""))

# ---------------------------------------------------------------- install
cells.append(code(f'''
#@markdown Installs the package if it is not already available. Safe to re-run — it never
#@markdown reinstalls over a working copy.
REPO = "{REPO}"

try:
    import clippr
    _msg = f"clippr {{clippr.__version__}} already available"
except ImportError:
    token = None
    try:
        from google.colab import userdata          # private-repo fallback
        token = userdata.get("GITHUB_TOKEN")
    except Exception:
        pass
    url = (f"git+https://{{token}}@github.com/{{REPO}}.git" if token
           else f"git+https://github.com/{{REPO}}.git")
    %pip install --quiet $url
    import clippr
    _msg = f"installed clippr {{clippr.__version__}}"

from IPython.display import HTML, display
display(HTML(
    f'<div style="border-left:3px solid #1a7f5a;padding:.5em .9em;'
    f'font-family:ui-monospace,monospace;font-size:13px;opacity:.85">{{_msg}}</div>'))
''', title="Setup — install CLIPPR"))

# ---------------------------------------------------------------- parameters
cells.append(code('''
#@markdown ### Target
#@markdown The RNA sequence the PPR should bind. **Its length sets the architecture** —
#@markdown 9, 14 or 19 bases give a 9S, 14S or 19S protein.
target_rna = "AAAAUGUGG"  #@param {type:"string"}

#@markdown ### Host
organism = "c_reinhardtii_nuclear"  #@param ["c_reinhardtii_nuclear"]

#@markdown ### Which enzyme sites must be absent
#@markdown `assembly` — this assembly's own chemistry (BsaI, BbsI)
#@markdown &nbsp;&nbsp;·&nbsp; `igem_rfc1000` — adds SapI, required by iGEM's Type IIS standard
#@markdown &nbsp;&nbsp;·&nbsp; `moclo_compat` — adds BsmBI to keep later MoClo levels open,
#@markdown a preference that can make some junctions infeasible
enzyme_profile = "igem_rfc1000"  #@param ["assembly", "igem_rfc1000", "moclo_compat"]

#@markdown ### Destination vector level
destination_level = "level0"  #@param ["level_minus1", "level0", "level1"]

#@markdown ### Fragments
#@markdown How many pieces to split the gene into. Leave at 0 to let the length decide —
#@markdown set it only if your vendor has an awkward limit.
n_fragments = 0  #@param {type:"integer"}

#@markdown ### Reproducibility
#@markdown The same seed always gives the same design.
seed = 42  #@param {type:"integer"}
write_files = True  #@param {type:"boolean"}
check_offtarget = True  #@param {type:"boolean"}
''', title="Design parameters — edit these"))

# ---------------------------------------------------------------- design + results
cells.append(code('''
#@markdown Picks cut positions and Golden Gate overhangs **first**, then codon-optimises with
#@markdown those positions locked — optimising first would let the optimiser rewrite the very
#@markdown bases the junctions depend on.
from clippr import DESTINATION_OVERHANGS, design_oneshot
from IPython.display import HTML, display

result = design_oneshot(
    target_rna,
    organism=organism,
    enzyme_profile=enzyme_profile,
    destination=DESTINATION_OVERHANGS[destination_level],
    n_fragments=n_fragments or None,
    seed=seed,
    check_offtarget=check_offtarget,
    outdir="clippr_output" if write_files else None,
)

qc = result["qc"]
TONE = {"PASS": "#1a7f5a", "WARNING": "#9a6b1f", "FAIL": "#a8402c"}
tone = TONE.get(qc["status"], "#6b7b75")


def _stat(label, value, hint=""):
    return (
        f'<div style="padding:.55em .9em .6em;border-left:1px solid rgba(128,145,138,.35)">'
        f'<div style="font-size:10.5px;letter-spacing:.09em;text-transform:uppercase;'
        f'opacity:.6">{label}</div>'
        f'<div style="font-size:19px;font-weight:600;font-variant-numeric:tabular-nums;'
        f'margin-top:.15em">{value}</div>'
        f'<div style="font-size:11.5px;opacity:.6">{hint}</div></div>')


ppr_code = result["ppr_code"]
if len(ppr_code) > 24:
    ppr_code = ppr_code[:24] + "…"

cards = "".join([
    _stat("architecture", result["architecture"], f'{len(result["protein"])} aa protein'),
    _stat("coding sequence", f'{len(result["cds"])} nt', ppr_code),
    _stat("fragments", len(result["oligos"]),
          "cuts at " + ", ".join(str(c) for c in result["cuts"])),
    _stat("fidelity", f'{result["fidelity"]:.3f}', "predicted, not measured"),
    _stat("GC", f'{qc["gc_pct"]:.1f}%',
          f'windows {qc["gc_window_min"]:.0f}–{qc["gc_window_max"]:.0f}%'),
    _stat("repeats", f'{qc["repeated_kmer_fraction"]:.1%}',
          f'longest {qc["longest_repeat"]} nt'),
])

warn = "".join(
    f'<div style="border-left:3px solid {TONE["WARNING"]};padding:.5em .9em;'
    f'margin-top:.7em;font-size:13px;line-height:1.5">{w}</div>'
    for w in result["warnings"])

RULE = "1px solid rgba(128,145,138,.35)"
card = (
    f'<div style="font-family:ui-sans-serif,system-ui,sans-serif;border:{RULE};'
    f'border-radius:5px;overflow:hidden;max-width:920px">'
    f'<div style="display:flex;align-items:center;gap:.8em;padding:.7em 1em;'
    f'border-bottom:{RULE}">'
    f'<span style="font-family:ui-monospace,monospace;font-size:16px;font-weight:600">'
    f'{result["target_rna"]}</span>'
    f'<span style="background:{tone};color:#fff;font-size:11px;font-weight:700;'
    f'letter-spacing:.07em;padding:.2em .7em;border-radius:99px">{qc["status"]}</span>'
    f'<span style="margin-left:auto;font-size:12.5px;opacity:.65">'
    f'{result["cost"]["total_eur"]:.2f} EUR list price · not a quote</span></div>'
    f'<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr))">'
    f'{cards}</div></div>{warn}')

display(HTML(card))
''', title="Design — run this"))

# ---------------------------------------------------------------- fragments
cells.append(md("""
## The fragments to order

One row per orderable piece. `oh5` and `oh3` are the four-base Golden Gate overhangs that
join each fragment to its neighbours.
"""))

cells.append(code('''
cols = {"fragment_id": "fragment", "assembly_order": "order", "aa_length": "residues",
        "oligo_length": "oligo nt", "oh5_coding_site_5to3": "oh5",
        "oh3_coding_site_5to3": "oh3"}
table = result["oligos"][list(cols)].rename(columns=cols)

table.style.hide(axis="index").set_properties(
    subset=["oh5", "oh3"], **{"font-family": "ui-monospace, monospace"}).set_table_styles([
        {"selector": "th", "props": [("text-align", "left"), ("font-size", "11px"),
                                     ("letter-spacing", ".07em"), ("text-transform", "uppercase"),
                                     ("opacity", ".65"), ("padding", ".4em .9em")]},
        {"selector": "td", "props": [("padding", ".4em .9em"),
                                     ("font-variant-numeric", "tabular-nums")]}])
''', title="Fragment table"))

# ---------------------------------------------------------------- audit
cells.append(md("""
## Why this design, and not another

Every junction records the overhangs it *could* have used and what became of each:

- **selected** — the one used
- **considered** — feasible, but another scored at least as well
- **rejected** — no synonymous codon arrangement could avoid an excluded enzyme site, so
  that junction is *impossible* under the active profile, not merely worse

`local realizations` counts the synonymous arrangements still available around a junction.
It is reported, never used to choose — but a junction with 2 is more fragile than one with
16, and that is worth seeing before you order.

**If a design looks surprising, read this rather than trusting it.**
"""))

cells.append(code('''
audit = result["audit"]
print(audit.report())

rejected = audit.rejected
print(f"\\n{len(rejected)} candidate overhang(s) ruled out entirely under "
      f"profile '{audit.enzyme_profile}'")
for d in rejected[:8]:
    print(f"    {d.sequence}  cut {d.junction_cut:>4d}   {d.reason}")
if len(rejected) > 8:
    print(f"    … and {len(rejected) - 8} more")
if not rejected:
    print("    (every achievable overhang was usable at every junction)")
''', title="Design audit"))

# ---------------------------------------------------------------- download
cells.append(md("""
## Take the files

Four artefacts: the order CSV, the oligos as FASTA, the assembled gene, and an annotated
GenBank record — every PPR repeat labelled with the base it reads — that opens directly in
Benchling or SnapGene.
"""))

cells.append(code('''
import os
from IPython.display import HTML, display

if not result["paths"]:
    display(HTML('<div style="opacity:.7">Set <code>write_files</code> to True in the '
                 'parameters cell and re-run.</div>'))
else:
    try:
        from google.colab import files as colab_files
    except ImportError:
        colab_files = None

    LABEL = {"oligo_csv": "Order sheet (CSV)", "oligo_fasta": "Oligos (FASTA)",
             "gene_fasta": "Assembled gene (FASTA)", "genbank": "Annotated GenBank"}
    rows = "".join(
        f'<tr><td style="padding:.35em .9em">{LABEL.get(k, k)}</td>'
        f'<td style="padding:.35em .9em;font-family:ui-monospace,monospace;font-size:12px;'
        f'opacity:.7">{os.path.basename(p)}</td>'
        f'<td style="padding:.35em .9em;text-align:right;font-variant-numeric:tabular-nums;'
        f'opacity:.7">{os.path.getsize(p):,} B</td></tr>'
        for k, p in result["paths"].items())
    display(HTML(f'<table style="font-family:ui-sans-serif,system-ui,sans-serif;'
                 f'font-size:13px;border-collapse:collapse">{rows}</table>'))

    if colab_files:
        for p in result["paths"].values():
            colab_files.download(p)
''', title="Download the design files"))

# ---------------------------------------------------------------- library
# ---------------------------------------------------------------- off-target
cells.append(md("""
## Does this target also exist in the chloroplast?

A PPR cannot tell which copy of a sequence you meant. If your target also occurs in an
endogenous chloroplast transcript, the protein binds there too and stops being specific to
your construct.

**This is a length problem, and the numbers are stark.** The *Chlamydomonas* chloroplast
genome is 203,828 bases and 34.5% GC. Measured over 200 random targets of each length:

| target length | occur somewhere in the host |
|---|---|
| **9 nt** | **97 of 200 — 48%** |
| 14 nt | 0 of 200 |
| 19 nt | 0 of 200 |

A nine-base sequence is simply not rare enough in a 204 kb genome. If you need a 9S design,
check it; if a target comes back flagged, lengthening it is the reliable fix.

Occurrence is a *necessary* condition for off-target binding, not a sufficient one — this
reports sequence, not affinity.
"""))

cells.append(code('''
from clippr.offtarget import architecture_advice, load_genome, report, scan

genome = load_genome()
gc = 100 * (genome.count("G") + genome.count("C")) / len(genome)
print(f"host: Chlamydomonas reinhardtii chloroplast, {len(genome):,} bp, {gc:.1f}% GC")
print()

print("expected occurrences by chance, for an average target:")
for n, e in architecture_advice(genome).items():
    print(f"  {n:>2}-nt target : {e:8.3f}")

print()
print(report([scan(target_rna, genome, max_mismatches=1)]))
''', title="Off-target check — does the host already contain this sequence?"))

cells.append(md("""
---

## Designing a whole library

For a set of regulators, what matters is **orthogonality**: PPRᵢ must bind UTRᵢ and not
UTRⱼ. The matrix below is the pairwise distance between targets — larger is better
separated. Targets that sit close together risk one PPR binding another's UTR.

> `orthogonal.py` is a **capability, not a validated result**. Every other part of this
> package is checked against a 200-design corpus; this one has unit tests only, because no
> ground truth for it exists.
"""))

cells.append(code('''
targets = "AAAAUGUGG, GCUAAAGAC, UUACACGUG"  #@param {type:"string"}

from clippr import design_library

target_list = [t.strip().upper() for t in targets.split(",") if t.strip()]
lib = design_library(target_list, codon_table=None, organism=organism,
                     enzyme_profile=enzyme_profile, seed=seed,
                     check_offtarget=check_offtarget,
                     outdir="clippr_library" if write_files else None,
                     on_progress=lambda i, n, t: print(f"  {i}/{n}  {t}", flush=True))

print()
print(lib.summary())
print()
print(lib.crosstalk())
''', title="Design the whole library"))

cells.append(code('''
lib.qc_table()
''', title="Library QC table"))

nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True, "name": "CLIPPR_designer.ipynb"},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}
DIAGRAM.write_text(pipeline_svg(), encoding="utf-8")
OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"wrote {DIAGRAM} ({DIAGRAM.stat().st_size} B)")
print(f"wrote {OUT}  ({len(cells)} cells: "
      f"{sum(1 for c in cells if c['cell_type']=='code')} code, "
      f"{sum(1 for c in cells if c['cell_type']=='markdown')} markdown)")
