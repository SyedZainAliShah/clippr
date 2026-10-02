"""Generate the CLIPPR Colab notebook.

Written as a generator rather than hand-edited JSON so cell order, titles and the form
metadata stay consistent, and so the whole notebook can be regenerated after an API change
instead of patched.
"""
from pathlib import Path

from notebook_kit import REPO, code, colab_url, md, setup_cell, write_notebook

OUT = Path("notebooks/CLIPPR_designer.ipynb")
NAME = "CLIPPR_designer.ipynb"
COLAB = colab_url(NAME)
INVENTORY_COLAB = colab_url("CLIPPR_inventory.ipynb")
DIAGRAM = Path("notebooks/pipeline.svg")
#: Measured, not recalled: the sum of `pytest --collect-only -q` on 2026-10-01.
#: Re-measure when tests are added or removed. A badge nothing recomputes goes stale,
#: and this one read 283 for long enough to be wrong by 644.
TESTS = 927
#: Served from the repository rather than inlined -- see pipeline_svg().
DIAGRAM_URL = f"https://raw.githubusercontent.com/{REPO}/main/{DIAGRAM.as_posix()}"


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

**Design a PPR protein that binds any RNA sequence you choose**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB})
[![License: MIT](https://img.shields.io/badge/License-MIT-1a7f5a.svg)](https://github.com/{REPO}/blob/main/LICENSE)
[![Tests](https://img.shields.io/badge/tests-{TESTS}%20passing-1a7f5a.svg)](https://github.com/{REPO})

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

**Four steps.** Run **0 · Setup** once, edit **1**, then run **2**, **3** and **4** in order.
Everything under *Optional extras* can be skipped — open one when you want that check.
Code is hidden behind the forms; use a cell's ⋮ menu → **Form → Show code** to read it.

> **Designing from the deposited 42-module kit instead?** That is a different workflow with a
> different set of answers, so it has its own notebook:
> [CLIPPR_inventory.ipynb]({INVENTORY_COLAB}).

> ##### Before you read any number
> **Predicted fidelity** comes from published ligation-count matrices (Pryor *et al.* 2020) —
> it is not a measured assembly efficiency in your hands. **QC** is a sequence-complexity
> check, not calibrated against vendor outcomes.
> Nothing here has been validated at the bench.
"""))

# ---------------------------------------------------------------- install
cells.append(md("""
## 0 · Setup

Run once — installs CLIPPR from GitHub, about a minute.
"""))

cells.append(setup_cell())

# ---------------------------------------------------------------- parameters
cells.append(md("""
## 1 · Your target and host

The only cell you have to edit. Everything else runs on what you set here.
"""))

cells.append(code('''
#@markdown # 1 · What should it bind?
#@markdown ---
#@markdown The RNA sequence your PPR will recognise. **Its length sets the architecture** —
#@markdown 9, 14 or 19 bases give a 9S, 14S or 19S protein.
#@markdown
#@markdown ⏱️ **Length also sets how long Step 2 takes, and the difference is large.** Measured
#@markdown once each on the project's own machine, screening off: **9 bases ≈ 0.5 s**,
#@markdown **14 bases ≈ 14 s**, **19 bases ≈ 3—4 minutes**. The codon optimiser has to route a
#@markdown longer coding sequence around more fixed junctions and more forbidden enzyme sites,
#@markdown so the cost climbs far faster than the length does. A 19-mer has not hung — it is
#@markdown working.
target_rna = "AAAAUGUGG"  #@param {type:"string"}

#@markdown # 2 · Where will it be expressed?
#@markdown ---
#@markdown Sets the codon usage and the genetic code (nuclear hosts use table 1,
#@markdown chloroplasts table 11).
#@markdown
#@markdown ⚠️ **The two Chlamydomonas entries are not interchangeable.** The nucleus is
#@markdown GC-rich and prefers Leu `CTG`; the chloroplast is AT-rich and prefers Leu `TTA`.
#@markdown Using one for the other produces DNA that looks fine and is wrong.
organism = "c_reinhardtii_nuclear"  #@param ["c_reinhardtii_nuclear", "c_reinhardtii_chloroplast", "e_coli", "s_cerevisiae", "a_thaliana_nuclear", "n_tabacum_chloroplast"]

#@markdown ### Your own codon usage — optional
#@markdown Leave all three blank to use the host above. Highest filled one wins.
#@markdown
#@markdown **`codon_table_file`** — a `codon,frequency` CSV or a CDS FASTA. Upload it with
#@markdown the folder icon in the sidebar, then put the filename here.
#@markdown
#@markdown **`kazusa_taxid`** — any NCBI taxonomy id Kazusa carries, e.g. `4577` for maize.
#@markdown Fetched and cached on first use.
#@markdown
#@markdown **`genetic_code_override`** — the genetic code for a table you supplied. Leave 0
#@markdown to inherit from the host; set 11 for anything organellar.
codon_table_file = ""  #@param {type:"string"}
kazusa_taxid = 0  #@param {type:"integer"}
genetic_code_override = 0  #@param {type:"integer"}

#@markdown # 3 · How will it be assembled?
#@markdown ---
#@markdown Which Type IIS enzyme cuts the fragments out, and which published mis-ligation
#@markdown table scores the junctions. `BsaI-HFv2` and `BbsI-HF` are the two measured in
#@markdown Pryor *et al.* 2020 at 25 °C over 18 h.
assembly_enzyme = "BsaI"  #@param ["BsaI", "BbsI", "BsmBI", "SapI"]
ligation_table = "BsaI-HFv2"  #@param ["BsaI-HFv2", "BbsI-HF"]

#@markdown ### Which enzyme sites must be absent
#@markdown `assembly` — this assembly's own chemistry only &nbsp;·&nbsp;
#@markdown `igem_rfc1000` — adds SapI, required by iGEM's Type IIS standard &nbsp;·&nbsp;
#@markdown `moclo_compat` — adds BsmBI to keep later MoClo levels open, a preference that can
#@markdown make some junctions infeasible.
#@markdown
#@markdown **`extra_blacklist`** — anything else this experiment needs kept clear,
#@markdown comma-separated. Any name Biopython knows, for example `EcoRI, BamHI, HindIII, NotI`.
enzyme_profile = "igem_rfc1000"  #@param ["assembly", "igem_rfc1000", "moclo_compat"]
extra_blacklist = ""  #@param {type:"string"}

#@markdown ### Destination vector
#@markdown Or type your own acceptor overhangs below as `5prime,3prime` coding sites, which
#@markdown override the level. The 3' entry is the **coding site**; the enzyme leaves its
#@markdown reverse complement.
destination_level = "level0"  #@param ["level_minus1", "level0", "level1"]
custom_destination = ""  #@param {type:"string"}

#@markdown # 4 · Anything else
#@markdown ---
#@markdown **`n_fragments`** — how many pieces to split the gene into. Leave at 0 to let the
#@markdown length decide; set it only if your vendor has an awkward limit.
#@markdown
#@markdown **`seed`** — the same seed always gives the same design.
#@markdown
#@markdown **`write_files`** writes the design files to disk, and **`check_offtarget`** checks
#@markdown the target against the host genome.
n_fragments = 0  #@param {type:"integer"}
seed = 42  #@param {type:"integer"}
write_files = True  #@param {type:"boolean"}
check_offtarget = True  #@param {type:"boolean"}
''', title="1 · Your target and host  (edit this)"))

# ---------------------------------------------------------------- design + results
cells.append(md("""
## 2 · Design it

Run it. Nothing to edit here.
"""))

cells.append(code('''
#@markdown Picks cut positions and Golden Gate overhangs **first**, then codon-optimises with
#@markdown those positions locked — optimising first would let the optimiser rewrite the very
#@markdown bases the junctions depend on.
from clippr import DESTINATION_OVERHANGS, design_oneshot, table_from_kazusa
from IPython.display import HTML, display

result = design_oneshot(
    target_rna,
    organism=organism,
    codon_table=(codon_table_file or
                 (table_from_kazusa(kazusa_taxid) if kazusa_taxid else None)),
    genetic_code=genetic_code_override or None,
    enzyme_profile=enzyme_profile,
    extra_blacklist=extra_blacklist,
    enzyme=assembly_enzyme,
    matrix=ligation_table,
    destination=(tuple(s.strip().upper() for s in custom_destination.split(",")[:2])
                 if custom_destination else DESTINATION_OVERHANGS[destination_level]),
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


# Labelled, because an unlabelled TNTNND... under "coding sequence" reads as the sequence.
ppr_code = result["ppr_code"]
if len(ppr_code) > 22:
    ppr_code = ppr_code[:22] + "…"
ppr_code = f"PPR code {ppr_code}"

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
    f'</div>'
    f'<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr))">'
    f'{cards}</div></div>{warn}')

display(HTML(card))
''', title="2 · Design it"))

# ---------------------------------------------------------------- fragments
cells.append(md("""
## 3 · The fragments to order

One row per orderable piece. `oh5` and `oh3` are the four-base Golden Gate overhangs that
join each fragment to its neighbours.
"""))

cells.append(code('''
cols = {"fragment_id": "fragment", "assembly_order": "order", "aa_length": "residues",
        "oligo_length": "oligo nt", "oh5_coding_site_5to3": "oh5",
        "oh3_coding_site_5to3": "oh3"}
table = result["oligos"][list(cols)].rename(columns=cols)

# pandas' .style needs jinja2, which Colab has but a bare local environment may not.
# Test that directly rather than catching the AttributeError pandas raises, which would also
# swallow real errors. Falls back to the plain frame: cosmetics should never break a cell.
try:
    import jinja2  # noqa: F401
    display(table.style.hide(axis="index").set_properties(
        subset=["oh5", "oh3"], **{"font-family": "ui-monospace, monospace"}).set_table_styles([
            {"selector": "th", "props": [("text-align", "left"), ("font-size", "11px"),
                                         ("letter-spacing", ".07em"),
                                         ("text-transform", "uppercase"),
                                         ("opacity", ".65"), ("padding", ".4em .9em")]},
            {"selector": "td", "props": [("padding", ".4em .9em"),
                                         ("font-variant-numeric", "tabular-nums")]}]))
except ImportError:
    display(table)
''', title="3 · The fragments to order"))

# ---------------------------------------------------------------- audit

# ---------------------------------------------------------------- download
cells.append(md("""
## 4 · Take the files

The order CSV, the oligos as FASTA, the assembled gene, and an annotated GenBank record —
every PPR repeat labelled with the base it reads — that opens in Benchling or SnapGene.
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
''', title="4 · Download the files"))

# ---------------------------------------------------------------- optional extras
cells.append(md("""
---

## Optional extras

*Nothing below is needed to place an order.* Use the arrow beside this heading to collapse
the whole section, or open any single item.
"""))

# ---------------------------------------------------------------- upload a codon table
cells.append(md("""
### Upload a codon table

Only if you want to supply your own codon usage from a file rather than typing a path.
"""))

cells.append(code('''
#@markdown Optional. Run this only if you want to upload a codon table from your computer
#@markdown rather than type a path. It puts the file in the runtime and fills in the
#@markdown filename for you — then re-run the Design cell.
try:
    from google.colab import files as _f
    _up = _f.upload()
    if _up:
        codon_table_file = list(_up)[0]
        print(f"using {codon_table_file}")
except ImportError:
    print("not running on Colab — put the file beside the notebook and give its path "
          "in codon_table_file instead.")
''', title="Upload a codon table"))

# ---------------------------------------------------------------- off-target
cells.append(md("""
### Does this target already occur in the host?

This looks for your target as an **exact match** in the *Chlamydomonas* chloroplast reference.
A match means the host already carries that sequence — a reason to look closer, **not** a
prediction that your PPR will bind there. Nothing here measures affinity.

**Only transcripts can matter.** A match in non-transcribed DNA is not an RNA target, and
neither is a reverse-complement match — the transcript from that locus carries the other
sequence. Both tiers are printed separately.

**What was searched.** 109 annotated features of one reference. Nine gene names appear
**twice** — `psbA` among them, as two identical copies at different places in the genome — so
a gene name does not identify a single location, and the CDS entries are coding spans rather
than UTR-inclusive transcripts. Your construct and the other members of your library are
**not** searched; for those, see the cross-talk cells below. Background rates and how they
were measured are in `CLIPPR_FOR_THE_WET_LAB.md`.

> The target this notebook ships with, `AAAAUGUGG`, **does** occur in the host, deliberately,
> so a first run shows you what a match looks like rather than a blank result.
"""))

cells.append(code('''
from clippr.offtarget import architecture_advice, load_genome, load_transcripts, report, scan

genome = load_genome()
transcripts = load_transcripts()

# The answer first. A background-rate table printed above the result made readers ask what
# the statistic was for before they had seen whether their own target was flagged.
print(report([scan(target_rna, genome, transcripts)]))

# Length is the lever that fixes a flagged target, so it belongs beside the result. Name the
# quantity: this is a GENOMIC, both-strand figure, while the verdict above is about
# transcripts. Using one to dismiss the other is a confusion a reviewer already hit here.
expected = architecture_advice(genome).get(len(target_rna))
print()
if expected is not None:
    print(f"For scale: an average {len(target_rna)}-base sequence occurs about "
          f"{expected:.1f} time(s) by chance")
    print("in this genome, counting both strands of the DNA. That is a different quantity")
    print("from the transcript verdict above and does not explain it away.")
    print()
    print("Length is what buys specificity: across 200 random trials a 14-base target")
    print("occurred nowhere in this genome, in either tier.")

print()
print(f"host: Chlamydomonas reinhardtii chloroplast, {len(genome):,} bp, "
      f"{len(transcripts)} annotated transcripts")
''', title="Does this target occur in the host?"))

cells.append(md("""
---

### Do you already own the parts?

Everything above designs DNA to be **synthesised**. The GRASP authors also deposited a
42-plasmid kit, and a lab holding it can assemble many PPRs from parts it already has.

| route | what you supply | what it constrains |
|---|---|---|
| **de novo synthesis** | new DNA | nothing — full synonymous freedom |
| **GRASP module kit** | parts you already hold | fixed to the deposited modules |

An *n*-base target needs *n*+1 modules, and each internal join consumes one linker pair. The
kit holds exactly three, so 19 bases is the longest it can build — the cell below says so
plainly when asked for more.

> Module identity, overhangs and plate positions come from **Dennis et al. 2025 Table S1**.
> This selects PPR modules only, not the acceptor plasmids.
"""))

cells.append(code('''
from clippr import parts_report, select_parts

plan = select_parts(target_rna)
print(parts_report(plan))
''', title="Build it from parts you own"))

cells.append(md("""
---

### Designing a whole library

For a set of regulators what matters is **orthogonality**: PPRᵢ must bind UTRᵢ and not UTRⱼ.
The matrix below is the pairwise distance between targets — larger is better separated.

> This section is a **capability, not a validated result**. Every other part of this package
> is checked against a 200-design corpus; this one has unit tests only, because no ground
> truth for it exists.
"""))

cells.append(code('''
#@markdown **These are three example targets, not derived from your target above.** Replace
#@markdown them with your own set. The point of this cell is the *pairwise* question — whether
#@markdown several regulators interfere with each other — which needs more than one target, so
#@markdown it cannot reuse the single target from the top of the notebook.
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
''', title="Design a whole library"))

cells.append(code('''
lib.qc_table()
''', title="Library QC table"))

cells.append(md("""
### DNA shared between members

Cross-talk asks whether two PPRs could bind each other's **target**. This asks whether two
**genes** share enough identical DNA to recombine — a different question. Every member
carries the same scaffold, so the answer is never trivially no.

`diversify_library` gives each member its own synonymous encoding of that scaffold. On five
9S designs it took the longest shared stretch from 107 nt to 47 and cleared every pair over
the 50 nt threshold; at six members it reaches only 77 nt and leaves 2 of 15 pairs above it.

A shared stretch is a *necessary* substrate for recombination, never a prediction that it
will happen, and 50 nt is a rule of thumb rather than a measured constant for this host.
"""))

cells.append(code('''
print(lib.homology())
''', title="DNA shared between members"))

cells.append(md("""
### Cross-talk, in two tiers

Sequence separation and *predicted binding* are different questions, so they stay apart:

| tier | what it is | status |
|---|---|---|
| **A — Hamming distance** | two targets differ in *k* of *n* positions | **the only thing that gates** |
| **B — predicted affinity** | the PPR designed for A, scored against B | an annotation; gates nothing |

Tier B needs a PPR specificity table. **CLIPPR does not ship one** — the available table
carries no licence, and it comes from **P-type** experiments while this scaffold is
**S-type**. A high score means *look*, never *fail*.

Leave the path blank and you get tier A alone, which is what gates anyway.
"""))

cells.append(code('''
ppr_score_table = ""  #@param {type:"string"}

from clippr import compare_tiers, load_ppr_scores

scores = load_ppr_scores(ppr_score_table) if ppr_score_table.strip() else None
print(lib.crosstalk(scores=scores))

if scores:
    cmp = compare_tiers(lib.targets, scores)
    print()
    print(f"pairs flagged by separation:      {len(cmp['hamming_flagged'])}")
    print(f"pairs flagged by predicted affinity: {len(cmp['affinity_flagged'])}")
    print(f"the model points somewhere separation does not: {cmp['tiers_disagree']}")
    print()
    print("A disagreement is a reading recommendation, not a failed design —")
    print("tier A alone decides what this library accepts.")
''', title="Cross-talk between targets"))


cells.append(md("""
---

## Going deeper — why this design, and not another

*You do not need this to order. Read it if a design looks surprising.*

The gene is cut into fragments, and each place two fragments rejoin is a **junction**. Every
junction needs a 4-base sticky end, the **overhang**, so the right two fragments anneal to
each other and to nothing else.

The constraint that makes this hard: the overhang is made of **your own coding sequence** at
that point, not a linker bolted on. So the only overhangs available are the ones synonymous
codons can spell there while leaving the protein unchanged.

Below, each junction lists what it could have used and what became of it — **selected**,
**considered** (feasible, but another scored at least as well), or **rejected** (no synonymous
arrangement avoids an excluded enzyme site, so that overhang is impossible here, not merely
worse).
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
''', title="Going deeper: the design audit"))

cells.append(md(f"""
---

## What each step does

| Step | What happens |
|---|---|
| **0 · Setup** | installs CLIPPR from GitHub |
| **1 · Target and host** | your target RNA, codon table, enzyme and destination |
| **2 · Design it** | target to protein, then cuts and overhangs chosen **before** the coding sequence is optimised, so the optimiser cannot rewrite the bases the junctions depend on |
| **3 · Fragments** | the oligo table you order |
| **4 · Files** | FASTA, GenBank and the oligo CSV |
| *Optional extras* | host occurrence check, the deposited-kit route, a whole library, cross-talk |
| *Going deeper* | every overhang considered at every junction, and why each was kept or ruled out |

**Predicted fidelity** is from published ligation-count matrices (Pryor *et al.* 2020), not a measured efficiency in your hands. **QC** is a sequence-complexity check, not calibrated against vendor outcomes. No sequence from this project has been synthesised.

Building from the deposited 42-module kit instead? [CLIPPR_inventory.ipynb]({INVENTORY_COLAB})
"""))

DIAGRAM.write_text(pipeline_svg(), encoding="utf-8")
print(f"wrote {DIAGRAM} ({DIAGRAM.stat().st_size} B)")
write_notebook(cells, OUT, NAME)
