"""Generate the PPR-opt reusable-inventory notebook.

Split out of the designer notebook. It is a second complete workflow for a different reader --
someone who already holds the deposited GRASP module kit -- and appending it underneath the
designer meant scrolling past a workflow you were not using to reach it. Reviewers reported
exactly that confusion.
"""
from pathlib import Path

from notebook_inventory_cells import cells as inventory_cells
from notebook_kit import (code, colab_url, glossary, md, setup_cell,
                          write_notebook)

OUT = Path("notebooks/PPROPT_inventory.ipynb")
NAME = "PPROPT_inventory.ipynb"
DESIGNER_COLAB = colab_url("PPROPT_designer.ipynb")

cells = [
    md(f"""
<div align="center">

# PPR-opt — reusable inventory

### Build PPRs from the deposited GRASP module kit instead of ordering new DNA

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({colab_url(NAME)})

**iGEM Marburg 2026**

</div>

---

**Use this notebook if you hold the kit.** If you want a PPR designed and synthesised from
scratch, use [PPROPT_designer.ipynb]({DESIGNER_COLAB}) instead — the other route, which
needs no kit.

Recode the kit once for your host, then *compile* targets from it: the DNA is ordered once
and reused. The kit can only spell the targets its modules cover, which is the trade you are
making against full synthesis.

**You need Supplementary Table S1** (`.xlsx`). It holds the module insert sequences and is not
distributed with this package — see `NOTICE.md`. Upload it at step 1. Every step after the
recode works from a *saved inventory* carrying its own sequences, so Table S1 is needed once.

**Run the cells top to bottom.** Each step depends on the one before it.
"""),
    setup_cell(),
    # Why this cell exists, kept here rather than in the form text above it: the inventory
    # workflow arrived after the designer route, so an older published build imports fine and
    # then raises ImportError twelve cells down. That is a maintainer's concern. The reader
    # wants to know whether they can carry on, and what to do if not.
    code('''
#@markdown Checks this notebook can run on the copy of PPR-opt you have installed.
#@markdown
#@markdown **Green?** Carry on. **Red?** Run the Setup cell above again, then
#@markdown **Runtime → Restart session**, then Setup once more.
import importlib.util

_needed = ["ppropt.inventories", "ppropt.workflow", "ppropt.ordering", "ppropt.substrates"]
_missing = [m for m in _needed if importlib.util.find_spec(m) is None]

from IPython.display import HTML, display

if _missing:
    display(HTML(
        '<div style="border-left:3px solid #a8402c;padding:.6em 1em;font-size:13px;'
        'line-height:1.55">The installed PPR-opt does not carry the reusable-inventory '
        'workflow.<br>Missing: <code>' + ", ".join(_missing) + '</code><br><br>'
        'Nothing below this cell will run until the published build includes it. The '
        'designer notebook is unaffected.</div>'))
else:
    display(HTML(
        '<div style="border-left:3px solid #1a7f5a;padding:.6em 1em;'
        'font-family:ui-monospace,monospace;font-size:13px;opacity:.85">'
        'inventory workflow available</div>'))
''', title="Setup check — does this build have the inventory workflow?"),
    code('''
#@markdown The host the kit is recoded for. The genetic code follows automatically —
#@markdown nuclear hosts use table 1, chloroplasts table 11.
organism = "c_reinhardtii_nuclear"  #@param ["c_reinhardtii_nuclear", "c_reinhardtii_chloroplast", "e_coli", "s_cerevisiae", "a_thaliana_nuclear", "n_tabacum_chloroplast"]

#@markdown Optional: your own codon usage as a `codon,frequency` CSV or a CDS FASTA.
#@markdown Upload it with the folder icon in the sidebar and put the filename here.
codon_table_file = ""  #@param {type:"string"}

#@markdown Optional: the genetic code for a table you supplied. Leave 0 to inherit from the
#@markdown host. Set 11 for anything organellar.
genetic_code_override = 0  #@param {type:"integer"}
''', title="Before you start — your host  (edit this)"),
]

cells += inventory_cells(md, code)
cells.append(glossary(md, [
    ("deposited kit", "The 42 published modules you physically hold. This notebook builds "
                      "from those rather than designing new DNA."),
    ("module", "One part in that kit: a stretch of DNA with fixed ends, ready to join to its "
               "neighbours."),
    ("interface", "The four-base ends by which modules join. They are held fixed while "
                  "recoding, so a rewritten module still drops into an assembly built from "
                  "untouched ones."),
    ("recode", "Rewrite a module's DNA for your host's preferred codons. The protein, the "
               "length and both interfaces stay exactly as they were."),
    ("inventory", "Your saved set of modules after recoding. Every later step works from "
                  "this rather than from Table S1 again."),
    ("pool", "Oligos ordered together in a single tube, which is cheaper than ordering each "
             "one separately."),
    ("front", "The set of candidates where nothing is strictly better than anything else — "
              "improving one property costs you another. A menu of honest trade-offs rather "
              "than a ranking, which is why you have to pick one yourself."),
    ("level 0", "The smallest assembly unit in this standard: a single part, wrapped so the "
                "enzyme releases it with the right ends. Level 0 parts combine into bigger "
                "constructs later."),
    ("digest check", "A simulated enzyme cut, run before export, confirming each ordered "
                     "piece really does release the fragment it claims to."),
    ("Supplementary Table S1", "The published spreadsheet holding the kit's module "
                               "sequences. It is not distributed with this package, so you "
                               "upload it yourself — see `NOTICE.md`."),
]))
write_notebook(cells, OUT, NAME)
