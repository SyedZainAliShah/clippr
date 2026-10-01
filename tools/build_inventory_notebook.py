"""Generate the CLIPPR reusable-inventory notebook.

Split out of the designer notebook. It is a second complete workflow for a different reader --
someone who already holds the deposited GRASP module kit -- and appending it underneath the
designer meant scrolling past a workflow you were not using to reach it. Reviewers reported
exactly that confusion.
"""
from pathlib import Path

from notebook_inventory_cells import cells as inventory_cells
from notebook_kit import code, colab_url, md, setup_cell, write_notebook

OUT = Path("notebooks/CLIPPR_inventory.ipynb")
NAME = "CLIPPR_inventory.ipynb"

cells = [
    md(f"""
<div align="center">

# CLIPPR — reusable inventory

### Build PPRs from the deposited GRASP module kit instead of ordering new DNA

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({colab_url(NAME)})

**iGEM Marburg 2026**

</div>

---

**Use this notebook if you hold the kit.** If you want a PPR designed and synthesised from
scratch, use `CLIPPR_designer.ipynb` instead — that is the other route, and it needs no kit.

Recode the kit once for your host, then *compile* targets from it: the DNA is ordered once
and reused. The kit can only spell the targets its modules cover, which is the trade you are
making against full synthesis.

**You need Supplementary Table S1** (`.xlsx`). It holds the module insert sequences and is not
distributed with this package — see `NOTICE.md`. Upload it at step 1. Every step after the
recode works from a *saved inventory* carrying its own sequences, so Table S1 is needed once.

**Run the cells top to bottom.** Each step depends on the one before it.
"""),
    setup_cell(),
    code('''
#@markdown Checks the installed build actually carries the inventory workflow. It arrived
#@markdown after the designer route, so an older published build will not have it, and an
#@markdown ImportError twelve cells down is a worse way to find that out.
import importlib.util

_needed = ["clippr.inventories", "clippr.workflow", "clippr.ordering", "clippr.substrates"]
_missing = [m for m in _needed if importlib.util.find_spec(m) is None]

from IPython.display import HTML, display

if _missing:
    display(HTML(
        '<div style="border-left:3px solid #a8402c;padding:.6em 1em;font-size:13px;'
        'line-height:1.55">The installed CLIPPR does not carry the reusable-inventory '
        'workflow.<br>Missing: <code>' + ", ".join(_missing) + '</code><br><br>'
        'Nothing below this cell will run until the published build includes it. The '
        'designer notebook is unaffected.</div>'))
else:
    display(HTML(
        '<div style="border-left:3px solid #1a7f5a;padding:.6em 1em;'
        'font-family:ui-monospace,monospace;font-size:13px;opacity:.85">'
        'inventory workflow available</div>'))
''', title="Check this build has the inventory workflow"),
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
''', title="Host — edit these"),
]

cells += inventory_cells(md, code)
write_notebook(cells, OUT, NAME)
