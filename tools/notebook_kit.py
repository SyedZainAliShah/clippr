"""Shared cell constructors for the PPR-opt notebooks.

Both notebooks are generated rather than hand-edited so cell order, titles and Colab form
metadata stay consistent, and so either can be regenerated after an API change instead of
patched. `build_notebook.py` builds the designer; `build_inventory_notebook.py` builds the
reusable-inventory route.
"""
import json
from pathlib import Path

REPO = "SyedZainAliShah/clippr"


def colab_url(name: str) -> str:
    return f"https://colab.research.google.com/github/{REPO}/blob/main/notebooks/{name}"


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text, title=None, form=True):
    meta = {"cellView": "form"} if (title and form) else {}
    body = (f'#@title {title} {{display-mode: "form"}}\n' if title else "") + text.strip("\n")
    return {"cell_type": "code", "execution_count": None, "metadata": meta,
            "outputs": [], "source": body.splitlines(True)}


#: The install cell. Identical in both notebooks, so it lives here: a fix to the private-repo
#: fallback or the version banner should never land in one notebook and not the other.
_SETUP = '''
#@markdown Installs PPR-opt or updates it when this notebook is newer than the copy you
#@markdown already have. Safe to re-run. Will take about a minute the first time, seconds after that.
#@markdown
#@markdown **If it fails:** run it again - a dropped network call is the usual cause. If it
#@markdown reports that it updated, use **Runtime → Restart session** and run it once more:
#@markdown the old copy stays loaded in memory until you do. Every cell below needs this one
#@markdown to have succeeded, so do not skip past a red error here.
REPO = "__REPO__"

import importlib.metadata as _md, json as _json, sys, urllib.request


def _installed():
    """(kind, commit) for the copy already on this machine, without importing it."""
    try:
        info = _json.loads(_md.distribution("ppropt").read_text("direct_url.json") or "{}")
    except Exception as exc:
        return ("absent", None) if isinstance(exc, _md.PackageNotFoundError) else ("other", None)
    if "vcs_info" in info:
        return "github", info["vcs_info"].get("commit_id")
    if "dir_info" in info or str(info.get("url", "")).startswith("file:"):
        return "local", None
    return "other", None


def _latest():
    """The current main commit, or None when GitHub cannot be reached."""
    try:
        req = urllib.request.Request(f"https://api.github.com/repos/{REPO}/commits/main",
                                     headers={"Accept": "application/vnd.github.sha"})
        with urllib.request.urlopen(req, timeout=10) as r:
            return r.read().decode().strip()
    except Exception:
        return None


# A bare `import ppropt` succeeding is not evidence that the copy is current. An older build
# imports fine and then fails deep in a later cell on a key it never returned, which reads as
# a bug in the notebook. The version string is static, so the commit is what distinguishes
# them. A local checkout is left alone: it is someone working on the package.
_kind, _have = _installed()
_want = _latest()
#: Whether a previous cell already imported it. Reinstalling cannot replace a module object
#: that is already in memory, and the commit check above reads the disk, so "current" says
#: nothing about what this kernel is actually running.
_already = "ppropt" in sys.modules

if _kind == "local":
    _msg = "local checkout in use, left alone"
elif _kind == "absent" or (_kind == "github" and _want and _have != _want):
    token = None
    try:
        from google.colab import userdata          # private-repo fallback
        token = userdata.get("GITHUB_TOKEN")
    except Exception:
        pass
    url = (f"git+https://{token}@github.com/{REPO}.git" if token
           else f"git+https://github.com/{REPO}.git")
    if _kind == "absent":
        %pip install --quiet $url
        import ppropt                               # fail here, not three cells later
        _msg = f"installed ppropt @ {(_want or '?')[:7]}"
    elif "ppropt" in sys.modules:
        # Already imported: reinstalling cannot replace the module object in memory, so the
        # only honest thing to report is that a restart is still required.
        %pip install --quiet --force-reinstall --no-deps $url
        _msg = (f"updated to {(_want or '?')[:7]} — now use Runtime > Restart session and "
                f"run this cell again; the old copy is still loaded until you do")
    else:
        # Dependencies are already present; a newly added one fails loudly on import.
        %pip install --quiet --force-reinstall --no-deps $url
        import ppropt
        _msg = f"updated to {(_want or '?')[:7]}"
elif _kind == "other" or not _have:
    _msg = "ppropt installed, but not from git, so it cannot be checked against main"
elif _want is None:
    _msg = f"ppropt @ {_have[:7]}, could not reach GitHub to check it is current"
else:
    _msg = f"ppropt @ {_have[:7]}, current with {REPO}@main"
    if _already:
        _msg += ("  —  already imported in this session; if it was updated since, "
                 "Runtime > Restart session")

from IPython.display import HTML, display
display(HTML(
    f'<div style="border-left:3px solid #1a7f5a;padding:.5em .9em;'
    f'font-family:ui-monospace,monospace;font-size:13px;opacity:.85">{_msg}</div>'))
'''


def setup_cell():
    return code(_SETUP.replace("__REPO__", REPO), title="Setup - install PPR-opt")


#: Terms both notebooks use. Written for someone who wants a PPR and has not met Golden Gate
#: assembly: no term is explained using another term from further down the table, and nothing
#: here is a claim about the design -- only about what the word means.
_SHARED_TERMS = [
    ("PPR protein", "What this tool designs. A chain of near-identical **repeats**, each one "
                    "gripping a single RNA base."),
    ("target", "The RNA sequence you want the protein to bind. You choose it."),
    ("repeat", "One unit of the protein. Two amino acids inside it decide which RNA base it "
               "grips, and that mapping is published, so the design is a lookup rather than "
               "a search."),
    ("architecture (9S / 14S / 19S)", "How many repeats the protein has. It follows from "
                                      "your target's length: 9 bases gives 9S, 14 gives 14S, "
                                      "19 gives 19S."),
    ("codon", "Three DNA letters spelling one amino acid. Most amino acids have several "
              "possible codons, so the same protein can be written many different ways."),
    ("codon optimisation", "Choosing among those spellings so the DNA is easier to "
                           "manufacture and easier for your host to read. The protein does "
                           "not change."),
    ("coding sequence (CDS)", "The DNA that spells the protein, start to finish."),
    # Not "an orderable piece": whether a fragment may be ordered is a verdict reached later,
    # and `tests/test_display_verdict.py` bans the phrase for exactly that reason.
    ("fragment", "One of the pieces a gene is split into. Long stretches of DNA are "
                 "unreliable to synthesise, so the gene is made in a few parts and joined."),
    ("overhang", "Four unpaired bases at the end of a fragment. Only its intended partner "
                 "matches it, so the pieces can join in one order and no other."),
    ("`oh5` / `oh3`", "The overhang at a fragment's start and at its end."),
    ("junction", "The point where two fragments meet, held together by their matching "
                 "overhangs."),
    ("Golden Gate assembly", "The method used to join the fragments: one enzyme cuts and the "
                             "pieces are stitched together in a single reaction."),
    ("BsaI", "That enzyme. It cuts *outside* the sequence it recognises, which is what lets "
             "the overhangs be any four bases you like."),
    ("wrapper", "Fixed DNA added to each end of a fragment so the enzyme has something to "
                "cut. It is not part of your gene, but it **is** part of what you order."),
    ("GC content", "The share of bases that are G or C. Very high or very low makes DNA hard "
                   "to synthesise."),
    ("oligo", "A piece of synthesised DNA. Here, one fragment."),
    ("fidelity", "The predicted chance the pieces assemble in the intended order. It comes "
                 "from published measurements of this enzyme, not from your bench."),
    ("coding-span QC", "Checks run on the gene itself."),
    ("final artifacts", "Checks run on the exact sequences you would send to a vendor, "
                        "wrapper included. This is the one that decides orderable."),
    ("orderable", "Passed every check this software declares. It is **not** vendor "
                  "acceptance, and **not** a claim that the protein works."),
]


def glossary(md, extra: list[tuple[str, str]] = ()):
    """A plain-language table of every term the notebook uses.

    Placed last on purpose. A reader meets these words while running the steps, and a
    glossary they have to scroll past to reach step 1 is one they will not read.
    """
    rows = "\n".join(f"| **{term}** | {meaning} |" for term, meaning in
                     [*_SHARED_TERMS, *extra])
    return md(f"""
---

## Words used here

Every term this notebook uses, in plain language. Nothing here is a claim about your
design — only what the word means.

| term | what it means |
|---|---|
{rows}
""")


def write_notebook(cells, out: Path, name: str) -> None:
    nb = {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True, "name": name},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }
    # newline="\n" or the generated file is platform-dependent. `.gitattributes` marks
    # *.ipynb as `-text`, so notebooks are stored byte-for-byte as committed and never
    # converted on checkout. Writing them from Windows without this put CRLF in the
    # repository, and `test_regenerating_would_not_change_it` then compared a CRLF checkout
    # against an LF regeneration on Linux and failed -- passing here and only here.
    out.write_text(json.dumps(nb, indent=1, ensure_ascii=False),
                   encoding="utf-8", newline="\n")
    print(f"wrote {out}  ({len(cells)} cells: "
          f"{sum(1 for c in cells if c['cell_type'] == 'code')} code, "
          f"{sum(1 for c in cells if c['cell_type'] == 'markdown')} markdown)")
