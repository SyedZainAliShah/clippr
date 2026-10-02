"""Shared cell constructors for the CLIPPR notebooks.

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
#@markdown Installs the package if it is not already available. Safe to re-run — it never
#@markdown reinstalls over a working copy. Takes about a minute the first time.
#@markdown
#@markdown **If it fails:** run it again — a dropped network call is the usual cause. If it
#@markdown still fails, use **Runtime → Restart session** and then run it once more, which
#@markdown clears a half-installed package. Every cell below needs this one to have
#@markdown succeeded, so do not skip past a red error here.
REPO = "__REPO__"

try:
    import clippr
    _msg = f"clippr {clippr.__version__} already available"
except ImportError:
    token = None
    try:
        from google.colab import userdata          # private-repo fallback
        token = userdata.get("GITHUB_TOKEN")
    except Exception:
        pass
    url = (f"git+https://{token}@github.com/{REPO}.git" if token
           else f"git+https://github.com/{REPO}.git")
    %pip install --quiet $url
    import clippr
    _msg = f"installed clippr {clippr.__version__}"

from IPython.display import HTML, display
display(HTML(
    f'<div style="border-left:3px solid #1a7f5a;padding:.5em .9em;'
    f'font-family:ui-monospace,monospace;font-size:13px;opacity:.85">{_msg}</div>'))
'''


def setup_cell():
    return code(_SETUP.replace("__REPO__", REPO), title="Setup — install CLIPPR")


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
    out.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {out}  ({len(cells)} cells: "
          f"{sum(1 for c in cells if c['cell_type'] == 'code')} code, "
          f"{sum(1 for c in cells if c['cell_type'] == 'markdown')} markdown)")
