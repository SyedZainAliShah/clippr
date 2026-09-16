"""Where downloaded reference data is cached.

This exists because `Path(__file__).resolve().parents[2]` means two different things. In a
source checkout it is the repository root, so `data/codon_tables` lands beside the code and
everything works. After `pip install` the same expression is the interpreter's `Lib`
directory, and a cache written there sits next to `site-packages` in a location no installer
owns -- verified on a clean install, which created `Lib/data/codon_tables` and
`Lib/data/genomes` and downloaded the Kazusa table and the whole chloroplast genome into
them. It is the same resolve-against-the-repository mistake `overhangs.py` records for the
ligation matrices, which ship inside the package precisely to avoid it.

The matrices are *data the package needs*, so they are packaged. These are *downloads the
package makes*, so they belong in a cache directory instead.
"""
from __future__ import annotations

import os
from pathlib import Path


def cache_dir(name: str) -> Path:
    """The directory for cached downloads called `name`, without creating it.

    A checkout keeps using its existing `data/` cache, so a developer's populated tables and
    genomes are still found and nothing re-downloads. Anything else -- an installed wheel
    above all -- gets the user's cache directory. `CLIPPR_CACHE_DIR` overrides both, which is
    how a clean-install check can point the cache somewhere disposable.

    A checkout is identified by `pyproject.toml` beside `src/`, not by the cache already
    existing, so a fresh clone behaves like a checkout on its very first run.
    """
    override = os.environ.get("CLIPPR_CACHE_DIR")
    if override:
        return Path(override) / name

    root = Path(__file__).resolve().parents[2]
    if (root / "pyproject.toml").is_file():
        return root / "data" / name

    base = (os.environ.get("LOCALAPPDATA")
            or os.environ.get("XDG_CACHE_HOME")
            or str(Path.home() / ".cache"))
    return Path(base) / "clippr" / name
