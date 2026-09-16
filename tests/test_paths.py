"""Where cached downloads land, in a checkout and after `pip install`.

The bug these pin: `Path(__file__).resolve().parents[2]` is the repository root in a
checkout but the interpreter's `Lib` directory in an installed wheel. A clean install
built from this package created `Lib/data/codon_tables` and `Lib/data/genomes` beside
`site-packages` and downloaded into them.
"""
from __future__ import annotations

from pathlib import Path

from clippr.paths import cache_dir


class TestCacheDir:
    def test_env_override_wins(self, monkeypatch, tmp_path):
        monkeypatch.setenv("CLIPPR_CACHE_DIR", str(tmp_path))
        assert cache_dir("genomes") == tmp_path / "genomes"

    def test_checkout_uses_the_repository_data_directory(self, monkeypatch):
        """So a developer's already-populated cache keeps being found."""
        monkeypatch.delenv("CLIPPR_CACHE_DIR", raising=False)
        root = Path(__file__).resolve().parents[1]
        assert (root / "pyproject.toml").is_file(), "test assumes it runs from a checkout"
        assert cache_dir("codon_tables") == root / "data" / "codon_tables"

    def test_installed_layout_never_writes_beside_site_packages(self, monkeypatch):
        """The actual defect: an installed package must not cache into the interpreter."""
        monkeypatch.delenv("CLIPPR_CACHE_DIR", raising=False)
        monkeypatch.setenv("LOCALAPPDATA", r"C:\Users\someone\AppData\Local")
        installed = Path(r"C:\env\Lib\site-packages\clippr\paths.py")
        monkeypatch.setattr("clippr.paths.__file__", str(installed))

        got = cache_dir("codon_tables")
        lib = installed.resolve().parents[2]            # the C:\env\Lib that caused the bug
        assert lib not in got.parents, f"{got} is inside the interpreter directory {lib}"
        assert got.parts[-2:] == ("clippr", "codon_tables")
