"""Design a whole library at once, and order it as one thing.

`design.design_oneshot` answers "what does this one target need?". A wet lab with fifty
regulators does not want fifty answers in fifty folders -- it wants **one order sheet**,
one pool, one cost, and one place to see which designs need attention.

Two things only become visible at library scale, and both are handled here:

**Pooling is where the money is.** An oligo pool is priced per pool, essentially flat
across the sizes this project produces. Ordering fifty designs as fifty pools costs fifty
times what ordering them as one pool costs, for the same DNA. `LibraryResult.cost`
compares the two so the difference is impossible to miss.

**Cross-talk is a property of the set, not of any member.** A design that is perfect alone
is useless if another target in the same library sits one base away from it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .design import design_oneshot
from .export import OPOOL_LIST_PRICE_EUR, PHOSPHORYLATION_EUR_PER_OLIGO
from .orthogonal import crosstalk_report, distance


@dataclass
class LibraryResult:
    """Every design in the library, plus what only makes sense across the whole set."""

    designs: dict[str, dict] = field(default_factory=dict)
    failed: dict[str, str] = field(default_factory=dict)

    @property
    def targets(self) -> list[str]:
        return list(self.designs)

    def oligos(self):
        """Every fragment in the library as one table, ready to send to a vendor."""
        import pandas as pd

        frames = []
        for target, r in self.designs.items():
            df = r["oligos"].copy()
            df.insert(0, "target_rna", target)
            df.insert(1, "architecture", r["architecture"])
            frames.append(df)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def qc_table(self):
        """One row per design: what passed, what needs looking at."""
        import pandas as pd

        rows = []
        for target, r in self.designs.items():
            qc = r["qc"]
            rows.append({
                "target_rna": target,
                "architecture": r["architecture"],
                "protein_aa": len(r["protein"]),
                "cds_nt": len(r["cds"]),
                "fragments": len(r["oligos"]),
                "qc": qc["status"],
                "gc_pct": qc["gc_pct"],
                "repeated_kmer_fraction": qc["repeated_kmer_fraction"],
                "longest_repeat": qc["longest_repeat"],
                "fidelity": r["fidelity"],
                "constraints_ok": r["constraints_ok"],
                "warnings": len(r["warnings"]),
            })
        return pd.DataFrame(rows)

    @property
    def cost(self) -> dict:
        """Pooled versus per-design cost. The comparison is the point.

        A list price, not a quote, and flat per pool across these sizes -- so the saving
        from pooling is a property of the price model, not of any design.
        """
        n_designs = len(self.designs)
        n_oligos = sum(len(r["oligos"]) for r in self.designs.values())
        total_bases = sum(int(r["oligos"]["oligo_length"].sum())
                          for r in self.designs.values())
        pooled = OPOOL_LIST_PRICE_EUR
        separate = OPOOL_LIST_PRICE_EUR * n_designs
        return {
            "n_designs": n_designs,
            "n_oligos": n_oligos,
            "total_bases": total_bases,
            "pooled_eur": pooled,
            "separate_pools_eur": separate,
            "saving_eur": separate - pooled,
            "phosphorylation_eur": round(PHOSPHORYLATION_EUR_PER_OLIGO * n_oligos, 2),
            "list_price_not_a_quote": True,
        }

    def crosstalk(self, metric: str = "uniform") -> str:
        """Pairwise separation between targets. Only defined for one length at a time."""
        by_len: dict[int, list[str]] = {}
        for t in self.targets:
            by_len.setdefault(len(t), []).append(t)
        blocks = []
        for length, group in sorted(by_len.items()):
            if len(group) < 2:
                continue
            blocks.append(f"{len(group)} targets of length {length}:")
            blocks.append(crosstalk_report(group, metric))
        return "\n\n".join(blocks) if blocks else "fewer than two targets of any one length"

    def closest_pairs(self, limit: int = 5) -> list[tuple[str, str, int]]:
        """The pairs most at risk of cross-reacting, closest first."""
        pairs = []
        for i, a in enumerate(self.targets):
            for b in self.targets[i + 1:]:
                if len(a) == len(b):
                    pairs.append((a, b, int(distance(a, b, [1.0] * len(a)))))
        return sorted(pairs, key=lambda p: p[2])[:limit]

    def summary(self) -> str:
        qc = self.qc_table()
        c = self.cost
        lines = [
            f"library of {len(self.designs)} designs"
            + (f", {len(self.failed)} failed" if self.failed else ""),
            f"  fragments      {c['n_oligos']} across all designs, "
            f"{c['total_bases']:,} bases",
            f"  pooled cost    {c['pooled_eur']:.2f} EUR as one pool, versus "
            f"{c['separate_pools_eur']:.2f} separately "
            f"— {c['saving_eur']:.2f} EUR saved (list price, not a quote)",
        ]
        if not qc.empty:
            counts = qc["qc"].value_counts().to_dict()
            lines.append("  QC             "
                         + ", ".join(f"{v} {k}" for k, v in counts.items()))
        close = self.closest_pairs(3)
        if close:
            lines.append("  closest targets "
                         + ", ".join(f"{a}/{b} differ at {d}" for a, b, d in close))
        for t, err in self.failed.items():
            lines.append(f"  FAILED {t}: {err}")
        return "\n".join(lines)


def design_library(targets, outdir: str | Path | None = None,
                   on_progress=None, **kwargs: Any) -> LibraryResult:
    """Design every target, collecting failures rather than stopping at the first one.

    One bad target in fifty should not cost you the other forty-nine, so a design that
    raises is recorded in `failed` and the run continues.
    """
    lib = LibraryResult()
    targets = [str(t).strip().upper() for t in targets if str(t).strip()]
    for i, target in enumerate(targets, 1):
        try:
            lib.designs[target] = design_oneshot(target, outdir=None, **kwargs)
        except Exception as exc:
            lib.failed[target] = f"{type(exc).__name__}: {exc}"
        if on_progress:
            on_progress(i, len(targets), target)

    if outdir is not None:
        write_library(lib, outdir)
    return lib


def write_library(lib: LibraryResult, outdir: str | Path) -> dict[str, str]:
    """Write the library as one order sheet, one QC table, and per-design GenBank."""
    from .export import write_genbank

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    oligos = lib.oligos()
    if not oligos.empty:
        p = out / "library_order_sheet.csv"
        oligos.to_csv(p, index=False)
        paths["order_sheet"] = str(p)

        # the minimal columns an oligo-pool vendor actually needs
        pool = oligos[["order_fragment_id", "oligo_sequence_5to3"]].rename(
            columns={"order_fragment_id": "Pool name", "oligo_sequence_5to3": "Sequence"})
        p = out / "library_opool.csv"
        pool.to_csv(p, index=False)
        paths["opool"] = str(p)

    qc = lib.qc_table()
    if not qc.empty:
        p = out / "library_qc.csv"
        qc.to_csv(p, index=False)
        paths["qc_table"] = str(p)

    gb = out / "genbank"
    gb.mkdir(exist_ok=True)
    for target, r in lib.designs.items():
        write_genbank(r["cds"], r["fragments"], target, gb / f"clippr_{target}.gb")
    paths["genbank_dir"] = str(gb)

    p = out / "library_summary.txt"
    p.write_text(lib.summary() + "\n\n" + lib.crosstalk() + "\n", encoding="utf-8")
    paths["summary"] = str(p)
    return paths
