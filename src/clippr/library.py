"""Design a whole library at once, and order it as one thing.

`design.design_oneshot` answers "what does this one target need?". A wet lab with fifty
regulators does not want fifty answers in fifty folders -- it wants **one order sheet**,
one pool, and one place to see which designs need attention.

Two things only become visible at library scale, and both are handled here:

**Cross-talk is a property of the set, not of any member.** A design that is perfect alone
is useless if another target in the same library sits one base away from it.

**So is DNA homology, and it is a different risk.** Cross-talk asks whether two PPRs could bind
each other's target; homology asks whether two *genes* share enough identical DNA to recombine.
Every member of a PPR library carries the same scaffold, so measured on five designs every pair
shared at least 59 nt before anything was done about it. `LibraryResult.homology` reports it and
`homology.diversify_library` reduces it.
"""
from __future__ import annotations

import json

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .crosstalk import report as two_tier_report
from .design import design_oneshot, synthesis_readiness
from .orthogonal import crosstalk_report, distance


@dataclass
class LibraryResult:
    """Every design in the library, plus what only makes sense across the whole set."""

    designs: dict[str, dict] = field(default_factory=dict)
    failed: dict[str, str] = field(default_factory=dict)

    @property
    def targets(self) -> list[str]:
        return list(self.designs)

    @staticmethod
    def _ready(result: dict) -> bool:
        """One authoritative verdict, shared with `design_oneshot`.

        This used to trust `synthesis_ready` when present and fall back to
        `status != "fail"` when not, which admitted `not_assessed` and any unrecognised
        status. The predicate now lives in one place and reads the evidence.
        """
        return synthesis_readiness(result)[0]

    @property
    def accepted(self) -> dict[str, dict]:
        """Designs that passed the final-artifact contract."""
        return {t: r for t, r in self.designs.items() if self._ready(r)}

    @property
    def rejected(self) -> dict[str, dict]:
        """Designs that were produced and must not be ordered."""
        return {t: r for t, r in self.designs.items() if not self._ready(r)}

    def rejection_reasons(self) -> dict[str, list[str]]:
        """Why each rejected design may not be ordered, from the shared predicate."""
        out = {}
        for target, r in self.rejected.items():
            why = synthesis_readiness(r)[1]
            out[target] = why or ["rejected without a recorded reason"]
        return out

    def oligos(self, include_rejected: bool = False):
        """Fragments ready to send to a vendor.

        Rejected designs are excluded by default. They used to be concatenated in with the
        rest, which is how a design carrying `synthesis_ready=False` reached a normal order
        sheet. Pass `include_rejected=True` only for a diagnostic export.
        """
        import pandas as pd

        source = self.designs if include_rejected else self.accepted
        frames = []
        for target, r in source.items():
            df = r["oligos"].copy()
            df.insert(0, "target_rna", target)
            df.insert(1, "architecture", r["architecture"])
            if include_rejected:
                df.insert(2, "synthesis_ready", self._ready(r))
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
                # `qc` scored the coding span only, so a rejected design showed "PASS"
                # with nothing beside it. The scope is now in the column name and the
                # final verdict travels in the same row.
                "cds_qc": qc["status"],
                "final_artifacts": (r.get("final_artifacts") or {}).get(
                    "status", "not_recorded"),
                "synthesis_ready": self._ready(r),
                "gc_pct": qc["gc_pct"],
                "repeated_kmer_fraction": qc["repeated_kmer_fraction"],
                "longest_repeat": qc["longest_repeat"],
                "fidelity": r["fidelity"],
                "constraints_ok": r["constraints_ok"],
                "warnings": len(r["warnings"]),
            })
        return pd.DataFrame(rows)

    @property
    def totals(self) -> dict:
        """How much DNA this library is, across every design."""
        return {
            "n_designs": len(self.designs),
            "n_oligos": sum(len(r["oligos"]) for r in self.designs.values()),
            "total_bases": sum(int(r["oligos"]["oligo_length"].sum())
                               for r in self.designs.values()),
        }

    def crosstalk(self, metric: str = "uniform", scores=None) -> str:
        """Pairwise separation between targets. Only defined for one length at a time.

        Without `scores`, a distance matrix: every pair, sequence separation only. Pass a
        PPR specificity table (see `clippr.crosstalk.load_ppr_scores`) to get the two-tier
        report instead -- pairs sorted by risk, annotated with predicted affinity. The
        annotation gates nothing; separation remains the only criterion either way.

        The two-tier report counts differing positions unweighted, so asking for both a
        weighted metric and a scoring table is refused rather than silently ignoring one.
        """
        if scores and metric != "uniform":
            raise ValueError(
                f"the two-tier report uses unweighted Hamming, so metric={metric!r} "
                f"would be ignored; call crosstalk(metric) and crosstalk(scores=...) "
                f"separately")
        by_len: dict[int, list[str]] = {}
        for t in self.targets:
            by_len.setdefault(len(t), []).append(t)
        blocks = []
        for length, group in sorted(by_len.items()):
            if len(group) < 2:
                continue
            blocks.append(f"{len(group)} targets of length {length}:")
            blocks.append(two_tier_report(group, scores) if scores
                          else crosstalk_report(group, metric))
        return "\n\n".join(blocks) if blocks else "fewer than two targets of any one length"

    def homology(self, limit: int = 10) -> str:
        """DNA shared between members — the recombination and synthesis risk, not RNA cross-talk.

        A separate question from `crosstalk`: that one asks whether two PPRs could bind each
        other's target, this one asks whether two *genes* share enough identical DNA to
        recombine. Members of a PPR library share a scaffold by construction, so the answer is
        never trivially no. See `homology.py` for what was measured.
        """
        from .homology import report as homology_report

        return homology_report({t: r["cds"] for t, r in self.designs.items()}, limit=limit)

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
        c = self.totals
        lines = [
            f"library of {len(self.designs)} designs"
            + (f", {len(self.failed)} failed" if self.failed else ""),
            f"  fragments      {c['n_oligos']} across all designs, "
            f"{c['total_bases']:,} bases",
        ]
        # The orderable count first. A summary that led with "QC 1 PASS" read as
        # success on a library whose every design was refused.
        lines.append(f"  synthesis      {len(self.accepted)} of {len(self.designs)} "
                     f"designs are orderable"
                     + (f", {len(self.rejected)} not orderable" if self.rejected else ""))
        if not qc.empty:
            counts = qc["cds_qc"].value_counts().to_dict()
            lines.append("  coding-span QC "
                         + ", ".join(f"{v} {k}" for k, v in counts.items())
                         + "  (the CDS only, not the ordered sequences)")
        close = self.closest_pairs(3)
        if close:
            lines.append("  closest targets "
                         + ", ".join(f"{a}/{b} differ at {d}" for a, b, d in close))
        for t, why in self.rejection_reasons().items():
            lines.append(f"  NOT ORDERABLE {t}: {why[0]}")
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

    # Accepted designs only. A normal order sheet must not contain a design the
    # final-artifact contract rejected.
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

    # Rejected designs are kept and named, never silently dropped.
    if lib.rejected:
        rejected = lib.oligos(include_rejected=True)
        rejected = rejected[~rejected["synthesis_ready"]]
        p = out / "library_REJECTED_oligos.csv"
        rejected.to_csv(p, index=False)
        paths["rejected_oligos"] = str(p)
        p = out / "library_REJECTED_reasons.json"
        p.write_text(json.dumps(lib.rejection_reasons(), indent=2, sort_keys=True) + "\n",
                     encoding="utf-8")
        paths["rejected_reasons"] = str(p)

    qc = lib.qc_table()
    if not qc.empty:
        p = out / "library_qc.csv"
        qc.to_csv(p, index=False)
        paths["qc_table"] = str(p)

    gb = out / "genbank"
    gb.mkdir(exist_ok=True)
    # A refused design must not be written under the name an orderable one gets.
    # `design_oneshot` marks its diagnostic files REJECTED precisely so a downstream step
    # cannot pick them up; writing the same record here as clippr_<target>.gb made a
    # sequence that must not be ordered indistinguishable from one that may be, on disk,
    # which is the one place the distinction has to survive.
    refused = set(lib.rejected)
    for target, r in lib.designs.items():
        stem = f"clippr_{target}_REJECTED" if target in refused else f"clippr_{target}"
        write_genbank(r["cds"], r["fragments"], target, gb / f"{stem}.gb",
                      genetic_code=r.get("genetic_code", 1))
    paths["genbank_dir"] = str(gb)

    p = out / "library_summary.txt"
    p.write_text(lib.summary() + "\n\n" + lib.crosstalk() + "\n\n"
                 + lib.homology() + "\n", encoding="utf-8")
    paths["summary"] = str(p)

    # The manifest ships with the package rather than being a separate step a caller can
    # forget. Without it the CSVs and GenBank files carry no record of the source, codon
    # table, dependencies or screening state that produced them, and a result package whose
    # provenance has to be reconstructed from memory is not reproducible.
    from . import manifest as _manifest

    failures = [{"target": t, "reason": why} for t, why in sorted(lib.failed.items())]
    # Coverage, not just crashes. `completion` counted `lib.failed` alone, so a library
    # whose every design was correctly refused still shipped a manifest saying "complete"
    # with an empty failure list.
    rejections = [{"target": t, "reasons": why}
                  for t, why in sorted(lib.rejection_reasons().items())]
    unassessed = sum(1 for r in lib.rejected.values()
                     if (r.get("final_artifacts") or {}).get("status")
                     not in ("fail", "pass", "warning"))
    coverage = {"requested": len(lib.designs) + len(lib.failed),
                "accepted": len(lib.accepted), "rejected": len(lib.rejected),
                "not_assessed": unassessed, "failed": len(lib.failed)}
    if lib.designs:
        record = _manifest.build(list(lib.designs.values()),
                                 completion=("complete" if not failures and not rejections
                                             else "partial"),
                                 failures=failures)
        record["synthesis_coverage"] = coverage
        record["rejections"] = rejections
    else:
        # Every target failed, so there is no codon table or enzyme profile that ran and
        # nothing to describe in `inputs`. The run still needs a provenance record: a
        # directory of failures with no statement of what produced them is the case the
        # manifest exists to prevent.
        record = {"schema_version": _manifest.SCHEMA_VERSION,
                  "route": "synthesis", "completion": "failed",
                  "software": _manifest._software(), "designs": [],
                  "failures": failures, "synthesis_coverage": coverage,
                  "rejections": rejections}
    paths["manifest"] = str(_manifest.write(record, out / "manifest.json"))
    return paths
