"""Why a design came out the way it did — the alternatives, and why they were rejected.

A design pipeline that returns only an answer has to be trusted. One that returns the
answer, the candidates it considered and the reason each was ruled out can be *checked*,
which is what a wet lab needs when a design looks surprising.

Scope is deliberately the scientific decision level -- inputs, constraints, candidates,
decisions, reasons. It is **not** an optimiser trace: how many iterations DNA Chisel ran or
which split was scored first tells a biologist nothing and would bury what does.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class OverhangDecision:
    """One candidate overhang at one junction, and what became of it."""

    junction_cut: int
    sequence: str
    status: Literal["selected", "rejected", "considered"]
    reason: str
    fidelity: float | None = None
    local_realizations: int | None = None

    def __str__(self) -> str:
        head = f"{self.sequence}  {self.status:10s} cut {self.junction_cut:>4d}"
        bits = []
        if self.fidelity is not None:
            bits.append(f"fidelity {self.fidelity:.3f}")
        if self.local_realizations is not None:
            bits.append(f"{self.local_realizations} local realizations")
        tail = f"  [{', '.join(bits)}]" if bits else ""
        return f"{head}  {self.reason}{tail}"


@dataclass(frozen=True)
class DesignAudit:
    """The decision record for one design."""

    target_rna: str
    architecture: str
    enzyme_profile: str
    enzymes_excluded: tuple[str, ...]
    genetic_code: int
    organism: str
    cuts: tuple[int, ...]
    overhang_decisions: tuple[OverhangDecision, ...] = ()
    fidelity: float | None = None
    fidelity_components: tuple[float, float] | None = None
    fidelity_ceiling: float | None = None
    constraints_satisfied: bool = True
    qc_status: str = ""
    findings: tuple[str, ...] = field(default_factory=tuple)

    @property
    def selected(self) -> tuple[OverhangDecision, ...]:
        return tuple(d for d in self.overhang_decisions if d.status == "selected")

    @property
    def rejected(self) -> tuple[OverhangDecision, ...]:
        return tuple(d for d in self.overhang_decisions if d.status == "rejected")

    def to_dataframe(self):
        """The decision record as a table, one row per candidate overhang.

        Every row carries the design-level context as well as the decision, so a single
        exported file is self-contained -- a lab notebook entry should not depend on
        remembering which profile was active when it was made.
        """
        import pandas as pd

        return pd.DataFrame([{
            "target_rna": self.target_rna,
            "architecture": self.architecture,
            "enzyme_profile": self.enzyme_profile,
            "enzymes_excluded": " ".join(self.enzymes_excluded),
            "genetic_code": self.genetic_code,
            "organism": self.organism,
            "junction_cut": d.junction_cut,
            "overhang": d.sequence,
            "status": d.status,
            "reason": d.reason,
            "local_realizations": d.local_realizations,
            "reaction_fidelity": self.fidelity,
            "qc_status": self.qc_status,
        } for d in self.overhang_decisions])

    def write_csv(self, path) -> str:
        """Write the decision record for a lab notebook or a submission appendix."""
        from pathlib import Path

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self.to_dataframe().to_csv(p, index=False)
        return str(p)

    def report(self) -> str:
        """A readable account of the design decisions, for a person."""
        lines = [
            f"design audit — {self.target_rna} ({self.architecture})",
            f"  host           {self.organism}, genetic code {self.genetic_code}",
            f"  enzyme profile {self.enzyme_profile} "
            f"({', '.join(self.enzymes_excluded) or 'nothing excluded'})",
            f"  cuts           {', '.join(str(c) for c in self.cuts)}",
        ]
        if self.fidelity is not None:
            line = f"  fidelity       {self.fidelity:.3f} predicted"
            if self.fidelity_ceiling is not None:
                line += (f"  (ceiling {self.fidelity_ceiling:.3f}, set by the "
                         f"destination pair)")
            lines.append(line)
        lines.append(f"  constraints    {'satisfied' if self.constraints_satisfied else 'NOT satisfied'}")
        lines.append(f"  QC             {self.qc_status}")

        if self.selected:
            lines.append("\n  selected overhangs")
            for d in self.selected:
                lines.append("    " + str(d))
        if self.rejected:
            lines.append(f"\n  rejected overhangs ({len(self.rejected)})")
            for d in self.rejected:
                lines.append("    " + str(d))
        if self.findings:
            lines.append("\n  findings")
            for f in self.findings:
                lines.append(f"    - {f}")
        return "\n".join(lines)
