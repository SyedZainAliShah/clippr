"""Demo: design an orthogonal CLIPPR regulator library and the matching PPR binders."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from clippr import crosstalk_report, describe, design_orthogonal_set  # noqa: E402

N = 8
LENGTH = 9

print("=" * 74)
print(f"CLIPPR regulator library: {N} orthogonal targets, {LENGTH} nt")
print("=" * 74)
print()

for metric in ("uniform", "weighted"):
    ts = design_orthogonal_set(N, length=LENGTH, metric=metric, seed=1)
    print(ts.report())
    print()
    print("  pairwise distances:")
    print(crosstalk_report(ts.targets, metric=metric))
    print()
    if metric == "uniform":
        chosen = ts
    print("-" * 74)
    print()

print("=" * 74)
print("PPR binders for the uniform-metric set")
print("=" * 74)
for i, t in enumerate(chosen.targets):
    d = describe(t)
    print(f"\n  [{i}] target 5'-{d['target_rna']}-3'   architecture {d['architecture']}")
    print(f"      PPR code   {d['ppr_code']}")
    print(f"      protein    {d['aa_length']} aa, {d['n_arelf']} ARELF motifs")
    print(f"                 {d['aa_sequence'][:64]}...")

print()
print("=" * 74)
print("Interpretation")
print("=" * 74)
print(f"""
  Minimum pairwise separation is {chosen.min_distance:.0f} of {LENGTH} positions, so the
  closest two targets differ at {chosen.min_distance:.0f} of {LENGTH} contacted bases. Every PPR in
  the library must discriminate its own target from {N - 1} others; the binding
  constraint is set by the closest pair, not the average.

  This maximises sequence separation only. It does not predict binding. Our
  measurement on the GRASP panel found the intended target was the dominant
  response for 16 of 30 engineered variants, so sequence separation is
  necessary but not sufficient - the library still needs validation.
""")
