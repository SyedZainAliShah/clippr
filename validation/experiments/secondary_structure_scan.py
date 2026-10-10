"""Secondary-structure scan of the ordered fragments, baseline against wrapper-aware.

IDT's own explanation of complexity screening says GC is a proxy: *"Regions with extreme GC
content may have higher propensity to form stable secondary structures."* CLIPPR enforces the
proxy and has never measured the quantity it stands for. This measures it.

It also answers a question the wrapper-aware work raised and did not address: that experiment
changed codons to lower GC on every ordered fragment. **Lowering GC could have made folding
worse, better or neither**, and nothing checked. A fix for one proxy that degrades the real
property would be a bad trade made invisibly.

Two figures per fragment:

  first-48 dG     the strongest fold in the first 48 nt, which is the span Twist names
                  (its stated concern is a hairpin with dG < -8 kcal/mol there)
  whole dG        minimum free energy of the entire ordered fragment

The whole-fragment figure has **no published threshold** and is reported as a distribution for
comparison between arms, not scored against a limit. A 257-mer folds substantially by nature;
a large negative number is not by itself a finding.

Thermodynamics are ViennaRNA with the Mathews 2004 DNA parameters. This is not IDT's
implementation and is not expected to reproduce its numbers exactly. Note that IDT's
OligoAnalyzer and primer3 both cap thermodynamic calculations at 60 bases, so neither can
evaluate a 257 nt fragment at all; that is why a folding engine is used here.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import RNA                                                            # noqa: E402

#: The span Twist names, and its stated limit there.
TWIST_WINDOW_NT = 48
TWIST_HAIRPIN_DG = -8.0


def _fold(seq: str) -> float:
    structure, mfe = RNA.fold_compound(seq).mfe()
    return float(mfe)


def _rows(path: Path) -> list[dict]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        if "__run__" in entry or "__run_identity__" in entry:
            continue
        out.append(entry)
    return out


def scan(rows: list[dict], label: str) -> dict:
    first48, whole, breaches = [], [], []
    for entry in rows:
        target = entry.get("target") or entry.get("target_rna")
        for o in entry.get("oligos") or []:
            seq = (o.get("oligo_sequence_5to3") or "").strip().upper()
            if not seq:
                continue
            head = _fold(seq[:TWIST_WINDOW_NT])
            full = _fold(seq)
            first48.append(head)
            whole.append(full)
            if head < TWIST_HAIRPIN_DG:
                breaches.append({"target": target,
                                 "fragment": o.get("order_fragment_id") or o.get("fragment_id"),
                                 "first48_dg": round(head, 2),
                                 "whole_dg": round(full, 2)})
    def dist(xs):
        xs = sorted(xs)
        return {"n": len(xs), "min": round(xs[0], 2),
                "p05": round(xs[int(0.05 * len(xs))], 2),
                "median": round(statistics.median(xs), 2),
                "max": round(xs[-1], 2), "mean": round(statistics.fmean(xs), 2)}

    return {"arm": label, "first48_dg": dist(first48), "whole_dg": dist(whole),
            "first48_breaching_twist_minus8": len(breaches),
            "breaches": breaches[:20]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", default="work/independent/corpus200_pinned.jsonl")
    ap.add_argument("--candidates", required=True,
                    help="the wrapper-aware candidate JSONL to compare against")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    RNA.params_load_DNA_Mathews2004()

    started = time.perf_counter()
    base = scan(_rows(ROOT / args.baseline), "baseline (shipped corpus)")
    cand = scan(_rows(Path(args.candidates)), "wrapper-aware candidates")

    # Did lowering GC change folding, and in which direction?
    delta_first48 = cand["first48_dg"]["median"] - base["first48_dg"]["median"]
    delta_whole = cand["whole_dg"]["median"] - base["whole_dg"]["median"]

    result = {
        "experiment": "secondary_structure_scan",
        "engine": f"ViennaRNA {RNA.__version__}, Mathews 2004 DNA parameters",
        "question": ("does the wrapper-aware GC fix change secondary structure, and does "
                     "any fragment breach the one published hairpin threshold we can find?"),
        "twist_rule": {"window_nt": TWIST_WINDOW_NT, "dg_limit": TWIST_HAIRPIN_DG,
                       "source": "Twist codon-optimisation FAQ: strong hairpins (dG < -8) "
                                 "in the first 48 base pairs"},
        "arms": [base, cand],
        "median_shift": {"first48_dg": round(delta_first48, 2),
                         "whole_dg": round(delta_whole, 2),
                         "reading": ("negative means the candidates fold more tightly than "
                                     "the baseline, which would be a cost; positive means "
                                     "they fold less tightly")},
        "scope": ("a folding prediction, not a measurement, and not IDT's implementation. "
                  "The whole-fragment figure has no published threshold and is reported for "
                  "comparison between arms only. Nothing here is a vendor acceptance."),
        "elapsed_seconds": round(time.perf_counter() - started, 2),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    for arm in (base, cand):
        print(f"=== {arm['arm']} ===")
        print(f"  first 48 nt dG : {arm['first48_dg']}")
        print(f"  whole fragment : {arm['whole_dg']}")
        print(f"  breaching Twist's -8 in the first 48: "
              f"{arm['first48_breaching_twist_minus8']} of {arm['first48_dg']['n']}")
    print(f"\nmedian shift, first 48 nt : {delta_first48:+.2f} kcal/mol")
    print(f"median shift, whole        : {delta_whole:+.2f} kcal/mol")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
