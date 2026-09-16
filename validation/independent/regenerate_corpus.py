"""Regenerate the 200-design corpus under the current source, with full run identity.

Every independent check so far (V5 digest-and-ligate, V6 constraints, V7 export semantics)
reads a corpus whose journal header predates content-complete run identity: it records the
configuration but not the source that produced the rows. That corpus is kept as historical
evidence and is not relabelled; this writes a new one that is pinned.

Long by design -- roughly three and a half hours, dominated by 19S -- so it runs through
`clippr.experiment`, which journals each row as it completes, checkpoints elapsed time per
item, quarantines an interrupted trailing write and refuses a resume whose identity differs.
An interrupted run resumes; it does not restart.

    python validation/independent/regenerate_corpus.py [--limit N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from integration import corpus                                    # noqa: E402

from clippr import experiment, manifest                           # noqa: E402
from clippr.design import design_oneshot                          # noqa: E402
from clippr.overhangs import cache_stats                          # noqa: E402

TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--wall-seconds", type=float, default=None)
    ap.add_argument("--out", default=str(ROOT / "work" / "independent"))
    args = ap.parse_args()

    table = json.loads(TABLE.read_text(encoding="utf-8"))
    identity = manifest.run_identity(
        table, 1, matrix="BsaI-HFv2", experiment="pinned 200-design corpus",
        codon_table_path=str(TABLE),
        codon_table_file_sha256=hashlib.sha256(TABLE.read_bytes()).hexdigest(),
        enzyme_profile="igem_rfc1000", destination="level0", seed=42,
        check_offtarget=False)

    targets = corpus()
    if args.limit:
        targets = targets[:args.limit]
    out = Path(args.out)
    print(f"{len(targets)} targets | source {identity['source_sha256'][:16]} | "
          f"scorer {identity['scorer_version']} | matrix {identity['matrix_fingerprint']}")

    started = time.perf_counter()

    def build(target: str) -> dict:
        d = design_oneshot(target, codon_table=table, genetic_code=1,
                           enzyme_profile="igem_rfc1000", matrix="BsaI-HFv2",
                           check_offtarget=False, seed=42)
        return {
            "target": target, "architecture": d["architecture"],
            "cds": d["cds"], "protein": d["protein"],
            "cds_sha256": hashlib.sha256(d["cds"].encode()).hexdigest(),
            "cuts": list(d["cuts"]), "junction_overhangs": list(d["junction_overhangs"]),
            "fidelity": d["fidelity"], "qc": d["qc"],
            "constraints_ok": bool(d["constraints_ok"]),
            "offtarget": d["offtarget"],
            "offtarget_status": d["audit"].offtarget_status,
            "oligos": d["oligos"].to_dict(orient="records"),
            "timings": d["timings"],
        }

    def progress(i, n, key):
        if i % 5 == 0:
            print(f"  [{i:>3}/{n}] {key}  elapsed {(time.perf_counter() - started) / 60:.1f} min",
                  flush=True)

    result = experiment.run(targets, build, path=out / "corpus200_pinned.jsonl",
                            identity=identity, wall_seconds=args.wall_seconds,
                            on_progress=progress)

    print(f"\nstatus {result['status']}  requested {result['requested']}  "
          f"new {result['completed_this_run']}  resumed {result['resumed']}  "
          f"rows {result['completed']}  failures {len(result['failures'])}")
    print(f"wall {result['wall_seconds_all_runs'] / 3600:.2f} h across all runs  "
          f"| budget exceeded {result['budget_exceeded']}  "
          f"| quarantined tail {result['quarantined_tail']}")
    print(f"planning cache: {cache_stats()}")

    designs = [r for r in result["rows"]]
    summary = {
        "identity": identity, "status": result["status"],
        "requested": result["requested"], "rows": result["completed"],
        "failures": result["failures"],
        "architectures": {a: sum(1 for d in designs if d["architecture"] == a)
                          for a in ("9S", "14S", "19S")},
        "qc_pass": sum(1 for d in designs if d["qc"]["status"] == "PASS"),
        "constraints_ok": sum(1 for d in designs if d["constraints_ok"]),
        "wall_seconds_all_runs": result["wall_seconds_all_runs"],
    }
    (out / "corpus200_pinned_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'corpus200_pinned.jsonl'} and its summary")
    return 0 if result["status"] == "complete" and not result["failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
