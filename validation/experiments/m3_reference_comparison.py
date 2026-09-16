"""M3 / E6 -- compare our recoded inventory with the reference's, on matched inputs.

Runs the reference implementation's own reusable-library workflow in its own environment and
compares it with our M3 pilot on the one objective both can be held to: codon adaptation of
the same modules under the same host table, scored here by one scorer applied to both.

**Which reference.** The pinned local checkout, v0.1.15 at commit 8882759a, not the v0.1.5
wheel every earlier oracle comparison in this project used. They are different references and
the older one cannot stand in for this workflow, which does not exist in it.

**What is held fixed.** The reference is given our Chlamydomonas codon table in its own
`codon,aa,frequency` format, so both systems optimise toward the same host. Module identity
joins directly once the naming conventions are normalised: our `pPR-1_1A_5N_AATG` is their
`1A_5N_AATG`, which their output suffixes `_v1`. Both inventories contain the `AATG` and
`AGGT` variants, so the config's `ppr_5prime_fusion_site` selects which variant an *assembly*
uses; it does not restrict this module-level comparison.

**Declared mismatches.** These are stated rather than aligned, and the comparison is limited
to common scope rather than presented as a head-to-head:

    GC band            reference 0.25-0.65 global, 0.15-0.85 windowed | ours 0.35-0.65 / 50 nt
    max homopolymer    reference 3                                    | ours 4
    repeat k           reference 16                                   | ours 20
    scope              reference redesigns all 42 modules from amino acids plus a coding mask
                       | ours recodes the 21 the six reference targets need, from DNA, with
                       the four-base interfaces frozen

A looser GC band and a stricter homopolymer limit push in opposite directions, so neither
system is simply the more constrained one. Runtime is reported per module for that reason,
and even then it compares different amounts of work on different inputs.

A reference failure is reported as an observed failure of that run under these inputs. It is
not evidence about the reference generally, and not evidence for us. The first attempt here
died on a `UnicodeEncodeError`; that was this harness inheriting cp1252, a fault of ours, and
it is not recorded as a reference result.

    python validation/experiments/m3_reference_comparison.py [--wall-seconds 3600]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

REFERENCE = ROOT.parent / "grasp-library-designer"
REF_PYTHON = REFERENCE / ".venv" / "Scripts" / "python.exe"
RUNNER = Path(__file__).resolve().parent / "m3_reference_runner.py"
TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"

#: Their ids carry a version suffix; ours carry a plasmid prefix. Strip both to compare.
OUR_PREFIX = "pPR-1_"
THEIR_SUFFIX = "_v1"

DECLARED_MISMATCHES = {
    "gc_band": "reference 0.25-0.65 global / 0.15-0.85 windowed; ours 0.35-0.65 over 50 nt",
    "max_homopolymer": "reference 3; ours 4",
    "repeat_k": "reference 16; ours 20",
    "scope": ("reference redesigns all 42 modules from amino acids plus a coding mask; ours "
              "recodes the 21 the six reference targets need, from DNA, with the four-base "
              "interfaces frozen"),
}


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_reference_codon_csv(destination: Path) -> int:
    """Our Chlamydomonas table in the reference's codon,aa,frequency form."""
    from Bio.Data import CodonTable

    forward = dict(CodonTable.unambiguous_dna_by_id[1].forward_table)
    table = json.loads(TABLE.read_text(encoding="utf-8"))
    lines = ["codon,aa,frequency"]
    for aa, codons in sorted(table.items()):
        for codon, value in sorted(codons.items()):
            dna = codon.replace("U", "T")
            if forward.get(dna) != aa and aa != "*":
                continue
            lines.append(f"{dna},{aa},{value * 1000:.4f}")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines) - 1


def cai(dna: str, table) -> float:
    """One scorer, applied to both systems' output. Neither system's own metric is used.

    Their `codon_score` and our `cai_after` are not the same statistic, so reading each
    system's self-reported number and putting the two in a table would compare definitions
    rather than designs.
    """
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    weights = {c.replace("U", "T"): v / max(d.values())
               for aa, d in table.items() for c, v in d.items() if max(d.values()) > 0}
    used = [weights[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3)
            if dna[i:i + 3] in weights and forward.get(dna[i:i + 3]) not in ("M", "W", None)]
    if not used or min(used) <= 0:
        return 0.0
    return math.exp(sum(math.log(w) for w in used) / len(used))


def translate(dna: str) -> str:
    from Bio.Seq import Seq
    return str(Seq(dna).translate())


def run_reference(work: Path, wall_seconds: float) -> tuple[dict | None, float]:
    """Execute the reference workflow in its own interpreter; None on an observed failure."""
    inp, outp = work / "input", work / "output"
    # The reference prints arrows and degree signs. A Windows subprocess inherits cp1252 and
    # dies encoding them -- a fault in this harness, not the reference, so force UTF-8 rather
    # than recording a crash that belongs to us.
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    started = time.perf_counter()
    try:
        proc = subprocess.run([str(REF_PYTHON), str(RUNNER), str(inp), str(outp)],
                              capture_output=True, text=True, encoding="utf-8",
                              errors="replace", timeout=wall_seconds,
                              cwd=str(REFERENCE), env=env)
    except subprocess.TimeoutExpired:
        print(f"OBSERVED FAILURE: reference run exceeded the {wall_seconds}s cap")
        return None, wall_seconds
    elapsed = time.perf_counter() - started

    if "JSON_BEGIN" not in proc.stdout:
        print(f"OBSERVED FAILURE: reference run produced no result in {elapsed:.1f}s "
              f"(exit {proc.returncode})")
        (work / "reference_stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
        print(f"  stderr written to {work / 'reference_stderr.txt'}")
        return None, elapsed
    return json.loads(proc.stdout.split("JSON_BEGIN")[1].split("JSON_END")[0]), elapsed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wall-seconds", type=float, default=3600.0)
    ap.add_argument("--table", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    ap.add_argument("--repeats", type=int, default=5,
                    help="reference runs to average; it does not reproduce under its own seed")
    args = ap.parse_args()

    from clippr.codons import complete_table
    from clippr.products import load_inserts

    if not REF_PYTHON.is_file():
        print(f"reference interpreter not found at {REF_PYTHON}")
        return 2

    # Archive rather than delete. A rerun used to remove the directory outright, which
    # destroyed the evidence a reviewer was in the middle of reading -- including the run
    # payloads that are the only replayable record of a non-reproducible reference.
    work = ROOT / "work" / "m3" / "reference"
    if work.exists():
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archived = work.with_name(f"reference_archived_{stamp}")
        shutil.move(str(work), str(archived))
        print(f"archived the previous comparison to {archived.name}")
    shutil.copytree(REFERENCE / "grasp_library_project" / "input", work / "input")
    codons = write_reference_codon_csv(work / "input" / "codon_usage.csv")

    revision = subprocess.run(["git", "-C", str(REFERENCE), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    print(f"reference: {REFERENCE.name} at {revision[:16]}")
    print(f"  interpreter {REF_PYTHON}")
    print(f"  supplied our Chlamydomonas table as {codons} codon rows  |  cap "
          f"{args.wall_seconds}s\n")

    # The reference does not reproduce under its own `seed=42`: five runs on identical
    # inputs gave five different inventories. Scoring one draw of a stochastic optimiser
    # against our deterministic result would report sampling noise as a difference between
    # the systems, so every run is kept and the spread is reported.
    payloads = []
    for attempt in range(1, args.repeats + 1):
        payload, _ = run_reference(work, args.wall_seconds)
        if payload is None:
            return 1
        payloads.append(payload)
        # Each run's complete inventory, kept as its own artefact. The reference is not
        # reproducible under its seed, so an aggregate alone cannot be replayed: without the
        # individual inventories no one can recheck which run produced which number.
        runs_dir = work / "runs"
        runs_dir.mkdir(parents=True, exist_ok=True)
        (runs_dir / f"run{attempt:02d}.json").write_text(
            json.dumps({**payload,
                        "inputs": {name: _sha256(path) for name, path in
                                   sorted((f.name, f) for f in (work / "input").iterdir()
                                          if f.is_file())},
                        "import_origin": payload.get("origin"),
                        "interpreter": str(REF_PYTHON),
                        "revision": revision},
                       indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"  reference run {attempt}/{args.repeats}: {payload['rows']} modules in "
              f"{payload['seconds']:.1f}s -> runs/run{attempt:02d}.json")

    payload = payloads[0]
    print(f"\nreference v{payload['version']}, config as run:")
    for key, value in sorted(payload["config"].items()):
        print(f"    {key}: {value}")

    pilot = json.loads((ROOT / "work" / "m3" / "inventory_pilot.json")
                       .read_text(encoding="utf-8"))
    ours = {k[len(OUR_PREFIX):] if k.startswith(OUR_PREFIX) else k: v
            for k, v in pilot["modules"].items()}
    originals_full, _ = load_inserts(args.table)
    originals = {k[len(OUR_PREFIX):] if k.startswith(OUR_PREFIX) else k: v
                 for k, v in originals_full.items()}
    table = complete_table(json.loads(TABLE.read_text(encoding="utf-8")), 1)

    runs = [{k[:-len(THEIR_SUFFIX)] if k.endswith(THEIR_SUFFIX) else k: v
             for k, v in p["modules"].items()} for p in payloads]
    theirs = runs[0]
    shared = sorted(set.intersection(*(set(r) for r in runs)) & set(ours))
    print(f"\nmodule identities: {len(theirs)} reference, {len(ours)} ours, "
          f"{len(shared)} shared after normalising names")
    if not shared:
        print("  no shared identities; the comparison cannot proceed on common scope")
        return 1

    distinct = sum(1 for pid in shared
                   if len({r[pid]["cds"] for r in runs}) > 1)
    print(f"  reference runs producing more than one sequence for the same module: "
          f"{distinct}/{len(shared)} -- it does not reproduce under its declared seed")

    # The two systems do not define a module's coding span the same way. Their CDS runs 1-2
    # codons longer because it includes the codon that spans the junction, completed by bases
    # the neighbouring module supplies; ours stops at the insert's own last whole codon. So
    # scoring each system's full output would compare different amino-acid content, and CAI
    # depends on that content.
    #
    # Our protein turns out to be a contiguous substring of theirs on every shared module, so
    # there is a real common region: the codons encoding our protein, in both. Every score
    # below is taken over exactly that region. A module where the alignment fails is reported
    # as unalignable rather than scored on a guess.
    rows, unalignable = [], []
    for pid in shared:
        record = ours[pid]
        frame = record["frame"]
        n = (len(record["sequence"]) - frame) // 3
        our_coding = record["sequence"][frame:frame + 3 * n]

        our_protein = translate(our_coding)
        scores, proteins_agree, extra = [], True, set()
        for run in runs:
            their_full = run[pid]["cds"]
            their_protein = translate(their_full)
            offset = their_protein.find(our_protein)
            if offset < 0:
                break
            their_coding = their_full[3 * offset:3 * (offset + len(our_protein))]
            scores.append(cai(their_coding, table))
            proteins_agree &= translate(their_coding) == our_protein
            extra.add(len(their_protein) - len(our_protein))
        if len(scores) != len(runs):
            unalignable.append(pid)
            continue

        original = originals.get(pid, "")
        orig_coding = original[frame:frame + 3 * n] if original else ""
        rows.append({
            "module": pid,
            "aligned_aa": len(our_protein),
            "reference_extra_codons": sorted(extra),
            "baseline_cai": round(cai(orig_coding, table), 4) if orig_coding else None,
            "ours_cai": round(cai(our_coding, table), 4),
            "reference_cai": round(sum(scores) / len(scores), 4),
            "reference_cai_runs": [round(s, 4) for s in scores],
            "reference_cai_min": round(min(scores), 4),
            "reference_cai_max": round(max(scores), 4),
            "ours_nt": len(our_coding),
            "reference_nt": 3 * len(our_protein),
            "same_protein": proteins_agree,
            "reference_qc": theirs[pid]["qc_status"],
            "reference_hard_constraints": theirs[pid]["hard_constraints_passed"],
        })

    ref_better = sum(1 for r in rows if r["reference_cai"] > r["ours_cai"] + 1e-9)
    we_better = sum(1 for r in rows if r["ours_cai"] > r["reference_cai"] + 1e-9)
    tied = len(rows) - ref_better - we_better
    ref_mean = sum(r["reference_cai"] for r in rows) / len(rows)
    our_mean = sum(r["ours_cai"] for r in rows) / len(rows)
    base = [r["baseline_cai"] for r in rows if r["baseline_cai"] is not None]
    same_length = sum(1 for r in rows if r["ours_nt"] == r["reference_nt"])
    same_protein = sum(1 for r in rows if r["same_protein"])

    extra = sorted({e for r in rows for e in r["reference_extra_codons"]})
    per_run_means = [sum(r["reference_cai_runs"][i] for r in rows) / len(rows)
                     for i in range(len(runs))]

    print(f"\naligned {len(rows)}/{len(shared)} shared modules onto a common coding region; "
          f"unalignable {unalignable or 'none'}")
    print(f"  the reference's span runs {min(extra)}-{max(extra)} codons longer, from the "
          f"junction-spanning codon\n  its neighbour completes; scores below are taken over "
          f"the common region only")
    print(f"\nCAI under our Chlamydomonas table, one scorer applied to both "
          f"({len(runs)} reference runs):")
    if base:
        print(f"  original deposited inventory  mean {sum(base) / len(base):.4f}")
    print(f"  ours                          mean {our_mean:.4f}  (deterministic)")
    print(f"  reference                     mean {ref_mean:.4f}  over runs "
          f"{[f'{m:.4f}' for m in per_run_means]}")
    print(f"  per module, against the reference's per-module mean: reference higher "
          f"{ref_better}, ours higher {we_better}, tied {tied}")
    print(f"  on the common region the two encode the same protein on "
          f"{same_protein}/{len(rows)} (equal length {same_length}/{len(rows)})")
    print(f"  reference QC PASS on its own terms "
          f"{sum(1 for r in rows if r['reference_qc'] == 'PASS')}/{len(rows)}, hard "
          f"constraints {sum(1 for r in rows if r['reference_hard_constraints'])}/{len(rows)}")

    print("\ncoverage and cost -- different amounts of work, so not a speed claim:")
    print(f"  reference {payload['rows']} modules in "
          f"{sum(p['seconds'] for p in payloads) / len(payloads):.1f}s mean "
          f"({sum(p['seconds'] for p in payloads) / len(payloads) / payload['rows']:.2f}s "
          f"each)")
    print(f"  ours      {len(pilot['modules'])} modules in "
          f"{pilot['budget']['elapsed_seconds']:.1f}s "
          f"({pilot['budget']['elapsed_seconds'] / len(pilot['modules']):.2f}s each)")

    print("\ndeclared mismatches (comparison limited to common scope):")
    for key, value in DECLARED_MISMATCHES.items():
        print(f"   {key}: {value}")
    print("   note: both inventories hold AATG and AGGT variants, so the config's\n"
          "         ppr_5prime_fusion_site selects an assembly variant and does not\n"
          "         restrict this module-level comparison.")

    (work / "comparison.json").write_text(json.dumps({
        "reference": {"version": payload["version"], "revision": revision,
                      "modules": payload["rows"], "runs": len(payloads),
                      "seconds_per_run": [round(p["seconds"], 1) for p in payloads],
                      "config_as_run": payload["config"]},
        "ours": {"modules": len(pilot["modules"]),
                 "seconds": pilot["budget"]["elapsed_seconds"]},
        "reproducibility": {
            "ours": "deterministic",
            "reference": (f"{distinct}/{len(shared)} modules received more than one sequence "
                          f"across {len(payloads)} runs at its declared seed=42"),
            "reference_run_means": [round(m, 4) for m in per_run_means],
        },
        "shared_modules": len(rows),
        "unalignable_modules": unalignable,
        "objective": ("CAI under data/codon_tables/kazusa_3055.json, computed here for both "
                      "systems over the aligned common coding region; neither system's "
                      "self-reported score is used"),
        "alignment": ("our protein is a contiguous substring of the reference's; the "
                      "reference span includes the junction-spanning codon its neighbour "
                      "completes, which is excluded from both scores"),
        "means": {"original": round(sum(base) / len(base), 4) if base else None,
                  "ours": round(our_mean, 4), "reference": round(ref_mean, 4)},
        "per_module_wins": {"reference": ref_better, "ours": we_better, "tied": tied},
        "protein_agreement": {"same_length": same_length, "same_protein": same_protein,
                              "of": len(rows)},
        "rows": rows,
        "declared_mismatches": DECLARED_MISMATCHES,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {work / 'comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
