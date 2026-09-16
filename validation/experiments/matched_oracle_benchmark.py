"""The matched oracle benchmark — both systems' *ordered sequences*, under both rule sets.

This is the acceptance requirement that stayed open through four review passes, and the
reason it stayed open is worth stating: until the substrate scope was settled, our deliverable
was a bare insert and theirs was a wrapped oligo, so any head-to-head would have compared two
different kinds of object. Both systems now emit a synthesis-ready sequence, so they can be
held to the same checker.

**What is compared.** One artefact per module per system: the sequence a vendor would be asked
to make. Theirs is `oligo_sequence_5to3` from `optimized_library.csv`; ours is
`substrates.build(...).sequence`. All 42 modules join by name.

**Neither system is graded only on its home rules.** The checker below is written here, imports
nothing from `clippr` or from the reference, and is run twice over both systems' output:

    CLIPPR regime      GC 0.35-0.65 over 50 nt, homopolymer <= 4, no BsaI/BbsI/SapI site
                       outside the two that release the fragment
    reference regime   global GC 0.25-0.65, windowed GC 0.15-0.85 over 50 nt,
                       homopolymer <= 3, no SapI/BsaI/BpiI site outside the intended pair

Reporting only "each system passes its own rules" would be two self-assessments printed side
by side. The cross terms are the benchmark.

**Intended sites are excluded by role, per system, and by site rather than by name.** Their
wrapper uses BsaI, ours uses BbsI. An occurrence is exempt only at the two positions that
release that system's fragment; a third anywhere is a breach for either. The exemption keys on
the recognition sequence because BpiI and BbsI are isoschizomers of GAAGAC -- an enzyme-name
test reported our own two wrapper sites as unintended BpiI in all 42 modules, which was an
artefact of the checker and not a finding about the design.

**Declared asymmetries**, stated rather than scored away:

  * Their oligos are one run of a stochastic search; ours are deterministic. Their spread
    across five preserved runs is reported alongside, on the CDS, which is the quantity their
    run records preserve.
  * They emit vendor-acceptance and QC flags we do not; we emit per-reaction ligation fidelity
    and pooled ordering they do not. Neither is counted for or against the other.
  * **The two wrappers are for different steps of the same workflow.** Digesting their oligos
    shows all 42 release `ACAT ... TTGT` -- the **level -1** destination, i.e. the synthetic
    insert for building the entry clone in pAGM1311. Ours is a level-0 BbsI cassette releasing
    the fragment an assembly consumes. A module plasmid is *made* at level -1 and *used* at
    level 0, so both artefacts are real and neither supersedes the other. Wrapper length and GC
    are therefore reported and never ranked.

    python validation/experiments/matched_oracle_benchmark.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

REFERENCE_CSV = ROOT / "work" / "m3" / "reference" / "output" / "optimized_library.csv"
REFERENCE_RUNS = ROOT / "work" / "m3" / "reference" / "runs"
OUR_INVENTORY = ROOT / "work" / "phaseb" / "inventory_recoded.json"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"

OUR_PREFIX = "pPR-1_"
THEIR_SUFFIX = "_v1"

#: Recognition sites, written out rather than imported, so agreement with the production
#: constants is corroboration and not an echo.
SITES = {"BsaI": "GGTCTC", "BbsI": "GAAGAC", "BpiI": "GAAGAC", "SapI": "GCTCTTC"}

REGIMES = {
    "clippr": {
        "window": 50, "window_gc": (0.35, 0.65), "global_gc": None,
        "max_homopolymer": 4, "enzymes": ("BsaI", "BbsI", "SapI"),
    },
    "reference": {
        "window": 50, "window_gc": (0.15, 0.85), "global_gc": (0.25, 0.65),
        "max_homopolymer": 3, "enzymes": ("SapI", "BsaI", "BpiI"),
    },
}

#: Which enzyme each system's wrapper legitimately carries twice.
WRAPPER_ENZYME = {"ours": "BbsI", "reference": "BsaI"}

#: Breach kinds that a *declared* threshold difference fully accounts for. A cross-regime
#: failure outside this set is not a policy difference; it is a defect.
DECLARED_DIFFERENCES = ("homopolymer longer than", "window GC", "global GC")


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def occurrences(text: str, needle: str) -> list[int]:
    found, at = [], text.find(needle)
    while at != -1:
        found.append(at)
        at = text.find(needle, at + 1)
    return found


def breaches(sequence: str, regime: dict, wrapper_enzyme: str) -> list[str]:
    """Every way `sequence` fails `regime`, as an ordered sequence.

    `wrapper_enzyme` names the enzyme whose site the construction legitimately carries; its
    first and last occurrence are exempt, and nothing else is. Exempting *all* occurrences of
    the wrapper enzyme would hide exactly the defect an internal site represents.
    """
    problems = []
    rules = regime

    # Exemption is by *site*, not by enzyme name. BpiI and BbsI are isoschizomers reading the
    # same GAAGAC, so an enzyme-name test flagged our own two wrapper sites as unintended BpiI
    # in all 42 modules -- a checker artefact, not a finding.
    wrapper_site = SITES[wrapper_enzyme]
    for name in rules["enzymes"]:
        site = SITES[name]
        hits = sorted(occurrences(sequence, site) + occurrences(sequence, rc(site)))
        if site == wrapper_site:
            if len(hits) < 2:
                problems.append(f"{name}: {len(hits)} sites, the wrapper needs 2")
            for at in hits[1:-1]:
                problems.append(f"internal {name} site at {at}")
        else:
            for at in hits:
                problems.append(f"unintended {name} site at {at}")

    if rules["global_gc"]:
        lo, hi = rules["global_gc"]
        gc = (sequence.count("G") + sequence.count("C")) / len(sequence)
        if not lo <= gc <= hi:
            problems.append(f"global GC {gc:.3f} outside {lo}-{hi}")

    lo, hi = rules["window_gc"]
    width = rules["window"]
    for start in range(0, max(1, len(sequence) - width + 1)):
        window = sequence[start:start + width]
        if len(window) < width:
            break
        gc = (window.count("G") + window.count("C")) / width
        if not lo <= gc <= hi:
            problems.append(f"window GC {gc:.3f} at {start} outside {lo}-{hi}")
            break

    run, previous = 1, ""
    for base in sequence:
        run = run + 1 if base == previous else 1
        previous = base
        if run > rules["max_homopolymer"]:
            problems.append(f"homopolymer longer than {rules['max_homopolymer']}")
            break
    return problems


def align_coding_spans(ours: str, theirs: str) -> dict:
    """Locate our translated span inside theirs, or refuse.

    **This is why the first version of this benchmark was wrong.** Only 4 of 42 modules have
    identical translated spans; in the other 38 the reference carries additional boundary
    amino acids. Scoring `cai(ours)` against `cai(theirs)` therefore compared codon adaptation
    over different regions and reported a mean of 0.6348 and 29 wins where the aligned figures
    are 0.6587 and 22.

    Neither representation is wrong -- they include different boundary context. Comparing them
    unaligned is what was wrong, and it is the same error as checking a constraint on a
    sub-span of what ships, moved into the comparison.

    An ambiguous or absent containment is refused rather than resolved by picking an offset.
    """
    from Bio.Seq import Seq

    ours_aa = str(Seq(ours).translate())
    theirs_aa = str(Seq(theirs).translate())
    offsets = [i for i in range(len(theirs_aa) - len(ours_aa) + 1)
               if theirs_aa[i:i + len(ours_aa)] == ours_aa]
    if len(offsets) != 1:
        return {"aligned": False,
                "reason": ("our protein is absent from theirs" if not offsets
                           else f"our protein occurs {len(offsets)} times in theirs"),
                "ours_aa": len(ours_aa), "reference_aa": len(theirs_aa)}
    start = offsets[0] * 3
    return {"aligned": True, "identical_span": ours_aa == theirs_aa,
            "ours_aa": len(ours_aa), "reference_aa": len(theirs_aa),
            "reference_nt_bounds": [start, start + len(ours) ],
            "reference_trimmed_aa": len(theirs_aa) - len(ours_aa),
            "reference_span": theirs[start:start + len(ours)]}


def cai(dna: str, table) -> float:
    """One scorer for both systems, over whole codons from the start of the given span."""
    import math

    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    weights = {c.replace("U", "T"): v / max(d.values())
               for aa, d in table.items() for c, v in d.items() if max(d.values()) > 0}
    used = [weights[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3)
            if dna[i:i + 3] in weights and forward.get(dna[i:i + 3]) not in ("M", "W", None)]
    if not used or min(used) <= 0:
        return 0.0
    return math.exp(sum(math.log(w) for w in used) / len(used))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "work" / "oracle" / "matched_benchmark.json"))
    args = ap.parse_args()

    for required in (REFERENCE_CSV, OUR_INVENTORY):
        if not required.is_file():
            print(f"missing {required}; run the reference comparison and Phase B first")
            return 2

    from clippr import inventories as inv
    from clippr import substrates as sub

    table = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
    ours_inv = inv.load(OUR_INVENTORY)

    theirs = {row["optimized_part_id"].removesuffix(THEIR_SUFFIX): row
              for row in csv.DictReader(REFERENCE_CSV.open(encoding="utf-8"))}
    ours = {s.module_id.removeprefix(OUR_PREFIX): s for s in sub.substrates_for(ours_inv)}

    shared = sorted(set(ours) & set(theirs))
    print(f"ours {len(ours)} | reference {len(theirs)} | joined {len(shared)}")
    if set(ours) ^ set(theirs):
        print(f"  unmatched: ours-only {sorted(set(ours) - set(theirs))}, "
              f"reference-only {sorted(set(theirs) - set(ours))}")

    rows, counts, unalignable = [], {}, []
    for system in ("ours", "reference"):
        for regime in REGIMES:
            counts[f"{system}/{regime}"] = 0

    for module in shared:
        our_seq = ours[module].sequence
        their_seq = theirs[module]["oligo_sequence_5to3"]
        record = ours_inv.modules[OUR_PREFIX + module]
        our_cds = record.dna[slice(*record.coding_interval)]
        their_cds = theirs[module]["optimized_cds"]
        alignment = align_coding_spans(our_cds, their_cds)
        row = {"module": module,
               "ours_nt": len(our_seq), "reference_nt": len(their_seq),
               "alignment": alignment,
               "ours_cai": round(cai(our_cds, table), 6),
               # Scored over the span both systems have in common. The whole-CDS figure is
               # kept beside it, labelled, so the difference the alignment makes is visible
               # rather than quietly corrected away.
               "reference_cai": (round(cai(alignment["reference_span"], table), 6)
                                 if alignment["aligned"] else None),
               "reference_cai_unaligned": round(cai(their_cds, table), 6)}
        if not alignment["aligned"]:
            unalignable.append({"module": module, **alignment})
        for system, sequence in (("ours", our_seq), ("reference", their_seq)):
            for regime, rules in REGIMES.items():
                found = breaches(sequence, rules, WRAPPER_ENZYME[system])
                row[f"{system}_{regime}"] = found
                if not found:
                    counts[f"{system}/{regime}"] += 1
        rows.append(row)

    print("\nclean ordered sequences, out of "
          f"{len(shared)} — each system under BOTH rule sets:")
    print(f"{'':12s}{'CLIPPR rules':>16s}{'reference rules':>18s}")
    for system in ("ours", "reference"):
        print(f"  {system:10s}{counts[f'{system}/clippr']:>16d}"
              f"{counts[f'{system}/reference']:>18d}")

    def offenders(system, regime):
        return {r["module"]: r[f"{system}_{regime}"][:2]
                for r in rows if r[f"{system}_{regime}"]}

    detail = {f"{s}/{g}": offenders(s, g) for s in ("ours", "reference") for g in REGIMES}
    for key, bad in detail.items():
        if bad:
            first = next(iter(bad.items()))
            print(f"  {key}: {len(bad)} failing, e.g. {first[0]} -> {first[1]}")

    # Every cross-regime failure is attributed to a declared threshold difference or it
    # is not. Requiring clean output under the *other* system's thresholds would make the
    # benchmark a test of whether we adopted their policy; an unexplained breach is a real
    # defect and is what this gate is for.
    unexplained: dict[str, dict[str, str]] = {}
    for system, other in (("ours", "reference"), ("reference", "clippr")):
        for row in rows:
            for problem in row[f"{system}_{other}"]:
                if not problem.startswith(DECLARED_DIFFERENCES):
                    unexplained.setdefault(f"{system}/{other}", {})[row["module"]] = problem
    print(chr(10) + "cross-regime failures traced to a declared threshold difference: "
          f"{'all' if not unexplained else 'NOT all'}")
    for key, bad in unexplained.items():
        print(f"  {key}: {len(bad)} unexplained, e.g. {next(iter(bad.items()))}")

    # --- codon adaptation, one scorer over the span both systems share ---
    scored = [r for r in rows if r["alignment"]["aligned"]]
    identical = sum(r["alignment"]["identical_span"] for r in scored)
    print(f"{chr(10)}coding-span alignment: {identical}/{len(rows)} identical as "
          f"translated, {len(scored) - identical} contained unambiguously, "
          f"{len(unalignable)} unalignable")
    for bad in unalignable:
        print("   REFUSED " + bad["module"] + ": " + bad["reason"])

    ours_cai = [r["ours_cai"] for r in scored]
    theirs_cai = [r["reference_cai"] for r in scored]
    theirs_unaligned = [r["reference_cai_unaligned"] for r in scored]
    wins = {"ours": sum(a > b for a, b in zip(ours_cai, theirs_cai)),
            "reference": sum(b > a for a, b in zip(ours_cai, theirs_cai)),
            "tied": sum(a == b for a, b in zip(ours_cai, theirs_cai))}
    print(f"CAI over the shared span: ours {statistics.mean(ours_cai):.6f}, "
          f"reference {statistics.mean(theirs_cai):.6f}")
    print(f"  per module: ours ahead on {wins['ours']}, reference ahead on "
          f"{wins['reference']}, tied {wins['tied']}")
    print(f"  (their whole CDS instead reads {statistics.mean(theirs_unaligned):.6f} "
          f"over a different region, which is not a comparison)")

    # --- the same alignment, over every preserved run rather than one saved CSV ---
    per_run = []
    for path in sorted(REFERENCE_RUNS.glob("run*.json")):
        modules = json.loads(path.read_text(encoding="utf-8"))["modules"]
        values, ahead = [], 0
        for row in scored:
            entry = modules.get(row["module"])
            if not entry:
                continue
            record = ours_inv.modules[OUR_PREFIX + row["module"]]
            got = align_coding_spans(record.dna[slice(*record.coding_interval)],
                                     entry["cds"])
            if not got["aligned"]:
                continue
            value = cai(got["reference_span"], table)
            values.append(value)
            ahead += row["ours_cai"] > value
        if values:
            per_run.append({"run": path.stem, "modules": len(values),
                            "reference_mean": round(statistics.mean(values), 6),
                            "ours_ahead": ahead,
                            "reference_ahead": len(values) - ahead})
    if per_run:
        print(f"{chr(10)}the same alignment over each preserved reference run:")
        for entry in per_run:
            print(f"   {entry['run']}  reference mean {entry['reference_mean']:.6f}  "
                  f"ours ahead on {entry['ours_ahead']} of {entry['modules']}")
        print("   The headline is one saved CSV; these are the runs behind it. A "
              "stochastic optimiser's" + chr(10) + "   single output is an observation, not its level.")

    # --- determinism, on the quantity the preserved runs record ---
    run_files = sorted(REFERENCE_RUNS.glob("run*.json"))
    spread = None
    if run_files:
        per_module: dict[str, set[str]] = {}
        for path in run_files:
            for name, record in json.loads(path.read_text(encoding="utf-8"))["modules"].items():
                per_module.setdefault(name, set()).add(record["cds"])
        varying = sum(1 for v in per_module.values() if len(v) > 1)
        spread = {"runs": len(run_files), "modules": len(per_module), "varying": varying}
        print(f"\ndeterminism: ours deterministic by construction; reference gave more than "
              f"one CDS for {varying}/{len(per_module)} modules across {len(run_files)} runs "
              f"at its declared seed")

    passed = counts["ours/clippr"] == len(shared) and not unexplained
    print(f"{chr(10)}BENCHMARK: {'PASS' if passed else 'FAIL'} "
          "(ours clean under its own declared rules, and every cross-regime failure on "
          "either side traced to a declared threshold difference)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "compared": "the ordered sequence each system emits per module",
        "modules": len(shared),
        "inputs": {
            "our_inventory_version": ours_inv.version,
            "reference_csv_sha256": sha256(REFERENCE_CSV),
            "codon_table_sha256": sha256(CODON_TABLE),
        },
        "regimes": REGIMES,
        "wrapper_enzyme": WRAPPER_ENZYME,
        "clean_counts": counts,
        "offenders": detail,
        "cai": {"scope": "scored over the translated span both systems share",
                "modules_scored": len(scored),
                "identical_spans": identical,
                "unalignable": unalignable,
                "ours_mean": round(statistics.mean(ours_cai), 6),
                "reference_mean": round(statistics.mean(theirs_cai), 6),
                "reference_mean_unaligned_do_not_quote":
                    round(statistics.mean(theirs_unaligned), 6),
                "per_module_wins": wins,
                "per_preserved_run": per_run},
        "reference_spread": spread,
        "unexplained_cross_regime_failures": unexplained,
        "declared_asymmetries": [
            "their oligos are one run of a stochastic search; ours are deterministic",
            "they emit vendor-acceptance and QC flags we do not; we emit per-reaction "
            "ligation fidelity and pooled ordering they do not",
            "their oligos all release ACAT...TTGT, the level -1 destination: they are the "
            "synthetic inserts for building the entry clone. Ours is a level-0 BbsI cassette "
            "releasing the fragment an assembly consumes. Different steps of one workflow, "
            "so wrapper length and GC are reported and never ranked",
        ],
        "rows": rows,
        "benchmark": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
