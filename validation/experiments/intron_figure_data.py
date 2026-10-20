"""Every number the intron wiki figures draw, generated on the enzyme the reaction uses.

One script so the figures cannot disagree with each other. Designs are cached to
`designs.json`, because 19S takes about nine minutes and nothing downstream of it is slow.

Scored on `introns.PART_MATRIX` (BbsI-HF). Figures built before 2026-10-16 used BsaI-HFv2,
which is the level 1 enzyme -- `assembly_spec.py:43` calls that "a reaction that does not
happen". Those figures read about 0.05 high.
"""
import json
import sys
import time
from pathlib import Path

from ppropt import DESTINATION_OVERHANGS, design_oneshot, introns
from ppropt.assembly import reverse_complement
from ppropt.overhangs import set_fidelity
from ppropt.synthesis_profile import STRICT_LEGACY as POLICY

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "work/introns")
OUT.mkdir(parents=True, exist_ok=True)
MATRIX = introns.PART_MATRIX
WINDOW = POLICY.window
LO, HI = POLICY.local_gc
WRAPPER_NT = 30  # BbsI wrapper, both ends together
OLIGO_CAP_NT = 350

TARGETS = {
    "9S": "AAAAUGUGG",
    "14S": "UUACACGUGCGUAC",
    "19S": "CUAUCACAUCACAUAAGCG",
}
HOSTS = [
    "c_reinhardtii_nuclear", "c_reinhardtii_chloroplast", "n_tabacum_chloroplast",
    "e_coli", "s_cerevisiae", "a_thaliana_nuclear",
]


def gc_trace(seq):
    """GC of every WINDOW-nt window, in order. What the band check actually looks at."""
    return [round((seq[i:i + WINDOW].count("G") + seq[i:i + WINDOW].count("C")) / WINDOW, 4)
            for i in range(len(seq) - WINDOW + 1)]


def longest_run(seq):
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def gc_of(seq):
    return round((seq.count("G") + seq.count("C")) / len(seq), 4)


def design(arch, seed=42):
    """Design once, cache forever. Keyed on arch and seed."""
    cache = OUT / "designs.json"
    have = json.loads(cache.read_text(encoding="utf-8")) if cache.is_file() else {}
    key = arch + ":" + str(seed)
    if key not in have:
        t0 = time.time()
        r = design_oneshot(
            TARGETS[arch], organism="c_reinhardtii_nuclear", codon_table=None,
            genetic_code=None, enzyme_profile="igem_rfc1000", enzyme="BsaI",
            matrix="BsaI-HFv2", destination=DESTINATION_OVERHANGS["level0"],
            seed=seed, check_offtarget=False, outdir=None)
        have[key] = {
            "cds": r["cds"],
            "protein": r["protein"],
            # The internal junctions the design chose, and the whole set the reaction sees.
            # `reaction_overhangs` already carries the acceptor as CTCA/CTCG -- the pair that
            # actually meets in the tube, not the CTCA/CGAG the spec stores.
            "junction_overhangs": list(r["junction_overhangs"]),
            "reaction_overhangs": list(r["reaction_overhangs"]),
            # `oligos` is a DataFrame, so iterating it would yield column names.
            "oligo_nt": [len(s) for s in r["oligos"]["oligo_sequence_5to3"]],
            "payload_nt": [len(s) for s in r["oligos"]["payload_5to3"]],
            "planned_fidelity": r["audit"].fidelity,
            "planned_ceiling": r["audit"].fidelity_ceiling,
            "seconds": round(time.time() - t0, 1),
        }
        cache.write_text(json.dumps(have, indent=1), encoding="utf-8")
        print("  designed " + key + " in " + str(have[key]["seconds"]) + "s", flush=True)
    return have[key]


# ---------------------------------------------------------------- figure 1: the two introns
def figure_introns():
    """rbcS2i1 against LHCBM1i2, on PPROPT's own band. No assembly, no matrix."""
    out = {"window_nt": WINDOW, "band": [LO, HI], "introns": {}}
    both = (
        ("rbcS2i1", introns.RBCS2I1, "X04472.1:1307..1451", 5.5),
        ("LHCBM1i2", introns.LHCBM1I2, "AF495473.1:527..779", 8.0),
    )
    for name, seq, accession, fold in both:
        trace = gc_trace(seq)
        problems = introns.sequence_problems(seq, POLICY)
        out["introns"][name] = {
            "length_nt": len(seq),
            "accession": accession,
            "fold_enhancement": fold,
            "overall_gc": gc_of(seq),
            "gc_trace": trace,
            "worst_window_gc": max(trace),
            "worst_window_at": trace.index(max(trace)),
            "best_window_gc": min(trace),
            "margin_to_ceiling": round(HI - max(trace), 4),
            "windows": len(trace),
            "windows_breaching": sum(1 for g in trace if g < LO or g > HI),
            "longest_homopolymer": longest_run(seq),
            "starts_gt": seq[:2] == "GT",
            "ends_ag": seq[-2:] == "AG",
            "passes": problems == [],
            "problems": problems,
            "frame_shift_if_unspliced": len(seq) % 3,
            "frames_with_stops": introns.frames_with_stops(seq),
        }
    return out


# ------------------------------------------------------- figure 2+3: sites, placement, cost
def figure_architectures():
    """Per architecture: where an intron can go, which are orderable, what it costs."""
    rows = {}
    for arch in TARGETS:
        d = design(arch)
        cds = d["cds"]
        report = introns.advise(cds, "c_reinhardtii_nuclear")
        # Baseline on the right matrix. The design was planned on BsaI-HFv2, so the
        # fidelity it reported for itself is not what this reaction gets.
        baseline = set_fidelity(d["reaction_overhangs"], MATRIX)
        exons = report["exon_lengths"]
        oligos = [e + WRAPPER_NT for e in exons]
        # Route A keeps the design's own fragments and grows one of them by the intron, so
        # its worst oligo is the largest baseline oligo plus 145 nt.
        route_a_worst = max(d["oligo_nt"]) + len(introns.RBCS2I1)
        rows[arch] = {
            "cds_nt": len(cds),
            "protein_aa": len(d["protein"]),
            "sites_gg": report["sites_available"],
            "site_positions": introns.insertion_sites(cds),
            "nt_per_site": round(len(cds) / max(report["sites_available"], 1), 1),
            "route_a_clean": report["route_a_embedded"]["sites_clean"],
            "route_a_positions": report["route_a_embedded"]["clean_positions"],
            "route_b_clean": report["route_b_own_part"]["sites_clean"],
            "route_b_positions": report["route_b_own_part"]["clean_positions"],
            "baseline_fidelity": round(baseline, 4),
            "baseline_overhangs": d["reaction_overhangs"],
            "baseline_oligo_nt": d["oligo_nt"],
            "planned_fidelity_wrong_matrix": round(d["planned_fidelity"], 4),
            "chosen": report["suggested"],
            "chosen_fidelity": report["predicted_fidelity"],
            "margin_to_floor": report["margin_to_floor"],
            "fidelity_by_count": report["fidelity_by_intron_count"],
            "exon_lengths": exons,
            "exon_oligo_nt": oligos,
            "intron_oligo_nt": len(introns.RBCS2I1) + 8 + WRAPPER_NT,
            "max_oligo_nt": max(oligos, default=0),
            "route_a_worst_oligo_nt": route_a_worst,
            "route_a_fits": route_a_worst <= OLIGO_CAP_NT,
            "warning": report.get("fidelity_warning"),
            "blocked_by": report.get("blocked_by"),
            "design_seconds": d["seconds"],
        }
        print("  " + arch + ": " + str(report["sites_available"]) + " sites, "
              + str(len(report["suggested"])) + " introns, "
              + str(report["predicted_fidelity"])
              + " (baseline " + format(baseline, ".4f") + ")", flush=True)
    return {"matrix": MATRIX, "floor": introns.FIDELITY_FLOOR,
            "oligo_cap_nt": OLIGO_CAP_NT, "wrapper_nt": WRAPPER_NT,
            "architectures": rows}


# --------------------------------------------------- figure 4: where the fidelity goes
def figure_decomposition():
    """Score the acceptor alone, the junctions alone, and both together.

    The point of the figure: the junctions are cheap and the backbone is not.
    """
    d = design("9S")
    cds = d["cds"]
    report = introns.advise(cds, "c_reinhardtii_nuclear")
    chosen = report["suggested"]
    junctions = [introns._overhang(cds, s) for s in chosen]
    frag = d["junction_overhangs"]

    rows = {}
    for level, pair in DESTINATION_OVERHANGS.items():
        # What meets in the tube is the 5' site and the *reverse complement* of the 3' one.
        # At level 0 that is CTCA against CTCG, one base apart, which is the whole finding.
        # Taking the stored pair literally scores a reaction that cannot occur.
        meeting = [pair[0], reverse_complement(pair[1])]
        rows[level] = {
            "stored": list(pair),
            "ends_that_meet": meeting,
            "acceptor_alone": round(set_fidelity(meeting, MATRIX), 4),
            # Route B cuts the exons *at* the intron positions, so the intron junctions
            # replace the design's fragment junctions rather than joining them. Adding both
            # sets double-counts, and on 9S it also repeats AATG, which reads as a collapse
            # to 0.14 that no real reaction would show.
            "with_introns": round(set_fidelity(meeting + junctions, MATRIX), 4),
            # The intron-free design, for the same acceptor, as the thing to beat.
            "intron_free": round(set_fidelity(meeting + frag, MATRIX), 4),
        }
    rows["none"] = {
        "stored": [],
        "ends_that_meet": [],
        "acceptor_alone": None,
        "with_introns": round(set_fidelity(junctions, MATRIX), 4),
        "intron_free": round(set_fidelity(frag, MATRIX), 4),
    }
    return {"matrix": MATRIX, "arch": "9S", "chosen": chosen,
            "junction_overhangs": junctions, "fragment_overhangs": frag,
            "vector": "pIM 046, confirmed by Jan Vollmer 2026-10-16",
            "note": ("`with_introns` is the Route B reaction: acceptor plus the intron "
                     "junctions, which are the only internal junctions it has. "
                     "`intron_free` is the same acceptor with the design's own fragment "
                     "junctions instead."),
            "levels": rows}


# ------------------------------------------------- figure 5: the seed is not the answer
def figure_seeds(seeds=(1, 7, 42, 99, 2024)):
    """Same target, same host, five seeds. How much of the headline is the draw."""
    rows = {}
    for seed in seeds:
        d = design("9S", seed=seed)
        report = introns.advise(d["cds"], "c_reinhardtii_nuclear")
        rows[str(seed)] = {
            "cds_gc": gc_of(d["cds"]),
            "sites_gg": report["sites_available"],
            "route_a_clean": report["route_a_embedded"]["sites_clean"],
            "route_b_clean": report["route_b_own_part"]["sites_clean"],
            "introns": len(report["suggested"]),
            "positions": report["suggested"],
            "fidelity": report["predicted_fidelity"],
            "margin_to_floor": report["margin_to_floor"],
            "fidelity_by_count": report["fidelity_by_intron_count"],
        }
        print("  seed " + str(seed) + ": " + str(len(report["suggested"]))
              + " introns, " + str(report["predicted_fidelity"]), flush=True)
    return {"matrix": MATRIX, "arch": "9S", "floor": introns.FIDELITY_FLOOR, "seeds": rows}


# --------------------------------------------------------------- figure 6: the host gate
def figure_hosts():
    """Which hosts an intron is even defined for, and what happens in the ones it is not."""
    i = introns.RBCS2I1
    rows = {}
    for host in HOSTS:
        rows[host] = {
            "spliceosome": introns.applicable(host),
            "evidence": introns.has_evidence(host),
            "why_not": introns.why_not(host),
        }
    return {
        "intron_nt": len(i),
        "frame_shift": len(i) % 3,
        "frames_with_stops": introns.frames_with_stops(i),
        "consequence": ("145 nt is not a multiple of three, and the intron carries a stop "
                        "in all three frames, so a host that cannot splice it does not make "
                        "less protein. It makes none."),
        "hosts": rows,
    }


if __name__ == "__main__":
    figures = {}
    steps = (
        ("introns", figure_introns),
        ("architectures", figure_architectures),
        ("decomposition", figure_decomposition),
        ("seeds", figure_seeds),
        ("hosts", figure_hosts),
    )
    for name, fn in steps:
        print("[" + name + "]", flush=True)
        figures[name] = fn()
    figures["generated"] = time.strftime("%Y-%m-%d")
    figures["matrix"] = MATRIX
    path = OUT / "figure_data.json"
    path.write_text(json.dumps(figures, indent=1), encoding="utf-8")
    print("\nwrote " + str(path))
