"""Where an intron could go, and whether putting one there would still pass our own rules.

**This module never changes a design.** It reads a coding sequence and reports candidate
positions. Nothing here inserts anything, and `design_oneshot` does not call it. A design
PPROPT produces contains no intron, and the advisory says so in every result it returns.

**Why the question exists.** PPROPT's default host is *Chlamydomonas reinhardtii*, where the
published work is blunt about codon optimisation alone being insufficient: Baier et al. 2018
built a codon-optimised, intron-free gene, got minimal expression, and recovered more than
100-fold the mRNA and about 8-fold the protein by inserting three copies of the *rbcS2*
intron 1. Codon optimisation is the whole of what PPROPT does for expression. See
`docs/intron_feasibility.md` for the assessment, including what it would cost to act on this
rather than merely report it.

**Three hosts are refused, not warned.** `e_coli` and both chloroplast tables have no nuclear
spliceosome, so the intron is never excised and is read as coding sequence. At 145 nt it
frameshifts everything downstream, and it carries a stop codon in all three frames. That
destroys the protein rather than merely failing to improve it, which is the one case where a
quiet default would be unforgivable.

Every figure quoted here was measured by other people in other constructs. Nothing in this
repository has tested whether any of it transfers to a designed PPR.
"""
from __future__ import annotations

import re

from .synthesis_profile import SynthesisProfile, resolve

#: *Chlamydomonas reinhardtii* rbcS2 intron 1, NCBI X04472.1 bases 1307..1451, 145 bp.
#: The canonical enhancing intron for this host and the one PPROPT's synthesis policy was
#: measured against: GC 0.579, no 50-nt window outside 0.35-0.65, longest homopolymer 4, and
#: no BsaI, BbsI or SapI site. It needs no exemption from our own rules.
RBCS2I1 = (
    "GTGAGTCGACGAGCAAGCCCGGCGGATCAGGCAGCGTGCTTGCAGATTTGACTTGCAACGCCCGCATTGT"
    "GTCGACGAAGGCTTTTGGCTCCTCTGTCGCTGTCTCAAGCAGCATCTAACCCTGCGTCGCCGTTTCCATT"
    "TGCAG"
)

#: Hosts with a nuclear spliceosome, where a spliceosomal intron is excised at all.
#: Named explicitly rather than derived from the genetic code: `e_coli` uses code 1 like the
#: nuclear eukaryotes but is a bacterium, so code alone would wrongly admit it.
SPLICEOSOMAL_HOSTS = frozenset({
    "c_reinhardtii_nuclear",
    "s_cerevisiae",
    "a_thaliana_nuclear",
})

#: The single host where intron-mediated enhancement is actually documented for this intron.
#: Having a spliceosome and being shown to benefit are different claims: Baier 2020 found
#: that most foreign introns failed in *Chlamydomonas* until their splice sites were adapted,
#: which is a caution against assuming introns travel between hosts.
EVIDENCE_HOSTS = frozenset({"c_reinhardtii_nuclear"})

#: The enzyme that wraps an intron ordered as its own part.
#:
#: BbsI, not BsaI, and the reason is arithmetic rather than taste. rbcS2i1's own worst 50-nt
#: window is 0.640, ten thousandths under this project's 0.65 ceiling. Wrapping it for
#: ordering shifts the window frame and some frame then includes the enzyme's recognition
#: site: BsaI's `GGTCTC` carries four G or C, BbsI's `GAAGAC` carries three. That one base
#: is 0.660 against 0.640 -- refused against accepted. An intron is not coding, so unlike a
#: fragment it cannot be recoded to compensate, and no choice of the arbitrary pads fixes it.
#:
#: It is also what the published system does: MoInClo flanks its level -1 intron modules
#: with BbsI, and the GRASP kit assembles level 0 with BpiI, which is the same site. The
#: idiomatic choice and the one that passes are the same choice.
PART_ENZYME: str = "BbsI"

#: Shortest exon worth designing. MoInClo's were 234-432 nt and Baier's about 240-583, so
#: this is a floor against arithmetic artefacts rather than a published rule: a search that
#: optimises fidelity alone will happily place two introns one base apart, which is not a
#: design.
MIN_EXON_NT = 100

#: Native *C. reinhardtii* average exon length, and the practical ceiling Baier 2018 gives.
#: The paper's own best construct has a 583 bp exon, so the ceiling is a guideline rather
#: than a threshold even to its authors.
NATIVE_EXON_NT = 240
MAX_EXON_NT = 500

#: What goes on every result, including the ones that succeed.
DISCLAIMER = (
    "This design does not take introns into consideration. The sequence PPROPT produced "
    "contains no intron and was neither optimised for nor checked against one. The "
    "positions reported here are a reading of published placement rules, not a "
    "modification of your design, and no expression effect has been measured here."
)

_STOPS = ("TAA", "TAG", "TGA")


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def part_problems(payload: str, enzyme: str,
                  profile: SynthesisProfile | None = None) -> list[str]:
    """Breaches on `payload` once it is wrapped for ordering, from the real validator.

    Routed through `substrates.synthesis_findings` with this route's own intended-site
    records, because a wrapper deliberately contains the enzyme's recognition site and a
    checker that does not know which sites are intended reports every ordered oligo as
    carrying a forbidden one. `sequence_problems` below is for bare DNA and will do exactly
    that if handed a wrapped oligo.
    """
    from . import substrates as _sub
    from .assembly import wrap_fragment, wrapper_intended_sites

    policy = profile or resolve(None)
    oligo = wrap_fragment(payload, enzyme)
    findings = _sub.synthesis_findings(
        oligo, policy, intended_sites=wrapper_intended_sites(payload, enzyme),
        enzymes=("BsaI", "BbsI", "SapI"))
    return list(findings["hard"])


def sequence_problems(seq: str, profile: SynthesisProfile | None = None) -> list[str]:
    """Every way a **bare** stretch of DNA breaks the declared synthesis rules.

    Not for a wrapped oligo: it has no notion of intended enzyme sites and would report the
    wrapper's own recognition site as a breach. Use `part_problems` for anything ordered.

    Deliberately re-derived here rather than routed through `substrates.synthesis_findings`:
    that function scores *ordered fragments*, which carry a wrapper and an entry overhang and
    expect a fragment record. This answers a narrower question about a bare stretch of DNA.
    """
    policy = profile or resolve(None)
    seq = (seq or "").upper()
    problems: list[str] = []
    if not seq:
        return ["empty sequence"]

    low, high = policy.local_gc
    window = policy.window
    if len(seq) >= window:
        worst = max((_gc(seq[i:i + window]), i) for i in range(len(seq) - window + 1))
        value, where = worst
        least = min(_gc(seq[i:i + window]) for i in range(len(seq) - window + 1))
        if value > high:
            problems.append(f"a {window}-nt window at {where} is {value:.3f} GC, over {high}")
        if least < low:
            problems.append(f"a {window}-nt window is {least:.3f} GC, under {low}")

    run = max((len(m.group()) for m in re.finditer(r"(.)\1*", seq)), default=0)
    if run > policy.max_homopolymer:
        problems.append(f"a homopolymer run of {run}, over {policy.max_homopolymer}")

    for name, site in (("BsaI", "GGTCTC"), ("BbsI", "GAAGAC"), ("SapI", "GCTCTTC")):
        rc = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if site in seq or rc in seq:
            problems.append(f"a {name} site")
    return problems


def junction_problems(cds: str, position: int, intron: str,
                      profile: SynthesisProfile | None = None) -> list[str]:
    """Problems the insertion itself creates, at the two new exon/intron junctions.

    Scoped to the windows that actually straddle the inserted block, and deliberately so.
    Asking instead whether the *whole* spliced sequence is clean conflates two questions: a
    coding sequence that already breaches the band somewhere else would then veto every
    candidate position, reporting "nowhere to put an intron" when the real finding is
    "this sequence has a problem that has nothing to do with introns".

    A real PPROPT coding sequence satisfies the band by construction, so in practice the two
    readings agree. They diverge exactly when something has already gone wrong, which is
    when a misleading answer would cost the most.
    """
    policy = profile or resolve(None)
    cds, intron = (cds or "").upper(), (intron or "").upper()
    spliced = cds[:position] + intron + cds[position:]
    window = policy.window
    low, high = policy.local_gc

    start = max(0, position - window + 1)
    end = min(len(spliced) - window, position + len(intron))
    problems: list[str] = []
    if end >= start:
        values = [(_gc(spliced[i:i + window]), i) for i in range(start, end + 1)]
        worst, where = max(values)
        least, least_at = min(values)

        # What the same neighbourhood looked like before the insertion. A breach the coding
        # sequence already had is not something inserting an intron caused, and reporting it
        # as one would blame the intron for the sequence's own problem.
        b_lo = max(0, position - window)
        b_hi = min(len(cds) - window, position + window)
        base = [_gc(cds[i:i + window]) for i in range(b_lo, b_hi + 1)] or [_gc(cds)]
        base_worst, base_least = max(base), min(base)

        if worst > high >= base_worst:
            problems.append(
                f"inserting here puts a {window}-nt window at {where} at {worst:.3f} GC, "
                f"over {high} (this stretch was {base_worst:.3f} before)")
        if least < low <= base_least:
            problems.append(
                f"inserting here puts a {window}-nt window at {least_at} at {least:.3f} GC, "
                f"under {low} (this stretch was {base_least:.3f} before)")

    # A homopolymer or an enzyme site can only be created across the two new seams.
    for seam, flank in ((position, cds[max(0, position - 8):position] + intron[:8]),
                        (position + len(intron),
                         intron[-8:] + cds[position:position + 8])):
        run = max((len(m.group()) for m in re.finditer(r"(.)\1*", flank)), default=0)
        if run > policy.max_homopolymer:
            problems.append(f"the seam at {seam} creates a homopolymer run of {run}")
        for name, site in (("BsaI", "GGTCTC"), ("BbsI", "GAAGAC"), ("SapI", "GCTCTTC")):
            rc = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            if site in flank or rc in flank:
                problems.append(f"the seam at {seam} creates a {name} site")
    return problems


def frames_with_stops(seq: str) -> list[int]:
    """Reading frames of `seq` that contain a stop codon. Used to show what unspliced costs."""
    seq = (seq or "").upper()
    out = []
    for frame in range(3):
        body = seq[frame:]
        body = body[:len(body) // 3 * 3]
        if any(body[i:i + 3] in _STOPS for i in range(0, len(body), 3)):
            out.append(frame)
    return out


def applicable(organism: str) -> bool:
    """Whether a spliceosomal intron is excised at all in this host."""
    return organism in SPLICEOSOMAL_HOSTS


def has_evidence(organism: str) -> bool:
    """Whether intron-mediated enhancement is documented in this host for this intron."""
    return organism in EVIDENCE_HOSTS


def why_not(organism: str) -> str | None:
    """The reason a host is refused, or None when it is not."""
    if applicable(organism):
        return None
    return (
        f"{organism} has no nuclear spliceosome, so the intron would never be excised. It "
        f"would be read as coding sequence: at {len(RBCS2I1)} nt it frameshifts everything "
        f"downstream, and it carries a stop codon in all three reading frames. Inserting "
        f"one here destroys the protein rather than failing to improve it."
    )


def insertion_sites(cds: str) -> list[int]:
    """Positions where an intron may be inserted, by Baier 2018's NG/GN rule.

    A position `i` means "between `cds[i-1]` and `cds[i]`", and qualifies when both are G.
    The rule is quantitative rather than stylistic: across all sixteen flanking combinations
    a G|G junction gave 132% of control colony yield and C|C gave 3%, with purine-flanked
    sites at 86-106% and pyrimidine-flanked at 15-16%.

    The first and last positions are excluded: an intron before the start codon or after the
    stop is not an intron in the coding sequence.
    """
    cds = (cds or "").upper()
    return [i for i in range(1, len(cds)) if cds[i - 1] == "G" and cds[i] == "G"]


def _overhang(cds: str, position: int) -> str:
    """The four bases a junction at `position` would expose."""
    return cds[position - 4:position]


def _select(cds: str, clean: list[int], spacing: int,
            intron: str) -> tuple[list[int], float | None]:
    """Pick the intron positions whose junction overhangs assemble best.

    Enumerates the sets that satisfy the exon-length rules and scores each with the same
    published mis-ligation data PPROPT uses for its own cut positions, so the number here is
    comparable with the fidelity the design already reports. Sets that repeat an overhang are
    discarded outright rather than scored: two identical four-base ends in one reaction is
    not a low-fidelity design, it is fragments that can swap.
    """
    import itertools

    from .assembly import reverse_complement
    from .overhangs import set_fidelity
    from .assembly_spec import DESTINATION_OVERHANGS

    if not clean:
        return [], None
    dest = DESTINATION_OVERHANGS["level0"]
    fixed = [dest[0], reverse_complement(dest[1])]
    want = max(1, round(len(cds) / spacing) - 1)

    def ok(combo) -> bool:
        bounds = [0] + list(combo) + [len(cds)]
        return all(MIN_EXON_NT <= b - a <= MAX_EXON_NT
                   for a, b in zip(bounds, bounds[1:]))

    best: tuple[float, tuple[int, ...]] | None = None
    for k in range(want, 0, -1):
        if len(clean) > 40 and k > 6:          # keep the enumeration bounded
            continue
        for combo in itertools.combinations(clean, k):
            if not ok(combo):
                continue
            ohs = [_overhang(cds, s) for s in combo]
            if len(set(ohs)) != len(ohs):
                continue
            score = set_fidelity(fixed + ohs, "BsaI-HFv2")
            if best is None or score > best[0]:
                best = (score, combo)
        if best is not None:
            break
    if best is None:
        return [], None
    return list(best[1]), round(best[0], 4)


def advise(cds: str, organism: str, *, intron: str = RBCS2I1,
           spacing: int = NATIVE_EXON_NT,
           profile: SynthesisProfile | None = None) -> dict:
    """Where introns could go in this coding sequence, for this host.

    Returns a report and nothing else. The coding sequence is not included in the result and
    is never modified -- a caller cannot accidentally take an intron-bearing sequence from
    here and order it.

    A candidate position is suggested only when inserting the intron there leaves the whole
    sequence inside the declared synthesis rules. That check matters: a 50-nt window spanning
    an exon/intron boundary sees GC-rich coding on one side, which is the same straddling
    window that refused every de novo design before the ordered fragment was constrained.
    """
    cds = (cds or "").upper()
    policy = profile or resolve(None)
    sites = insertion_sites(cds)

    report = {
        "organism": organism,
        "applicable": applicable(organism),
        "why_not": why_not(organism),
        "evidence_in_this_host": has_evidence(organism),
        "sites_available": len(sites),
        "suggested": [],
        "exon_lengths": [],
        "spacing_nt": spacing,
        "rejected_for_synthesis": [],
        "intron": {
            "name": "rbcS2i1",
            "length_nt": len(intron),
            "source": "C. reinhardtii rbcS2 intron 1, NCBI X04472.1:1307..1451",
            "passes_our_rules": sequence_problems(intron, policy) == [],
        },
        "evidence": (
            "In Chlamydomonas reinhardtii, three copies of rbcS2 intron 1 in a "
            "codon-optimised 1,662 bp gene gave >100-fold mRNA and ~8-fold protein over the "
            "intron-free control (Baier 2018); one copy gives ~5.5-fold (Baier 2020). The "
            "intron-free control was itself codon-optimised, so codon optimisation alone was "
            "not sufficient. Nothing comparable has been measured for a designed PPR."
        ),
        "disclaimer": DISCLAIMER,
    }

    if not report["applicable"]:
        return report

    if not has_evidence(organism):
        report["evidence"] += (
            f" {organism} has a spliceosome, so an intron would be excised, but rbcS2i1 has "
            f"not been shown to enhance expression there. Foreign introns frequently fail: "
            f"Baier 2020 found most needed their splice sites adapted to the host consensus "
            f"first, and the Arabidopsis UBQ10 intron did not enhance in Chlamydomonas at all."
        )

    # Two routes, evaluated separately, because they fail for different reasons and only
    # one of them works. Route A puts the intron inside an exon's oligo, so a window
    # straddling the seam is synthesised as one piece and must satisfy the band. Route B
    # orders the intron as its own part, so that window never exists at synthesis time --
    # what must pass instead is the intron's own wrapped oligo.
    embedded, rejected = [], []
    for position in sites:
        problems = junction_problems(cds, position, intron, policy)
        if problems:
            rejected.append({"position": position, "problems": problems})
        else:
            embedded.append(position)
    report["route_a_embedded"] = {
        "clean_positions": embedded,
        "sites_clean": len(embedded),
        "rejected": rejected[:10],
    }

    clean = []
    for position in sites:
        if position < 4 or position + 4 > len(cds):
            continue
        payload = cds[position - 4:position] + intron + cds[position:position + 4]
        if not part_problems(payload, PART_ENZYME, policy):
            clean.append(position)
    report["route_b_own_part"] = {
        "enzyme": PART_ENZYME,
        "clean_positions": clean,
        "sites_clean": len(clean),
    }
    report["sites_clean"] = len(clean)
    report["clean_positions"] = clean
    report["rejected_for_synthesis"] = rejected[:10]

    # Choose the set, not the first thing that fits.
    #
    # Taking the first clean site past each spacing target is the obvious walk and it is
    # badly wrong: the overhang at an intron junction is the four bases before it, every
    # candidate sits between two Gs, and so the naive walk repeatedly picks the *same*
    # four-mer. Identical overhangs in one Golden Gate reaction let fragments swap places.
    # Measured on the 9S design, that walk gives a predicted fidelity of 0.18 against a
    # baseline of 0.83; on 19S it gives 0.003. Choosing the set for fidelity instead
    # recovers 0.80 on 9S.
    chosen, chosen_fidelity = _select(cds, clean, spacing, intron)
    report["suggested"] = chosen
    report["predicted_fidelity"] = chosen_fidelity

    # Say when the placement is not worth having. Measured across the three architectures
    # with fixed codons: 9S reaches 0.80 against a 0.83 baseline, 14S 0.68, 19S 0.48. The
    # longer the protein the more introns the spacing rule wants, the more junctions compete
    # in one reaction, and the smaller the pool of usable overhangs they must come from.
    if chosen_fidelity is not None:
        if chosen_fidelity < 0.5:
            report["fidelity_warning"] = (
                f"predicted assembly fidelity falls to {chosen_fidelity:.3f}. About half the "
                f"assembled plasmids would carry the wrong arrangement. This placement is "
                f"not usable as it stands.")
        elif chosen_fidelity < 0.75:
            report["fidelity_warning"] = (
                f"predicted assembly fidelity falls to {chosen_fidelity:.3f}, against "
                f"roughly 0.83 without introns. Expect to screen more colonies, and treat "
                f"this as a trade-off to decide rather than a plan to follow.")
    bounds = [0] + chosen + [len(cds)]
    report["exon_lengths"] = [b - a for a, b in zip(bounds, bounds[1:])]

    if not chosen:
        if not clean:
            report["blocked_by"] = (
                f"none of the {len(sites)} candidate positions survives the synthesis rules: "
                f"the intron and the coding sequence each satisfy the GC band alone, but "
                f"every junction between them puts a {policy.window}-nt window over "
                f"{policy.local_gc[1]}. This is the straddling-window problem that refused "
                f"every de novo design before the ordered fragment was constrained, and the "
                f"same remedy applies: the codons would have to be chosen with the intron in "
                f"place, not after it."
            )
        else:
            report["blocked_by"] = (
                f"{len(clean)} position(s) survive the synthesis rules -- {clean} -- but none "
                f"leaves every exon at or under {MAX_EXON_NT} nt, so no placement satisfies "
                f"both rules at once. Re-optimising the codons with the intron in place "
                f"would open more positions; choosing from these anyway means accepting an "
                f"exon longer than the published guideline."
            )
    return report
