"""One call from a target RNA to an orderable, QC'd assembly plan.

The pipeline is `ppr` -> `arelf` -> `overhangs` -> `codons` -> `assembly` -> `qc` ->
`export`, and the ordering is forced by a dependency that is easy to get backwards: the
Golden Gate overhangs must be chosen *before* the coding sequence is optimised, because
they are then locked into it. Optimise first and the codon optimiser is free to change the
very bases the junctions depend on.

Cut selection is bounded rather than exhaustive. There are thousands of legal ways to
split a PPR into fragments, but predicted fidelity reaches its ceiling of 1.000 easily, so
the search tries the most balanced splits first and stops at the first proven optimum.
Balanced splits are tried first because fragment cost and synthesis risk both rise with
the longest fragment, not the average.
"""
from __future__ import annotations

import hashlib as _hashlib


from collections.abc import Mapping
from pathlib import Path

from . import constants as C
from .arelf import (achievable_overhangs, balanced_cuts, realization_count,
                    safe_overhangs)
from .audit import DesignAudit, OverhangDecision
from .assembly import build_oligos, split_cds
from .codons import (optimize_cds, table_from_cds_fasta, table_from_csv,
                     table_from_genome, table_from_kazusa)
from .experiment import stage as _stage
from .export import (write_fasta, write_gene_fasta, write_genbank,
                     write_oligo_csv)
from .overhangs import (best_set, enumerate_candidates, fidelity_components,
                        reaction_overhangs, set_fidelity)
from .ppr import describe
from .qc import synthesis_qc

#: Host and fragment-length choices live in `policy.py`, not here -- they are decisions
#: this project made, not facts about GRASP, and the split exists so that stays visible.
#:
#: `ORGANISMS` holds only the verified nuclear entry: taxid 3055 returns a
#: 21-amino-acid Chlamydomonas nuclear table with Leu CUG at 0.73. No Kazusa taxid for the
#: Chlamydomonas *chloroplast* table has been confirmed, and guessing one would silently
#: produce wrong DNA -- the whole failure this package is built to avoid. For chloroplast
#: work, pass `codon_table=` a CDS FASTA or CSV and `genetic_code=11`.
ORGANISMS = C.ORGANISMS
MAX_FRAGMENT_AA = C.MAX_FRAGMENT_AA
MIN_FRAGMENT_AA = C.MIN_FRAGMENT_AA

#: How many candidate splits to score for fidelity before settling.
SPLIT_BUDGET: int = 40

#: How many of those splits to run through codon optimisation before giving up. Each
#: attempt costs a full optimisation. The usual reason for a second attempt -- an overhang
#: that forces a blacklisted site -- is now prevented up front by `arelf.safe_overhangs`,
#: so this is a backstop for genuinely hard proteins, not the common path.
OPTIMISE_ATTEMPTS: int = 3


def _resolve_codon_table(codon_table, organism: str, genetic_code: int):
    """Accept a dict, a path to a CSV or CDS FASTA, or fall back to the named organism."""
    if isinstance(codon_table, Mapping):
        return dict(codon_table)
    if codon_table is not None:
        path = Path(codon_table)
        if not path.exists():
            raise FileNotFoundError(f"codon table not found: {path}")
        if path.suffix.lower() == ".csv":
            return table_from_csv(path)
        return table_from_cds_fasta(path, genetic_code=genetic_code)
    if organism not in ORGANISMS:
        raise ValueError(
            f"unknown organism {organism!r}; known: {sorted(ORGANISMS)}. For any other "
            f"host, pass codon_table= a CDS FASTA or CSV together with the right "
            f"genetic_code."
        )
    source, _ = ORGANISMS[organism]
    # An int is a Kazusa lookup; a str is a genome accession whose annotated genes are
    # counted. Organelles take the second path because Kazusa has no entry for them.
    return (table_from_genome(source, genetic_code) if isinstance(source, str)
            else table_from_kazusa(source))


def _plan_fragments(protein: str, n_fragments: int, destination, matrix: str,
                    enzyme_profile=C.DEFAULT_ENZYME_PROFILE,
                    want: int | None = None):
    """Pick cut positions and their overhangs together, best balance first.

    Returns (cuts, overhangs, fidelity). Stops at the first split whose reaction reaches
    the metric's ceiling; otherwise returns the best seen within `SPLIT_BUDGET`.

    `want` is how many ceiling-scoring plans to collect before stopping, defaulting to
    `OPTIMISE_ATTEMPTS`. `design_oneshot` leaves it alone so its behaviour is unchanged;
    `search.explore` raises it, because three plans is barely a search.
    """
    want = OPTIMISE_ATTEMPTS if want is None else want
    splits = balanced_cuts(protein, n_fragments, MIN_FRAGMENT_AA, MAX_FRAGMENT_AA,
                           limit=SPLIT_BUDGET)
    if not splits:
        raise ValueError(
            f"no legal split of {len(protein)} residues into {n_fragments} fragments of "
            f"{MIN_FRAGMENT_AA}-{MAX_FRAGMENT_AA} aa")

    fixed = reaction_overhangs([], "level0") if destination is None else [
        destination[0], _rc(destination[1])]

    # Adding overhangs to a reaction can only take probability away, so the destination
    # pair scored alone is the ceiling any split can reach -- and it is reachable, when
    # the junctions cross-talk with nothing. Stopping there rather than at 1.000 matters:
    # the level 0 pair alone scores about 0.83, so a search waiting for 1.000 would never
    # exit early and would score every candidate split for nothing.
    ceiling = set_fidelity(fixed, matrix)

    out = []
    for cuts in splits:
        # safe_overhangs, not achievable_overhangs: an overhang that forces a blacklisted
        # site hands the codon optimiser an unsatisfiable problem, and no retry recovers.
        try:
            options = enumerate_candidates(
                [",".join(safe_overhangs(protein, c, enzymes=enzyme_profile))
                 for c in cuts])
        except ValueError:
            continue
        try:
            chosen, score = best_set(options, fixed=fixed, matrix=matrix)
        except ValueError:
            continue
        out.append((list(cuts), chosen, score))
        # Enough plans already at the ceiling: nothing later can score better, and each
        # extra split costs a full overhang enumeration. Scoring all of them regressed a
        # 19S design from 6 s to 190 s.
        if sum(1 for _, _, s in out if s >= ceiling - 1e-12) >= want:
            break
    if not out:
        raise ValueError("no split produced a valid overhang set")
    # best fidelity first; ties keep the balanced order balanced_cuts produced
    out.sort(key=lambda x: -x[2])
    return out, ceiling


def _rc(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _parts_plan(target_rna: str):
    """Can the deposited GRASP kit realise this target? Never raises.

    A design should still be produced when the kit cannot help, so a failure here is recorded
    as an unavailable route rather than allowed to stop the de novo pipeline.
    """
    from .parts import select

    return select(target_rna)


def plan_candidates(protein: str, n_fragments: int, destination, matrix: str,
                    enzyme_profile=C.DEFAULT_ENZYME_PROFILE, want: int | None = None):
    """Every assembly plan this design would consider, best predicted fidelity first.

    `design_oneshot` walks this list and keeps the first plan whose coding sequence
    satisfies every constraint, so everything after that point is discarded unexamined.
    Exposed because those discarded plans are exactly what `search.py` needs in order to
    answer "is the chosen design good relative to what was available?" — and reusing the
    same enumeration guarantees the search explores the real alternatives rather than a
    separately generated set that might differ.

    Returns `(plans, ceiling)` where each plan is `(cuts, overhangs, fidelity)` and
    `ceiling` is the best fidelity the destination pair permits. `want` widens the
    enumeration past the three plans `design_oneshot` needs.
    """
    return _plan_fragments(protein, n_fragments, destination, matrix, enzyme_profile, want)


def _overhang_decisions(protein, cuts, chosen, profile) -> list[OverhangDecision]:
    """Every candidate at every junction, and why it was or was not used.

    Three outcomes, and the distinction matters: *rejected* means no synonymous
    realization exists under the active profile -- the junction is impossible, not merely
    worse. *Considered* means it was viable and another was preferred.
    """
    out: list[OverhangDecision] = []
    for cut, pick in zip(cuts, chosen):
        feasible = set(safe_overhangs(protein, cut, enzymes=profile))
        for oh in achievable_overhangs(protein, cut):
            if oh not in feasible:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="rejected",
                    reason=f"no locally feasible synonymous realization under {profile}",
                    local_realizations=0))
            elif oh == pick:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="selected",
                    reason="highest predicted fidelity among feasible candidates",
                    local_realizations=realization_count(protein, cut, oh, profile)))
            else:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="considered",
                    reason="feasible, but another candidate scored at least as well",
                    local_realizations=realization_count(protein, cut, oh, profile)))
    return out


FINAL_ARTIFACT_SCHEMA = 1


def final_artifact_report(oligos, enzyme: str = "BsaI", profile=None,
                          enzymes=None) -> dict:
    """Judge the sequences that will actually be ordered, not the coding span.

    Scope is `wrapped_oligos`: each `oligo_sequence_5to3` after wrapping, checked through
    `substrates.synthesis_findings` with this route's own intended-site records. Translation
    is not rechecked here -- it is a property of the coding span, which `qc` covers.

    Status is `pass`, `warning` or `fail`. A missing check is reported as `not_assessed`
    rather than counted as a pass.
    """
    from . import substrates as _sub
    from .assembly import wrapper_intended_sites
    from .synthesis_profile import resolve

    from .substrates import rc as _rc
    from .assembly import enzyme_geometry as _geometry

    policy = resolve(profile)
    rows = oligos.to_dict("records") if oligos is not None and len(oligos) else []
    fragments, hard, advisory, unassessed = [], [], [], []

    if not rows:
        # An empty table is not a clean table. This used to return "pass".
        return {
            "schema": FINAL_ARTIFACT_SCHEMA, "scope": "wrapped_oligos", "enzyme": enzyme,
            "policy": policy.name, "policy_rules": _policy_rules(policy, enzymes),
            "status": "not_assessed", "hard": [], "advisory": [],
            "not_assessed": ["no oligo rows were supplied, so nothing was checked"],
            "fragments": [],
            "means": "nothing was assessed; this is not a pass",
        }

    for row in rows:
        sequence = (row.get("oligo_sequence_5to3") or "").strip()
        payload = (row.get("payload_5to3") or "").strip()
        row_enzyme = row.get("wrap_enzyme") or enzyme
        ident = row.get("order_fragment_id")
        missing = []
        if not sequence:
            missing.append("no oligo_sequence_5to3")
        if not payload:
            missing.append("no payload_5to3, so the construction sites cannot be located")
        elif payload not in sequence:
            missing.append("the declared payload is not inside the ordered sequence")
        if missing:
            entry = {"order_fragment_id": ident, "length_nt": len(sequence),
                     "hard": [], "advisory": [], "intended_sites_located": False,
                     "not_assessed": missing}
            fragments.append(entry)
            unassessed += [f"{ident}: {m}" for m in missing]
            continue

        intended = wrapper_intended_sites(payload, row_enzyme)
        # Derived offsets are where the sites *should* be. Check they are there: a
        # fabricated wrapper with no recognition site at all reported them as located.
        site, _spacer, _oh = _geometry(row_enzyme)
        absent = [f"{nm} {strand} site missing at {at}"
                  for nm, strand, at in intended
                  if sequence[at:at + len(site)] != (site if strand == "forward"
                                                     else _rc(site))]
        findings = _sub.synthesis_findings(sequence, policy, intended_sites=intended,
                                           enzymes=enzymes)
        entry = {
            "order_fragment_id": ident,
            "length_nt": len(sequence),
            "sha256": _hashlib.sha256(sequence.encode()).hexdigest(),
            "hard": findings["hard"] + absent,
            "advisory": findings["target"],
            "intended_sites_located": not absent,
        }
        fragments.append(entry)
        hard += [f"{ident}: {m}" for m in entry["hard"]]
        advisory += [f"{ident}: {m}" for m in findings["target"]]

    status = ("fail" if hard else
              "not_assessed" if unassessed else
              "warning" if advisory else "pass")
    return {
        "schema": FINAL_ARTIFACT_SCHEMA,
        "scope": "wrapped_oligos",
        "enzyme": enzyme,
        "policy": policy.name,
        "policy_rules": _policy_rules(policy, enzymes),
        "status": status,
        "hard": hard,
        "advisory": advisory,
        "not_assessed": unassessed,
        "fragments": fragments,
        "means": ("the ordered sequences satisfy the declared software checks under this "
                  "policy; not vendor acceptance and not biological safety"),
    }


#: The only final-artifact verdicts that may be ordered. Named once and read everywhere,
#: because the project previously carried two predicates: `status != "fail"` in the library
#: and a stricter conjunction in `design_oneshot`. `not_assessed` sat between them, so one
#: path refused it and the other wrote it a normal order sheet.
ORDERABLE_STATUSES = ("pass", "warning")


def synthesis_readiness(result: dict) -> tuple[bool, list[str]]:
    """The one authoritative readiness verdict, with its reasons.

    Judged from the recorded evidence rather than from a `synthesis_ready` flag, which a
    caller may have copied, edited or carried over from an earlier verdict: a flag that
    disagrees with the record is itself a reason to refuse. Absent, unrecognised and
    self-contradictory information are all refusals, never defaults to acceptance.
    """
    why: list[str] = []
    report = result.get("final_artifacts")
    if not isinstance(report, dict) or "status" not in report:
        why.append("no final-artifact assessment was recorded")
    else:
        status = report.get("status")
        if status not in ORDERABLE_STATUSES:
            why.append(f"the final-artifact status is {status!r}, not one of "
                       f"{', '.join(ORDERABLE_STATUSES)}")
        if report.get("hard"):
            why.append(f"{len(report['hard'])} hard policy breach(es) on the ordered "
                       f"sequences: {report['hard'][0]}")

    if "constraints_ok" not in result:
        why.append("no optimiser constraint result was recorded")
    elif not result["constraints_ok"]:
        why.append("the codon optimiser could not satisfy every constraint")

    ready = not why
    flag = result.get("synthesis_ready")
    if flag is not None and bool(flag) != ready:
        why.append(f"the recorded synthesis_ready={bool(flag)} disagrees with the "
                   f"assessment on file")
        ready = False
    return ready, why


def _policy_rules(policy, enzymes=None) -> dict:
    """The numeric rules actually applied. A policy name does not identify a custom one."""
    from . import constants as _C

    return {"local_gc": list(policy.local_gc), "window": policy.window,
            "global_gc": list(policy.global_gc) if policy.global_gc else None,
            "max_homopolymer": policy.max_homopolymer,
            "forbidden_sites": sorted(enzymes) if enzymes
            else sorted(_C.enzymes_for(_C.DEFAULT_ENZYME_PROFILE))}


def design_oneshot(
    target_rna: str,
    organism: str = "c_reinhardtii_nuclear",
    genetic_code: int | None = None,
    codon_table=None,
    n_fragments: int | None = None,
    destination: tuple[str, str] | None = None,
    enzyme: str = "BsaI",
    enzyme_profile: str | tuple[str, ...] = C.DEFAULT_ENZYME_PROFILE,
    extra_blacklist: tuple[str, ...] | str = (),
    matrix: str = "BsaI-HFv2",
    seed: int = 42,
    check_offtarget: bool = True,
    avoid_sequences: tuple[str, ...] = (),
    lock_prefix: str = "",
    outdir: str | Path | None = None,
) -> dict:
    """Design a complete PPR binder for `target_rna`, ready to order.

    `architecture` is not a parameter: it follows from the target's length, since a PPR
    needs one repeat per base. Pass a 9, 14 or 19 base target to get 9S, 14S or 19S.

    Returns the protein, the coding sequence, the oligo table, the QC verdict, the
    predicted ligation fidelity and -- when `outdir` is given -- the paths written. No cost
    estimate: pricing was removed in 44446fe because we cannot obtain quotes.

    `avoid_sequences` forbids specific DNA verbatim, and `lock_prefix` fixes the first bases
    of the coding sequence to a given synonymous encoding. Both exist so
    `homology.diversify_library` can stop library members from sharing their scaffold
    sequence, and both default to nothing, so a design that does not ask for them is
    unaffected.
    """
    d = describe(target_rna)
    protein, architecture = d["aa_sequence"], d["architecture"]

    if genetic_code is None:
        genetic_code = ORGANISMS.get(organism, (None, 1))[1]
    table = _resolve_codon_table(codon_table, organism, genetic_code)

    if n_fragments is None:
        n_fragments = max(1, -(-len(protein) // MAX_FRAGMENT_AA))

    # A profile says which sites this project always excludes. `extra_blacklist` is for
    # the site a particular experiment needs clear -- an EcoRI in a downstream vector, say
    # -- which is a property of that experiment, not of PPROPT's policy, so it is added
    # here rather than becoming a new profile.
    if extra_blacklist:
        extra = (tuple(e.strip() for e in extra_blacklist.split(",") if e.strip())
                 if isinstance(extra_blacklist, str) else tuple(extra_blacklist))
        enzyme_profile = tuple(dict.fromkeys(C.enzymes_for(enzyme_profile) + extra))

    timings: dict[str, float] = {}
    with _stage(timings, "planning"):
        plans, _ceiling = _plan_fragments(protein, n_fragments, destination, matrix,
                                          enzyme_profile)

    # A locked overhang can create a blacklisted enzyme site that no synonymous change
    # can remove, because EnforceSequence has frozen those bases -- DNA Chisel then
    # reports NoSolutionError and returns a sequence that satisfies nothing. Measured on
    # the corpus, this hits 13 of 50 19S targets. The cure is a different split, so try
    # the next one rather than hand back a coding sequence that failed its constraints.
    attempt = fallback = None
    with _stage(timings, "optimization"):
        for cuts, overhangs, fidelity in plans[:OPTIMISE_ATTEMPTS]:
            locked = {3 * c - 4: o for c, o in zip(cuts, overhangs)}
            if lock_prefix:
                prefix = str(lock_prefix).upper().replace("U", "T")
                # Overlapping a junction overhang would give DNA Chisel two EnforceSequence
                # constraints on the same bases; if they disagreed the design would fail with a
                # message about neither of them. Refuse up front instead.
                clash = [s for s in locked if s < len(prefix)]
                if clash:
                    raise ValueError(
                        f"lock_prefix covers {len(prefix)} nt but a junction overhang is locked "
                        f"at {min(clash)}; they would conflict")
                locked[0] = prefix
            opt = optimize_cds(protein, locked_sites=locked, codon_table=table,
                               genetic_code=genetic_code, seed=seed,
                               enzymes=enzyme_profile, avoid_sequences=avoid_sequences)
            if opt["constraints_ok"]:
                attempt = (cuts, overhangs, fidelity, opt)
                break
            if fallback is None:
                fallback = (cuts, overhangs, fidelity, opt)
        if attempt is None:
            if fallback is None:
                raise ValueError(f"no viable design for {target_rna}")
            attempt = fallback
        cuts, overhangs, fidelity, opt = attempt
        cds = opt["cds"]

    decisions = _overhang_decisions(protein, cuts, overhangs, enzyme_profile)

    # split_cds re-derives the overhangs from the CDS and raises if the lock slipped
    with _stage(timings, "assembly"):
        fragments = split_cds(cds, cuts, overhangs, destination, enzyme)
        oligos = build_oligos(cds, cuts, overhangs, destination, enzyme, prefix=target_rna)
    with _stage(timings, "qc"):
        qc = synthesis_qc(cds)

    dest = destination or C.DESTINATION_OVERHANGS["level0"]
    reaction = reaction_overhangs(overhangs, "level0") if destination is None else [
        dest[0], *overhangs, _rc(dest[1])]

    paths: dict[str, str] = {}
    # Judged on the sequences that will be ordered, before any of them is written. An
    # insert can pass its own GC band while the appended junction and wrapper push the
    # ordered oligo outside it.
    final_artifacts = final_artifact_report(oligos, enzyme=enzyme,
                                            enzymes=C.enzymes_for(enzyme_profile))
    # The export route and the reported verdict come from one predicate. Routing on
    # `status == "fail"` alone sent `not_assessed` down the normal path while the
    # conjunction below reported it unready, so the files and the flag disagreed.
    ready, not_ready_because = synthesis_readiness(
        {"final_artifacts": final_artifacts, "constraints_ok": opt["constraints_ok"]})

    if outdir is not None and not ready:
        # A hard breach must not leave a full set of normal order files behind. The
        # diagnostic route keeps the sequences for inspection and names them as failed.
        # Everything is still written, because a rejected design is worth inspecting -- but
        # every filename and key says REJECTED, and no normal order-file key exists, so a
        # downstream step cannot pick these up as orderable.
        with _stage(timings, "export"):
            out = Path(outdir)
            out.mkdir(parents=True, exist_ok=True)
            stem = f"ppropt_{target_rna}_REJECTED"
            paths = {
                "diagnostic_oligo_csv": str(write_oligo_csv(oligos,
                                                            out / f"{stem}_oligos.csv")),
                "diagnostic_oligo_fasta": str(write_fasta(oligos,
                                                          out / f"{stem}_oligos.fasta")),
                "diagnostic_gene_fasta": str(write_gene_fasta(
                    cds, target_rna, out / f"{stem}_gene.fasta", architecture)),
                "diagnostic_genbank": str(write_genbank(
                    cds, fragments, target_rna, out / f"{stem}.gb", genetic_code, enzyme,
                    destination)),
            }
    elif outdir is not None:
        with _stage(timings, "export"):
            out = Path(outdir)
            stem = f"ppropt_{target_rna}"
            paths = {
                "oligo_csv": str(write_oligo_csv(oligos, out / f"{stem}_oligos.csv")),
                "oligo_fasta": str(write_fasta(oligos, out / f"{stem}_oligos.fasta")),
                "gene_fasta": str(write_gene_fasta(cds, target_rna,
                                                   out / f"{stem}_gene.fasta", architecture)),
                "genbank": str(write_genbank(cds, fragments, target_rna, out / f"{stem}.gb",
                                             genetic_code, enzyme, destination)),
            }

    # Does the target also occur where a PPR could bind it? A PPR binds RNA, so the tier
    # that matters is an annotated transcript in the sense orientation, not genomic DNA on
    # either strand. Measured over 200 random targets: 16% of 9-nt targets occur in a
    # transcript (48% occur somewhere in the DNA, which overstates it threefold); 14-nt
    # and 19-nt targets occur in neither.
    offtarget = None
    if check_offtarget:
        try:
            from .offtarget import scan as _scan
            with _stage(timings, "screening"):
                offtarget = _scan(d["target_rna"])
        except (OSError, ValueError, ImportError) as exc:
            # Only genuine unavailability -- no network, no cached genome -- is tolerated.
            # A broad `except` here previously hid a signature mismatch behind a cheerful
            # "not checked", which is exactly how a silent no-op ships.
            offtarget = {"verdict": "not checked",
                         "error": f"{type(exc).__name__}: {exc}",
                         "n_transcript": None, "n_genomic": None, "genes": []}

    # Screening that did not run must be visible in the audit. The matches themselves are
    # reported through `warnings` below, so this line only records whether it happened.
    if not check_offtarget:
        offtarget_status = "not assessed — screening disabled"
    elif offtarget.get("verdict") == "not checked":
        offtarget_status = f"not assessed — {offtarget['error']}"
    else:
        offtarget_status = "screened against host transcripts and genomic DNA"

    warnings = list(qc["warnings"])
    if not opt["constraints_ok"]:
        warnings.append("codon optimiser could not satisfy every constraint")
    if offtarget and offtarget.get("n_transcript"):
        genes = ", ".join(offtarget["genes"][:3])
        warnings.append(
            f"this target occurs inside {offtarget['n_transcript']} annotated host "
            f"transcript(s) ({genes}), in the sense orientation — RNA a PPR could bind "
            f"as well as your construct. A longer target often reduces exact matches, "
            f"but screen the actual candidate against the chosen reference rather than "
            f"assuming absence; longer targets are not guaranteed to be absent. This "
            f"counts exact sequence matches only: it predicts no binding affinity, and "
            f"mismatch-tolerant interactions are not assessed, so it does not bound the "
            f"set of possible off-targets in either direction."
        )
    elif offtarget and offtarget.get("n_genomic"):
        warnings.append(
            f"this target occurs {offtarget['n_genomic']}x in host genomic DNA and in no "
            f"annotated transcript of the supplied reference. That is an absence of "
            f"evidence in one annotation, not an absence of RNA: unannotated transcripts "
            f"and other interaction mechanisms are not assessed."
        )
    backbone = set_fidelity([dest[0], _rc(dest[1])], matrix)
    if backbone < 0.95:
        warnings.append(
            f"the destination pair {dest[0]}/{_rc(dest[1])} scores {backbone:.3f} on its "
            f"own, and the whole reaction scores {fidelity:.3f}. The backbone is the "
            f"limit here, not the junctions we chose."
        )

    # The final verdict is part of the one-line summary. "QC PASS" alone was true of a
    # design the final-artifact check rejects, and it was the whole printed outcome.
    summary = (
        f"{target_rna} ({architecture}) -> {len(protein)} aa, {len(cds)} nt, "
        f"{len(oligos)} fragments. Fidelity {fidelity:.3f}. "
        f"Coding-span QC {qc['status']}. "
        + ("Orderable under the declared software checks."
           if ready else f"NOT ORDERABLE (final artifacts: {final_artifacts['status']})."))

    audit = DesignAudit(
        target_rna=d["target_rna"],
        architecture=architecture,
        enzyme_profile=enzyme_profile if isinstance(enzyme_profile, str) else "custom",
        enzymes_excluded=tuple(C.enzymes_for(enzyme_profile)),
        genetic_code=genetic_code,
        organism=organism,
        cuts=tuple(cuts),
        overhang_decisions=tuple(decisions),
        fidelity=fidelity,
        fidelity_components=fidelity_components(reaction, matrix),
        fidelity_ceiling=backbone,
        constraints_satisfied=bool(opt["constraints_ok"]),
        qc_status=qc["status"],
        final_status=final_artifacts["status"],
        synthesis_ready=ready,
        not_ready_because=tuple(not_ready_because),
        offtarget_status=offtarget_status,
        findings=tuple(warnings),
    )

    return {
        "summary": summary,
        "audit": audit,
        "target_rna": d["target_rna"],
        "architecture": architecture,
        "ppr_code": d["ppr_code"],
        "protein": protein,
        "cds": cds,
        "cuts": cuts,
        "junction_overhangs": overhangs,
        "reaction_overhangs": reaction,
        "fragments": fragments,
        "oligos": oligos,
        "qc": qc,
        # The coding span's complexity check and the final ordered sequence's policy check
        # are separate results with different scopes; neither substitutes for the other.
        "final_artifacts": final_artifacts,
        # The same verdict that chose the export route above, so a file set and a flag
        # cannot disagree. Named so it cannot be read as more than it is.
        "synthesis_ready": ready,
        "not_ready_because": not_ready_because,
        "synthesis_ready_means": ("the declared software checks passed on the ordered "
                                  "sequences; not vendor acceptance, not biological safety"),
        "offtarget": offtarget,
        "fidelity": fidelity,
        "constraints_ok": opt["constraints_ok"],
        "warnings": warnings,
        "paths": paths,
        #: Wall seconds per pipeline stage. Nondeterministic by nature, so nothing that
        #: fingerprints a design may include it -- see `validation/check_regression.py`.
        "timings": dict(timings),
        # The resolved inputs, not the requested ones. `organism` names a table but
        # `codon_table` may override it, and `extra_blacklist` is merged into the profile,
        # so a caller cannot reconstruct either from its own arguments. Reported so a design
        # can be reproduced exactly, and so `search.py` can re-optimise under identical
        # settings rather than re-deriving them.
        # Whether the deposited GRASP kit could build this target without synthesising
        # anything. Always reported, because "you already own these parts" is an answer the
        # de novo route cannot give and a lab holding the kit would want first.
        "parts_plan": _parts_plan(d["target_rna"]),
        "codon_table": table,
        "genetic_code": genetic_code,
        "enzyme_profile_effective": tuple(C.enzymes_for(enzyme_profile)),
        "n_fragments": n_fragments,
    }
