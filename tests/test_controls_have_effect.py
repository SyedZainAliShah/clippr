"""Every parameter offered must actually change the outcome.

A control that changes nothing is worse than no control: it implies a decision was made
when none was. This project has already shipped one -- `enzyme_profile` sat in the notebook
form for a while without being passed through, so switching it did nothing at all.

Each test varies one parameter and asserts the result differs in the way that parameter
should cause. Where a parameter can legitimately be a no-op, the test is built so it cannot
be: the EagI case below picks a site that *is* present, because blacklisting a site the
sequence never had proves nothing.

**That bar is the point of the file, and it is the rule for anything added here.** Every
precondition is measured rather than assumed. Two concrete cases from building this sweep:
the BASE fixture at seed 42 happens to contain no SapI site even unconstrained, so a test
asserting one would have passed vacuously; and 2 of the first 3 scaffold encodings are
infeasible for this target, so `lock_prefix` has to pick a feasible one rather than index 0.

**Coverage.** `design_oneshot` offers 16 parameters and the notebook form offers 23
controls. Every one of them has a test here, except `target_rna`, which is the subject of
every test in the suite and needs no separate control. A control that only *appears* in the
notebook is not covered -- `TestNotebookControlsAreWired` requires it to be acted on, and
`TestTheWiringGuardItself` proves that guard rejects a control that is merely displayed.

Nothing here keys on the package's brand name or on a notebook's filename: the notebook is
found by content and output files are checked by their path *keys*, so a rename cannot make
this sweep fail for a reason that has nothing to do with the controls.
"""
from __future__ import annotations

import ast
import contextlib
import io
import json
import re
from pathlib import Path

import pytest
from Bio import Restriction

from ppropt import design_oneshot
from ppropt.biology import N_TERMINAL
from ppropt.codons import table_from_kazusa
from ppropt.design import _resolve_codon_table
from ppropt.homology import scaffold_encoding
from ppropt.policy import ORGANISMS

from test_crosstalk import TABLE as PPR_SCORE_TABLE
from test_design import TOY_TABLE
from test_notebook import as_python

BASE = dict(target_rna="AAAAUGUGG", seed=42, check_offtarget=False,
            codon_table=TOY_TABLE)


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.fixture(scope="module")
def ref():
    return design_oneshot(**BASE)


@pytest.fixture(scope="module")
def sapi_pair():
    """(assembly, igem_rfc1000) at a seed whose unconstrained sequence really has SapI."""
    kw = {**BASE, "seed": 0}
    return (design_oneshot(**{**kw, "enzyme_profile": "assembly"}),
            design_oneshot(**kw))


class TestControlsHaveEffect:
    def test_enzyme_changes_the_ordered_oligo(self, ref):
        alt = design_oneshot(**{**BASE, "enzyme": "BbsI"})
        a = ref["oligos"].iloc[0]["oligo_sequence_5to3"]
        b = alt["oligos"].iloc[0]["oligo_sequence_5to3"]
        assert a != b
        assert str(Restriction.BsaI.site) in a
        assert str(Restriction.BbsI.site) in b

    def test_matrix_changes_predicted_fidelity(self, ref):
        alt = design_oneshot(**{**BASE, "matrix": "BbsI-HF"})
        assert abs(alt["fidelity"] - ref["fidelity"]) > 1e-9

    def test_enzyme_profile_changes_what_is_excluded(self, ref):
        alt = design_oneshot(**{**BASE, "enzyme_profile": "assembly"})
        assert alt["audit"].enzymes_excluded != ref["audit"].enzymes_excluded

    def test_enzyme_profile_changes_the_dna(self, sapi_pair):
        """The site is genuinely there without the constraint, and gone with it.

        The test above asserts only that the *audit record* differs, which a parameter
        that reached nothing but the record would also satisfy. This one asserts the DNA.
        SapI avoidance is load-bearing rather than an edge case: GCTCTTC falls out of the
        repeat's own E-L-F-D motif (GAG CTC TTC GAC), which recurs about ten times per
        design, and under `assembly` 24 of 30 seeds at 9S carry at least one site with the
        real Chlamydomonas table. Seed 0 is used because the file's seed 42 happens to
        produce a clean sequence even unconstrained, where the assertion would prove
        nothing -- 10 of the first 16 toy-table seeds carry a site, 42 is not one.
        """
        assembly, default = sapi_pair
        site = str(Restriction.SapI.site)
        assert site == "GCTCTTC"

        before = assembly["cds"].count(site) + assembly["cds"].count(rc(site))
        assert before > 0, "pick a seed whose unconstrained sequence really has the site"
        assert default["cds"].count(site) + default["cds"].count(rc(site)) == 0

        assert assembly["protein"] == default["protein"]
        assert assembly["constraints_ok"] and default["constraints_ok"]

    def test_enzyme_profile_reaches_the_qc_verdict(self, sapi_pair):
        """QC must judge the profile the caller declared, not the module default.

        It did not: `qc.BLACKLIST_ENZYMES` was applied unconditionally, so `assembly`
        failed a design for the SapI site it had deliberately allowed while the
        ordered-oligo gate -- which does receive the caller's profile -- passed the same
        design. The reverse case is worse: `moclo_compat` exists only to exclude BsmBI and
        the verdict never checked for it.
        """
        assembly, default = sapi_pair
        assert "SapI" in default["qc"]["blacklist_enzymes"]
        assert "SapI" not in assembly["qc"]["blacklist_enzymes"]
        assert assembly["qc"]["blacklist_enzymes"] != default["qc"]["blacklist_enzymes"]

        # The site is present, so the old module-default gate reported it here.
        assert not any("SapI" in h for h in assembly["qc"]["blacklist_hits"])
        assert not any("SapI" in f for f in assembly["qc"]["failures"])

    def test_extra_blacklist_removes_a_site_that_was_there(self, ref):
        """Uses a site the default sequence actually contains.

        Blacklisting an absent site is a legitimate no-op and proves nothing, which is
        how an earlier version of this test wrongly reported the parameter as ignored.
        """
        site = str(Restriction.EagI.site)
        before = ref["cds"].count(site) + ref["cds"].count(rc(site))
        assert before > 0, "pick a site the unconstrained sequence really contains"

        alt = design_oneshot(**{**BASE, "extra_blacklist": "EagI"})
        after = alt["cds"].count(site) + alt["cds"].count(rc(site))
        assert after == 0
        assert "EagI" in alt["audit"].enzymes_excluded

    def test_destination_changes_the_reaction(self, ref):
        alt = design_oneshot(**{**BASE, "destination": ("GGAG", "CGCT")})
        assert alt["reaction_overhangs"][0] == "GGAG"
        assert alt["reaction_overhangs"] != ref["reaction_overhangs"]

    def test_n_fragments_is_obeyed(self, ref):
        alt = design_oneshot(**{**BASE, "n_fragments": 5})
        assert len(alt["oligos"]) == 5 != len(ref["oligos"])

    def test_seed_reproduces_exactly(self, ref):
        assert design_oneshot(**BASE)["cds"] == ref["cds"]

    def test_check_offtarget_toggles_the_report(self, ref):
        assert ref["offtarget"] is None
        on = design_oneshot(**{**BASE, "check_offtarget": True})
        assert on["offtarget"] is not None

    def test_codon_table_changes_the_sequence(self, ref):
        sparse = {aa: {next(iter(c)): 1.0} for aa, c in TOY_TABLE.items()}
        alt = design_oneshot(**{**BASE, "codon_table": sparse})
        assert alt["cds"] != ref["cds"]

    def test_genetic_code_reaches_the_audit(self, ref):
        assert ref["audit"].genetic_code == 1
        alt = design_oneshot(**{**BASE, "genetic_code": 11})
        assert alt["audit"].genetic_code == 11

    def test_avoid_sequences_removes_dna_that_was_there(self, ref):
        """The motif is a slice of the unconstrained design, so it is genuinely present.

        `avoid_sequences` had no test in any file before this one -- the single real hole
        the sweep found. An arbitrary motif would prove nothing for the same reason the
        EagI case above picks a site that exists: forbidding DNA the design never carried
        is a legitimate no-op. Measured here, the slice occurs twice before and never
        after, and only the codons move.
        """
        motif = ref["cds"][90:102]
        assert len(motif) == 12
        before = ref["cds"].count(motif)
        assert before > 0, "the motif must come from the unconstrained design"

        alt = design_oneshot(**{**BASE, "avoid_sequences": (motif,)})
        assert alt["cds"].count(motif) == 0
        assert alt["cds"] != ref["cds"]
        assert alt["protein"] == ref["protein"], "a synonymous change only"
        assert alt["constraints_ok"]

    def test_lock_prefix_fixes_the_start_of_the_coding_sequence(self, ref):
        """A locked scaffold encoding appears verbatim, and the free design does not have it.

        The encoding is searched for rather than taken at index 0 because **most of them
        are infeasible for this target**: measured over the first three, only index 1
        satisfies its constraints, the other two fall back. That is the documented effect
        of a locked prefix freezing bases a blacklisted site needs -- see `design.py` on
        `NoSolutionError` -- so asserting `constraints_ok` on a fixed index would be a
        flaky assertion about DNA Chisel rather than about this parameter.
        """
        encodings = [scaffold_encoding(i, N_TERMINAL, ref["codon_table"],
                                       genetic_code=ref["genetic_code"],
                                       enzymes=ref["enzyme_profile_effective"])
                     for i in range(6)]
        assert len({*encodings}) > 1, "the encodings must differ to begin with"

        for prefix in encodings:
            assert not ref["cds"].startswith(prefix), (
                "the free design must not already start with the prefix, or locking it "
                "proves nothing")
            alt = design_oneshot(**{**BASE, "lock_prefix": prefix})
            assert alt["cds"].startswith(prefix), "the lock must hold regardless"
            assert alt["protein"] == ref["protein"]
            if alt["constraints_ok"]:
                assert alt["cds"] != ref["cds"]
                return
        pytest.fail("no scaffold encoding gave a feasible locked design")

    def test_wrapper_aware_changes_the_sequence_that_would_be_ordered(self, ref):
        """On by default since 2026-10-12, and this target is the case that moved.

        Not a no-op either way: `applied` is asserted true for the default, so the second
        optimisation pass really ran, and switching it off really does hand back a design
        the final-artifact check refuses. That refusal is the documented reason the switch
        exists, and `test_display_verdict.py` builds its rejected fixture the same way.
        """
        assert ref["wrapper_aware"] == {"requested": True, "applied": True}
        assert ref["final_artifacts"]["status"] == "pass"
        assert ref["synthesis_ready"] is True

        off = design_oneshot(**{**BASE, "wrapper_aware": False})
        assert off["wrapper_aware"]["requested"] is False
        assert off["wrapper_aware"]["applied"] is False
        assert off["cds"] != ref["cds"]
        assert off["protein"] == ref["protein"]
        assert off["final_artifacts"]["status"] == "fail"
        assert off["synthesis_ready"] is False

    def test_outdir_writes_the_order_files(self, ref, tmp_path):
        """Keyed on the `paths` keys and the file count, never on the filenames.

        The generated names carry the package name, which is being renamed, so asserting
        on them would make this test fail for a reason that is not about this parameter.
        """
        assert ref["paths"] == {}, "no directory means nothing is written"
        assert not any(tmp_path.iterdir())

        alt = design_oneshot(**{**BASE, "outdir": tmp_path})
        assert set(alt["paths"]) == {"genbank", "gene_fasta", "oligo_csv", "oligo_fasta"}
        written = sorted(tmp_path.iterdir())
        assert len(written) == 4
        assert all(p.stat().st_size > 0 for p in written), "an empty file is not an export"

    def test_organism_changes_the_coding_sequence(self, host_codon_table):
        """BASE's explicit table is dropped here, because with it `organism` is inert.

        That is correct behaviour, not a defect: `_resolve_codon_table` takes an explicit
        mapping first, so a caller who supplies a table has already answered the question
        the dropdown asks. The parameter is therefore exercised the only way it can have
        an effect -- with no table passed -- and both hosts below go through Kazusa, so
        neither is the slow genome-derived route.
        """
        a = design_oneshot("AAAAUGUGG", organism="c_reinhardtii_nuclear", seed=42,
                           check_offtarget=False)
        b = design_oneshot("AAAAUGUGG", organism="e_coli", seed=42, check_offtarget=False)
        assert a["cds"] != b["cds"]
        assert a["protein"] == b["protein"], "the host changes the spelling, not the protein"
        assert a["codon_table"] != b["codon_table"]

    def test_organism_picks_the_genetic_code_and_the_table(self, host_codon_table):
        """The two-context trap in CLAUDE.md, as a control: the tables are near mirrors.

        Checked at `_resolve_codon_table` rather than through a design because the
        chloroplast table is derived by counting a genome's own CDS features, and a full
        chloroplast design measured at 233 s -- too slow to be worth it when this is the
        one function the parameter reaches.
        """
        assert ORGANISMS["c_reinhardtii_nuclear"][1] == 1
        assert ORGANISMS["c_reinhardtii_chloroplast"][1] == 11

        nuclear = _resolve_codon_table(None, "c_reinhardtii_nuclear", 1)
        plastid = _resolve_codon_table(None, "c_reinhardtii_chloroplast", 11)
        assert max(nuclear["L"], key=nuclear["L"].get) == "CTG"
        assert max(plastid["L"], key=plastid["L"].get) == "TTA"
        assert nuclear != plastid


# --------------------------------------------------------------------------------------
# The notebook form, exercised as a form
#
# The notebook is the only interface most users touch, and eleven of its controls reach a
# function `design_oneshot` never sees. Those are tested by running the generated cell's
# own source with the form fields set -- which is what typing into a Colab form does --
# rather than by calling the library function the cell happens to call, because the wiring
# between the two is exactly what has broken twice before.
#
# `test_display_verdict.py` established this pattern for the result card; these helpers
# generalise it.
# --------------------------------------------------------------------------------------

NOTEBOOKS = Path(__file__).resolve().parents[1] / "notebooks"


def designer_notebook() -> Path:
    """The designer notebook, found by content rather than by name.

    The filename carries the package's brand, which has been renamed once already. Keying
    on it would make every test below fail for a reason that has nothing to do with form
    controls.

    **Deliberately not a skip**, on either branch. A skip disarms these tests exactly when
    it matters: move or rename the notebook and every control test below goes green while
    nothing at all is checked. The notebook is in the repository, so its absence is a
    failure and not a missing optional input. (External inputs are different: the cached
    annotation and the Kazusa tables below do skip, because those genuinely may not be
    present -- which is why 194 tests in this suite skip today.)
    """
    hits = [p for p in sorted(NOTEBOOKS.glob("*.ipynb"))
            if "design_oneshot(" in p.read_text(encoding="utf-8")]
    assert hits, (
        f"no notebook in {NOTEBOOKS} calls design_oneshot, so nothing below was checked")
    assert len(hits) == 1, (
        f"expected one designer notebook, found {[p.name for p in hits]}; "
        f"this helper picks the notebook to test and cannot choose between two")
    return hits[0]


@pytest.fixture(scope="module")
def cells() -> list[str]:
    """Every code cell, with IPython magics blanked so `ast` and `exec` both accept them."""
    nb = json.loads(designer_notebook().read_text(encoding="utf-8"))
    return [as_python("".join(c["source"])) for c in nb["cells"]
            if c["cell_type"] == "code"]


@pytest.fixture(scope="module")
def notebook_source(cells) -> str:
    return "\n".join(cells)


@pytest.fixture(scope="module")
def declared(notebook_source) -> set[str]:
    names = {line.split("=")[0].strip()
             for line in notebook_source.splitlines() if "#@param" in line}
    assert names, "no form controls found"
    return names


@pytest.fixture(scope="module")
def host_annotation():
    """The cached chloroplast annotation the window picker screens against."""
    from ppropt.offtarget import load_transcripts
    try:
        return load_transcripts()
    except Exception as exc:                        # pragma: no cover - needs the cache
        pytest.skip(f"chloroplast annotation unavailable: {exc}")


@pytest.fixture(scope="module")
def host_codon_table():
    """The default host's Kazusa table, which the organism dropdown resolves to."""
    try:
        return table_from_kazusa(ORGANISMS["c_reinhardtii_nuclear"][0])
    except Exception as exc:                        # pragma: no cover - needs the cache
        pytest.skip(f"host codon table unavailable: {exc}")


def cell_declaring(cells: list[str], name: str) -> str:
    """The one cell whose form declares `name`."""
    hits = [c for c in cells if re.search(rf"^{re.escape(name)}\s*=.*#@param", c, re.M)]
    assert len(hits) == 1, f"{name}: expected one declaring cell, found {len(hits)}"
    return hits[0]


def with_form_values(src: str, **values) -> str:
    """`src` with form fields set -- the one thing a Colab user can actually do."""
    for name, value in values.items():
        pattern = rf"^{re.escape(name)}\s*=.*#@param.*$"
        assert re.search(pattern, src, re.M), f"{name} is not a form control in this cell"
        # A function replacement, not a string: a Windows path's backslashes would
        # otherwise be read as group references and corrupt the value silently.
        src = re.sub(pattern, lambda _m, n=name, v=value: f"{n} = {v!r}",
                     src, count=1, flags=re.M)
    return src


def run_cell(src: str, env: dict | None = None) -> tuple[dict, str]:
    """Execute a cell, returning its namespace and everything it printed."""
    env = dict(env or {})
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        exec(compile(src, "<notebook cell>", "exec"), env)
    return env, out.getvalue()


#: Calls whose arguments are shown to the reader rather than acted on.
DISPLAY_CALLS = {"print", "display"}


def acted_on(source: str, name: str) -> bool:
    """Is `name` used for anything other than being shown to the reader?

    The predicate this replaces asked whether the name appeared on some line that was not
    its own `#@param` declaration. A control that is printed in a summary and never used
    satisfies that -- which is within one step of the defect this whole file exists for,
    so the guard could not have caught a second instance of it.

    A control counts as acted on when it, or a variable carrying its value, either
    **reaches a call that is not `print`/`display`**, or **decides control flow**. Value
    flow through plain assignments is followed, so a control formatted into a message that
    is then printed is still display-only, however many hops it takes to get there.

    Known limit, stated rather than hidden: a control whose only non-display call is into
    a function that itself merely formats text would pass. Nothing in the notebook does
    that today -- the result card takes already-computed values, not raw controls.
    """
    tree = ast.parse(source)
    parents = {child: node for node in ast.walk(tree)
               for child in ast.iter_child_nodes(node)}

    def loads(node) -> set[str]:
        return {n.id for n in ast.walk(node)
                if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}

    # Everything that could be carrying this control's value.
    tainted = {name}
    while True:
        grown = False
        for node in ast.walk(tree):
            if not isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                continue
            if node.value is None or not loads(node.value) & tainted:
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                for n in ast.walk(target):
                    if isinstance(n, ast.Name) and n.id not in tainted:
                        tainted.add(n.id)
                        grown = True
        if not grown:
            break

    def reaches_work(node) -> bool:
        current = node
        while current in parents:
            parent = parents[current]
            if isinstance(parent, ast.Call):
                func = parent.func
                called = (func.id if isinstance(func, ast.Name)
                          else getattr(func, "attr", ""))
                return called not in DISPLAY_CALLS
            if (isinstance(parent, (ast.If, ast.While, ast.IfExp))
                    and current is parent.test):
                return True
            if isinstance(parent, ast.comprehension) and current in parent.ifs:
                return True
            current = parent
        return False

    return any(isinstance(node, ast.Name) and node.id in tainted
               and isinstance(node.ctx, ast.Load) and reaches_work(node)
               for node in ast.walk(tree))


class TestTheWiringGuardItself:
    """A guard that cannot fail is not a guard, so these are its controls.

    Written before the guard was pointed at the real notebook. The last test pins the gap
    that motivated the rewrite, so nobody weakens the predicate back to a substring match
    without a test going red.
    """

    DISPLAY_ONLY = 'foo = "x"  #@param {type:"string"}\nprint(f"foo is {foo}")\n'

    def test_a_control_only_printed_is_not_wired(self):
        assert not acted_on(self.DISPLAY_ONLY, "foo")

    def test_a_control_formatted_into_a_printed_message_is_not_wired(self):
        """The subtler form, and the one a summary string actually looks like."""
        src = ('foo = "x"  #@param {type:"string"}\n'
               'msg = f"foo is {foo}"\n'
               'print(msg)\n')
        assert not acted_on(src, "foo")

    def test_a_control_assigned_but_never_read_is_not_wired(self):
        src = 'foo = "x"  #@param {type:"string"}\nbar = foo\n'
        assert not acted_on(src, "foo")

    def test_a_control_passed_to_a_function_is_wired(self):
        src = 'foo = "x"  #@param {type:"string"}\ndesign(target=foo)\n'
        assert acted_on(src, "foo")

    def test_a_control_that_decides_control_flow_is_wired(self):
        """`apply_to_target` is exactly this shape: it reaches no call, only an `if`."""
        src = 'foo = True  #@param {type:"boolean"}\nif foo:\n    go()\n'
        assert acted_on(src, "foo")

    def test_a_control_reaching_a_call_through_a_variable_is_wired(self):
        """`codon_table_file` is this shape, so the taint closure is load-bearing."""
        src = ('foo = ""  #@param {type:"string"}\n'
               'chosen = foo or None\n'
               'design(codon_table=chosen)\n')
        assert acted_on(src, "foo")

    def test_the_predicate_this_replaced_accepted_the_display_only_case(self):
        """The old guard asked only whether the name appeared on another line."""
        old = [line for line in self.DISPLAY_ONLY.splitlines()
               if "foo" in line and "#@param" not in line]
        assert old, "the old predicate accepted a display-only control"
        assert not acted_on(self.DISPLAY_ONLY, "foo"), "the new one must refuse it"


class TestNotebookControlsAreWired:
    """Every #@param in the notebook must be acted on, not merely displayed."""

    def test_every_form_control_is_acted_on(self, notebook_source, declared):
        inert = sorted(n for n in declared if not acted_on(notebook_source, n))
        assert not inert, (
            f"declared in a form but only displayed: {inert}. A control that changes "
            f"nothing implies a decision was made when none was.")

    def test_every_control_this_file_tests_is_still_a_control(self, declared):
        """Ties the two halves of the file together.

        The effect tests below rewrite named `#@param` lines. If a control were renamed
        or dropped, `with_form_values` would raise inside each affected test -- but a
        control removed from the form while an effect test still *passed* against a stale
        copy of the cell is the quieter failure, so the names are checked against the form
        once, here. Names rather than a count: a diff should read as a rename.
        """
        tested = {
            # Passed straight to a `design_oneshot` keyword, so the effect is the
            # corresponding control in `TestControlsHaveEffect` above.
            "target_rna", "organism", "genetic_code_override", "assembly_enzyme",
            "ligation_table", "enzyme_profile", "extra_blacklist", "destination_level",
            "custom_destination", "n_fragments", "seed", "check_offtarget",
            # Reached only from the notebook, so tested by running their own cells.
            "codon_table_file", "kazusa_taxid", "search_sequence", "sequence_origin",
            "window_lengths", "intended_gene", "screen_against_host", "apply_to_target",
            "targets", "ppr_score_table", "write_files",
        }
        assert not tested - declared, (
            f"this file tests controls the form no longer offers: "
            f"{sorted(tested - declared)}")
        assert not declared - tested, (
            f"the form offers controls with no effect test: {sorted(declared - tested)}")


#: The window picker's form as the notebook ships it.
PICKER_FORM = dict(search_sequence="", sequence_origin="native", window_lengths="9",
                   intended_gene="rbcL", screen_against_host=True,
                   apply_to_target=False)


@pytest.fixture(scope="module")
def picker(cells, host_annotation):
    return cell_declaring(cells, "apply_to_target")


def pick(picker: str, **overrides) -> tuple[dict, str]:
    """Run the window-picker cell with its form set, from the shipped defaults."""
    src = with_form_values(picker, **{**PICKER_FORM, **overrides})
    return run_cell(src, {"target_rna": "AAAAUGUGG"})


class TestTheWindowPickerControls:
    """The six controls on the window picker, none of which `design_oneshot` ever sees.

    Each is run through the generated cell itself. Two of them -- `sequence_origin` and
    `apply_to_target` -- have no effect anywhere in `src/` at all: their whole behaviour is
    a predicate in this cell's own source. Calling `select_window` directly would therefore
    report them as working no matter what the notebook did with them.
    """

    def test_window_lengths_chooses_that_architecture(self, picker):
        """Each option returns a window of its own length, so the field is not decorative."""
        chosen = {}
        for field in ("9", "14", "19"):
            env, _ = pick(picker, window_lengths=field)
            selection = env["selection"]
            assert selection.found, f"{field}: no window, so the length proves nothing"
            chosen[field] = selection.chosen.sequence
            assert len(chosen[field]) == int(field)
            assert selection.lengths == (int(field),)
        assert len(set(chosen.values())) == 3, "three lengths must give three windows"

    def test_screen_against_host_changes_what_is_rejected(self, picker):
        """Off searches the supplied region only, and says so instead of naming a reference.

        Non-vacuous in both directions: screening on genuinely rejects windows that
        screening off accepts, and the window finally chosen differs.
        """
        on, _ = pick(picker, screen_against_host=True)
        off, _ = pick(picker, screen_against_host=False)
        a, b = on["selection"], off["selection"]

        assert a.found and b.found
        assert len(a.rejected) > len(b.rejected) > 0, (
            "the host screen must reject windows the region-only screen accepts, and the "
            "region-only screen must reject some too, or neither half is measuring anything")
        assert a.chosen.sequence != b.chosen.sequence

        assert a.reference and a.reference != "none -- supplied sequence only"
        assert a.reference_binding is not None
        assert b.reference == "none -- supplied sequence only"
        assert b.reference_binding is None, "an unscreened reference must never be named"

    def test_sequence_origin_decides_whether_intent_is_declared(self, picker):
        """`origin` reaches no library decision; this cell's `native` predicate is its effect.

        `SearchSequence.origin` is validated on construction and then read nowhere in
        `src/` -- `select_window` never consults it. What the control actually does is
        decide, here, whether the named gene's own occurrences are declared as the intended
        site. The consequence is stark and measured: declaring intent leaves a usable
        window in rbcL, and withholding it rejects **every** candidate, because a native
        gene's own sequence is in the reference and counts against every window of itself.

        So this test has to run the cell. Asserting against `select_window` would have
        reported the control as working even if the notebook ignored it entirely.
        """
        native, _ = pick(picker, sequence_origin="native")
        designed, _ = pick(picker, sequence_origin="designed")
        unknown, _ = pick(picker, sequence_origin="unknown")

        assert native["selection"].intended_declared is True
        assert native["selection"].found, "intent declared must leave a window"

        for name, env in (("designed", designed), ("unknown", unknown)):
            selection = env["selection"]
            assert selection.intended_declared is False, name
            assert not selection.found, (
                f"{name}: withholding intent must cost the window, or the control is inert")
            assert len(selection.rejected) > len(native["selection"].rejected)

    def test_intended_gene_selects_which_gene_is_exempt(self, picker):
        """Two genes, two different exemptions, two different chosen windows.

        Blank is its own branch: nothing to search and nothing declared, which the cell
        says rather than silently rejecting every window in the user's own gene.
        """
        rbcl, _ = pick(picker, intended_gene="rbcL")
        psba, _ = pick(picker, intended_gene="psbA")

        assert rbcl["selection"].found and psba["selection"].found
        assert rbcl["selection"].intended_declared
        assert psba["selection"].intended_declared
        assert rbcl["selection"].chosen.sequence != psba["selection"].chosen.sequence
        assert len(rbcl["selection"].rejected) != len(psba["selection"].rejected)

        blank, printed = pick(picker, intended_gene="")
        assert "selection" not in blank, "with no gene and no region there is nothing to pick"
        assert "Nothing to search" in printed

    def test_search_sequence_is_what_gets_searched(self, picker):
        """A pasted region is searched; blank falls back to the named gene's annotated span.

        Both branches are asserted, because the fallback is the one most users will hit and
        it is the branch a careless change would drop.
        """
        region = "AAAAUGUGGCCCCCCCCCCCCGGGGGGGGG"
        pasted, _ = pick(picker, search_sequence=region)
        blank, printed = pick(picker, search_sequence="")

        assert pasted["selection"].found and blank["selection"].found
        assert pasted["selection"].chosen.source_id == "supplied"
        # Both are RNA already: windows come back as RNA, which is what a PPR reads.
        assert pasted["selection"].chosen.sequence in region

        assert blank["selection"].chosen.source_id == "rbcL"
        assert "annotated span of rbcL" in printed
        assert pasted["selection"].chosen.sequence != blank["selection"].chosen.sequence

    def test_apply_to_target_decides_whether_the_target_changes(self, picker):
        """Left off the cell reports and changes nothing; ticked it sets `target_rna`.

        The chosen window is asserted *different* from the starting target first. Without
        that, "unchanged" would be true for the uninteresting reason that there was nothing
        to change, and a cell that always overwrote the target would pass.
        """
        off, off_printed = pick(picker, apply_to_target=False)
        window = off["selection"].chosen.sequence
        assert window != "AAAAUGUGG", (
            "the picker must choose something other than the starting target, or neither "
            "branch of this control is observable")

        assert off["target_rna"] == "AAAAUGUGG", "a Run all must not redesign silently"
        assert "is still AAAAUGUGG" in off_printed

        on, on_printed = pick(picker, apply_to_target=True)
        assert on["selection"].chosen.sequence == window
        assert on["target_rna"] == window
        assert f"target_rna is now {window}" in on_printed


@pytest.fixture(scope="module")
def table_resolution(cells):
    """The design cell's codon-table resolution, runnable with the settings form set.

    The settings cell is pure assignment, and the design cell's head -- everything before
    the `design_oneshot` call -- is the three-way precedence the form documents as "highest
    filled one wins". Sliced rather than reimplemented, so this cannot agree with a copy of
    the rule while production does something else.
    """
    settings = cell_declaring(cells, "kazusa_taxid")
    design = next(c for c in cells if "chosen_codon_table" in c)
    head = design[:design.index("result = design_oneshot(")]

    def resolve(**form):
        env, _ = run_cell(with_form_values(settings, **form) + "\n" + head)
        return env

    return resolve


class TestTheCodonTableControls:
    """`codon_table_file` and `kazusa_taxid`, the two ways to override the host dropdown.

    Both are notebook-only: `design_oneshot` takes one `codon_table`, and which of the
    three sources fills it is decided in the notebook. That decision is the control.
    """

    def test_an_empty_form_inherits_the_host_dropdown(self, table_resolution):
        """The baseline the other two cases must differ from."""
        env = table_resolution(codon_table_file="", kazusa_taxid=0)
        assert env["chosen_codon_table"] is None, "None means 'use the organism dropdown'"
        assert env["chosen_table_source"] == env["organism"]

    def test_kazusa_taxid_replaces_the_host_table(self, table_resolution,
                                                  host_codon_table):
        """A typed taxid must give that organism's table, not the dropdown's.

        This is the defect the design cell's own comment records: the library cell once
        passed `codon_table=None` and fell back to the dropdown, so typing a taxid gave
        that host's design and a Chlamydomonas library two cells later.
        """
        taxid = ORGANISMS["e_coli"][0]
        try:
            expected = table_from_kazusa(taxid)
        except Exception as exc:                    # pragma: no cover - needs the cache
            pytest.skip(f"Kazusa table for {taxid} unavailable: {exc}")

        assert expected != host_codon_table, "the two tables must differ to begin with"

        env = table_resolution(codon_table_file="", kazusa_taxid=taxid)
        assert env["chosen_codon_table"] == expected
        assert str(taxid) in env["chosen_table_source"]

    def test_codon_table_file_wins_over_both(self, table_resolution):
        """The documented precedence: a filename beats a taxid, which beats the dropdown."""
        env = table_resolution(codon_table_file="usage.csv",
                               kazusa_taxid=ORGANISMS["e_coli"][0])
        assert env["chosen_codon_table"] == "usage.csv"
        assert "usage.csv" in env["chosen_table_source"]

    def test_a_supplied_table_file_changes_the_codons_chosen(self, ref, tmp_path):
        """End to end, because the resolution above only proves which path was taken.

        The file is the toy table with Leu forced to `TTA`. Measured: the free design
        spells Leu `CTG` five times and never `TTA`, and the design from this file spells
        all nine `TTA`. Picking a table that merely differed somewhere would have left the
        assertion unable to say the file was read rather than ignored.
        """
        rows = ["codon,frequency"]
        for aa, codons in TOY_TABLE.items():
            for codon in sorted(codons):
                dna = codon.replace("U", "T")
                rows.append(f"{dna},"
                            f"{(1.0 if dna == 'TTA' else 0.0) if aa == 'L' else codons[codon]}")
        csv = tmp_path / "usage.csv"
        csv.write_text("\n".join(rows), encoding="utf-8")

        def leucines(result):
            return [result["cds"][3 * i:3 * i + 3]
                    for i, aa in enumerate(result["protein"]) if aa == "L"]

        before = leucines(ref)
        assert before.count("TTA") == 0 and before.count("CTG") > 0, (
            "the free design must not already prefer TTA, or the file proves nothing")

        alt = design_oneshot(**{**BASE, "codon_table": str(csv)})
        after = leucines(alt)
        assert len(after) == len(before) > 0
        assert after.count("TTA") == len(after), "every Leu must take the file's codon"
        assert alt["constraints_ok"]
        assert alt["protein"] == ref["protein"]


@pytest.fixture(scope="module")
def library_cell(cells):
    return cell_declaring(cells, "targets")


#: What the library cell inherits from the cells above it.
LIBRARY_CONTEXT = dict(chosen_codon_table=TOY_TABLE, organism="c_reinhardtii_nuclear",
                       enzyme_profile="igem_rfc1000", seed=42, check_offtarget=False,
                       write_files=False)


def run_library(library_cell, targets: str) -> tuple[dict, str]:
    return run_cell(with_form_values(library_cell, targets=targets),
                    dict(LIBRARY_CONTEXT))


@pytest.fixture(scope="module")
def two_member_library(library_cell):
    env, _ = run_library(library_cell, "AAAAUGUGG, GCUAAAGAC")
    assert not env["lib"].failed
    return env["lib"]


class TestTheLibraryControls:
    """`targets` and `ppr_score_table`, both reached only from the library cells."""

    def test_targets_decides_which_members_the_library_has(self, library_cell):
        """Two target sets, two libraries. Run through the cell, so the parsing counts.

        `targets` is a comma-separated string that the cell splits, strips and upper-cases
        before `design_library` ever sees it. That parsing is part of the control.
        """
        a, _ = run_library(library_cell, "AAAAUGUGG, GCUAAAGAC")
        b, _ = run_library(library_cell, " aaaaugugg ,UUACACGUG, ")

        assert a["target_list"] == ["AAAAUGUGG", "GCUAAAGAC"]
        assert b["target_list"] == ["AAAAUGUGG", "UUACACGUG"], (
            "whitespace and case are normalised and a trailing comma is not a target")

        assert list(a["lib"].targets) != list(b["lib"].targets)
        assert not a["lib"].failed and not b["lib"].failed
        assert a["lib"].crosstalk() != b["lib"].crosstalk(), (
            "a different target set must give a different pairwise answer")

    def test_ppr_score_table_adds_the_affinity_tier(self, cells, two_member_library,
                                                    tmp_path):
        """Blank gives tier A alone; a path adds tier B as an annotation.

        **The real table is deliberately not vendored** -- PPRmatcher carries no licence,
        and it was derived from P-type motifs while this scaffold is S-type (see CLAUDE.md
        on `crosstalk.py`). So the stand-in is the minimal table `test_crosstalk.py`
        already keeps inline, in the documented PPRmatcher layout. What is under test is
        the control's effect on the output, not the table's biology: with no path the
        affinity tier is absent, and with one it appears, carrying its own caveat that it
        gates nothing.
        """
        cell = cell_declaring(cells, "ppr_score_table")
        scores = tmp_path / "scores.tsv"
        scores.write_text(PPR_SCORE_TABLE, encoding="utf-8")

        blank, tier_a = run_cell(with_form_values(cell, ppr_score_table=""),
                                 {"lib": two_member_library})
        given, both = run_cell(with_form_values(cell, ppr_score_table=str(scores)),
                               {"lib": two_member_library})

        assert blank["scores"] is None
        assert given["scores"] is not None and set(given["scores"]) == {"TN", "NN",
                                                                        "TD", "ND"}
        assert both != tier_a

        for marker in ("A->B", "flagged by predicted affinity", "gates nothing"):
            assert marker not in tier_a, f"{marker!r} must not appear without a table"
            assert marker in both, f"{marker!r} must appear once a table is given"


class TestWriteFilesSelectsTheOutputDirectory:
    """`write_files` is the notebook's name for `outdir`, and it must decide every call.

    Paired with `test_outdir_writes_the_order_files`, which proves the two branches differ
    in output. This one proves `write_files` is what chooses between them -- structurally,
    because the alternative is a design run per notebook cell for a parameter whose effect
    is already measured.
    """

    @staticmethod
    def _outdir_arguments(cells):
        out = []
        for cell in cells:
            for node in ast.walk(ast.parse(cell)):
                if not isinstance(node, ast.Call):
                    continue
                func = node.func
                name = func.id if isinstance(func, ast.Name) else getattr(func, "attr", "")
                if name not in {"design_oneshot", "design_library"}:
                    continue
                out.append((name, next((kw.value for kw in node.keywords
                                        if kw.arg == "outdir"), None)))
        return out

    def test_every_design_call_takes_its_outdir_from_write_files(self, cells):
        calls = self._outdir_arguments(cells)
        assert calls, "no design call found in the notebook"

        for name, outdir in calls:
            assert outdir is not None, f"{name} is called without an outdir argument"
            assert isinstance(outdir, ast.IfExp), (
                f"{name} passes an unconditional outdir, so write_files cannot turn "
                f"writing off")
            assert isinstance(outdir.test, ast.Name), (
                f"{name}: outdir is not switched on a plain form control")
            assert outdir.test.id == "write_files", (
                f"{name} switches its outdir on {outdir.test.id!r}, not on write_files")
            assert isinstance(outdir.orelse, ast.Constant) and outdir.orelse.value is None, (
                f"{name}: write_files off must pass None, which is what stops the write")
