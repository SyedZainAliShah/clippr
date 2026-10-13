"""Static checks on the inventory Colab notebook.

**Why these exist.** `tests/test_workflow.py` passes `codon_table=` to every workflow task
explicitly, so it can never catch the notebook omitting it. That is exactly what happened:
Step 6 called `select_interface` without the host table, the binding guard correctly refused
for any organism other than the default, and the cell then raised `KeyError: 'inventory'`
reaching into the artefacts of a result it never checked. Two defects, both invisible to a
suite that only exercises the library.

Static like `test_notebook.py` -- parse and cross-reference, no execution.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "PPROPT_inventory.ipynb"

#: Workflow tasks whose results the notebook measures under one host codon table. A front,
#: a recode and a search are all bound to the table they were computed with; passing a
#: different one (or none, which resolves the cached default) reports one host's numbers for
#: another host's DNA.
TABLE_BOUND_TASKS = {
    "recode_for_host", "optimise_collection", "explore_interfaces", "select_interface",
}


@pytest.fixture(scope="module")
def cells():
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def _calls(source):
    """Every `w.<task>(...)` call in a cell, as (task, ast.Call)."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []
    out = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "w"):
            out.append((node.func.attr, node))
    return out


class TestTheHostTableReachesEveryTaskBoundToIt:
    """The defect: Step 6 resolved the default table instead of the chosen host's."""

    def test_every_table_bound_task_is_given_the_host_table(self, cells):
        """Not passing it is silent on the default organism and a refusal on every other.

        `select_interface` shipped without `codon_table=`. Against `c_reinhardtii_nuclear`
        the resolved default happened to match the front's binding, so the step worked and
        looked correct. Against `e_coli` the sha differed, `binding_mismatches` reported it,
        and the notebook died on the refusal it had asked for.
        """
        missing = []
        for source in cells:
            for task, call in _calls(source):
                if task not in TABLE_BOUND_TASKS:
                    continue
                passed = {kw.arg for kw in call.keywords}
                positional = any(isinstance(a, ast.Name) and a.id == "host_codon_table"
                                 for a in call.args)
                if "host_codon_table" not in passed and not positional:
                    if not any(isinstance(kw.value, ast.Name)
                               and kw.value.id == "host_codon_table"
                               for kw in call.keywords):
                        missing.append(task)
        assert not missing, (
            f"these tasks are bound to the host codon table but the notebook does not pass "
            f"it: {sorted(set(missing))}. The default resolves to the cached Chlamydomonas "
            f"table, so the cell works for one organism and refuses for the rest.")


class TestARefusalIsPrintedRatherThanIndexed:
    """The second defect: `.artefacts[...]` on a result nobody checked."""

    def test_every_result_whose_artefacts_are_read_passed_through_done(self, cells):
        """`KeyError: 'inventory'` is what a refusal looked like to the user.

        Every workflow task can refuse, and a refused `WorkflowResult` carries no artefacts.
        The notebook's contract is that a name holding a task result is either `None` or
        succeeded, which `done()` establishes: it returns the result on success and prints
        the failures and returns `None` otherwise. A name assigned straight from `w.<task>()`
        and then indexed has skipped that, so a designed refusal surfaces as a traceback.
        """
        from_task, cleared, indexed = set(), set(), {}
        for source in cells:
            try:
                tree = ast.parse(source)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign):
                    names = [t.id for t in node.targets if isinstance(t, ast.Name)]
                    value = node.value
                    if (isinstance(value, ast.Call) and isinstance(value.func, ast.Attribute)
                            and isinstance(value.func.value, ast.Name)
                            and value.func.value.id == "w"):
                        from_task.update(names)
                    elif (isinstance(value, ast.Call) and isinstance(value.func, ast.Name)
                            and value.func.id == "done"):
                        cleared.update(names)
                # `"key" in x.artefacts` makes the KeyError impossible by construction,
                # which is the right guard where ok=False still writes an artefact: one
                # unspellable target does not invalidate the targets that compiled.
                if (isinstance(node, ast.Compare)
                        and any(isinstance(o, ast.In) for o in node.ops)):
                    for side in node.comparators:
                        if (isinstance(side, ast.Attribute) and side.attr == "artefacts"
                                and isinstance(side.value, ast.Name)):
                            cleared.add(side.value.id)
                if (isinstance(node, ast.Subscript)
                        and isinstance(node.value, ast.Attribute)
                        and node.value.attr == "artefacts"
                        and isinstance(node.value.value, ast.Name)):
                    indexed.setdefault(node.value.value.id, ast.unparse(node))

        offenders = {name: expr for name, expr in indexed.items()
                     if name in from_task and name not in cleared}
        assert not offenders, (
            "these results are indexed for artefacts without passing through done() or "
            "checking the key is present, so a refusal surfaces as a KeyError instead of "
            f"its reason: {offenders}")

    def test_the_helper_the_cells_depend_on_is_defined(self, cells):
        """The guards above are worthless if `done` is not in the notebook's namespace."""
        assert any("def done(result):" in source for source in cells), (
            "no cell defines done(); every refusal guard in this notebook calls it")
