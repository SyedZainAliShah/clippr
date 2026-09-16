"""Run the reference's reusable-library workflow and print its result as JSON.

Executed by `m3_reference_comparison.py` under the *reference's own* interpreter, never ours.
It lives in its own file rather than inside a string so that what runs against the reference
is readable, diffable and reviewable on its own terms.

    <reference venv>/python m3_reference_runner.py <input_dir> <output_dir>
"""
import json
import sys
import time
from pathlib import Path

import pandas as pd
import grasp_library as g

inp, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

config = g.build_default_config(inp)
codon_data = g.load_codon_usage(config["codon_usage_file"])[1]
parts = pd.read_csv(inp / "parts.csv")

started = time.perf_counter()
result = g.run_library_redesign_and_anneal(parts=parts, codon_data=codon_data, config=config,
                                           input_dir=inp, output_dir=out, seed=42,
                                           log=lambda *a, **k: None)
elapsed = time.perf_counter() - started

library = result["optimized_library"]
modules = {}
for _, row in library.iterrows():
    modules[str(row["part_id"])] = {
        "cds": str(row["optimized_cds"]),
        "aa": str(row["aa_sequence"]),
        "objective_score": float(row["objective_score"]),
        "codon_score": float(row["codon_score"]),
        "qc_status": str(row["qc_status"]),
        "hard_constraints_passed": bool(row["hard_constraints_passed"]),
    }

reported = ("site_blacklist", "ppr_5prime_fusion_site", "genetic_code",
            "synthesis", "optimizer", "pareto")
payload = {
    "version": getattr(g, "__version__", "?"),
    "origin": g.__file__,
    "seconds": elapsed,
    "rows": int(library.shape[0]),
    "modules": modules,
    "config": {k: str(config[k]) for k in reported if k in config},
}
print("JSON_BEGIN" + json.dumps(payload) + "JSON_END")
