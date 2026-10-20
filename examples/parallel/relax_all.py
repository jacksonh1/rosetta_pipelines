#!/usr/bin/env python3
"""Relax every PDB in a directory in parallel, each worker writing its own output.

Usage:
    python relax_all.py <input_dir> <output_dir> [--workers N]

Demonstrates the two `spawn` requirements: the worker function (`relax_one`) is defined at module
level, and the launch is guarded by `if __name__ == "__main__":`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from rosetta_pipelines.pool import run_pool
from rosetta_pipelines.relax import RelaxParams, relax_structure

PARAMS = RelaxParams(nstruct=5, relax_repeats=5)

# Set at launch and read inside the worker. Under `spawn` the module is re-imported in each worker,
# so a module-level global defined before the guarded launch is available to the worker.
OUTPUT_DIR = Path("relaxed")


def relax_one(pdb_path: str) -> str:
    """Relax one structure and return the output path. Module-level so it is picklable by `spawn`."""
    out = OUTPUT_DIR / (Path(pdb_path).stem + ".relaxed.pdb")
    result = relax_structure(pdb_path, str(out), params=PARAMS, focus_chain=None)
    return f"{result.out_pdb}  (best={result.best_score:.1f} REU)"


def main() -> None:
    global OUTPUT_DIR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_dir", help="Directory of input .pdb files")
    ap.add_argument("output_dir", help="Directory for relaxed outputs")
    ap.add_argument("--workers", type=int, default=4, help="Number of worker processes")
    args = ap.parse_args()

    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(str(p) for p in Path(args.input_dir).glob("*.pdb"))
    if not pdbs:
        raise SystemExit(f"No .pdb files found in {args.input_dir}")

    for done in run_pool(pdbs, relax_one, max_workers=args.workers):
        print("relaxed ->", done)


if __name__ == "__main__":  # required under spawn
    main()
