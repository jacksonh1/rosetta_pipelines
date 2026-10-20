#!/usr/bin/env python3
"""Relaxation demo: coordinate-restrained FastRelax of a structure.

A worked example of ``rosetta_pipelines.relax_structure`` — the "structure relaxation" operation.
Relaxation repacks sidechains and runs FastRelax while holding the backbone near the input with
coordinate restraints, so the structure settles into a low-energy conformation without drifting away
from the template.

Here the input is ``W.pdb`` (receptor chain A + a 33-residue peptide on chain C). The demo relaxes
just the peptide chain C plus a 6 Å neighbor shell (everything else frozen), which keeps the run
fast, and writes the result to ``out/W_relaxed.pdb``.

Run it (from this folder, with the env active):

    python relax_demo.py
"""

from pathlib import Path

# The whole public relaxation API is these two names.
from rosetta_pipelines import RelaxParams, relax_structure

HERE = Path(__file__).parent
INPUT_PDB = HERE / "W.pdb"
OUT_DIR = HERE / "out"


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    out_pdb = OUT_DIR / "W_relaxed.pdb"

    # RelaxParams holds the FastRelax engine knobs. The coordinate restraint (default weight 0.5)
    # keeps the backbone near the input; nstruct independent trajectories are run and the best kept.
    params = RelaxParams(
        nstruct=3,          # independent relax trajectories; the lowest-scoring is kept
        relax_repeats=5,    # FastRelax fa_rep ramp depth per trajectory
        neighbor_dist=6.0,  # radius (Å) of the shell relaxed around the focus chain
    )

    result = relax_structure(
        INPUT_PDB,
        out_pdb,
        params=params,
        focus_chain="C",    # relax chain C + its 6 Å shell, freeze the rest.
                            # Use focus_chain=None to relax the WHOLE structure (slower).
    )

    # `result` is a RelaxResult dataclass — everything the run measured.
    print(f"\nRelaxed chain C (+ shell) of {INPUT_PDB.name}")
    print(f"  best score (REU)    : {result.best_score:.2f}")
    print(f"  mean / worst (REU)  : {result.mean_score:.2f} / {result.worst_score:.2f}")
    print(f"  trajectory spread   : {result.score_std:.2f} (std across {len(result.trajectory_scores)})")
    print(f"  Cα drift, whole/mob : {result.ca_rmsd_drift:.3f} / {result.ca_rmsd_drift_mobile:.3f} Å")
    print(f"  restraint energy    : {result.coord_cst_energy:.2f} REU")
    print(f"  wrote relaxed PDB   : {result.out_pdb}")

    # --- Other ways to call it (not run here) -----------------------------------------------------
    # Whole-structure relax (every residue mobile):
    #     relax_structure(INPUT_PDB, out_pdb, params=params, focus_chain=None)
    #
    # Free relax (let the backbone move without a template tether):
    #     relax_structure(INPUT_PDB, out_pdb, params=RelaxParams(coord_cst_weight=0.0))
    #
    # Form disulfides before relaxing (PDB numbering, ((chain, resnum), (chain, resnum))):
    #     relax_structure(INPUT_PDB, out_pdb, focus_chain="A",
    #                     disulfide_pairs=[(("A", 261), ("A", 321)), (("A", 367), ("A", 425))])


if __name__ == "__main__":
    main()
