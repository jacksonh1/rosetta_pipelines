#!/usr/bin/env python3
"""Threading demo: put a poly-lysine sequence onto one chain of a complex and relax it.

This is a worked example of ``rosetta_pipelines.thread_and_relax`` — the "sequence threading"
operation. Threading replaces the amino-acid identities on a single chain with a new sequence of the
*same length*, keeping that chain's backbone and rebuilding its sidechains, then runs a
coordinate-restrained FastRelax to settle the result.

Here the target is chain C of ``W.pdb`` (a 33-residue peptide bound to receptor chain A). We thread a
poly-lysine sequence (lysine at every position) and relax chain C plus a 6 Å neighbor shell, leaving
the rest of the complex frozen. The relaxed structure is written to ``out/W_polyK.pdb``.

Run it (from anywhere, with the env active):

    python thread_polylysine.py

The only hard rule for threading: ``new_seq`` must be exactly as long as the target chain. Chain C
here is 33 residues, so the poly-lysine sequence is 33 lysines.
"""

from pathlib import Path

# The whole public threading API is these two names.
from rosetta_pipelines import RelaxParams, thread_and_relax

HERE = Path(__file__).parent
INPUT_PDB = HERE / "W.pdb"
OUT_DIR = HERE / "out"

# --- what to thread -------------------------------------------------------------------------------
TARGET_CHAIN = "C"
CHAIN_C_LENGTH = 33  # number of residues on chain C of W.pdb
NEW_SEQUENCE = "K" * CHAIN_C_LENGTH  # poly-lysine, same length as the chain


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    out_pdb = OUT_DIR / "W_polyK.pdb"

    # RelaxParams holds the FastRelax engine knobs (how thorough the relax is). Defaults are sensible;
    # here we keep the coordinate restraint on (default weight 0.5) so the backbone stays near the
    # input, and run 3 independent trajectories, keeping the lowest-scoring.
    params = RelaxParams(
        nstruct=3,          # independent relax trajectories; the best is kept
        relax_repeats=5,    # FastRelax fa_rep ramp depth per trajectory
        neighbor_dist=6.0,  # radius (Å) of the shell relaxed around chain C
    )

    # thread_and_relax: mutate chain C to the new sequence, then relax.
    result = thread_and_relax(
        INPUT_PDB,                   # template complex
        chain_id=TARGET_CHAIN,       # which chain to re-sequence
        new_seq=NEW_SEQUENCE,        # the new sequence (length must match the chain)
        out_pdb=out_pdb,             # where to write the relaxed structure
        params=params,
        thread_scope="all",          # "all": rebuild every position. ("changed": only positions that
                                     #   differ from the template — pointless here, we replace all 33.)
        relax_scope="neighborhood",  # relax chain C + its 6 Å shell. ("full": whole complex;
                                     #   "mutations": only the changed positions + shell.)
    )

    # `result` is a RelaxResult dataclass — everything the run measured.
    print(f"\nThreaded {CHAIN_C_LENGTH}-residue poly-lysine onto chain {TARGET_CHAIN} of {INPUT_PDB.name}")
    print(f"  positions mutated  : {result.n_threaded}")
    print(f"  best score (REU)   : {result.best_score:.2f}")
    print(f"  trajectory spread  : {result.score_std:.2f} (std across {len(result.trajectory_scores)})")
    print(f"  Cα drift (mobile)  : {result.ca_rmsd_drift_mobile:.3f} Å")
    print(f"  wrote relaxed PDB  : {result.out_pdb}")


if __name__ == "__main__":
    main()
