#!/usr/bin/env python3
"""Terminus-extension demo: add residues to a chain end and relax the addition.

A worked example of ``rosetta_pipelines.extend_terminus`` — add a short residue run to the N- or
C-terminus of a chain and run a coordinate-restrained FastRelax so the addition settles while the
pre-existing structure is held near the template.

Here the input is ``W.pdb`` (receptor chain A + a 33-residue peptide on chain C). The demo appends a
``GGY`` tag to the **C-terminus** of chain C — e.g. a glycine linker plus a tyrosine for A280
quantitation — and relaxes the whole peptide chain, writing the result to ``out/W_Cext.pdb``.

Run it (from this folder, with the env active):

    python extend_terminus.py
"""

from pathlib import Path

# The whole public terminus API is these two names.
from rosetta_pipelines import RelaxParams, extend_terminus

HERE = Path(__file__).parent
INPUT_PDB = HERE / "W.pdb"
OUT_DIR = HERE / "out"

TARGET_CHAIN = "C"
ADD_SEQUENCE = "GGY"   # residues to add (here: Gly-Gly linker + Tyr quantitation tag)
END = "C"              # "C" appends (numbering preserved); "N" prepends (chain renumbered from 1)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    out_pdb = OUT_DIR / "W_Cext.pdb"

    params = RelaxParams(nstruct=3, relax_repeats=5)

    result = extend_terminus(
        INPUT_PDB,
        chain_id=TARGET_CHAIN,
        add_seq=ADD_SEQUENCE,
        end=END,
        out_pdb=out_pdb,
        params=params,
        relax_scope="chain",   # relax the whole peptide chain. ("added": only the new residues —
                              #   but that can't relieve the new terminal bond's strain on the
                              #   neighbouring existing residue, so "chain" is the default.)
    )

    # `result` is a RelaxResult dataclass.
    print(f"\nAppended '{ADD_SEQUENCE}' to the {END}-terminus of chain {TARGET_CHAIN} of {INPUT_PDB.name}")
    print(f"  best score (REU)   : {result.best_score:.2f}")
    print(f"  trajectory spread  : {result.score_std:.2f} (std across {len(result.trajectory_scores)})")
    print(f"  Cα drift (mobile)  : {result.ca_rmsd_drift_mobile:.3f} Å")
    print(f"  wrote extended PDB : {result.out_pdb}")

    # --- Other ways to call it (not run here) -----------------------------------------------------
    # N-terminal extension (prepend; the chain is renumbered from 1):
    #     extend_terminus(INPUT_PDB, "C", "W", end="N", out_pdb=out_pdb)   # e.g. an N-terminal Trp tag
    #
    # Settle only the newly-added residues instead of the whole chain:
    #     extend_terminus(INPUT_PDB, "C", "GGY", end="C", out_pdb=out_pdb, relax_scope="added")


if __name__ == "__main__":
    main()
