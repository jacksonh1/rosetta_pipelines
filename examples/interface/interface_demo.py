#!/usr/bin/env python3
"""Interface-analysis demo: pack-separated binding energy + interface descriptors.

A worked example of ``rosetta_pipelines.run_interface_analyzer`` — the "interface analysis" operation.
It pulls the two partners of a complex apart, repacks their sidechains in the unbound state, and
rescores, giving the physically standard binding ΔG (not a rigid frozen-coordinate estimate). It also
reports geometry/quality descriptors (buried area, shape complementarity, packing, buried-unsat
H-bonds).

Here the input is ``W.pdb``: receptor chain **A** vs peptide chain **C**.

Run it (from this folder, with the env active):

    python interface_demo.py

Or use the installed command-line entry point (no script needed):

    rosetta-interface-analyzer --pdb W.pdb --receptor A --peptide C
"""

from pathlib import Path

# The whole public interface API is this one function.
from rosetta_pipelines import run_interface_analyzer

HERE = Path(__file__).parent
INPUT_PDB = HERE / "W.pdb"


def main() -> None:
    result = run_interface_analyzer(
        INPUT_PDB,
        receptor_chains=["A"],   # side 1 (one or more chain IDs)
        peptide_chains=["C"],    # side 2
        pack_separated=True,     # repack sidechains in the unbound state -> standard binding ΔG.
                                # pack_separated=False gives a rigid interaction-energy estimate.
        pack_input=False,        # leave the bound complex as-is (it is already relaxed here).
    )

    # `result` is an InterfaceAnalyzerResult dataclass. Energies are REU; areas are Å².
    print(f"\nInterface {result.interface} of {INPUT_PDB.name}")
    print("  --- energetics ---")
    print(f"  dG_separated          : {result.dG_separated:.2f} REU   (binding energy; lower = tighter)")
    print(f"  dG_dSASA_ratio        : {result.dG_dSASA_ratio:.4f}      (energy per 100 Å² buried)")
    print(f"  crossterm_interface_E : {result.crossterm_interface_energy:.2f} REU   (rigid, no repack)")
    print("  --- geometry / quality ---")
    print(f"  dSASA (buried area)   : {result.dSASA:.1f} Å²")
    print(f"  n_interface_residues  : {result.n_interface_residues}")
    print(f"  sc_value (shape comp) : {result.sc_value:.3f}          (0-1; higher = tighter fit)")
    print(f"  packstat              : {result.packstat:.3f}          (0-1; higher = better packed)")
    print(f"  delta_unsat_hbonds    : {result.delta_unsat_hbonds:.1f}            (buried unsatisfied; lower = better)")
    print(f"  interface_hbonds      : {result.interface_hbonds}  (E={result.total_hbond_E:.2f} REU)")


if __name__ == "__main__":
    main()
