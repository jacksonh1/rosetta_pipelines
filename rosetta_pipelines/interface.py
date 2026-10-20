#!/usr/bin/env python3
"""Analyze a two-partner (e.g. domain–peptide) interface with PyRosetta's InterfaceAnalyzerMover.

Unlike a rigid interaction-energy estimate (bound vs. unbound scored at frozen coordinates), this
reports the **pack-separated** binding energy: the partners are pulled apart AND their sidechains are
repacked in the unbound state before rescoring. That removes the affinity a frozen-unbound estimate
wrongly credits to interface rotamers that were optimized for the bound state — the relevant
correction when scoring threaded/relaxed designs.

It also returns the interface descriptors InterfaceAnalyzerMover computes for free: buried surface
area (dSASA), shape complementarity (sc_value), packing quality (packstat), buried-unsatisfied
H-bonds, interface H-bond count/energy, and the energy-per-buried-area ratio.

Energies are in Rosetta Energy Units (REU); areas in Å².

Command-line use:
  python -m rosetta_pipelines.interface --pdb complex.pdb --receptor A --peptide B
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

from pyrosetta import rosetta

from rosetta_pipelines import pose_tools, rosetta_init


@dataclass
class InterfaceAnalyzerResult:
    pdb: str
    interface: str
    """Interface spec passed to the mover, e.g. 'A_B' or 'AB_C'."""
    weights: str
    pack_separated: bool
    pack_input: bool

    # --- energetics ---
    dG_separated: float
    """Binding energy: E_complex - E_unbound. If pack_separated, the unbound state is repacked first,
    so this is the physically standard binding ΔG."""
    dG_dSASA_ratio: float
    """dG_separated per buried area (Rosetta convention: per 100 Å²). Energy density — a
    size-normalized ranking metric."""
    crossterm_interface_energy: float
    """Rigid cross-interface interaction energy (no unbound repacking), kept for comparison against
    dG_separated: the gap is the sidechain-reorganization cost."""

    # --- geometry / quality descriptors ---
    dSASA: float
    """Total buried interface surface area (Å²)."""
    n_interface_residues: int
    sc_value: float
    """Lawrence & Colman shape complementarity (0–1; higher = tighter fit)."""
    packstat: float
    """Interface packing quality (0–1; higher = better packed)."""
    delta_unsat_hbonds: float
    """Buried polar atoms that lose an H-bond partner on complex formation (lower = better; each is
    an energetic liability)."""
    interface_hbonds: int
    total_hbond_E: float


def _assert_chains_present(pose: "rosetta.core.pose.Pose", chains: list[str]) -> None:
    pdbi = pose.pdb_info()
    if pdbi is None:
        raise ValueError("Pose has no PDBInfo; cannot select chains by PDB chain ID.")
    present = {pdbi.chain(i) for i in range(1, pose.total_residue() + 1)}
    missing = sorted(set(chains) - present)
    if missing:
        raise ValueError(
            f"Requested chain(s) {missing} not found in PDB. Present chains: {sorted(present)}"
        )


def run_interface_analyzer(
    pdb: str | Path,
    receptor_chains: list[str],
    peptide_chains: list[str],
    weights: str = "ref2015",
    pack_separated: bool = True,
    pack_input: bool = False,
    extra_res_fa: list[str] | None = None,
    mute: str = "all",
) -> InterfaceAnalyzerResult:
    """Run InterfaceAnalyzerMover on a complex PDB and return interface metrics.

    Args:
        pdb: Input complex PDB.
        receptor_chains: Chain ID(s) on side 1, e.g. ['A'] or ['A', 'B'].
        peptide_chains: Chain ID(s) on side 2, e.g. ['B'].
        weights: Rosetta weights set (e.g. 'ref2015', 'beta_nov16').
        pack_separated: Repack sidechains in the unbound state before scoring. True gives the
            physically standard binding ΔG; False gives a rigid interaction-energy estimate.
        pack_input: Repack the bound complex before scoring. Leave False for already-relaxed inputs.
        extra_res_fa: Optional .params files for noncanonical residues.
        mute: PyRosetta mute level ('all', 'protocols', 'none').
    """
    rosetta_init.ensure_init(extra_res_fa=extra_res_fa or (), mute=mute)

    pdb = str(pdb)
    # Keep zero-occupancy atoms (Rosetta drops them by default) and tolerate spaces in the path.
    pose = pose_tools.load_pose(pdb)
    sfxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(weights)

    _assert_chains_present(pose, receptor_chains)
    _assert_chains_present(pose, peptide_chains)

    overlap = set(receptor_chains) & set(peptide_chains)
    if overlap:
        raise ValueError(
            f"Receptor and peptide chain selections overlap ({sorted(overlap)}); "
            "check receptor/peptide chains."
        )

    # Interface spec is "side1_side2", e.g. "AB_C". Set it via the setter (default-construct +
    # set_interface) rather than the string constructor: the string constructor was removed in newer
    # PyRosetta builds (2026.39+), but set_interface accepts the same string across versions.
    interface = "".join(receptor_chains) + "_" + "".join(peptide_chains)

    iam = rosetta.protocols.analysis.InterfaceAnalyzerMover()
    iam.set_interface(interface)
    iam.set_scorefunction(sfxn)
    iam.set_pack_separated(pack_separated)
    iam.set_pack_input(pack_input)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_sc(True)
    iam.set_compute_interface_delta_hbond_unsat(True)
    iam.set_use_tracer(False)
    iam.apply(pose)

    # Several InterfaceData fields are vector1 (index 1 = whole interface, 2/3 = per side); the scalar
    # getters return the whole-interface value directly.
    data = iam.get_all_data()

    return InterfaceAnalyzerResult(
        pdb=pdb,
        interface=interface,
        weights=weights,
        pack_separated=pack_separated,
        pack_input=pack_input,
        dG_separated=float(iam.get_interface_dG()),
        dG_dSASA_ratio=float(data.dG_dSASA_ratio),
        crossterm_interface_energy=float(iam.get_crossterm_interface_energy()),
        dSASA=float(iam.get_interface_delta_sasa()),
        n_interface_residues=int(iam.get_num_interface_residues()),
        sc_value=float(data.sc_value),
        packstat=float(iam.get_interface_packstat()),
        delta_unsat_hbonds=float(iam.get_interface_delta_hbond_unsat()),
        interface_hbonds=int(data.interface_hbonds),
        total_hbond_E=float(iam.get_total_Hbond_E()),
    )


def _parse_chain_string(s: str) -> list[str]:
    """Accept 'A', 'A,B', 'A B', or 'AB' → ['A', 'B']."""
    s = s.replace(",", " ").strip()
    chains = [c for c in s.split() if c]
    if len(chains) == 1 and len(chains[0]) > 1 and chains[0].isalpha():
        chains = list(chains[0])
    return chains


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True, help="Input complex PDB")
    ap.add_argument("--receptor", required=True, help="Receptor chain ID(s), e.g. A or A,B")
    ap.add_argument("--peptide", required=True, help="Peptide chain ID(s), e.g. B")
    ap.add_argument("--weights", default="ref2015", help="Rosetta weights (e.g. ref2015)")
    ap.add_argument(
        "--no_pack_separated",
        action="store_true",
        help="Skip unbound-state repacking (rigid dG)",
    )
    ap.add_argument(
        "--pack_input",
        action="store_true",
        help="Repack the bound complex before scoring",
    )
    ap.add_argument(
        "--extra_res_fa",
        nargs="*",
        default=[],
        help="Optional .params files for noncanonical residues (space-separated)",
    )
    ap.add_argument("--mute", default="all", help='PyRosetta mute level ("all"/"protocols"/"none")')
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    result = run_interface_analyzer(
        pdb=args.pdb,
        receptor_chains=_parse_chain_string(args.receptor),
        peptide_chains=_parse_chain_string(args.peptide),
        weights=args.weights,
        pack_separated=not args.no_pack_separated,
        pack_input=args.pack_input,
        extra_res_fa=args.extra_res_fa,
        mute=args.mute,
    )
    for field_name, value in result.__dict__.items():
        print(f"{field_name:28s} {value}")


if __name__ == "__main__":
    main()
