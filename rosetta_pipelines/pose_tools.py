#!/usr/bin/env python3
"""Pose surgery helpers: load a pose, map chains, form disulfides, and thread a new sequence onto a
chain.

Threading replaces the residue identities on one chain of a template PDB with a new sequence of the
same length, keeping the backbone and rebuilding sidechains from ideal internal coordinates. The
downstream relax (see :mod:`rosetta_pipelines.fastrelax`) then settles the result.
"""

import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Literal

from pyrosetta import rosetta

# No init at import time (deliberate): importing this module must not lock in PyRosetta flags before
# a caller that needs -extra_res_fa can set them. Every pose-creating entry point calls
# ensure_init() first, so ordinary work still "just runs" with -mute all.


# ---- mapping helpers ----

ONE_TO_THREE = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


@contextmanager
def _rosetta_readable_path(pdb_path: str | Path):
    """Yield a path Rosetta's file reader can open.

    Rosetta's ``pose_from_file`` fails to open a path containing whitespace (it mis-tokenizes the
    name), which is easy to hit when a directory in the path has a space in it. When the path has a
    space, copy it into a space-free temp directory and yield that copy (cleaned up afterwards);
    otherwise yield the path unchanged. Writing (``dump_pdb``) is unaffected, so only reads need
    this."""
    pdb_path = str(pdb_path)
    if " " not in pdb_path:
        yield pdb_path
        return
    tmp_dir = tempfile.mkdtemp(prefix="rosetta_pipelines_")  # temp dirs have no spaces
    tmp_pdb = os.path.join(tmp_dir, os.path.basename(pdb_path).replace(" ", "_"))
    shutil.copy(pdb_path, tmp_pdb)
    try:
        yield tmp_pdb
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def load_pose(pdb_path, *, ignore_zero_occupancy: bool = False) -> rosetta.core.pose.Pose:
    """Load a pose from a PDB with explicit control of zero-occupancy atoms.

    ``ignore_zero_occupancy=False`` (default) **keeps** zero-occupancy atoms. Rosetta's own default
    drops them silently on load, which can delete an *internal* residue (e.g. a receptor loop the
    depositor left unmodelled) and introduce a chain break the caller never asked for — a silent
    structure edit. Pass ``True`` to reproduce Rosetta's drop-on-load. Per-load
    (``ImportPoseOptions``), so it is a normal knob, not a process-global init flag.

    Paths containing spaces are loaded via a space-free temp copy (see
    :func:`_rosetta_readable_path`), since Rosetta's reader cannot open them directly."""
    opts = rosetta.core.import_pose.ImportPoseOptions()
    opts.set_ignore_zero_occupancy(ignore_zero_occupancy)
    with _rosetta_readable_path(pdb_path) as safe_path:
        return rosetta.core.import_pose.pose_from_file(
            safe_path, opts, False, rosetta.core.import_pose.FileType.PDB_file
        )


def chain_indices_from_pdbinfo(pose: rosetta.core.pose.Pose, chain_id: str) -> list[int]:
    """Pose indices (1-based) of every residue on ``chain_id``, via the PDBInfo chain IDs."""
    info = pose.pdb_info()
    if info is None:
        raise ValueError("Pose lacks PDBInfo; cannot map chain IDs")
    indices = [i for i in range(1, pose.total_residue() + 1) if info.chain(i) == chain_id]
    if not indices:
        raise ValueError(f"Chain {chain_id} not found in pose")
    return indices


def add_disulfides_from_pairs(
    pose: rosetta.core.pose.Pose,
    pairs: list[tuple[int, int]],
):
    """Add disulfide bonds for specified pose-index pairs. Returns the pairs actually added."""
    conf = pose.conformation()

    def is_cys(idx: int) -> bool:
        res = pose.residue(idx)
        return res.name3() == "CYS" and res.has("SG")

    pairs_added: list[tuple[int, int]] = []
    for i, j in pairs:
        if i == j:
            raise ValueError("Disulfide pair indices must differ")
        if not is_cys(i) or not is_cys(j):
            raise ValueError(
                f"Residues {i} and/or {j} are not cysteines with SG present; cannot form disulfide"
            )
        if (
            pose.residue(i).type().is_disulfide_bonded()
            or pose.residue(j).type().is_disulfide_bonded()
        ):
            continue  # skip already bonded cysteines
        conf.add_disulfide_bond(i, j)
        pairs_added.append((i, j))

    return pairs_added


def add_disulfides_from_chain_positions(
    pose: rosetta.core.pose.Pose,
    pairs: list[tuple[tuple[str, int], tuple[str, int]]],
):
    """Add disulfides from chain/resseq pairs using PDB numbering.

    Args:
        pairs: list of ((chain_id, resseq_i), (chain_id, resseq_j)).
    """
    info = pose.pdb_info()
    if info is None:
        raise ValueError("Pose lacks PDBInfo; cannot map chain/resseq to pose indices")

    mapped: list[tuple[int, int]] = []
    mapping_dict = {}
    for (chain_i, res_i), (chain_j, res_j) in pairs:
        # icode defaults to a blank string; passing an int breaks newer PyRosetta bindings
        idx_i = info.pdb2pose(chain_i, res_i)
        idx_j = info.pdb2pose(chain_j, res_j)
        if idx_i <= 0 or idx_j <= 0:
            raise ValueError(
                f"Could not map pair ({chain_i},{res_i})-({chain_j},{res_j}) to pose indices"
            )
        mapped.append((idx_i, idx_j))
        mapping_dict[(idx_i, idx_j)] = ((chain_i, res_i), (chain_j, res_j))

    idx_pairs = add_disulfides_from_pairs(pose, mapped)
    return [mapping_dict[pair] for pair in idx_pairs]


# ---- in-pose threading (no disk round-trip, no residue-name string editing) ----


def _disulfide_partner(pose: rosetta.core.pose.Pose, i: int) -> int | None:
    """Pose index of residue ``i``'s disulfide partner, or None if it is not disulfide-bonded."""
    res = pose.residue(i)
    rt = res.type()
    if not rt.is_disulfide_bonded():
        return None
    sg = res.atom_index("SG")
    for c in range(1, rt.n_possible_residue_connections() + 1):
        if rt.residue_connect_atom_index(c) == sg:
            partner = res.residue_connection_partner(c)
            if partner != 0:
                return partner
    return None


def break_chain_disulfides(pose: rosetta.core.pose.Pose, chain_indices: list[int]):
    """Break every disulfide with an endpoint in ``chain_indices``, leaving other chains' bonds
    intact. Threading a new sequence can remove or move the chain's cysteines, so any bond Rosetta
    auto-detected on load must be cleared first — otherwise a mutated residue is left with a dangling
    disulfide connection."""
    conf = pose.conformation()
    broken: set[tuple[int, int]] = set()
    for i in chain_indices:
        partner = _disulfide_partner(pose, i)
        if partner is None:
            continue
        pair = (min(i, partner), max(i, partner))
        if pair in broken:
            continue
        rosetta.core.conformation.break_disulfide(conf, pair[0], pair[1])
        broken.add(pair)


def thread_sequence_onto_chain(
    pose: rosetta.core.pose.Pose,
    chain_id: str,
    new_seq: str,
    thread_scope: Literal["all", "changed"] = "all",
) -> list[int]:
    """Mutate ``chain_id`` in-place to match ``new_seq`` and return the pose indices actually
    mutated (backbone preserved, sidechain rebuilt from ideal internal coordinates).

    ``thread_scope="changed"`` mutates only positions whose current residue already differs from the
    target — unchanged positions keep their input rotamer as the relax seed instead of being rebuilt
    from ideal coords; ``"all"`` rebuilds every position. ASSUMES any disulfides on the chain were
    already broken (see :func:`break_chain_disulfides`)."""
    if thread_scope not in {"all", "changed"}:
        raise ValueError(f"thread_scope must be 'all' or 'changed', got {thread_scope!r}")
    indices = chain_indices_from_pdbinfo(pose, chain_id)
    new_seq = new_seq.upper()
    if len(indices) != len(new_seq):
        raise ValueError(
            f"Sequence length ({len(new_seq)}) must equal chain {chain_id} "
            f"residue count ({len(indices)})."
        )
    changed: list[int] = []
    for idx, aa in zip(indices, new_seq, strict=True):
        if aa not in ONE_TO_THREE:
            raise ValueError(f"Unsupported residue code '{aa}'")
        if thread_scope == "changed" and pose.residue(idx).name1() == aa:
            continue  # already this residue — keep the input rotamer as the seed
        rosetta.protocols.simple_moves.MutateResidue(idx, ONE_TO_THREE[aa]).apply(pose)
        changed.append(idx)
    return changed


# ---- terminus extension (add residues to a chain end) ----


def add_nterminus_residues(pose, chain_id: str, seq: str) -> tuple[int, int]:
    """Prepend ``seq`` to the N-terminus of ``chain_id`` in-place; rebuild the PDBInfo; return the
    (start, end) pose indices of the added residues.

    RENUMBERS the target chain continuously from 1 (the first prepended residue becomes resnum 1), so
    the caller must re-derive sequence/disulfides from the resulting pose rather than carry the
    parent's numbering. Pure pose surgery — no relax, no write."""
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("prepend sequence must be non-empty")
    info = pose.pdb_info()
    original_chain_map = {}  # pose_idx -> (chain, resnum, icode), captured before mutation
    if info:
        for i in range(1, pose.total_residue() + 1):
            original_chain_map[i] = (info.chain(i), info.number(i), info.icode(i))

    chain_indices = chain_indices_from_pdbinfo(pose, chain_id)
    chain_start = chain_indices[0]
    target_pose_chain = pose.chain(chain_start)

    rts = pose.residue_type_set_for_pose()
    for aa in reversed(seq):
        if aa not in ONE_TO_THREE:
            raise ValueError(f"Unsupported residue code '{aa}'")
        restype = rts.name_map(ONE_TO_THREE[aa])
        new_res = rosetta.core.conformation.ResidueFactory.create_residue(restype)
        pose.prepend_polymer_residue_before_seqpos(new_res, chain_start, True)

    # Prepending corrupts the PDBInfo; rebuild it. i<chain_start: unchanged (chains before target);
    # chain_start..+n-1: the new residues (target chain, renumbered from 1); rest: shifted by +n, with
    # the target chain renumbered continuously and other chains restored from the original map.
    n_added = len(seq)
    info = pose.pdb_info()
    if info is not None:
        for i in range(1, pose.total_residue() + 1):
            if i < chain_start:
                orig_chain, orig_resnum, orig_icode = original_chain_map[i]
                info.set_resinfo(i, orig_chain, orig_resnum, orig_icode)
            elif i < chain_start + n_added:
                info.set_resinfo(i, chain_id, i - chain_start + 1, " ")
            else:
                orig_idx = i - n_added
                orig_chain, orig_resnum, orig_icode = original_chain_map[orig_idx]
                if orig_chain == chain_id:
                    info.set_resinfo(i, chain_id, i - chain_start + 1, " ")
                else:
                    info.set_resinfo(i, orig_chain, orig_resnum, orig_icode)

    target_chain_indices = [
        i for i in range(1, pose.total_residue() + 1) if pose.chain(i) == target_pose_chain
    ]
    return target_chain_indices[0], target_chain_indices[n_added - 1]


def add_cterminus_residues(pose, chain_id: str, seq: str) -> tuple[int, int]:
    """Append ``seq`` to the C-terminus of ``chain_id`` in-place; rebuild the PDBInfo; return the
    (start, end) pose indices of the added residues.

    Numbering is PRESERVED (added residues continue past the last resnum), so a caller can carry the
    parent's per-residue annotations forward. Pure pose surgery — no relax, no write."""
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("append sequence must be non-empty")
    info = pose.pdb_info()
    original_chain_map = {}  # pose_idx -> (chain, resnum, icode), captured before mutation
    if info:
        for i in range(1, pose.total_residue() + 1):
            original_chain_map[i] = (info.chain(i), info.number(i), info.icode(i))

    chain_indices = chain_indices_from_pdbinfo(pose, chain_id)
    original_chain_end = chain_indices[-1]
    original_chain_last_resnum = original_chain_map[original_chain_end][1]

    rts = pose.residue_type_set_for_pose()
    current_end = original_chain_end  # track the growing C-terminus
    for aa in seq:
        if aa not in ONE_TO_THREE:
            raise ValueError(f"Unsupported residue code '{aa}'")
        restype = rts.name_map(ONE_TO_THREE[aa])
        new_res = rosetta.core.conformation.ResidueFactory.create_residue(restype)
        pose.append_polymer_residue_after_seqpos(new_res, current_end, True)
        current_end += 1

    # <=original_chain_end: unchanged; +1..+n: the new residues (numbering continues from the last
    # original resnum); rest: post-target chains shifted by +n, restored from the original map.
    n_added = len(seq)
    info = pose.pdb_info()
    if info is not None:
        for i in range(1, pose.total_residue() + 1):
            if i <= original_chain_end:
                orig_chain, orig_resnum, orig_icode = original_chain_map[i]
                info.set_resinfo(i, orig_chain, orig_resnum, orig_icode)
            elif i <= original_chain_end + n_added:
                info.set_resinfo(
                    i, chain_id, original_chain_last_resnum + (i - original_chain_end), " "
                )
            else:
                orig_idx = i - n_added
                orig_chain, orig_resnum, orig_icode = original_chain_map[orig_idx]
                info.set_resinfo(i, orig_chain, orig_resnum, orig_icode)

    return original_chain_end + 1, original_chain_end + n_added


def write_pose_preserving_pdbinfo(pose, out_pdb: str | Path) -> None:
    """Write ``pose`` to ``out_pdb`` with its PDBInfo chain IDs / residue numbers preserved.

    ``pose.dump_pdb`` ignores PDBInfo and relabels chains A,B,C… sequentially, so the dump is re-read
    and each atom's chain/resnum rewritten from the PDBInfo (mapped by pose index). Needed after a
    terminus edit, which leaves non-sequential chain/residue numbering on the PDBInfo."""
    temp_pdb = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".pdb")
    temp_pdb.close()
    pose.dump_pdb(temp_pdb.name)
    info = pose.pdb_info()
    # (dump_chain, dump_resnum) -> pose index, in first-seen order (pose order)
    dump_residues = []
    with open(temp_pdb.name, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                key = (line[21], int(line[22:26].strip()))
                if key not in [r[0] for r in dump_residues]:
                    dump_residues.append((key, len(dump_residues) + 1))
    residue_map = {
        key: (info.chain(pose_idx), info.number(pose_idx)) for key, pose_idx in dump_residues
    }
    with open(temp_pdb.name, "r") as f_in:
        lines = f_in.readlines()
    with open(str(out_pdb), "w") as f_out:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                key = (line[21], int(line[22:26].strip()))
                correct_chain, correct_resnum = residue_map.get(key, key)
                f_out.write(line[:21] + correct_chain + f"{correct_resnum:4d}" + line[26:])
            else:
                f_out.write(line)
    os.unlink(temp_pdb.name)
