"""High-level relax entry points: relax a structure as-is, or thread a new sequence onto one chain
and relax the result.

Both wrap the restrained FastRelax primitive in :mod:`rosetta_pipelines.fastrelax`, choosing which
residues are mobile and assembling a :class:`RelaxResult` of every measurement the run produced. The
deliverable is the written ``out_pdb`` (the kept, lowest-scoring pose).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import pyrosetta

from rosetta_pipelines import fastrelax as fr
from rosetta_pipelines import pose_tools as pt
from rosetta_pipelines import rosetta_init

# ((chain_a, resnum_a), (chain_b, resnum_b)) in source PDB numbering — a disulfide bond.
DisulfidePair = tuple[tuple[str, int], tuple[str, int]]


@dataclass
class RelaxParams:
    """Engine knobs of a coordinate-restrained FastRelax — the mechanics shared by every relax
    operation, with NO scope choice (which residues are mobile is chosen by the entry point: a plain
    structure relax scopes by ``focus_chain``, threading by ``relax_scope``).

    Rosetta mapping for the passthroughs:

    ==========================  ===============================================================
    RelaxParams                 Rosetta name / how it is applied
    ==========================  ===============================================================
    ``coord_cst_weight``        the ``coordinate_constraint`` score term's weight. NB there is no
                                ``-relax:coord_cst_weight`` CLI flag — the build rejects it.
    ``coord_cst_sd``            ``AtomCoordinateCstMover.cst_sd`` — the ``HarmonicFunc`` width (Å).
    ``constrain_sidechains``    ``AtomCoordinateCstMover.cst_sidechain`` (backbone-only when False).
    ``relax_repeats``           ``FastRelax`` repeats (Rosetta ``-relax:default_repeats``).
    ``nstruct``                 independent trajectories; the lowest-scoring is kept.
    ``min_type``                ``FastRelax.min_type`` (Rosetta ``-run:min_type``).
    ``cartesian``               ``FastRelax.cartesian`` (Rosetta ``-relax:cartesian``).
    ``weights``                 score-function name (Rosetta ``-score:weights``, e.g. ``ref2015``).
    ``base_seed``               RNG base; the relax is seeded once from it (reproducible run-to-run).
    ``extra_rotamers``          bundles ``ex1`` + ``ex2`` extra-rotamer sampling.
    ==========================  ===============================================================
    """

    neighbor_dist: float = 6.0
    """Shell radius (Å) for a scoped (``neighborhood`` / ``mutations``) relax; ignored by a full
    (whole-pose) relax."""
    relax_repeats: int = 5
    """FastRelax fa_rep ramp depth within one trajectory."""
    nstruct: int = 5
    """Independent trajectories; the lowest-scoring is kept."""
    coord_cst_weight: float = 0.5
    """``coordinate_constraint`` scorefxn weight. ``0.0`` disables the restraint (free drift)."""
    coord_cst_sd: float = 0.5
    """Harmonic width (Å) of the coordinate restraint — how far an atom drifts before the penalty
    bites. Set explicitly, so behavior never depends on the PyRosetta build default."""
    constrain_sidechains: bool = False
    """Also tether sidechain heavy atoms to the template (default: backbone/CA only, so sidechains
    stay free to repack)."""
    extra_rotamers: bool = True
    """ex1 + ex2 extra-rotamer sampling during repack."""
    min_type: str = "lbfgs_armijo_nonmonotone"
    """Minimizer type. ``lbfgs_armijo_nonmonotone`` is Rosetta's (and FastRelax's) own default; an
    ``lbfgs_*`` type is also required with ``cartesian``."""
    cartesian: bool = False
    """Cartesian-space relax. Requires a cart-enabled ``weights`` set (e.g. ``ref2015_cart``); the
    entry point fails loud otherwise."""
    base_seed: int = 1
    """RNG base; the relax is seeded once from it at the start, so it is reproducible run-to-run while
    the ``nstruct`` trajectories stay independent."""
    weights: str = "ref2015"
    """Score-function name."""
    ignore_zero_occupancy: bool = False
    """Zero-occupancy atom handling at load. ``False`` (default) KEEPS them — Rosetta's own default
    drops them silently, which can delete an internal residue and introduce a chain break."""

    def traj_kwargs(self) -> dict:
        """Unpack the trajectory knobs into ``run_fastrelax_trajectories`` primitives (``neighbor_dist``
        is consumed upstream when building the movemap, not here)."""
        return dict(
            coord_cst_weight=self.coord_cst_weight,
            coord_cst_sd=self.coord_cst_sd,
            constrain_sidechains=self.constrain_sidechains,
            relax_repeats=self.relax_repeats,
            nstruct=self.nstruct,
            base_seed=self.base_seed,
            extra_rotamers=self.extra_rotamers,
            min_type=self.min_type,
            cartesian=self.cartesian,
            weights=self.weights,
        )


@dataclass
class RelaxResult:
    """Everything one relax run produced. ``out_pdb`` (the kept, lowest-scoring pose written to disk)
    is the deliverable; the rest is the measurement."""

    out_pdb: str
    params: RelaxParams
    trajectory_scores: list[float]
    chosen_index: int
    best_score: float
    mean_score: float
    worst_score: float
    score_std: float
    coord_cst_energy: float
    ca_rmsd_drift: float
    ca_rmsd_drift_mobile: float
    formed_disulfides: list[DisulfidePair] | None
    score_terms: dict[str, float]
    wall_seconds: float
    n_threaded: int | None = None


def _assert_cartesian_consistent(params: RelaxParams) -> None:
    """A cartesian relax silently produces garbage bond geometry unless the weights set carries a
    non-zero ``cart_bonded`` term and the minimizer is an ``lbfgs_*`` type. ``ref2015`` has
    ``cart_bonded`` 0; ``ref2015_cart`` has 0.5 — running cartesian against plain ``ref2015`` leaves
    bond lengths/angles unrestrained with nothing complaining. Checked against the actual score
    function rather than a name whitelist."""
    if not params.cartesian:
        return
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(
        params.weights
    )
    cart_weight = scorefxn.get_weight(pyrosetta.rosetta.core.scoring.cart_bonded)
    if cart_weight == 0.0:
        raise ValueError(
            f"cartesian=True needs a cart_bonded-enabled weights set, but "
            f"{params.weights!r} has cart_bonded=0 (bond geometry would go unrestrained). "
            f"Use e.g. 'ref2015_cart'."
        )
    if not params.min_type.startswith("lbfgs"):
        raise ValueError(f"cartesian=True needs an 'lbfgs_*' minimizer, got {params.min_type!r}.")


def _finish(pose, out_pdb, *, params, metrics, formed, n_threaded, t0) -> RelaxResult:
    """Write the kept pose and assemble the result (shared tail of both entry points)."""
    Path(out_pdb).parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(out_pdb))
    return RelaxResult(
        out_pdb=str(out_pdb),
        params=params,
        formed_disulfides=formed,
        n_threaded=n_threaded,
        wall_seconds=time.time() - t0,
        **metrics,
    )


def relax_structure(
    pdb_path: str | Path,
    out_pdb: str | Path,
    *,
    params: RelaxParams | None = None,
    focus_chain: str | None = None,
    disulfide_pairs: list[DisulfidePair] | None = None,
) -> RelaxResult:
    """Repack + coordinate-restrained FastRelax of a structure (no threading).

    Scope is the one choice: ``focus_chain=None`` relaxes the whole structure; a chain id relaxes
    that chain + its ``params.neighbor_dist`` shell, everything else frozen. Forms ``disulfide_pairs``
    before packing (rotamers built as CYD), runs ``params.nstruct`` trajectories, and writes the
    lowest-scoring to ``out_pdb``."""
    params = params or RelaxParams()
    rosetta_init.ensure_init()
    _assert_cartesian_consistent(params)
    t0 = time.time()
    pose = pt.load_pose(pdb_path, ignore_zero_occupancy=params.ignore_zero_occupancy)
    formed = None
    if disulfide_pairs is not None:
        formed = pt.add_disulfides_from_chain_positions(pose, disulfide_pairs)
    sel = pyrosetta.rosetta.core.select.residue_selector
    focus = None if focus_chain is None else sel.ChainSelector(focus_chain)
    mmf, prevent = fr.movemap_for_mobile(focus, params.neighbor_dist)
    best, metrics = fr.run_fastrelax_trajectories(
        pose, movemap_factory=mmf, prevent_selector=prevent, **params.traj_kwargs()
    )
    return _finish(
        best, out_pdb, params=params, metrics=metrics, formed=formed, n_threaded=None, t0=t0
    )


def thread_and_relax(
    pdb_path: str | Path,
    chain_id: str,
    new_seq: str,
    out_pdb: str | Path,
    *,
    params: RelaxParams | None = None,
    thread_scope: str = "changed",
    relax_scope: str = "neighborhood",
    disulfide_pairs: list[DisulfidePair] | None = None,
) -> RelaxResult:
    """Thread ``new_seq`` onto ``chain_id`` then relax the result.

    In-pose (no disk round-trip): load once, break the target chain's auto-detected disulfides,
    thread the sequence (``thread_scope``: ``"changed"`` mutates only differing positions, keeping the
    rest's input rotamers as the relax seed; ``"all"`` rebuilds every position), (re)form
    ``disulfide_pairs`` before packing, scope the mobile region per ``relax_scope`` (``"full"`` → whole
    pose; ``"neighborhood"`` → ``chain_id`` + shell; ``"mutations"`` → the just-threaded positions +
    shell), run the trajectories, and write the lowest-scoring to ``out_pdb``."""
    params = params or RelaxParams()
    rosetta_init.ensure_init()
    _assert_cartesian_consistent(params)
    t0 = time.time()
    pose = pt.load_pose(pdb_path, ignore_zero_occupancy=params.ignore_zero_occupancy)
    chain_indices = pt.chain_indices_from_pdbinfo(pose, chain_id)
    pt.break_chain_disulfides(pose, chain_indices)
    changed = pt.thread_sequence_onto_chain(pose, chain_id, new_seq, thread_scope)
    formed = None
    if disulfide_pairs is not None:
        formed = pt.add_disulfides_from_chain_positions(pose, disulfide_pairs)
    sel = pyrosetta.rosetta.core.select.residue_selector
    if relax_scope == "full":
        focus = None
    elif relax_scope == "neighborhood":
        focus = sel.ChainSelector(chain_id)
    elif relax_scope == "mutations":
        if not changed:
            raise ValueError(
                "relax_scope='mutations' but no positions changed (nothing to relax around); "
                "use 'neighborhood' or 'full', or thread_scope='all'."
            )
        idx = sel.ResidueIndexSelector()
        for i in changed:
            idx.append_index(i)
        focus = idx
    else:
        raise ValueError(
            f"relax_scope must be 'full', 'neighborhood', or 'mutations', got {relax_scope!r}"
        )
    mmf, prevent = fr.movemap_for_mobile(focus, params.neighbor_dist)
    best, metrics = fr.run_fastrelax_trajectories(
        pose, movemap_factory=mmf, prevent_selector=prevent, **params.traj_kwargs()
    )
    return _finish(
        best, out_pdb, params=params, metrics=metrics, formed=formed, n_threaded=len(changed), t0=t0
    )
