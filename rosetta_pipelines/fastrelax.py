"""The Rosetta FastRelax primitive and its sub-helpers.

Field-standard, template-preserving relax:
  * coordinate constraints to the input keep the template backbone from drifting;
  * the SAME repack task scopes packing that the movemap scopes minimization (otherwise FastRelax
    silently repacks the whole pose regardless of the movemap);
  * disulfides are formed BEFORE packing (rotamers built as CYD) — callers do this;
  * nstruct independent trajectories are run and the lowest-scoring kept, because one FastRelax
    trajectory only finds one local minimum.

``run_fastrelax_trajectories`` is *the* relaxation primitive: it takes an in-memory pose + a movemap
scope and runs the reseeded restrained FastRelax loop. It does no ``pyrosetta.init`` (the entry
points that build the pose call :func:`rosetta_pipelines.rosetta_init.ensure_init` first).
``movemap_for_mobile`` turns a focus residue selector into the (MoveMapFactory, prevent_selector)
pair the primitive consumes.
"""

from pyrosetta import rosetta

# ref2015 ships with coordinate_constraint weight 0, so a coordinate constraint contributes nothing
# until the weight is set.
#
# 0.5 IS Rosetta's own value, reached by a different route than one might expect. There is no
# `-relax:coord_cst_weight` option (the build rejects it: "OptionKey with id relax:coord_cst_weight
# not found") -- `coord_cst_weight` is a *relax-script* command, a per-stage scale factor. But the
# weight RelaxProtocolBase installs when constraining to start coords is 0.5, measured directly: init
# with `-relax:constrain_relax_to_start_coords`, run FastRelax, and the scorefxn's
# coordinate_constraint weight reads 0.5 afterwards (41 constraints on a 10-mer, N/CA/C/O + OXT,
# 0.28 A CA drift). So 0.5 with sd 0.5 A reproduces the relax application's own tether.
#
# Only `weight / sd**2` matters for a harmonic restraint, so weight and `coord_cst_sd` are
# redundant knobs: tune ONE. We parameterise the weight and leave `coord_cst_sd` at Rosetta's
# `-relax:coord_cst_stdev` default of 0.5 A (set explicitly, not inherited).
_DEFAULT_COORD_CST_WEIGHT = 0.5

# The relax script pinned by name rather than inherited. FastRelax(scorefxn, repeats) with no script
# argument silently takes whatever the Rosetta database's default script is; that script defines the
# entire fa_rep ramp and repack/min schedule, so a change to the database default would change this
# protocol's output with no error and no log line (we run muted). MonomerRelax2019 IS the current
# default in current builds -- pinning it is a no-op today and insurance later.
#
# Do NOT "fix" constraint holding by switching to the database's `no_cst_ramping` script: that one is
# the *legacy* ramp_repack_min ladder, so it would silently downgrade the fa_rep schedule to the
# pre-2019 protocol. Constraint holding is controlled by ramp_down_constraints (see below).
_RELAX_SCRIPT = "MonomerRelax2019"

# FastRelax is NOT deterministic per starting pose (three clones of one packed pose relaxed to
# 413.220 / 412.215 / 413.480 REU with no reseeding — the packer draws from the global RNG, so
# trajectories differ by construction). So the nstruct trajectories are already independent. We seed
# the RNG ONCE, from base_seed, at the start of the relax: the trajectories then draw sequentially
# from one stream (Rosetta's own idiom), which keeps the whole relax reproducible run-to-run without
# falsely assuming per-trajectory determinism.
_RELAX_BASE_SEED = 1

# Disulfide SG-SG restraint held during relax. dslf_fa13 (the ref2015 disulfide term) is a SOFT
# knowledge-based potential, so near a flexible region (e.g. a freshly-prepended terminus) the backbone
# coordinate constraint can out-pull it and leave the bond STRETCHED — Rosetta still calls it a
# disulfide (its detection cutoff is ~2.5 A), but GROMACS pdb2gmx builds a disulfide only within ~10%
# of the ideal (~2.24 A max for amber99sb-ildn) and otherwise SILENTLY drops it, so an over-long bond
# vanishes from the MD topology. A tight harmonic SG-SG restraint to the ideal length keeps every
# declared disulfide at ~2.03 A regardless of backbone flexibility. Added to the RESTRAINED scorefxn
# only (see run_fastrelax_trajectories); the constraint-free ranking/reporting scorefxn never sees it.
_DISULFIDE_SG_IDEAL = 2.03  # A, equilibrium SG-SG of a disulfide
_DISULFIDE_SG_SD = 0.1  # A, harmonic width
_DISULFIDE_CST_WEIGHT = 1.0  # atom_pair_constraint weight on the restrained scorefxn


def _make_scorefxn(coord_cst_weight: float, weights: str = "ref2015", atom_pair_cst_weight: float = 0.0):
    """``weights`` score function with explicit (possibly-zero) coordinate_constraint and
    atom_pair_constraint weights. ``atom_pair_cst_weight`` > 0 enables the disulfide SG-SG restraint
    (see :data:`_DISULFIDE_CST_WEIGHT`); the constraint-free ranking scorefxn leaves both at 0."""
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(weights)
    scorefxn.set_weight(rosetta.core.scoring.coordinate_constraint, coord_cst_weight)
    scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, atom_pair_cst_weight)
    return scorefxn


def _disulfide_sg_pairs(pose) -> list[tuple[int, int]]:
    """The (i, j) pose-index pairs of every declared disulfide in ``pose`` (each once, i < j).

    Reads the pose's own disulfide bonds (whatever load-time auto-detection or an explicit
    ``add_disulfides_*`` call declared), so it needs no external knowledge."""
    pairs: set[tuple[int, int]] = set()
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        rt = res.type()
        if not rt.is_disulfide_bonded() or not res.has("SG"):
            continue
        sg = res.atom_index("SG")
        for c in range(1, rt.n_possible_residue_connections() + 1):
            if rt.residue_connect_atom_index(c) == sg:
                j = res.residue_connection_partner(c)
                if j != 0:
                    pairs.add((min(i, j), max(i, j)))
    return sorted(pairs)


def _repack_task_factory(extra_rotamers: bool, prevent_selector=None):
    """Repack-only (no design) TaskFactory. ``extra_rotamers`` adds ex1/ex2 + the input rotamer for
    better sidechain sampling; ``prevent_selector`` (a residue selector) freezes everything it picks.
    Extra rotamers are set via task operations, not init flags, so they apply regardless of how/when
    ``pyrosetta.init`` was called elsewhere."""
    op = rosetta.core.pack.task.operation
    tf = rosetta.core.pack.task.TaskFactory()
    tf.push_back(op.InitializeFromCommandline())
    tf.push_back(op.RestrictToRepacking())
    tf.push_back(op.IncludeCurrent())  # keep the input sidechain as a candidate rotamer
    if extra_rotamers:
        ex = op.ExtraRotamersGeneric()
        ex.ex1(True)
        ex.ex2(True)
        tf.push_back(ex)
    if prevent_selector is not None:
        tf.push_back(op.OperateOnResidueSubset(op.PreventRepackingRLT(), prevent_selector))
    return tf


def _fast_relax_once(
    pose,
    scorefxn,
    task_factory,
    movemap_factory,
    relax_repeats,
    *,
    coord_cst_sd: float | None = None,
    constrain_sidechains: bool = False,
    min_type: str = "lbfgs_armijo_nonmonotone",
    cartesian: bool = False,
):
    """One FastRelax trajectory with explicit backbone coordinate constraints to the pose's current
    (template) coordinates, so the backbone does not drift.

    NB: FastRelax.constrain_relax_to_start_coords() was verified to add NO constraints in some
    PyRosetta builds (pose.constraint_set stays empty and the structure drifts as if unconstrained),
    so constraints are added directly with AtomCoordinateCstMover. ``scorefxn`` must carry a non-zero
    coordinate_constraint weight (see :func:`_make_scorefxn`). The task factory scopes repacking to
    match the movemap's minimization scope.

    ``coord_cst_sd`` (Å) sets the harmonic width of the restraint; ``None`` keeps the mover's build
    default (0.5 Å — callers should pass it explicitly rather than inherit it). ``min_type`` /
    ``cartesian`` are the FastRelax minimizer knobs (a cartesian relax needs an ``lbfgs_*`` min_type
    and a cart-enabled scorefxn — the caller's responsibility).

    ``constrain_sidechains`` also tethers sidechain heavy atoms; the default (False) tethers backbone
    heavy atoms only — N, CA, C, O, 4 per residue, plus OXT at the C-terminus. This is a DELIBERATE
    departure from the canonical all-heavy-atom restrained relax (Nivon, Moretti & Baker 2013, PLoS
    ONE 8:e59004), which also constrains sidechains via ``-relax:coord_constrain_sidechains``. That
    published protocol assumes sidechains from an experimental structure; on a threading path the
    sidechains were just rebuilt from ideal internal coordinates, so tethering them would restrain
    them to arbitrary geometry."""
    cst = rosetta.protocols.relax.AtomCoordinateCstMover()  # constrain to current coords
    if coord_cst_sd is not None:
        cst.cst_sd(coord_cst_sd)
    cst.cst_sidechain(constrain_sidechains)
    cst.apply(pose)
    # Pin every declared disulfide's SG-SG at ideal length so a soft dslf_fa13 term cannot let the
    # bond stretch under the backbone tether (which would make pdb2gmx silently drop it downstream —
    # see _DISULFIDE_SG_IDEAL). ``scorefxn`` must carry a non-zero atom_pair_constraint weight.
    func = rosetta.core.scoring.func
    APC = rosetta.core.scoring.constraints.AtomPairConstraint
    AtomID = rosetta.core.id.AtomID
    for i, j in _disulfide_sg_pairs(pose):
        ai = AtomID(pose.residue(i).atom_index("SG"), i)
        aj = AtomID(pose.residue(j).atom_index("SG"), j)
        pose.add_constraint(APC(ai, aj, func.HarmonicFunc(_DISULFIDE_SG_IDEAL, _DISULFIDE_SG_SD)))
    relax = rosetta.protocols.relax.FastRelax(scorefxn, relax_repeats, _RELAX_SCRIPT)
    # Hold the coordinate restraints for the WHOLE trajectory. The relax script carries a
    # coord_cst_weight 1.0 -> 0.5 -> 0.0 ladder that releases them partway through every repeat
    # (measured: 1.18 A CA drift with ramping on vs 0.0001 A with it off).
    #
    # This is NOT redundant with the option default. `-relax:ramp_constraints` does default to false,
    # but RelaxProtocolBase flips ramping ON whenever `-relax:constrain_relax_to_start_coords` is
    # set -- measured: with that one flag at init and nothing else, the option still reads false
    # while `FastRelax.ramp_down_constraints()` reads True, and drift goes 0.28 A -> 0.94 A. We stay
    # unramped only because nothing here sets that flag. Adding it to any init flags, for any reason,
    # would silently un-restrain this relax. Hence the explicit call.
    relax.ramp_down_constraints(False)
    relax.set_task_factory(task_factory)
    relax.set_movemap_factory(movemap_factory)
    relax.cartesian(cartesian)
    relax.min_type(min_type)
    relax.apply(pose)


def _ca_rmsd_nosuper(pose_a, pose_b, indices: list[int] | None = None) -> float:
    """No-superposition Cα RMSD (Å) between two poses of identical residue order/count.

    ``indices`` (1-based pose residue numbers) restricts the average to a residue subset; ``None``
    averages over the whole pose. Averaging over the whole pose understates the mobile region's real
    movement by ~``sqrt(N_mobile / N_total)`` when most of the pose is frozen — pass the mobile set to
    get that region's honest displacement. Fails loud on a residue-count mismatch or if no Cα are
    comparable."""
    n = pose_a.total_residue()
    if n != pose_b.total_residue():
        raise ValueError(f"pose residue counts differ ({n} vs {pose_b.total_residue()})")
    which = range(1, n + 1) if indices is None else indices
    sq = 0.0
    counted = 0
    for i in which:
        ra, rb = pose_a.residue(i), pose_b.residue(i)
        if not (ra.has("CA") and rb.has("CA")):
            continue
        sq += (ra.xyz("CA") - rb.xyz("CA")).length_squared()
        counted += 1
    if counted == 0:
        raise ValueError("no Cα atoms to compare for drift")
    return (sq / counted) ** 0.5


def _weighted_score_terms(scorefxn, pose) -> dict[str, float]:
    """Per-term *weighted* energy of ``pose`` for every nonzero-weight score type."""
    scorefxn(pose)  # ensure the energy cache is current
    energies = pose.energies().total_energies()
    weights = scorefxn.weights()
    return {
        str(st).rsplit(".", 1)[-1]: float(weights[st] * energies[st])
        for st in scorefxn.get_nonzero_weighted_scoretypes()
    }


def run_fastrelax_trajectories(
    pose,
    *,
    movemap_factory,
    prevent_selector=None,
    coord_cst_weight: float = _DEFAULT_COORD_CST_WEIGHT,
    coord_cst_sd: float | None = None,
    constrain_sidechains: bool = False,
    relax_repeats: int = 5,
    nstruct: int = 5,
    base_seed: int = _RELAX_BASE_SEED,
    extra_rotamers: bool = True,
    min_type: str = "lbfgs_armijo_nonmonotone",
    cartesian: bool = False,
    weights: str = "ref2015",
) -> tuple["rosetta.core.pose.Pose", dict]:
    """The relax core: rebuild sidechains (one repack) then run ``nstruct`` restrained FastRelax
    trajectories from that packed start; return ``(best_pose, metrics)``.

    ``movemap_factory`` scopes which DOFs minimize; ``prevent_selector`` (``None`` = repack the whole
    pose) freezes residues outside the mobile region — the same repack task drives both the initial
    pack and every trajectory, so packing and minimization stay scoped together. ``nstruct`` is
    sampling breadth (independent trajectories guarding against one bad local min); ``relax_repeats``
    is the fa_rep ramp depth *within* a trajectory. Disulfides must already be formed on ``pose``.

    Two score functions are used deliberately. The *restrained* one (non-zero
    ``coordinate_constraint`` weight) RUNS the relax, so each trajectory feels the tether holding it
    near the template. The *constraint-free* one (``coordinate_constraint`` = 0) RANKS the
    trajectories, fills the reported scores, and computes ``score_terms`` — so no reported number
    carries the spring energy of the restraint (which scales with atom count and drift and does not
    cancel across structures). ``coord_cst_energy`` keeps the restraint magnitude of the kept pose
    visible. The kept pose is returned with its constraint set cleared, so a later in-process rescore
    with a constraint-enabled scorefxn cannot silently resurrect the restraint energy.

    ``metrics`` carries the per-trajectory scores, the chosen index, best/mean/worst/std, the kept
    pose's coordinate-restraint energy, the no-superposition Cα drift of the kept pose from the input
    template (whole-pose and mobile-only), and its per-term weighted score breakdown."""
    assert nstruct >= 1, "nstruct must be >= 1"
    template = pose.clone()  # Cα reference for drift (backbone is untouched by the repack below)
    # restrained: RUNS the relax (backbone coord cst + disulfide SG-SG cst)
    scorefxn = _make_scorefxn(coord_cst_weight, weights, atom_pair_cst_weight=_DISULFIDE_CST_WEIGHT)
    free_scorefxn = _make_scorefxn(0.0, weights)  # constraint-free: RANKS + reports (both csts 0)
    tf = _repack_task_factory(extra_rotamers, prevent_selector=prevent_selector)

    # Seed ONCE, before the first stochastic op (the initial pack, which also draws from the global
    # RNG). Everything downstream — the pack and all nstruct trajectories — then draws sequentially
    # from one stream, so the whole relax is reproducible run-to-run while the trajectories stay
    # independent of each other. (Seeding only inside the loop would leave the pack keyed to leftover
    # global RNG state and break reproducibility across processes/reruns.)
    rosetta.numeric.random.rg().set_seed(base_seed)
    packer = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)

    best_pose = None
    best_score = None
    chosen = -1
    scores: list[float] = []
    for k in range(nstruct):
        trial = pose.clone()
        _fast_relax_once(
            trial,
            scorefxn,
            tf,
            movemap_factory,
            relax_repeats,
            coord_cst_sd=coord_cst_sd,
            constrain_sidechains=constrain_sidechains,
            min_type=min_type,
            cartesian=cartesian,
        )
        # Rank on the constraint-free score. The coordinate constraints are still attached to `trial`,
        # but free_scorefxn's coordinate_constraint weight is 0, so they contribute nothing.
        score = float(free_scorefxn(trial))
        scores.append(score)
        if best_score is None or score < best_score:
            best_pose, best_score, chosen = trial, score, k

    mean = sum(scores) / len(scores)
    std = (sum((s - mean) ** 2 for s in scores) / len(scores)) ** 0.5
    # Restraint energy of the kept pose, read BEFORE clearing its constraints (the restrained scorefxn
    # scores the coordinate_constraint term; the constraints are still on best_pose here).
    scorefxn(best_pose)
    cst_st = rosetta.core.scoring.coordinate_constraint
    coord_cst_energy = float(scorefxn.weights()[cst_st] * best_pose.energies().total_energies()[cst_st])
    # Return an honestly constraint-free pose.
    best_pose.remove_constraints()
    drift = _ca_rmsd_nosuper(template, best_pose)  # whole-pose: how much relax moved the structure
    # Drift over the residues that were ALLOWED to move (complement of prevent_selector). Whole-pose
    # drift understates the mobile region when most of the pose is frozen; report both.
    if prevent_selector is None:
        drift_mobile = drift  # full scope — every residue is mobile
    else:
        prevented = prevent_selector.apply(best_pose)  # 1-indexed vector<bool>
        mobile_idx = [i for i in range(1, best_pose.total_residue() + 1) if not prevented[i]]
        drift_mobile = _ca_rmsd_nosuper(template, best_pose, indices=mobile_idx)
    metrics = {
        "trajectory_scores": scores,
        "chosen_index": chosen,
        "best_score": best_score,
        "mean_score": mean,
        "worst_score": max(scores),
        "score_std": std,
        "coord_cst_energy": coord_cst_energy,
        "ca_rmsd_drift": drift,
        "ca_rmsd_drift_mobile": drift_mobile,
        "score_terms": _weighted_score_terms(free_scorefxn, best_pose),
    }
    # Report the spread (converged vs under-sampled) and the honest drift — "is it relaxed?".
    print(
        f"[relax] nstruct={nstruct} repeats={relax_repeats} "
        f"score best={min(scores):.1f} mean={mean:.1f} worst={max(scores):.1f} "
        f"cst={coord_cst_energy:.1f} caΔ={drift:.2f}A caΔmob={drift_mobile:.2f}A"
    )
    return best_pose, metrics


def movemap_for_mobile(focus_selector=None, neighbor_dist: float = 0.0, *, bb: bool = True, chi: bool = True):
    """(MoveMapFactory, prevent_selector) for a mobile region — the one mobile-set policy, scope-free.

    ``focus_selector`` (a residue selector) is the set of residues to relax; ``None`` makes the whole
    pose mobile (no ``prevent_selector``). ``neighbor_dist > 0`` expands the focus by that radius into a
    ``NeighborhoodResidueSelector`` (any chain); ``0`` relaxes exactly the focus with no shell. ``bb`` /
    ``chi`` pick backbone and/or sidechain minimization on the mobile set. ``prevent_selector`` (the
    complement of the mobile set) is what the packer must not repack.

    This is the single translator between a caller's scope and the primitive: each caller builds its
    own focus selector (a chain, the just-threaded positions, the just-added residues) and calls this.
    The coordinate constraint (not the movemap) holds the mobile set near the template, so
    ``coord_cst_weight=0`` frees the whole mobile set and ``>0`` tethers it uniformly."""
    sel = rosetta.core.select.residue_selector
    mm = rosetta.core.select.movemap
    mmf = mm.MoveMapFactory()
    if focus_selector is None:
        mmf.all_bb(bb)
        mmf.all_chi(chi)
        return mmf, None
    mobile = focus_selector
    if neighbor_dist and neighbor_dist > 0:
        mobile = sel.NeighborhoodResidueSelector(
            focus_selector, distance=neighbor_dist, include_focus_in_subset=True
        )
    act = mm.move_map_action.mm_enable
    mmf.all_bb(False)
    mmf.all_chi(False)
    if bb:
        mmf.add_bb_action(act, mobile)
    if chi:
        mmf.add_chi_action(act, mobile)
    return mmf, sel.NotResidueSelector(mobile)
