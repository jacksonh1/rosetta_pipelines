"""Process-pool helper for running PyRosetta over many independent inputs in parallel.

Use this to process **many PDBs at once** — one structure per task — not to parallelize the
PyRosetta ops *within* a single structure. Each task is expected to be an independent, self-contained
unit of work (load a structure, relax/score it, write its own output).

Three PyRosetta multiprocessing hazards are silent and easy to forget. They are all encoded here so
callers can't hit them by accident:

1. ``pyrosetta.init()`` is process-global and one-shot. Parallelism must therefore be *process*-based,
   never threads: each worker process gets its own single init. This pool inits PyRosetta exactly once
   per worker, in the ``initializer``, before any task runs.

2. ``fork()`` copies an already-initialized PyRosetta — *including its RNG state* — into every child,
   so forked workers produce correlated/identical stochastic trajectories with nothing complaining.
   This pool forces the ``"spawn"`` start method: every worker is a fresh interpreter that inits with
   its own seed. **Do not change this to fork.**

3. ``-extra_res_fa`` params paths are resolved against the worker's cwd at init time, which under
   ``spawn`` need not match the launch dir. This pool resolves them to absolute paths for you
   (``strict=True`` so a missing file fails loud at launch), so callers pass ordinary (possibly
   relative) paths.

Constraints on the caller (both are ``spawn`` requirements, not style choices):

- ``worker_fn`` and every item must be picklable — ``spawn`` pickles ``worker_fn`` by reference and
  each item by value. So ``worker_fn`` must be a **module-level** function (not a lambda or a closure),
  and items should be small (prefer file paths over whole in-memory objects).
- The launching module must be guarded by ``if __name__ == "__main__":`` — ``spawn`` re-imports it in
  every worker, and an unguarded launch would recurse.

Example — relax every PDB in a directory in parallel, each worker writing its own output::

    # driver module, e.g. relax_all.py
    from pathlib import Path
    from rosetta_pipelines.pool import run_pool
    from rosetta_pipelines.relax import relax_structure, RelaxParams

    PARAMS = RelaxParams(nstruct=5)

    def relax_one(pdb_path: str) -> str:          # module-level: picklable by reference
        out = str(Path(pdb_path).with_suffix(".relaxed.pdb"))
        relax_structure(pdb_path, out, params=PARAMS, focus_chain=None)
        return out

    def main() -> None:
        pdbs = [str(p) for p in Path("inputs").glob("*.pdb")]
        for out in run_pool(pdbs, relax_one, max_workers=8):
            print("relaxed ->", out)

    if __name__ == "__main__":                     # required under spawn
        main()

When you need noncanonical residues (e.g. covalent warheads), pass their ``.params`` files once via
``extra_res_fa=[...]``; the pool resolves them to absolute and inits each worker with them.
"""

from __future__ import annotations

import multiprocessing as mp
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def _worker_init(mute: str, extra_res_fa: tuple[str, ...]) -> None:
    """Runs once per worker process, before any task — the first and only PyRosetta init in this
    process, so its flags actually bind. The heavy import is deferred to here (not module top) so this
    module is import-safe in a driver that has not yet spawned the pool.

    ``extra_res_fa`` is already absolute (resolved in :func:`run_pool`)."""
    from rosetta_pipelines import rosetta_init

    rosetta_init.ensure_init(extra_res_fa=extra_res_fa, mute=mute)


def run_pool(
    items: Iterable[T],
    worker_fn: Callable[[T], R],
    *,
    max_workers: int,
    extra_res_fa: Iterable[str | Path] = (),
    mute: str = "all",
    max_tasks_per_child: int = 50,
) -> Iterator[R]:
    """Run ``worker_fn(item)`` for every item across ``max_workers`` spawned worker processes, each
    with its own single PyRosetta init.

    Yields results as they complete (order **not** preserved). A task that raises re-raises here on
    its ``yield`` — fail loud, no swallowed worker exceptions.

    Args:
        items: one input per task; each must be picklable (prefer small paths/ids over big objects).
        worker_fn: a **module-level** function run once per item in a worker.
        max_workers: number of worker processes.
        extra_res_fa: ``.params`` files for noncanonical residues, applied at each worker's init.
            Resolved to absolute here; a missing file fails loud immediately.
        mute: PyRosetta mute level ('all', 'protocols', 'none').
        max_tasks_per_child: recycle a worker after this many tasks to bound PyRosetta's memory growth
            over a long run.
    """
    resolved = tuple(str(Path(p).resolve(strict=True)) for p in extra_res_fa)

    # spawn (not fork): each worker inits PyRosetta fresh with its own RNG seed. See module docstring.
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_worker_init,
        initargs=(mute, resolved),
        max_tasks_per_child=max_tasks_per_child,
    ) as ex:
        futures = [ex.submit(worker_fn, item) for item in items]
        for fut in as_completed(futures):
            yield fut.result()
