"""Single source of truth for initializing PyRosetta in this process.

``pyrosetta.init`` is process-global and one-shot: the first call sets every flag, and every later
call is a silent no-op. So init must happen exactly once, before the first pose op, carrying the
flags the run actually needs. This module centralizes that in one function:

* :func:`ensure_init` — **lazy default.** Every entry point that creates a pose (or touches the
  Rosetta database, e.g. builds a score function) calls it first. It inits with ``-mute all`` if
  nothing has inited yet, else no-ops. This is what lets ordinary work "just run" with no init call
  in user scripts.

* To use ``-extra_res_fa`` (noncanonical residues: covalent warheads, ligands, NCAAs), those flags
  must ride the **first** init. Since :func:`ensure_init` only inits when nothing has yet, call
  ``ensure_init(extra_res_fa=[...])`` before the first pose op. If PyRosetta is already initialized
  when ``extra_res_fa`` is requested, this fails loud — the flags cannot bind after the fact, and
  silently dropping them would leave the residues out of the residue type set.

There is deliberately no import-time init: importing a module must not lock in flags before a caller
that needs ``extra_res_fa`` gets a chance to set them.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pyrosetta
from pyrosetta import rosetta

# The output-changing init flags this process actually initialized with, recorded on the first init so
# a later ensure_init can fail loud on a *mismatch* (a different flag set that would be silently
# ignored) while allowing a matching re-request to no-op. ``None`` = init has not gone through this
# function (either not yet, or PyRosetta was inited by some other path we cannot vouch for).
_INIT_EXTRA_RES_FA: tuple[str, ...] | None = None


def ensure_init(extra_res_fa: Iterable[str | Path] = (), mute: str = "all") -> None:
    """Initialize PyRosetta once for this process (idempotent via ``was_init_called``).

    ``pyrosetta.init`` binds every flag on the first call and silently no-ops thereafter, so a flag
    requested *after* init would be dropped without error. This guards the only output-changing init
    flag — ``extra_res_fa`` — by recording what the first init used and failing loud if a later call
    requests a **different** set (a matching re-request is a harmless no-op). Other init flags (``mute``)
    do not change output and are not guarded.

    Args:
        extra_res_fa: ``.params`` files for noncanonical residues. Must ride the first init; requesting
            a set that differs from the first init's fails loud.
        mute: PyRosetta mute level ('all', 'protocols', 'none').
    """
    global _INIT_EXTRA_RES_FA
    extra_res_fa = tuple(str(p) for p in extra_res_fa)
    if rosetta.basic.was_init_called():
        active = _INIT_EXTRA_RES_FA or ()
        # An empty request is always a fine no-op ("just ensure inited"). Only a non-empty request
        # that differs from what actually bound is the silent-drop hazard.
        if extra_res_fa and sorted(extra_res_fa) != sorted(active):
            raise RuntimeError(
                "PyRosetta is already initialized with "
                f"extra_res_fa={active!r}, so the requested extra_res_fa={extra_res_fa!r} cannot bind "
                "and would be silently ignored (the noncanonical residues would be absent from the "
                "residue type set). One flag set per process: call ensure_init(extra_res_fa=...) with "
                "the full set before the first pose op."
            )
        return
    flags: list[str] = []
    if mute and mute.lower() != "none":
        flags += ["-mute", mute]
    for params in extra_res_fa:
        flags += ["-extra_res_fa", params]
    pyrosetta.init(" ".join(flags))
    _INIT_EXTRA_RES_FA = extra_res_fa
