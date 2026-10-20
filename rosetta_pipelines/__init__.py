"""rosetta_pipelines — small, self-contained PyRosetta pipelines.

Three operations on protein structures, each a plain function that takes a PDB and returns a typed
result:

* :func:`relax_structure` — coordinate-restrained FastRelax of a structure (optionally scoped to one
  chain + a neighbor shell).
* :func:`thread_and_relax` — thread a new sequence onto one chain, then relax the result.
* :func:`run_interface_analyzer` — pack-separated binding ΔG + interface descriptors for a two-partner
  complex.

PyRosetta is initialized lazily on the first pose op (see :mod:`rosetta_pipelines.rosetta_init`), so
ordinary scripts need no ``pyrosetta.init`` call. To run across many structures in parallel, use
:func:`rosetta_pipelines.pool.run_pool`.
"""

from rosetta_pipelines.interface import InterfaceAnalyzerResult, run_interface_analyzer
from rosetta_pipelines.relax import (
    RelaxParams,
    RelaxResult,
    relax_structure,
    thread_and_relax,
)
from rosetta_pipelines.rosetta_init import ensure_init

__all__ = [
    "RelaxParams",
    "RelaxResult",
    "relax_structure",
    "thread_and_relax",
    "InterfaceAnalyzerResult",
    "run_interface_analyzer",
    "ensure_init",
]

# Single source of truth for the package version (SemVer: MAJOR.MINOR.PATCH). pyproject.toml reads
# this via [tool.setuptools.dynamic], so this is the only line to edit when cutting a release.
__version__ = "0.1.0"
