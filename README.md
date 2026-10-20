# rosetta_pipelines

Small, self-contained PyRosetta pipelines for working with protein structures. Three operations, each
a plain Python function that takes a PDB file and returns a typed result:

| Operation | Function | What it does |
|-----------|----------|--------------|
| **Structure relaxation** | `relax_structure` | Coordinate-restrained FastRelax of a structure, optionally scoped to one chain plus a neighbor shell. |
| **Sequence threading** | `thread_and_relax` | Thread a new sequence onto one chain (backbone preserved, sidechains rebuilt), then relax. |
| **Terminus extension** | `extend_terminus` | Add residues to a chain's N- or C-terminus, then relax the addition. |
| **Interface analysis** | `run_interface_analyzer` | Pack-separated binding ΔG and interface descriptors for a two-partner complex. |

A process-pool helper (`rosetta_pipelines.pool.run_pool`) runs any of these across many structures in
parallel, each in its own worker process.

The runtime code depends only on the Python standard library and **PyRosetta**.

---

## Installation

PyRosetta is not on the standard PyPI index and requires a license (free for academic use — see
<https://www.pyrosetta.org/home/licensing-pyrosetta>). The quickest path is the provided conda
environment, which installs PyRosetta via RosettaCommons' official `pyrosetta-installer` and then
installs this package:

```bash
make env
conda activate rosetta_pipelines
```

### Installing into an existing environment

If you already have an environment with PyRosetta, just install the package:

```bash
pip install .
```

If that environment does **not** yet have PyRosetta:

```bash
make install-pyrosetta   # pip install pyrosetta-installer + run it
pip install .
```

`pip install .` itself pulls no dependencies (the code is stdlib + PyRosetta only), so it will not try
to resolve PyRosetta from PyPI and fail — PyRosetta is installed separately by the steps above.

### Verify

```bash
python -c "import rosetta_pipelines; print(rosetta_pipelines.__version__)"
python -c "import pyrosetta; print('pyrosetta ok')"
```

---

## Usage

PyRosetta is initialized lazily on the first pose operation, so you never need to call
`pyrosetta.init()` yourself.

### Relax a structure

```python
from rosetta_pipelines import relax_structure, RelaxParams

result = relax_structure(
    "input.pdb",
    "relaxed.pdb",
    params=RelaxParams(nstruct=5, relax_repeats=5),
    focus_chain=None,          # None = relax the whole structure
)
print(result.best_score, result.ca_rmsd_drift)
```

Scope the relax to one chain plus its neighbor shell (everything else frozen):

```python
result = relax_structure(
    "complex.pdb", "relaxed.pdb",
    params=RelaxParams(neighbor_dist=6.0),
    focus_chain="B",           # relax chain B + residues within 6 Å, freeze the rest
)
```

A complete, runnable relaxation demo lives in [`examples/relaxation/`](examples/relaxation/).

Form disulfides before relaxing (PDB numbering, `((chain, resnum), (chain, resnum))`):

```python
result = relax_structure(
    "peptide.pdb", "relaxed.pdb",
    disulfide_pairs=[(("B", 3), ("B", 20))],
)
```

### Thread a sequence onto a chain and relax

```python
from rosetta_pipelines import thread_and_relax, RelaxParams

result = thread_and_relax(
    "template.pdb",
    chain_id="B",
    new_seq="ACDEFGHIKLMNPQRSTVWY",   # length must equal chain B's residue count
    out_pdb="threaded_relaxed.pdb",
    params=RelaxParams(nstruct=5),
    thread_scope="changed",            # mutate only positions that differ from the template
    relax_scope="neighborhood",        # relax chain B + shell; "full" or "mutations" also valid
)
print("positions mutated:", result.n_threaded)
```

A complete, runnable threading demo (poly-lysine onto a peptide chain, with a PDB included) lives in
[`examples/threading/`](examples/threading/).

### Extend a terminus (add residues to a chain end)

```python
from rosetta_pipelines import extend_terminus, RelaxParams

result = extend_terminus(
    "peptide.pdb",
    chain_id="C",
    add_seq="GGY",             # residues to add
    end="C",                   # "C" appends (numbering preserved); "N" prepends (chain renumbered from 1)
    out_pdb="extended.pdb",
    params=RelaxParams(nstruct=3),
    relax_scope="chain",       # "chain" (default) relaxes the whole chain; "added" only the new residues
)
```

There is no neighbor shell (the addition sits at a chain end), and the output keeps its PDBInfo chain
IDs / residue numbers. A complete, runnable demo lives in [`examples/terminus/`](examples/terminus/).

### Analyze an interface

```python
from rosetta_pipelines import run_interface_analyzer

res = run_interface_analyzer(
    "complex.pdb",
    receptor_chains=["A"],
    peptide_chains=["B"],
    pack_separated=True,     # repack sidechains in the unbound state (standard binding ΔG)
)
print(res.dG_separated, res.dSASA, res.sc_value, res.packstat)
```

Or from the command line (installed as a console script):

```bash
rosetta-interface-analyzer --pdb complex.pdb --receptor A --peptide B
# equivalently:
python -m rosetta_pipelines.interface --pdb complex.pdb --receptor A --peptide B
```

A complete, runnable interface demo lives in [`examples/interface/`](examples/interface/).

### Run across many structures in parallel

See [`examples/parallel/`](examples/parallel/) for a complete, runnable driver. The worker function
must be defined at module level and the launch guarded by `if __name__ == "__main__":` (both are
`spawn` requirements).

```python
from rosetta_pipelines.pool import run_pool

for out in run_pool(pdb_paths, relax_one, max_workers=8):
    print("done:", out)
```

---

## Notes on behavior

These defaults encode some non-obvious PyRosetta details:

- **Template-preserving relax.** Relaxation applies explicit backbone coordinate constraints to the
  input (weight 0.5, σ 0.5 Å), so the backbone stays near the template instead of drifting. Set
  `coord_cst_weight=0.0` for a free relax.
- **Reported scores are constraint-free.** The restraint drives the relax but is excluded from the
  reported `best_score`/`mean_score`/etc.; the restraint magnitude is reported separately as
  `coord_cst_energy`. The returned pose has its constraints cleared.
- **`nstruct` independent trajectories** are run and the lowest-scoring is kept (one FastRelax
  trajectory finds only one local minimum). The run is seeded once from `base_seed`, so it is
  reproducible run-to-run while the trajectories stay independent.
- **Disulfides are pinned** at their ideal SG–SG length (2.03 Å) during relax, so a soft disulfide
  term cannot let a bond stretch under the backbone tether.
- **Cartesian relax** requires a cart-enabled weights set (e.g. `ref2015_cart`) and an `lbfgs_*`
  minimizer; the code fails loud otherwise rather than silently producing distorted bond geometry.
- **Zero-occupancy atoms are kept** on load by default (Rosetta's own default drops them silently,
  which can delete an unmodelled internal residue and introduce a chain break).
- **Parallelism uses `spawn`, never `fork`** — a forked worker would inherit PyRosetta's RNG state and
  produce correlated trajectories.

---

## Development

```bash
make test       # run the test suite
make format     # black, line length 99
```

---

## License

This project's code is released under the [MIT License](LICENSE).

**PyRosetta is not included or redistributed here** — it is a separate, required dependency with its
own license. PyRosetta is distributed by RosettaCommons and is free for academic/non-commercial use;
commercial use requires a license from the University of Washington. You must obtain PyRosetta and
agree to its terms separately: <https://www.pyrosetta.org/home/licensing-pyrosetta>. The MIT license
on this code grants no rights to PyRosetta.
