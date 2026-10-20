# Parallel batch demo

A worked example of **running across many structures in parallel** with
`rosetta_pipelines.pool.run_pool`: process one structure per worker process, each with its own
PyRosetta init. Use this to relax/score a whole directory of PDBs at once — not to parallelize the
work *within* a single structure.

## Files

- `relax_all.py` — relaxes every `*.pdb` in an input directory across N worker processes, each
  writing its own relaxed output.

## Run it

With the environment active (`conda activate rosetta_pipelines`):

```bash
cd examples/parallel

# Make a small input set (reuse the example complex under a couple of names):
mkdir -p inputs
cp ../interface/W.pdb inputs/a.pdb
cp ../interface/W.pdb inputs/b.pdb

python relax_all.py inputs outputs --workers 2
```

Each input is relaxed independently; results land in `outputs/`, printed as they finish (order is
not preserved).

## The core call

```python
from rosetta_pipelines.pool import run_pool

for result in run_pool(pdb_paths, relax_one, max_workers=8):
    print("done:", result)
```

## Two hard rules (both are `spawn` requirements)

`run_pool` forces the `spawn` start method — a `fork`ed worker would inherit PyRosetta's RNG state and
produce correlated trajectories. Because of `spawn`:

1. **The worker function must be defined at module level** (not a lambda or a closure), so it is
   picklable by reference. See `relax_one` in `relax_all.py`.
2. **The launch must be guarded** by `if __name__ == "__main__":`, because `spawn` re-imports the
   launching module in every worker.

Pass each worker small, picklable inputs (file paths or ids), not large in-memory objects. For
noncanonical residues, pass their `.params` files once via `run_pool(..., extra_res_fa=[...])`.
