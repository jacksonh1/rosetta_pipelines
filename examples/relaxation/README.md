# Relaxation demo

A worked example of **structure relaxation** with `rosetta_pipelines.relax_structure`: repack
sidechains and run coordinate-restrained FastRelax so a structure settles into a low-energy
conformation without drifting from the template.

## Files

- `W.pdb` — input complex: receptor chain **A** + a 33-residue peptide on chain **C**.
- `relax_demo.py` — relaxes chain C plus a 6 Å shell (rest frozen) and writes `out/W_relaxed.pdb`.

## Run it

With the environment active (`conda activate rosetta_pipelines`):

```bash
cd examples/relaxation
python relax_demo.py
```

## The core call

```python
from rosetta_pipelines import RelaxParams, relax_structure

result = relax_structure(
    "W.pdb",
    "out/W_relaxed.pdb",
    params=RelaxParams(nstruct=3, neighbor_dist=6.0),
    focus_chain="C",        # relax chain C + 6 Å shell; None = whole structure
)
print(result.best_score, result.ca_rmsd_drift_mobile)
```

## Scope options

- `focus_chain="C"` — relax one chain plus its `neighbor_dist` shell, freeze everything else (fast).
- `focus_chain=None` — relax the whole structure (slower; every residue mobile).
- `RelaxParams(coord_cst_weight=0.0)` — free relax, no template tether (the backbone moves freely).
- `disulfide_pairs=[(("A", 261), ("A", 321)), ...]` — form disulfides (PDB numbering) before relaxing.

## What the result tells you

`relax_structure` returns a `RelaxResult` with, among other fields:

- `best_score` / `mean_score` / `worst_score` — total energy (REU) across the kept trajectories.
- `score_std` — spread across trajectories (small = converged, large = under-sampled).
- `ca_rmsd_drift` / `ca_rmsd_drift_mobile` — how far the structure moved from the input (whole-pose vs
  just the mobile region). A scoped relax freezes most residues, so the mobile value is the honest one.
- `coord_cst_energy` — magnitude of the restraint holding the pose to the template (reported
  separately; it is *excluded* from the scores above).
