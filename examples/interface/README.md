# Interface-analysis demo

A worked example of **interface analysis** with `rosetta_pipelines.run_interface_analyzer`: compute
the pack-separated binding ΔG and interface descriptors for a two-partner complex.

"Pack-separated" means the partners are pulled apart **and** their sidechains are repacked in the
unbound state before rescoring — the physically standard binding ΔG, rather than a rigid
frozen-coordinate estimate.

## Files

- `W.pdb` — input complex: receptor chain **A** + peptide chain **C**.
- `interface_demo.py` — analyzes the A vs C interface and prints the metrics.

## Run it

With the environment active (`conda activate rosetta_pipelines`):

```bash
cd examples/interface
python interface_demo.py
```

Or use the installed command-line entry point — no script needed:

```bash
rosetta-interface-analyzer --pdb W.pdb --receptor A --peptide C
```

## The core call

```python
from rosetta_pipelines import run_interface_analyzer

result = run_interface_analyzer(
    "W.pdb",
    receptor_chains=["A"],   # side 1
    peptide_chains=["C"],    # side 2
    pack_separated=True,     # repack unbound state -> standard binding ΔG
)
print(result.dG_separated, result.dSASA, result.sc_value, result.packstat)
```

## What the result tells you

`run_interface_analyzer` returns an `InterfaceAnalyzerResult` (energies in REU, areas in Å²):

| Field | Meaning |
|-------|---------|
| `dG_separated` | Binding energy E_complex − E_unbound (lower = tighter binding). |
| `dG_dSASA_ratio` | Binding energy per 100 Å² buried — a size-normalized ranking metric. |
| `crossterm_interface_energy` | Rigid cross-interface energy (no unbound repack); the gap vs `dG_separated` is the sidechain-reorganization cost. |
| `dSASA` | Total buried interface surface area. |
| `sc_value` | Lawrence & Colman shape complementarity (0–1; higher = tighter fit). |
| `packstat` | Interface packing quality (0–1; higher = better packed). |
| `delta_unsat_hbonds` | Buried polar atoms that lose an H-bond partner on binding (lower = better). |
| `interface_hbonds` / `total_hbond_E` | Interface H-bond count and total energy. |

Multiple chains per side are allowed, e.g. `receptor_chains=["A", "B"]`, `peptide_chains=["C"]`.
