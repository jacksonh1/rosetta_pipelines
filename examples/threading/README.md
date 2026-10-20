# Threading demo

A worked example of **sequence threading** with `rosetta_pipelines.thread_and_relax`: replace the
sequence on one chain of a complex with a new sequence of the same length, then relax.

## Files

- `W.pdb` — input complex: receptor chain **A** + a 33-residue peptide on chain **C**.
- `thread_polylysine.py` — threads a poly-lysine sequence (33 lysines) onto chain C and relaxes
  chain C plus a 6 Å shell, writing `out/W_polyK.pdb`.

## Run it

With the environment active (`conda activate rosetta_pipelines`):

```bash
cd examples/threading
python thread_polylysine.py
```

Expected output (scores will vary slightly by build):

```
Threaded 33-residue poly-lysine onto chain C of W.pdb
  positions mutated  : 33
  best score (REU)   : ...
  trajectory spread  : ...
  Cα drift (mobile)  : ...
  wrote relaxed PDB  : out/W_polyK.pdb
```

## The core call

```python
from rosetta_pipelines import RelaxParams, thread_and_relax

result = thread_and_relax(
    "W.pdb",
    chain_id="C",
    new_seq="K" * 33,          # same length as chain C
    out_pdb="out/W_polyK.pdb",
    params=RelaxParams(nstruct=3, neighbor_dist=6.0),
    thread_scope="all",         # rebuild every position
    relax_scope="neighborhood", # relax chain C + 6 Å shell, freeze the rest
)
```

The one hard rule: `new_seq` must be exactly as long as the target chain (chain C is 33 residues, so
the poly-lysine sequence is 33 lysines). To thread a *designed* sequence instead of poly-lysine, just
pass it as `new_seq`.
```
