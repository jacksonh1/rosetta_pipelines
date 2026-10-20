# Terminus-extension demo

A worked example of **terminus extension** with `rosetta_pipelines.extend_terminus`: add a short
residue run to a chain's N- or C-terminus and relax the addition, holding the pre-existing structure
near the template.

## Files

- `W.pdb` — input complex: receptor chain **A** + a 33-residue peptide on chain **C**.
- `extend_terminus.py` — appends a `GGY` tag to the C-terminus of chain C and relaxes the chain,
  writing `out/W_Cext.pdb`.

## Run it

With the environment active (`conda activate rosetta_pipelines`):

```bash
cd examples/terminus
python extend_terminus.py
```

## The core call

```python
from rosetta_pipelines import RelaxParams, extend_terminus

result = extend_terminus(
    "W.pdb",
    chain_id="C",
    add_seq="GGY",          # residues to add
    end="C",                # "C" appends; "N" prepends
    out_pdb="out/W_Cext.pdb",
    params=RelaxParams(nstruct=3),
    relax_scope="chain",    # "chain" (default) or "added"
)
```

## Behavior notes

- **`end="C"`** appends — the chain's existing numbering is preserved and the new residues continue
  past the last residue number.
- **`end="N"`** prepends — the chain is renumbered continuously from 1 (so re-derive any
  residue-number-based annotations from the output).
- **`relax_scope="chain"`** (default) relaxes the whole peptide chain. `relax_scope="added"` relaxes
  only the new residues — but that can't relieve the strain of the new terminal peptide bond, which
  sits on the *neighbouring* existing residue, so `"chain"` is preferred.
- There is **no neighbor shell** (the addition is at a chain end); `params.neighbor_dist` is ignored.
- The output is written with its PDBInfo chain IDs / residue numbers preserved (a terminus edit
  leaves non-sequential numbering that a plain dump would relabel).
- `disulfide_pairs=[...]` (PDB numbering) can be passed to (re)form disulfides before relaxing.
