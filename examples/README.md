# Examples

Each subfolder is a self-contained, runnable demo of one operation (with its own input PDB and a
README). Activate the environment first (`conda activate rosetta_pipelines`).

| Folder | Operation | Script |
|--------|-----------|--------|
| [`relaxation/`](relaxation/) | Structure relaxation (`relax_structure`) | `relax_demo.py` |
| [`threading/`](threading/) | Sequence threading (`thread_and_relax`) | `thread_polylysine.py` |
| [`terminus/`](terminus/) | Terminus extension (`extend_terminus`) | `extend_terminus.py` |
| [`interface/`](interface/) | Interface analysis (`run_interface_analyzer`) | `interface_demo.py` |
| [`parallel/`](parallel/) | Batch over many structures (`run_pool`) | `relax_all.py` |

Quick run:

```bash
conda activate rosetta_pipelines
cd examples/relaxation && python relax_demo.py
cd ../threading     && python thread_polylysine.py
cd ../interface     && python interface_demo.py
```

Each demo writes its outputs to an `out/` folder (gitignored). The relax-based demos run real
FastRelax and take a couple of minutes; the interface demo is faster.
