"""Lightweight tests for the pure-Python logic in rosetta_pipelines.

The modules import PyRosetta at import time, so the whole suite is skipped when PyRosetta is not
installed. Nothing here runs an actual Rosetta protocol (that needs structure inputs and is slow);
these cover the argument plumbing and param packing that is easy to get wrong.
"""

import pytest

pytest.importorskip("pyrosetta")


def test_public_api_exports():
    import rosetta_pipelines as rp

    for name in (
        "relax_structure",
        "thread_and_relax",
        "run_interface_analyzer",
        "RelaxParams",
        "RelaxResult",
        "ensure_init",
    ):
        assert hasattr(rp, name)


def test_relax_params_traj_kwargs():
    from rosetta_pipelines import RelaxParams

    kw = RelaxParams(nstruct=3, weights="ref2015").traj_kwargs()
    # neighbor_dist is consumed when building the movemap, not passed to the trajectory primitive.
    assert "neighbor_dist" not in kw
    assert kw["nstruct"] == 3
    assert kw["weights"] == "ref2015"
    assert set(kw) == {
        "coord_cst_weight",
        "coord_cst_sd",
        "constrain_sidechains",
        "relax_repeats",
        "nstruct",
        "base_seed",
        "extra_rotamers",
        "min_type",
        "cartesian",
        "weights",
    }


@pytest.mark.parametrize(
    "arg,expected",
    [
        ("A", ["A"]),
        ("A,B", ["A", "B"]),
        ("A B", ["A", "B"]),
        ("AB", ["A", "B"]),
    ],
)
def test_parse_chain_string(arg, expected):
    from rosetta_pipelines.interface import _parse_chain_string

    assert _parse_chain_string(arg) == expected
