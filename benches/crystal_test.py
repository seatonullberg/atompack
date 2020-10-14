import numpy as np

from atompack.crystal import Basis, LatticeParameters, UnitCell, Crystal
from atompack.spacegroup import Spacegroup


def get_cubic_unit_cell():
    basis = Basis.primitive("X")
    lattparams = LatticeParameters.cubic(10)
    spg = Spacegroup(225)
    return (basis, lattparams, spg)


def bench_unit_cell(basis, lattice_parameters, spacegroup):
    return UnitCell(basis, lattice_parameters, spacegroup)


def test_unit_cell_cubic(benchmark):
    args = get_cubic_unit_cell()
    res = benchmark.pedantic(bench_unit_cell, args, rounds=10, iterations=100)
    assert len(res.atoms) == 4
