import numpy as np

from atompack.crystal import Basis, LatticeParameters, UnitCell, Crystal
from atompack.symmetry import Spacegroup

###############
#    Setup    #
###############


def get_cubic_unit_cell():
    basis = Basis.primitive("X")
    lattparams = LatticeParameters.cubic(10)
    spg = Spacegroup(225)
    return (basis, lattparams, spg)


###################################
#    Benchmark Implementations    #
###################################


def bench_unit_cell(basis, lattice_parameters, spacegroup):
    return UnitCell(basis, lattice_parameters, spacegroup)


def bench_crystal_supercell(unit_cell, extent):
    return Crystal(unit_cell).supercell(extent).finish()


############################
#    Benchmark Wrappers    #
############################


def test_unit_cell_cubic(benchmark):
    args = get_cubic_unit_cell()
    res = benchmark.pedantic(bench_unit_cell, args, rounds=10, iterations=100)
    assert len(res.atoms) == 4


def test_crystal_cubic_supercell_2x2x2(benchmark):
    args = get_cubic_unit_cell()
    args = (UnitCell(*args), (2, 2, 2))
    res = benchmark.pedantic(bench_crystal_supercell, args, rounds=10, iterations=100)


def test_crystal_cubic_supercell_3x3x3(benchmark):
    args = get_cubic_unit_cell()
    args = (UnitCell(*args), (3, 3, 3))
    res = benchmark.pedantic(bench_crystal_supercell, args, rounds=10, iterations=100)


def test_crystal_cubic_supercell_4x4x4(benchmark):
    args = get_cubic_unit_cell()
    args = (UnitCell(*args), (4, 4, 4))
    res = benchmark.pedantic(bench_crystal_supercell, args, rounds=10, iterations=100)
