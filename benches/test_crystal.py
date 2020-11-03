import numpy as np

from atompack.crystal import (Basis, Crystal, LatticeParameters, Transform, UnitCell)
from atompack.symmetry import Spacegroup

###############
#    Setup    #
###############


def get_cubic_unit_cell():
    basis = Basis.primitive("X")
    lattparams = LatticeParameters.cubic(10)
    spg = Spacegroup(225)
    return UnitCell(basis, lattparams, spg)


###################################
#    Benchmark Implementations    #
###################################


def bench_crystal_supercell(crystal, transform):
    return transform.apply(crystal)


############################
#    Benchmark Wrappers    #
############################


def test_crystal_supercell_cubic_3x3x3(benchmark):
    unit_cell = get_cubic_unit_cell()
    crystal = Crystal(unit_cell)
    transform = Transform().supercell((3, 3, 3))
    args = (crystal, transform)
    res = benchmark.pedantic(bench_crystal_supercell, args)
