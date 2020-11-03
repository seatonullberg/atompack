import copy

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
    # copy here to prevent exponential growth
    return transform.apply(copy.deepcopy(crystal))


def bench_crystal_from_json(json_data):
    return Crystal.from_json(json_data)


def bench_crystal_to_json(crystal):
    return crystal.to_json()


def bench_crystal_deepcopy(crystal):
    return copy.deepcopy(crystal)


############################
#    Benchmark Wrappers    #
############################


def test_crystal_supercell(benchmark):
    unit_cell = get_cubic_unit_cell()
    crystal = Crystal(unit_cell)
    transform = Transform().supercell((5, 5, 5))
    res = benchmark.pedantic(
        bench_crystal_supercell,
        (crystal, transform),
        rounds=10,
        iterations=100,
    )
    assert len(res.atoms) == 500


def test_crystal_to_json(benchmark):
    unit_cell = get_cubic_unit_cell()
    crystal = Crystal(unit_cell)
    res = benchmark.pedantic(
        bench_crystal_to_json,
        (crystal,),
        rounds=100,
        iterations=1000,
    )
    new_crystal = Crystal.from_json(res)
    assert len(new_crystal.atoms) == len(crystal.atoms)


def test_crystal_from_json(benchmark):
    unit_cell = get_cubic_unit_cell()
    crystal = Crystal(unit_cell)
    json_data = crystal.to_json()
    res = benchmark.pedantic(
        bench_crystal_from_json,
        (json_data,),
        rounds=100,
        iterations=1000,
    )
    assert len(crystal.atoms) == len(res.atoms)


def test_crystal_deepcopy(benchmark):
    unit_cell = get_cubic_unit_cell()
    crystal = Crystal(unit_cell)
    res = benchmark.pedantic(
        bench_crystal_deepcopy,
        (crystal,),
        rounds=100,
        iterations=1000,
    )
    assert len(res.atoms) == 4
