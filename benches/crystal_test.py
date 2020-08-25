import numpy as np

from atompack.atom import Atom
from atompack.crystal import Crystal, UnitCell


def get_cubic_unit_cell():
    a, b, c = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    atoms = [
        Atom(np.array([0, 0, 0])),
        Atom(np.array([1.425, 1.425, 1.425])),
    ]
    return UnitCell(atoms, a, b, c, alpha, beta, gamma)


def bench_crystal(unit_cell, scale, orientation, rotation, tolerance):
    return Crystal(unit_cell, scale, orientation, rotation, tolerance)


def test_bench_crystal_cubic_100_orientation(benchmark):
    unit_cell = get_cubic_unit_cell()
    scale = (1, 1, 1)
    orientation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tolerance = 1.0e-6
    args = (unit_cell, scale, orientation, rotation, tolerance)
    benchmark.pedantic(bench_crystal, args=args, rounds=10, iterations=100)


def test_bench_crystal_cubic_110_orientation(benchmark):
    unit_cell = get_cubic_unit_cell()
    scale = (1, 1, 1)
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tolerance = 1.0e-6
    args = (unit_cell, scale, orientation, rotation, tolerance)
    benchmark.pedantic(bench_crystal, args=args, rounds=10, iterations=100)


def test_bench_crystal_cubic_111_orientation(benchmark):
    unit_cell = get_cubic_unit_cell()
    scale = (1, 1, 1)
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tolerance = 1.0e-6
    args = (unit_cell, scale, orientation, rotation, tolerance)
    benchmark.pedantic(bench_crystal, args=args, rounds=10, iterations=100)


def test_bench_crystal_cubic_2x2x2_supercell(benchmark):
    unit_cell = get_cubic_unit_cell()
    scale = (2, 2, 2)
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tolerance = 1.0e-6
    args = (unit_cell, scale, orientation, rotation, tolerance)
    benchmark.pedantic(bench_crystal, args=args, rounds=10, iterations=100)


def test_bench_crystal_cubic_3x3x3_supercell(benchmark):
    unit_cell = get_cubic_unit_cell()
    scale = (3, 3, 3)
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tolerance = 1.0e-6
    args = (unit_cell, scale, orientation, rotation, tolerance)
    benchmark.pedantic(bench_crystal, args=args, rounds=10, iterations=100)
