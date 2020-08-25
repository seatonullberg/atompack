import numpy as np

from _cell import cell_contains


def bench_cell_contains(cell, position, tolerance):
    return cell_contains(cell, position, tolerance)


def test_bench_cell_contains_inside(benchmark):
    cell = np.identity(3)
    position = np.array([0.5, 0.5, 0.5])
    tolerance = 1.0e-6
    args = (cell, position, tolerance)
    benchmark.pedantic(bench_cell_contains, args=args, rounds=10, iterations=100)


def test_bench_cell_contains_surface(benchmark):
    cell = np.identity(3)
    position = np.array([1.0, 1.0, 0.5])
    tolerance = 1.0e-6
    args = (cell, position, tolerance)
    benchmark.pedantic(bench_cell_contains, args=args, rounds=10, iterations=100)


def test_bench_cell_contains_outside(benchmark):
    cell = np.identity(3)
    position = np.array([1.0, 1.0, 1.5])
    tolerance = 1.0e-6
    args = (cell, position, tolerance)
    benchmark.pedantic(bench_cell_contains, args=args, rounds=10, iterations=100)
