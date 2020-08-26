import numpy as np

from atompack._cell import cell_contains, cell_enforce


def bench_cell_contains(cell, position, tolerance):
    return cell_contains(cell, position, tolerance)


def bench_cell_enforce(cell, position, tolerance):
    return cell_enforce(cell, position, tolerance)


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


def test_bench_cell_enforce_inside(benchmark):
    cell = np.identity(3)
    position = np.array([0.75, 0.75, 0.25])
    tolerance = 1.0e-6
    args = (cell, position, tolerance)
    benchmark.pedantic(bench_cell_enforce, args=args, rounds=10, iterations=100)


def test_bench_cell_enforce_outside(benchmark):
    cell = np.identity(3)
    position = np.array([1.75, 1.75, -0.75])
    tolerance = 1.0e-6
    args = (cell, position, tolerance)
    benchmark.pedantic(bench_cell_enforce, args=args, rounds=10, iterations=100)
