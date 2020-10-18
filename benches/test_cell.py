import numpy as np

from atompack.cell import is_point_in_cell, wrap_point_into_cell

###################################
#    Benchmark Implementations    #
###################################


def bench_is_point_in_cell(cell, point):
    return is_point_in_cell(cell, point)


def bench_wrap_point_into_cell(cell, point):
    return wrap_point_into_cell(cell, point)


############################
#    Benchmark Wrappers    #
############################


def test_is_point_in_cell_inside(benchmark):
    cell = np.identity(3)
    point = np.array([0.5, 0.5, 0.5])
    args = (cell, point)
    res = benchmark.pedantic(bench_is_point_in_cell, args, rounds=10, iterations=100)
    assert res


def test_is_point_in_cell_outside(benchmark):
    cell = np.identity(3)
    point = np.array([1.5, 1.5, 1.5])
    args = (cell, point)
    res = benchmark.pedantic(bench_is_point_in_cell, args, rounds=10, iterations=100)
    assert not res


def test_wrap_point_into_cell_inside(benchmark):
    cell = np.identity(3)
    point = np.array([0.5, 0.5, 0.5])
    target = np.array([0.5, 0.5, 0.5])
    args = (cell, point)
    res = benchmark.pedantic(bench_wrap_point_into_cell, args, rounds=10, iterations=100)
    assert np.array_equal(point, target)


def test_wrap_point_into_cell_outside(benchmark):
    cell = np.identity(3)
    point = np.array([1.5, 1.5, 1.5])
    target = np.array([0.5, 0.5, 0.5])
    args = (cell, point)
    res = benchmark.pedantic(bench_wrap_point_into_cell, args, rounds=10, iterations=100)
    assert np.array_equal(point, target)
