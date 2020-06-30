import numpy as np

from atompack.bindings import (cell_contains, load_libatompack, metric_tensor, nearest_neighbor)


def test_cell_contains_inside():
    lib = load_libatompack()
    position = np.array([0.5, 0.5, 0.5])
    cell = np.identity(3)
    tolerance = 1.0e-6
    res = cell_contains(lib, cell, position, tolerance)
    assert res == True


def test_cell_contains_surface():
    lib = load_libatompack()
    position = np.array([0.5, 0.5, 1.0])
    cell = np.identity(3)
    tolerance = 1.0e-6
    res = cell_contains(lib, cell, position, tolerance)
    assert res == True


def test_cell_contains_outside():
    lib = load_libatompack()
    position = np.array([1.5, 1.5, 1.5])
    cell = np.identity(3)
    tolerance = 1.0e-6
    res = cell_contains(lib, cell, position, tolerance)
    assert res == False


def test_metric_tensor():
    lib = load_libatompack()
    a, b, c = 2.0, 2.0, 2.0
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    target_res = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ])
    res = metric_tensor(lib, a, b, c, alpha, beta, gamma)
    assert np.allclose(res, target_res)


def test_nearest_neighbor():
    lib = load_libatompack()
    position = np.array([0, 0, 0])
    positions = np.array([
        [0.2, 0.2, 0.2],
        [0.5, 0.5, 0.5],
        [0.9, 0.9, 0.9],
    ])
    cell = np.identity(3)
    periodicity = np.array([1, 1, 1])
    tolerance = 1.0e-6
    index, distance = nearest_neighbor(lib, position, positions, cell, periodicity, tolerance)
    assert index == 2
    assert distance - 0.173205081 < tolerance
