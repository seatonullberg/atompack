import numpy as np

from _cell import cell_contains


def test_cell_contains_inside():
    cell = np.identity(3)
    position = np.array([0.5, 0.5, 0.5])
    tolerance = 1.0e-6
    assert cell_contains(cell, position, tolerance)


def test_cell_contains_surface():
    cell = np.identity(3)
    position = np.array([1.0, 1.0, 0.5])
    tolerance = 1.0e-6
    assert cell_contains(cell, position, tolerance)


def test_cell_contains_outside():
    cell = np.identity(3)
    position = np.array([1.0, 1.0, 1.5])
    tolerance = 1.0e-6
    assert not cell_contains(cell, position, tolerance)
