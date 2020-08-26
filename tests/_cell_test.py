import numpy as np

from _cell import cell_contains, cell_enforce


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


def test_cell_enforce_inside():
    cell = np.identity(3)
    position = np.array([0.5, 0.5, 0.5])
    tolerance = 1.0e-6
    cell_enforce(cell, position, tolerance)
    target_position = np.array([0.5, 0.5, 0.5])
    assert np.array_equal(position, target_position)


def test_cell_enforce_surface():
    cell = np.identity(3)
    position = np.array([1.0, 1.0, 0.5])
    tolerance = 1.0e-6
    cell_enforce(cell, position, tolerance)
    # surfaces are reduced
    target_position = np.array([0.0, 0.0, 0.5])
    assert np.array_equal(position, target_position)


def test_cell_enforce_above():
    cell = np.identity(3)
    position = np.array([0.5, 0.5, 1.75])
    tolerance = 1.0e-6
    cell_enforce(cell, position, tolerance)
    target_position = np.array([0.5, 0.5, 0.75])
    assert np.array_equal(position, target_position)


def test_cell_enforce_below():
    cell = np.identity(3)
    position = np.array([0.5, 0.5, -1.75])
    tolerance = 1.0e-6
    cell_enforce(cell, position, tolerance)
    target_position = np.array([0.5, 0.5, 0.25])
    assert np.array_equal(position, target_position)
