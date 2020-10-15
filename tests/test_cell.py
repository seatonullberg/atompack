import numpy as np

from atompack.cell import is_point_in_cell, wrap_point_into_cell


def test_is_point_in_cell():
    cell = np.identity(3)
    # edge point
    point = np.array([0.0, 0.0, 0.0])
    assert is_point_in_cell(cell, point)
    # surface point
    point = np.array([0.0, 0.5, 0.0])
    assert is_point_in_cell(cell, point)
    # interior point
    point = np.array([0.5, 0.5, 0.5])
    assert is_point_in_cell(cell, point)
    # exterior point (positive)
    point = np.array([1.5, 1.5, 1.5])
    assert not is_point_in_cell(cell, point)
    # exterior point (negative)
    point = np.array([-0.5, -0.5, -0.5])
    assert not is_point_in_cell(cell, point)


def test_wrap_point_into_cell():
    cell = np.identity(3)
    # edge point
    point = np.array([0.0, 0.0, 0.0])
    target = np.array([0.0, 0.0, 0.0])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # surface point
    point = np.array([0.0, 0.5, 0.0])
    target = np.array([0.0, 0.5, 0.0])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # interior point
    point = np.array([0.5, 0.5, 0.5])
    target = np.array([0.5, 0.5, 0.5])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # exterior point (near) (positive)
    point = np.array([1.1, 1.1, 1.1])
    target = np.array([0.1, 0.1, 0.1])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # exterior point (far) (positive)
    point = np.array([5.9, 5.9, 5.9])
    target = np.array([0.9, 0.9, 0.9])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # exterior point (near) (negative)
    point = np.array([-0.1, -0.1, -0.1])
    target = np.array([0.9, 0.9, 0.9])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
    # exterior point (far) (positive)
    point = np.array([-5.9, -5.9, -5.9])
    target = np.array([0.1, 0.1, 0.1])
    res = wrap_point_into_cell(cell, point)
    assert np.allclose(res, target)
