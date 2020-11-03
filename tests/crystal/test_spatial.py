import numpy as np
import pytest

from atompack.crystal.spatial import MillerIndex, Orientation, Plane

###########################
#    MillerIndex Tests    #
###########################


@pytest.mark.parametrize("test_input,expectation", [
    (np.array([1 / 2, 2 / 3, 1]), MillerIndex((4, 3, 2))),
    (np.array([-1 / 2, -2 / 3, -1]), MillerIndex((-4, -3, -2))),
    (np.array([1, np.inf, np.inf]), MillerIndex((1, 0, 0))),
    (np.array([-1, np.inf, np.inf]), MillerIndex((-1, 0, 0))),
    (np.array([1 / 2, np.inf, np.inf]), MillerIndex((2, 0, 0))),
    (np.array([-1 / 2, np.inf, np.inf]), MillerIndex((-2, 0, 0))),
    (np.array([1, 1, 1]), MillerIndex((1, 1, 1))),
    (np.array([-1, -1, -1]), MillerIndex((-1, -1, -1))),
])
def test_miller_index_from_intercepts(test_input, expectation):
    res = MillerIndex.from_intercepts(test_input)
    assert res == expectation


@pytest.mark.parametrize("test_input,expectation", [
    (MillerIndex((4, 3, 2)), np.array([1 / 2, 2 / 3, 1])),
    (MillerIndex((-4, -3, -2)), np.array([-1 / 2, -2 / 3, -1])),
    (MillerIndex((1, 0, 0)), np.array([1, np.inf, np.inf])),
    (MillerIndex((-1, 0, 0)), np.array([-1, np.inf, np.inf])),
    (MillerIndex((2, 0, 0)), np.array([1 / 2, np.inf, np.inf])),
    (MillerIndex((-2, 0, 0)), np.array([-1 / 2, np.inf, np.inf])),
    (MillerIndex((1, 1, 1)), np.array([1, 1, 1])),
    (MillerIndex((-1, -1, -1)), np.array([-1, -1, -1])),
])
def test_miller_index_intercepts(test_input, expectation):
    res = test_input.intercepts
    assert np.allclose(res, expectation)


def test_miller_index_equality():
    hkl = (1, 2, 3)
    uvw = (3, 2, 1)
    # valid equality
    assert MillerIndex(hkl) == MillerIndex(hkl)
    # valid inequality
    assert MillerIndex(hkl) != MillerIndex(uvw)
    # invalid inequality
    assert MillerIndex(hkl) != hkl


###########################
#    Orientation Tests    #
###########################


def test_orientation_miller_indices():
    plane = MillerIndex((1, 0, 0))
    direction = MillerIndex((1, 2, 0))
    orientation = Orientation.from_miller_indices(plane, direction)
    res_plane, res_direction = orientation.as_miller_indices()
    assert res_plane == plane
    assert res_direction == direction


#####################
#    Plane Tests    #
#####################


@pytest.mark.parametrize("test_input,expectation", [
    (MillerIndex((4, 3, 2)), np.array([[1 / 2, 0, 0], [0, 2 / 3, 0], [0, 0, 1]])),
    (MillerIndex((1, 0, 0)), np.array([[1, 0, 0], [1, 1, 0], [1, 0, 1]])),
])
def test_plane_from_miller_index(test_input, expectation):
    res = Plane.from_miller_index(test_input)
    assert np.allclose(res.coplanar_points, expectation)


def test_plane_coefficients():
    plane = Plane(np.array([
        [1, 2, 3],
        [4, 6, 9],
        [12, 11, 9],
    ]))
    assert np.allclose(plane.coefficients, np.array([30, -48, 17, -15]))
