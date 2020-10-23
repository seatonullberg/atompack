import numpy as np

from atompack.crystal.spatial import MillerIndex, Orientation

###########################
#    MillerIndex Tests    #
###########################

def test_miller_index_reciprocal():
    miller = MillerIndex((3, 2, 0))
    assert np.allclose(miller.reciprocal, np.array([1.0, 2 / 3, np.inf]))

def test_miller_indices_equality():
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
