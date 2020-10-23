import numpy as np

from atompack.crystal.spatial import MillerIndex

###########################
#    MillerIndex Tests    #
###########################

def test_miller_index_reciprocal():
    miller = MillerIndex((3, 2, 0))
    assert np.allclose(miller.reciprocal, np.array([1.0, 2 / 3, np.inf]))
