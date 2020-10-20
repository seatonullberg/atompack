import numpy as np

from atompack.crystal.spatial import LatticeVectors, MillerIndex


def test_miller_index_reciprocal():
    miller = MillerIndex((3, 2, 0))
    assert np.allclose(miller.reciprocal, np.array([1.0, 2 / 3, np.inf]))


def test_lattice_vectors_contain():
    vectors = LatticeVectors(np.identity(3))
    # edge point
    point = np.array([0.0, 0.0, 0.0])
    assert vectors.contain(point)
    # surface point
    point = np.array([0.0, 0.5, 0.0])
    assert vectors.contain(point)
    # interior point
    point = np.array([0.5, 0.5, 0.5])
    assert vectors.contain(point)
    # exterior point (positive)
    point = np.array([1.5, 1.5, 1.5])
    assert not vectors.contain(point)
    # exterior point (negative)
    point = np.array([-0.5, -0.5, -0.5])
    assert not vectors.contain(point)


def test_lattice_vectors_wrap():
    vectors = LatticeVectors(np.identity(3))
    # edge point
    point = np.array([0.0, 0.0, 0.0])
    target = np.array([0.0, 0.0, 0.0])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # surface point
    point = np.array([0.0, 0.5, 0.0])
    target = np.array([0.0, 0.5, 0.0])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # interior point
    point = np.array([0.5, 0.5, 0.5])
    target = np.array([0.5, 0.5, 0.5])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # exterior point (near) (positive)
    point = np.array([1.1, 1.1, 1.1])
    target = np.array([0.1, 0.1, 0.1])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # exterior point (far) (positive)
    point = np.array([5.9, 5.9, 5.9])
    target = np.array([0.9, 0.9, 0.9])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # exterior point (near) (negative)
    point = np.array([-0.1, -0.1, -0.1])
    target = np.array([0.9, 0.9, 0.9])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
    # exterior point (far) (positive)
    point = np.array([-5.9, -5.9, -5.9])
    target = np.array([0.1, 0.1, 0.1])
    res = vectors.wrap(point)
    assert np.allclose(res, target)
