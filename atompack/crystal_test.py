from atompack.atom import Atom
from atompack.crystal import Crystal

import numpy as np


def test_crystal_cubic_100_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])),
                 (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.identity(3)
    rotation = np.identity(3)
    size = (1, 1, 1)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation,
                      rotation, size)
    assert np.allclose(crystal.basis, np.array([[1, 0, 0], [0, 1, 0], [0, 0,
                                                                       1]]))
    assert len(crystal) == 2


def test_crystal_cubic_100_super_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])),
                 (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.identity(3)
    rotation = np.identity(3)
    size = (2, 2, 2)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation,
                      rotation, size)
    assert np.allclose(crystal.basis, np.array([[2, 0, 0], [0, 2, 0], [0, 0,
                                                                       2]]))
    assert len(crystal) == 16


def test_crystal_cubic_110_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])),
                 (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    rotation = np.identity(3)
    size = (1, 1, 1)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation,
                      rotation, size)
    assert np.allclose(
        crystal.basis,
        np.array([[np.sqrt(2), 0, 0], [0, np.sqrt(2), 0], [0, 0, 1]]))
    assert len(crystal) == 4
    assert 1 == 2
