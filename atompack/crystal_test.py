from atompack.atom import Atom
from atompack.crystal import Crystal

import numpy as np


def test_crystal_cubic_100():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])),
                 (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.identity(3)
    rotation = np.identity(3)
    size = (1, 1, 1)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation,
                      rotation, size)
    assert len(crystal) == 2
    assert np.allclose(crystal.basis, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
