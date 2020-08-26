from atompack.atom import Atom
from atompack.bond import Bond

import numpy as np


def test_bond_vector():
    a = Atom(np.zeros(3))
    b = Atom(np.array([1, 1, 1]))
    bond = Bond(a, b)
    res = bond.vector
    target = np.array([1, 1, 1])
    assert np.array_equal(res, target)
    a.position = np.array([0.5, 0.5, 0.5])
    b.position = np.array([0.75, 0.75, 0.75])
    res = bond.vector
    target = np.array([0.25, 0.25, 0.25])
    assert np.array_equal(res, target)
