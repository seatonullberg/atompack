from atompack.atom import Atom
from atompack.internal import search_for_atom

import numpy as np


def test_search_for_atom_occupied():
    atoms = [Atom(position=np.array([0, 0, 0]))]
    position = np.array([0, 0, 0])
    res = search_for_atom(atoms, position, 1e-6)
    assert res == 0


def test_search_for_atom_unoccupied():
    atoms = []
    position = np.array([0, 0, 0])
    res = search_for_atom(atoms, position, 1e-6)
    assert res is None
