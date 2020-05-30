import numpy as np

from atompack.atom import Atom
from atompack.crystal import Crystal


def test_crystal_cubic_100_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell)
    assert np.allclose(
        crystal.basis,
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        rtol=1e-6,
        atol=1e-6,
    )
    assert len(crystal) == 2


def test_crystal_cubic_100_super_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    size = (2, 2, 2)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, size=size)
    assert np.allclose(
        crystal.basis,
        np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
        rtol=1e-6,
        atol=1e-6,
    )
    assert len(crystal) == 16


def test_crystal_cubic_110_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation=orientation)
    assert np.allclose(
        crystal.basis,
        np.array([[np.sqrt(2), 0, 0], [0, np.sqrt(2), 0], [0, 0, 1]]),
        rtol=1e-6,
        atol=1e-6,
    )
    assert len(crystal) == 4


def test_crystal_cubic_111_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation=orientation)
    print(crystal.basis[0][0])
    print(crystal.basis[1][1])
    print(crystal.basis[2][2])
    for atom in crystal.atoms:
        print(atom.position)
    assert np.allclose(
        crystal.basis,
        np.array([[np.sqrt(2), 0, 0], [0, np.sqrt(6), 0], [0, 0, np.sqrt(3) / 2]]),
        rtol=1e-6,
        atol=1e-6,
    )
    assert len(crystal) == 4
