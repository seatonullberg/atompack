import numpy as np

from atompack.atom import Atom
from atompack.crystal import Crystal


def test_crystal_cubic_100_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell)
    target_basis = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert np.allclose(crystal.basis, target_basis, atol=1e-6)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    assert np.allclose(positions, target_positions, atol=1e-6)


def test_crystal_cubic_100_super_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    size = (2, 2, 2)
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, size=size)
    assert np.allclose(crystal.basis, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]), rtol=1e-6, atol=1e-6)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.0, 0.0, 1.0],
        [0.5, 0.5, 1.5],
        [0.0, 1.0, 0.0],
        [0.5, 1.5, 0.5],
        [0.0, 1.0, 1.0],
        [0.5, 1.5, 1.5],
        [1.0, 0.0, 0.0],
        [1.5, 0.5, 0.5],
        [1.0, 0.0, 1.0],
        [1.5, 0.5, 1.5],
        [1.0, 1.0, 0.0],
        [1.5, 1.5, 0.5],
        [1.0, 1.0, 1.0],
        [1.5, 1.5, 1.5],
    ])
    assert np.allclose(positions, target_positions, atol=1e-6)


def test_crystal_cubic_110_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation=orientation)
    target_basis = np.array([[np.sqrt(2), 0, 0], [0, 1, 0], [0, 0, np.sqrt(2)]])
    assert np.allclose(crystal.basis, target_basis, atol=1e-6)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.70710678],
        [0.70710678, 0.0, 0.70710678],
        [0.70710678, 0.5, 0.0],
    ])
    assert np.allclose(positions, target_positions, atol=1e-6)


# TODO: verify all positions
def test_crystal_cubic_111_unit_cell():
    a, b, c = 1, 1, 1
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    unit_cell = [(Atom(symbol="Fe"), np.array([0, 0, 0])), (Atom(symbol="Fe"), np.array([0.5, 0.5, 0.5]))]
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    crystal = Crystal(a, b, c, alpha, beta, gamma, unit_cell, orientation=orientation)
    target_basis = np.array([[np.sqrt(2), 0, 0], [0, np.sqrt(6), 0], [0, 0, np.sqrt(3)]])
    assert np.allclose(crystal.basis, target_basis, atol=1e-6)
    positions = np.array([atom.position for atom in crystal.atoms])
    assert len(crystal) == 12
