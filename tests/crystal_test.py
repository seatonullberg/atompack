import numpy as np
from numpy.testing import assert_array_almost_equal

from atompack.atom import Atom
from atompack.crystal import Crystal


def test_crystal_cubic_100_unit_cell():
    a, b, c = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    crystal = Crystal(lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma)
    target_basis = np.array([[2.85, 0, 0], [0, 2.85, 0], [0, 0, 2.85]])
    assert_array_almost_equal(crystal.basis, target_basis)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([[0.0, 0.0, 0.0], [1.425, 1.425, 1.425]])
    assert_array_almost_equal(positions, target_positions)


def test_crystal_cubic_100_super_cell():
    a, b, c = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    duplicates = (2, 2, 2)
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    crystal = Crystal(lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, duplicates=duplicates)
    target_basis = np.array([[5.7, 0, 0], [0, 5.7, 0], [0, 0, 5.7]])
    assert_array_almost_equal(crystal.basis, target_basis)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([
        [0.00000000, 0.00000000, 0.00000000],
        [1.42500000, 1.42500000, 1.42500000],
        [0.00000000, 0.00000000, 2.85000000],
        [1.42500000, 1.42500000, 4.27500000],
        [0.00000000, 2.85000000, 0.00000000],
        [1.42500000, 4.27500000, 1.42500000],
        [0.00000000, 2.85000000, 2.85000000],
        [1.42500000, 4.27500000, 4.27500000],
        [2.85000000, 0.00000000, 0.00000000],
        [4.27500000, 1.42500000, 1.42500000],
        [2.85000000, 0.00000000, 2.85000000],
        [4.27500000, 1.42500000, 4.27500000],
        [2.85000000, 2.85000000, 0.00000000],
        [4.27500000, 4.27500000, 1.42500000],
        [2.85000000, 2.85000000, 2.85000000],
        [4.27500000, 4.27500000, 4.27500000],
    ])
    assert_array_almost_equal(positions, target_positions)


def test_crystal_cubic_110_unit_cell():
    a, b, c = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    crystal = Crystal(lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, orientation=orientation)
    target_basis = np.array([[np.sqrt(2) * 2.85, 0, 0], [0, 2.85, 0], [0, 0, np.sqrt(2) * 2.85]])
    assert_array_almost_equal(crystal.basis, target_basis)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.425, 2.015254],
        [2.015254, 0.0, 2.015254],
        [2.015254, 1.425, 0.0],
    ])
    assert_array_almost_equal(positions, target_positions)


def test_crystal_cubic_111_unit_cell():
    a, b, c = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    crystal = Crystal(lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, orientation=orientation)
    target_basis = np.array([[np.sqrt(2) * 2.85, 0, 0], [0, np.sqrt(6) * 2.85, 0], [0, 0, np.sqrt(3) * 2.85]])
    assert_array_almost_equal(crystal.basis, target_basis)
    positions = np.array([atom.position for atom in crystal.atoms])
    target_positions = np.array([
        [0.00000000, 0.00000000, 0.00000000],
        [0.00000000, 0.00000000, 2.46817240],
        [0.00000000, 4.65403051, 1.64544827],
        [0.00000000, 4.65403051, 4.11362067],
        [2.01525433, 1.16350763, 1.64544827],
        [2.01525433, 1.16350763, 4.11362067],
        [2.01525433, 5.81753814, 3.29089653],
        [2.01525433, 5.81753814, 0.82272413],
        [0.00000000, 2.32701526, 3.29089653],
        [0.00000000, 2.32701526, 0.82272413],
        [2.01525433, 3.49052288, 0.00000000],
        [2.01525433, 3.49052288, 2.46817240],
    ])
    assert_array_almost_equal(positions, target_positions)
