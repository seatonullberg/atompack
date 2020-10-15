import numpy as np
import pytest

from atompack.crystal import Basis, Crystal, LatticeParameters, UnitCell
from atompack.spacegroup import Spacegroup


def test_basis_not_fractional():
    specie = "Fe"
    site = np.array([1.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(site, specie)])
    site = np.array([0, 0, -1.5])
    with pytest.raises(ValueError):
        _ = Basis([(site, specie)])


def test_basis_apply_spacegroup_primitive_cubic():
    # primitive basis of aluminum
    basis = Basis.primitive("Al")
    # FCC spacegroup from Hermann Mauguin symbol
    spg = Spacegroup("F m -3 m")
    new_basis = basis.apply_spacegroup(spg)
    assert len(new_basis) == 4
    target = [
        (np.array([0.0, 0.0, 0.0]), "Al"),
        (np.array([0.0, 0.5, 0.5]), "Al"),
        (np.array([0.5, 0.0, 0.5]), "Al"),
        (np.array([0.5, 0.5, 0.0]), "Al"),
    ]
    for i, (_site, _specie) in enumerate(target):
        assert np.allclose(new_basis[i][0], _site)
        assert new_basis[i][1] == _specie


def test_basis_apply_spacegroup_complex_monoclinic():
    # complex basis of arbitrary species
    basis = Basis([
        (np.array([0.0, 0.0, 0.0]), "X"),
        (np.array([0.1, 0.2, 0.3]), "Y"),
    ])
    # Monoclinic spacegroup from international number
    spg = Spacegroup(3)
    new_basis = basis.apply_spacegroup(spg)
    assert len(new_basis) == 3
    target = [
        (np.array([0.0, 0.0, 0.0]), "X"),
        (np.array([0.1, 0.2, 0.3]), "Y"),
        (np.array([0.9, 0.2, 0.7]), "Y"),
    ]
    for i, (_site, _specie) in enumerate(target):
        assert np.allclose(new_basis[i][0], _site)
        assert new_basis[i][1] == _specie


def test_unit_cell_cubic():
    # primitive basis of iron
    basis = Basis.primitive("Fe")
    # cubic lattice parameters
    lattparams = LatticeParameters.cubic(2.85)
    # BCC spacegroup
    spg = Spacegroup("I m -3 m")
    # build the unit cell
    unit_cell = UnitCell(basis, lattparams, spg)
    assert len(unit_cell.atoms) == 2
    target = [
        (np.array([0.0, 0.0, 0.0]), "Fe"),
        (np.array([1.425, 1.425, 1.425]), "Fe"),
    ]
    for i, (_site, _specie) in enumerate(target):
        assert np.allclose(unit_cell.atoms[i].position, _site)
        assert unit_cell.atoms[i].specie == _specie


def test_unit_cell_hexagonal():
    # primitive basis of arbitrary species
    basis = Basis.primitive("X")
    # cubic lattice parameters
    lattparams = LatticeParameters.hexagonal(4.0, 2.5)
    # hexagonal spacegroup
    spg = Spacegroup("P 6/m c c")
    # build the unit cell
    unit_cell = UnitCell(basis, lattparams, spg)
    assert len(unit_cell.atoms) == 2
    target = [
        (np.array([0.0, 0.0, 0.0]), "X"),
        (np.array([0.0, 0.0, 1.25]), "X"),
    ]
    for i, (_site, _specie) in enumerate(target):
        assert np.allclose(unit_cell.atoms[i].position, _site)
        assert unit_cell.atoms[i].specie == _specie


def test_crystal_supercell_cubic():
    # primitive basis of iron
    basis = Basis.primitive("Fe")
    # cubic lattice parameters
    lattparams = LatticeParameters.cubic(2.85)
    # BCC spacegroup
    spg = Spacegroup("I m -3 m")
    # build the crystal
    crystal = Crystal(basis, lattparams, spg)
    # generate supercell
    crystal.supercell((1, 2, 3)).finish()
    assert len(crystal.atoms) == 12
    assert np.allclose(crystal.lattice_vectors, crystal.unit_cell.lattice_vectors * np.array([1, 2, 3]))
    # generate another supercell
    crystal.supercell((2, 2, 2)).finish()
    assert len(crystal.atoms) == 96
    assert np.allclose(crystal.lattice_vectors, crystal.unit_cell.lattice_vectors * np.array([2, 4, 6]))
    # reset the supercell
    crystal.reset()
    assert len(crystal.atoms) == 2
    assert np.allclose(crystal.lattice_vectors, crystal.unit_cell.lattice_vectors)
