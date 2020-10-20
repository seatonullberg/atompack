import numpy as np

from atompack.crystal.components import Basis, LatticeParameters
from atompack.crystal.crystal import Crystal, UnitCell
from atompack.crystal.spatial import MillerIndex
from atompack.symmetry import Spacegroup


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
        ("Fe", np.array([0.0, 0.0, 0.0])),
        ("Fe", np.array([1.425, 1.425, 1.425])),
    ]
    for i, (_specie, _site) in enumerate(target):
        assert unit_cell.atoms[i].specie == _specie
        assert np.allclose(unit_cell.atoms[i].position, _site)


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
        ("X", np.array([0.0, 0.0, 0.0])),
        ("X", np.array([0.0, 0.0, 1.25])),
    ]
    for i, (_specie, _site) in enumerate(target):
        assert unit_cell.atoms[i].specie == _specie
        assert np.allclose(unit_cell.atoms[i].position, _site)


def test_unit_cell_to_from_json():
    basis = Basis.primitive("X")
    params = LatticeParameters.cubic(10)
    spg = Spacegroup("F m -3 m")
    unit_cell = UnitCell(basis, params, spg)
    json_data = unit_cell.to_json()
    new_unit_cell = UnitCell.from_json(json_data)
    # test basis
    assert new_unit_cell.basis[0][0] == unit_cell.basis[0][0]
    assert np.array_equal(new_unit_cell.basis[0][1], unit_cell.basis[0][1])
    # test lattice parameters
    assert new_unit_cell.lattice_parameters.a == unit_cell.lattice_parameters.a
    assert new_unit_cell.lattice_parameters.alpha == unit_cell.lattice_parameters.alpha
    # test spacegroup
    assert new_unit_cell.spacegroup == unit_cell.spacegroup
    # test atoms/bonds
    assert len(new_unit_cell.atoms) == len(unit_cell.atoms) == 4
    assert len(new_unit_cell.bonds) == len(unit_cell.bonds) == 0
    assert np.array_equal(
        new_unit_cell.atoms[0].position,
        unit_cell.atoms[0].position,
    )
    assert new_unit_cell.atoms[0].specie == unit_cell.atoms[0].specie


def test_crystal_supercell_cubic():
    # primitive basis of iron
    basis = Basis.primitive("Fe")
    # cubic lattice parameters
    lattparams = LatticeParameters.cubic(2.85)
    # BCC spacegroup
    spg = Spacegroup("I m -3 m")
    # build the crystal
    unit_cell = UnitCell(basis, lattparams, spg)
    crystal = Crystal(unit_cell)
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


# def test_crystal_cut_cubic():
#     # primitive basis of iron
#     basis = Basis.primitive("Fe")
#     # cubic lattice parameters
#     lattparams = LatticeParameters.cubic(2.85)
#     # BCC spacegroup
#     spg = Spacegroup("I m -3 m")
#     # build the crystal
#     unit_cell = UnitCell(basis, lattparams, spg)
#     crystal = Crystal(unit_cell)
#     # define a plane to cut along
#     plane = MillerIndex((1, 1, 1))
#     # cut along plane
#     crystal.cut(plane).finish()
#     assert len(crystal.atoms) == 2


def test_crystal_to_from_json():
    basis = Basis.primitive("Fe")
    lattparams = LatticeParameters.cubic(2.85)
    spg = Spacegroup("I m -3 m")
    unit_cell = UnitCell(basis, lattparams, spg)
    crystal = Crystal(unit_cell)
    json_data = crystal.to_json()
    new_crystal = Crystal.from_json(json_data)
    assert np.array_equal(
        new_crystal.lattice_vectors,
        crystal.lattice_vectors,
    )
    assert len(new_crystal.atoms) == len(crystal.atoms) == 2
    assert np.array_equal(
        new_crystal.unit_cell.lattice_vectors,
        crystal.unit_cell.lattice_vectors,
    )
