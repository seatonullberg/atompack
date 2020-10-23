import numpy as np

from atompack.crystal.components import Basis, LatticeParameters
from atompack.crystal.crystal import Crystal, UnitCell
from atompack.crystal.spatial import MillerIndex
from atompack.symmetry import Spacegroup

########################
#    UnitCell Tests    #
########################

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
    for i, (specie, site) in enumerate(target):
        assert unit_cell.atoms[i].specie == specie
        assert np.allclose(unit_cell.atoms[i].position, site)


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
    for i, (specie, site) in enumerate(target):
        assert unit_cell.atoms[i].specie == specie
        assert np.allclose(unit_cell.atoms[i].position, site)


def test_unit_cell_to_from_json():
    basis = Basis.primitive("X")
    params = LatticeParameters.cubic(10)
    spg = Spacegroup("F m -3 m")
    unit_cell = UnitCell(basis, params, spg)
    json_data = unit_cell.to_json()
    res = UnitCell.from_json(json_data)
    # test basis
    assert res.basis[0][0] == unit_cell.basis[0][0]
    assert np.array_equal(res.basis[0][1], unit_cell.basis[0][1])
    # test lattice parameters
    assert res.lattice_parameters.a == unit_cell.lattice_parameters.a
    assert res.lattice_parameters.alpha == unit_cell.lattice_parameters.alpha
    # test spacegroup
    assert res.spacegroup == unit_cell.spacegroup
    # test atoms/bonds
    assert len(res.atoms) == len(unit_cell.atoms) == 4
    assert len(res.bonds) == len(unit_cell.bonds) == 0
    assert np.array_equal(
        res.atoms[0].position,
        unit_cell.atoms[0].position,
    )
    assert res.atoms[0].specie == unit_cell.atoms[0].specie

#######################
#    Crystal Tests    #
#######################

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
    target_vectors = np.array([
        [2.85, 0.0, 0.0],
        [0.0, 2.85, 0.0],
        [0.0, 0.0, 2.85],
    ])
    # generate supercell
    extent0 = (1, 2, 3)
    crystal.supercell(extent0).finish()
    assert len(crystal.atoms) == 12
    assert np.allclose(crystal.lattice_vectors.vectors, target_vectors * np.array(extent0), atol=1E-6)
    # generate another supercell
    extent1 = (2, 2, 2)
    crystal.supercell(extent1).finish()
    assert len(crystal.atoms) == 96
    assert np.allclose(crystal.lattice_vectors.vectors,
                       target_vectors * np.array(extent0) * np.array(extent1),
                       atol=1E-6)
    # reset the supercell
    crystal.reset()
    assert len(crystal.atoms) == 2
    assert np.allclose(crystal.lattice_vectors.vectors, target_vectors, atol=1E-6)


def test_crystal_to_from_json():
    basis = Basis.primitive("Fe")
    lattparams = LatticeParameters.cubic(2.85)
    spg = Spacegroup("I m -3 m")
    unit_cell = UnitCell(basis, lattparams, spg)
    crystal = Crystal(unit_cell)
    json_data = crystal.to_json()
    res = Crystal.from_json(json_data)
    assert np.array_equal(
        res.lattice_vectors.vectors,
        crystal.lattice_vectors.vectors,
    )
    assert len(res.atoms) == len(crystal.atoms) == 2
