import numpy as np

from atompack.crystal import (Basis, Crystal, LatticeParameters, Transform, UnitCell)
from atompack.symmetry import Spacegroup


def test_transform_supercell():
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
    # generate supercell transform
    supercell_size = (1, 2, 3)
    transform = Transform().supercell(supercell_size)
    # apply the transform
    crystal = transform.apply(crystal)
    assert len(crystal.atoms) == 12
    assert np.allclose(crystal.lattice_vectors.vectors, target_vectors * np.array(supercell_size), atol=1E-6)
    # reapply the transform
    transform.apply(crystal)
    assert len(crystal.atoms) == 72
    assert np.allclose(crystal.lattice_vectors.vectors, target_vectors * np.array(supercell_size)**2, atol=1E-6)
