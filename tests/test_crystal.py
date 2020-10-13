import numpy as np
import pytest

from atompack.crystal import Basis
from atompack.spacegroup import Spacegroup


def test_basis_not_fractional():
    specie = "Fe"
    site = np.array([1.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(site, specie)])
    site = np.array([-0.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(site, specie)])


def test_basis_apply_spacegroup():
    # primitive basis of aluminum
    basis = Basis.primitive("Al")
    # FCC spacegroup from Hermann Mauguin symbol
    spg = Spacegroup("F m -3 m")
    atoms = basis.apply_spacegroup(spg)
    assert len(atoms) == 4
