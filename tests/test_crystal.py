import numpy as np
import pytest

from atompack.crystal import Basis


def test_basis_not_fractional():
    specie = "Fe"
    site = np.array([1.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(specie, site)])
    site = np.array([-0.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(specie, site)])
