import numpy as np
import pytest

from atompack.crystal.components import Basis, LatticeParameters
from atompack.symmetry import Spacegroup


def test_basis_not_fractional():
    specie = "Fe"
    site = np.array([1.5, 0, 0])
    with pytest.raises(ValueError):
        _ = Basis([(specie, site)])
    site = np.array([0, 0, -1.5])
    with pytest.raises(ValueError):
        _ = Basis([(specie, site)])


def test_basis_apply_spacegroup_primitive_cubic():
    # primitive basis of aluminum
    basis = Basis.primitive("Al")
    # FCC spacegroup from Hermann Mauguin symbol
    spg = Spacegroup("F m -3 m")
    new_basis = basis.apply_spacegroup(spg)
    assert len(new_basis) == 4
    target = [
        ("Al", np.array([0.0, 0.0, 0.0])),
        ("Al", np.array([0.0, 0.5, 0.5])),
        ("Al", np.array([0.5, 0.0, 0.5])),
        ("Al", np.array([0.5, 0.5, 0.0])),
    ]
    for i, (
            _specie,
            _site,
    ) in enumerate(target):
        assert new_basis[i][0] == _specie
        assert np.allclose(new_basis[i][1], _site)


def test_basis_apply_spacegroup_complex_monoclinic():
    # complex basis of arbitrary species
    basis = Basis([
        ("X", np.array([0.0, 0.0, 0.0])),
        ("Y", np.array([0.1, 0.2, 0.3])),
    ])
    # Monoclinic spacegroup from international number
    spg = Spacegroup(3)
    new_basis = basis.apply_spacegroup(spg)
    assert len(new_basis) == 3
    target = [
        ("X", np.array([0.0, 0.0, 0.0])),
        ("Y", np.array([0.1, 0.2, 0.3])),
        ("Y", np.array([0.9, 0.2, 0.7])),
    ]
    for i, (_specie, _site) in enumerate(target):
        assert new_basis[i][0] == _specie
        assert np.allclose(new_basis[i][1], _site)


def test_basis_to_from_json():
    basis = Basis.primitive("X")
    json_data = basis.to_json()
    new_basis = Basis.from_json(json_data)
    assert len(basis) == len(new_basis) == 1
    assert new_basis[0][0] == basis[0][0]
    assert np.array_equal(new_basis[0][1], basis[0][1])


def test_lattice_parameters_to_from_json():
    params = LatticeParameters.cubic(10)
    json_data = params.to_json()
    new_params = LatticeParameters.from_json(json_data)
    assert new_params.a == params.a
    assert new_params.alpha == params.alpha
