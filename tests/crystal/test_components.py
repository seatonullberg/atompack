import numpy as np
import pytest

from atompack.crystal.components import (Basis, LatticeParameters, LatticeVectors)
from atompack.symmetry import Spacegroup

#####################
#    Basis Tests    #
#####################


def test_basis_init_not_fractional():
    specie = "Fe"
    positive_site = np.array([1.5, 0, 0])
    negative_site = np.array([-1.5, 0, 0])
    # test a positive value
    with pytest.raises(ValueError):
        _ = Basis([(specie, positive_site)])
    # test a negative value
    with pytest.raises(ValueError):
        _ = Basis([(specie, negative_site)])


def test_basis_apply_spacegroup_primitive_cubic():
    # primitive basis of aluminum
    basis = Basis.primitive("Al")
    # FCC spacegroup from Hermann Mauguin symbol
    spg = Spacegroup("F m -3 m")
    res = basis.apply_spacegroup(spg)
    assert len(res) == 4
    target = [
        ("Al", np.array([0.0, 0.0, 0.0])),
        ("Al", np.array([0.0, 0.5, 0.5])),
        ("Al", np.array([0.5, 0.0, 0.5])),
        ("Al", np.array([0.5, 0.5, 0.0])),
    ]
    for (target_specie, target_site), (res_specie, res_site) in zip(target, res):
        assert target_specie == res_specie
        assert np.allclose(target_site, res_site)


def test_basis_apply_spacegroup_complex_monoclinic():
    # complex basis of arbitrary species
    basis = Basis([
        ("X", np.array([0.0, 0.0, 0.0])),
        ("Y", np.array([0.1, 0.2, 0.3])),
    ])
    # Monoclinic spacegroup from international number
    spg = Spacegroup(3)
    res = basis.apply_spacegroup(spg)
    assert len(res) == 3
    target = [
        ("X", np.array([0.0, 0.0, 0.0])),
        ("Y", np.array([0.1, 0.2, 0.3])),
        ("Y", np.array([0.9, 0.2, 0.7])),
    ]
    for (target_specie, target_site), (res_specie, res_site) in zip(target, res):
        assert target_specie == res_specie
        assert np.allclose(target_site, res_site)


def test_basis_to_from_json():
    basis = Basis.primitive("X")
    json_data = basis.to_json()
    res = Basis.from_json(json_data)
    assert len(basis) == len(res) == 1
    assert res[0][0] == basis[0][0]
    assert np.array_equal(res[0][1], basis[0][1])


#################################
#    LatticeParameters Tests    #
#################################


def test_lattice_parameters_to_from_json():
    params = LatticeParameters.cubic(10)
    json_data = params.to_json()
    res = LatticeParameters.from_json(json_data)
    assert res.a == params.a
    assert res.alpha == params.alpha


def test_lattice_parameters_unused_kwargs():
    kwargs = {"a": 5.0, "b": 0.0, "alpha": 0.0}
    params = LatticeParameters.cubic(**kwargs)
    assert params.a == params.b == params.c == 5.0
    assert params.alpha == params.beta == params.gamma == np.pi / 2


##############################
#    LatticeVectors Tests    #
##############################


@pytest.mark.parametrize("test_input,expectation", [
    (np.array([0.0, 0.0, 0.0]), True),
    (np.array([0.0, 0.5, 0.5]), True),
    (np.array([0.5, 0.5, 0.5]), True),
    (np.array([1.5, 1.5, 1.5]), False),
    (np.array([-0.5, -0.5, -0.5]), False),
])
def test_lattice_vectors_contain(test_input, expectation):
    vectors = LatticeVectors(np.identity(3))
    assert vectors.contain(test_input) is expectation


@pytest.mark.parametrize("test_input,expectation", [
    (np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])),
    (np.array([0.0, 0.5, 0.0]), np.array([0.0, 0.5, 0.0])),
    (np.array([0.5, 0.5, 0.5]), np.array([0.5, 0.5, 0.5])),
    (np.array([1.1, 1.1, 1.1]), np.array([0.1, 0.1, 0.1])),
    (np.array([5.9, 5.9, 5.9]), np.array([0.9, 0.9, 0.9])),
    (np.array([-0.1, -0.1, -0.1]), np.array([0.9, 0.9, 0.9])),
    (np.array([-5.9, -5.9, -5.9]), np.array([0.1, 0.1, 0.1])),
])
def test_lattice_vectors_wrap(test_input, expectation):
    vectors = LatticeVectors(np.identity(3))
    res = vectors.wrap(test_input)
    assert np.allclose(res, expectation)


def test_lattice_vectors_to_from_json():
    vectors = LatticeVectors(np.identity(3))
    json_data = vectors.to_json()
    res = LatticeVectors.from_json(json_data)
    assert np.allclose(res.vectors, vectors.vectors)
