import pytest

from atompack.symmetry import Spacegroup


def test_spacegroup_init_valid():
    international_number = 1
    hermann_mauguin = "P 1"
    # init from number
    assert Spacegroup(international_number).hermann_mauguin == hermann_mauguin
    # init from symbol
    assert Spacegroup(hermann_mauguin).international_number == international_number


def test_spacegroup_init_invalid():
    international_number = 231
    hermann_mauguin = "x"
    # init from number
    with pytest.raises(ValueError):
        _ = Spacegroup(international_number)
    # init from symbol
    with pytest.raises(ValueError):
        _ = Spacegroup(hermann_mauguin)
    # init from invalid type
    with pytest.raises(TypeError):
        _ = Spacegroup([])


def test_spacegroup_equality():
    international_number = 1
    hermann_mauguin = "P 1"
    # valid equality
    assert Spacegroup(international_number) == Spacegroup(hermann_mauguin)
    # valid inequality
    assert Spacegroup(international_number) != Spacegroup(international_number + 1)
    # invalid comparison
    assert Spacegroup(international_number) != international_number
