import numpy as np

from atompack.atom import Atom
from atompack.bond import Bond


def test_bond_to_from_json():
    bond = Bond((0, 1), test="test")
    json_data = bond.to_json()
    new_bond = Bond.from_json(json_data)
    assert new_bond.indices == bond.indices
    assert new_bond["test"] == bond["test"]
