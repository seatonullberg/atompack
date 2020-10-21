import numpy as np

from atompack.bond import Bond


def test_bond_to_from_json():
    indices = (0, 1)
    test_value = "test"
    bond = Bond(indices, test_value=test_value)
    json_data = bond.to_json()
    new_bond = Bond.from_json(json_data)
    assert new_bond.indices == bond.indices
    assert new_bond["test_value"] == bond["test_value"]
