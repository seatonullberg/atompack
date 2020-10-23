import numpy as np

from atompack.bond import Bond

####################
#    Bond Tests    #
####################

def test_bond_to_from_json():
    indices = (0, 1)
    test_value = "test"
    bond = Bond(indices, test_value=test_value)
    json_data = bond.to_json()
    res = Bond.from_json(json_data)
    assert res.indices == bond.indices
    assert res["test_value"] == bond["test_value"]
