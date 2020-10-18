import numpy as np

from atompack.atom import Atom
from atompack.bond import Bond


def test_bond_to_from_json():
    atom0 = Atom(np.zeros(3), "X", test="test")
    atom1 = Atom(np.ones(3), "Y", test="test")
    bond = Bond((atom0, atom1), test="test")
    json_data = bond.to_json()
    new_bond = Bond.from_json(json_data)
    new_atom0, new_atom1 = new_bond.endpoints
    for old, new in zip([atom0, atom1], [new_atom0, new_atom1]):
        assert np.array_equal(new.position, old.position)
        assert new.specie == old.specie
        assert new["test"] == old["test"]
    assert new_bond["test"] == bond["test"]
