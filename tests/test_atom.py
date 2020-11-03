import numpy as np

from atompack.atom import Atom

####################
#    Atom Tests    #
####################


def test_atom_to_from_json():
    specie = "X"
    position = np.zeros(3)
    test_value = "test"
    atom = Atom(specie, position, test_value=test_value)
    json_data = atom.to_json()
    res = Atom.from_json(json_data)
    assert res.specie == atom.specie
    assert np.allclose(res.position, atom.position)
    assert res["test_value"] == atom["test_value"]
