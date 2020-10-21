import numpy as np

from atompack.atom import Atom


def test_atom_to_from_json():
    specie = "X"
    position = np.zeros(3)
    test_value = "test"
    atom = Atom(specie, position, test_value=test_value)
    json_data = atom.to_json()
    new_atom = Atom.from_json(json_data)
    assert new_atom.specie == atom.specie
    assert np.allclose(new_atom.position, atom.position)
    assert new_atom["test_value"] == atom["test_value"]
