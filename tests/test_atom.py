import numpy as np

from atompack.atom import Atom


def test_atom_to_from_json():
    position = np.zeros(3)
    specie = "X"
    atom = Atom(position, specie, test="test")
    json_data = atom.to_json()
    new_atom = Atom.from_json(json_data)
    assert np.array_equal(new_atom.position, atom.position)
    assert new_atom.specie == atom.specie == specie
    assert new_atom["test"] == atom["test"] == "test"
