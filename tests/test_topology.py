import numpy as np

from atompack.atom import Atom
from atompack.bond import Bond
from atompack.topology import Topology


def test_topology_remove_atom_removes_bonds():
    top = Topology()
    top.insert_atom(Atom("X", np.zeros(3), test="test"))
    top.insert_atom(Atom("X", np.zeros(3), test="test"))
    top.insert_bond(Bond((0, 1), test="test"))
    assert len(top.atoms) == 2
    assert len(top.bonds) == 1
    top.remove_atom(0)
    assert len(top.atoms) == 1
    assert len(top.bonds) == 0


def test_topology_to_from_json():
    top = Topology()
    top.insert_atom(Atom("X", np.zeros(3), test="test"))
    top.insert_atom(Atom("X", np.zeros(3), test="test"))
    top.insert_bond(Bond((0, 1), test="test"))
    json_data = top.to_json()
    new_top = Topology.from_json(json_data)
    assert new_top.atoms[0].specie == top.atoms[0].specie
    assert np.array_equal(new_top.atoms[0].position, top.atoms[0].position)
    assert new_top.atoms[0]["test"] == top.atoms[0]["test"]
    assert new_top.bonds[0].indices == top.bonds[0].indices
    assert new_top.bonds[0]["test"] == top.bonds[0]["test"]
