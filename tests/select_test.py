from atompack.atom import Atom
from atompack.select import Random
from atompack.topology import Topology
import numpy as np


def test_random_selection_from_number():
    t = Topology()
    t.insert(Atom(np.zeros(3), symbol="H"))
    t.insert(Atom(np.zeros(3), symbol="H"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    rand = Random.from_number(1, lambda x: x.symbol == "He")
    indices = rand(t)
    assert len(indices) == 1
    assert t.atoms[indices[0]].symbol == "He"


def test_random_selection_from_percentage():
    t = Topology()
    t.insert(Atom(np.zeros(3), symbol="H"))
    t.insert(Atom(np.zeros(3), symbol="H"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    t.insert(Atom(np.zeros(3), symbol="He"))
    rand = Random.from_percentage(0.2, lambda x: x.symbol == "He")
    indices = rand(t)
    assert len(indices) == 1
    assert t.atoms[indices[0]].symbol == "He"
