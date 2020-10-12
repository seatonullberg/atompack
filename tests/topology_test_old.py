# import numpy as np

# from atompack.atom import Atom
# from atompack.bond import Bond
# from atompack.topology import Topology

# def test_topology_insert():
#     t = Topology()
#     t.insert(Atom(np.zeros(3)))
#     assert len(t._graph.vs) == 1

# def test_topology_connect():
#     t = Topology()
#     a = Atom(np.zeros(3))
#     b = Atom(np.zeros(3))
#     ai = t.insert(a)
#     bi = t.insert(b)
#     bond = Bond(kind="sigma")
#     t.connect(ai, bi, bond=bond)
#     assert len(t._graph.es) == 1
#     assert t.bonds[0][0] == 0
#     assert t.bonds[0][1] == 1
#     assert t.bonds[0][2].kind == "sigma"

# def test_topology_disconnect():
#     t = Topology()
#     a = t.insert(Atom(np.zeros(3)))
#     b = t.insert(Atom(np.zeros(3)))
#     t.connect(a, b)
#     t.disconnect(a, b)
#     assert len(t._graph.es) == 0

# def test_topology_nearest():
#     t = Topology()
#     position = np.zeros(3)
#     t.insert(Atom(np.array([1.0, 0.0, 0.0])))
#     t.insert(Atom(np.array([0.0, 1.0, 0.0])))
#     t.insert(Atom(np.array([0.0, 0.0, 1.0])))
#     t.insert(Atom(np.array([0.5, 0.5, 0.5])))
#     assert t.nearest(position) == 3

# def test_topology_translate():
#     t = Topology()
#     t.insert(Atom(np.array([0.5, 0.5, 0.5])))
#     t.insert(Atom(np.array([1.0, 1.0, 1.0])))
#     t.translate(np.array([1.0, 2.0, 3.0]))
#     assert np.array_equal(t.atoms[0].position, np.array([1.5, 2.5, 3.5]))
#     assert np.array_equal(t.atoms[1].position, np.array([2.0, 3.0, 4.0]))
#     t.translate(np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]))
#     assert np.array_equal(t.atoms[0].position, np.array([2.5, 3.5, 4.5]))
#     assert np.array_equal(t.atoms[1].position, np.array([4.0, 5.0, 6.0]))

# def test_topology_merge():
#     t = Topology()
#     t.insert(Atom(np.zeros(3)))
#     new_t = Topology()
#     new_t.insert(Atom(np.zeros(3)))
#     t.merge(new_t)
#     assert len(t._graph.vs) == 2
#     assert len(new_t._graph.vs) == 1

# def test_topology_remove():
#     t = Topology()
#     t.insert(Atom(np.zeros(3)))
#     t.remove(0)
#     assert len(t._graph.vs) == 0

# def test_topology_atoms():
#     t = Topology()
#     atom = Atom(np.zeros(3), symbol="H")
#     t.insert(atom)
#     assert t.atoms[0].symbol == "H"

# def test_topology_bonds():
#     t = Topology()
#     a = Atom(np.zeros(3))
#     b = Atom(np.zeros(3))
#     t.insert(a)
#     t.insert(b)
#     t.connect(0, 1, Bond(kind="sigma"))
#     assert t.bonds[0][0] == 0
#     assert t.bonds[0][1] == 1
#     assert t.bonds[0][2].kind == "sigma"
