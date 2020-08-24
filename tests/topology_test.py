import numpy as np

from atompack.topology import Topology


def test_topology_insert():
    t = Topology()
    t.insert(np.zeros(3))
    assert len(t._graph.vs) == 1


def test_topology_connect():
    t = Topology()
    a = t.insert(np.zeros(3))
    b = t.insert(np.zeros(3))
    t.connect(a, b)
    assert len(t._graph.es) == 1


def test_topology_disconnect():
    t = Topology()
    a = t.insert(np.zeros(3))
    b = t.insert(np.zeros(3))
    t.connect(a, b)
    t.disconnect(a, b)
    assert len(t._graph.es) == 0


def test_topology_nearest():
    t = Topology()
    position = np.zeros(3)
    t.insert(np.array([1.0, 0.0, 0.0]))
    t.insert(np.array([0.0, 1.0, 0.0]))
    t.insert(np.array([0.0, 0.0, 1.0]))
    t.insert(np.array([0.5, 0.5, 0.5]))
    assert t.nearest(position) == 3


def test_topology_merge():
    t = Topology()
    t.insert(np.zeros(3))
    new_t = Topology()
    new_t.insert(np.zeros(3))
    t.merge(new_t)
    assert len(t._graph.vs) == 2
    assert len(new_t._graph.vs) == 1


def test_topology_remove():
    t = Topology()
    t.insert(np.zeros(3))
    t.remove(0)
    assert len(t._graph.vs) == 0
