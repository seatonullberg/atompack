from typing import List

import numpy as np
from igraph import Graph, Vertex

from atompack.atom import Atom


class Topology(object):
    """An undirected graph of atoms."""

    def __init__(self) -> None:
        self._graph = Graph()

    def insert(self, atom: Atom) -> int:
        """Insert an atom into the topology and return its index."""
        self._graph.add_vertices(1)
        self._graph.vs[-1]["atom"] = atom
        return len(self._graph.vs) - 1

    def connect(self, a: int, b: int) -> None:
        """Creates an edge between indices `a` and `b`."""
        self._graph.add_edges([(a, b)])

    def disconnect(self, a: int, b: int) -> None:
        """Destroys the edge between indices `a` and `b`."""
        eid = self._graph.get_eid(a, b)
        self._graph.delete_edges(eid)

    def nearest(self, position: np.ndarray) -> int:
        """Returns the index of the atom located nearest to `position`."""
        index = 0
        distance = np.inf
        for vertex in self._graph.vs:
            d = np.linalg.norm(position - vertex["atom"].position)
            if d < distance:
                index = vertex.index
                distance = d
        return index

    def merge(self, other: 'Topology') -> List[int]:
        """Combines two topologies and returns the indices of the new atoms."""
        return [self.insert(vertex["atom"]) for vertex in other._graph.vs]

    def remove(self, index: int) -> None:
        """Removes an atom from the topology by index.
        
        Note:
            This operation invalidates existing indices.
            Indices are always continuous: deleting a point 
            will cause all points added chronologically 
            afterwards to be renumbered.
        """
        self._graph.delete_vertices(index)

    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all atoms in the topology."""
        return [vertex["atom"] for vertex in self._graph.vs]
