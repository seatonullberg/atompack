from typing import List, Optional, Tuple

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

    def translate(self, translation: np.ndarray) -> None:
        """Translates all atoms in the topology.
        
        Args:
            translation: The translation to apply.
                If `translation` is a 1x3 vector that vector will be applied to each atom.
                If `translation` is a Nx3 matrix each row will be mapped to each atom.  
        """
        shape = translation.shape
        n_atoms = len(self._graph.vs)
        if shape == (3,):
            for vertex in self._graph.vs:
                vertex["atom"]._position += translation
        elif shape == (n_atoms, 3):
            for i, vertex in enumerate(self._graph.vs):
                vertex["atom"]._position += translation[i]
        else:
            raise ValueError

    def merge(self, other: 'Topology', edges: Optional[List[Tuple[int, int]]] = None) -> List[int]:
        """Combines two topologies and returns the indices of the merged atoms.
        
        Args:
            other: The topology to merge in.
            edges: List of edges to create upon merge.
                In each tuple, the first index should be from `self` 
                and the second index should be from `other`.
        """
        if edges is None:
            edges = []
        indices = [self.insert(vertex["atom"]) for vertex in other._graph.vs]
        for edge in edges:
            self.connect(edge[0], indices[edge[1]])
        return indices

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
