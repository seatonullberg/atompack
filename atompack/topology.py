from typing import List

import numpy as np
from igraph import Graph, Vertex


class Topology(object):
    """A specialized graph with spatial data."""

    def __init__(self) -> None:
        self._graph = Graph()

    def insert(self, position: np.ndarray) -> int:
        """Insert a point into the topology and return its index."""
        self._graph.add_vertices(1)
        self._graph.vs[-1]["position"] = position
        return len(self._graph.vs) - 1

    def connect(self, a: int, b: int) -> None:
        """Creates an edge between indices `a` and `b`."""
        self._graph.add_edges([(a, b)])

    def disconnect(self, a: int, b: int) -> None:
        """Destroys the edge between indices `a` and `b`."""
        eid = self._graph.get_eid(a, b)
        self._graph.delete_edges(eid)

    def nearest(self, position: np.ndarray) -> int:
        """Returns the index of the point located nearest to `position`."""
        index = 0
        distance = np.inf
        for i in range(len(self._graph.vs)):
            point = self._graph.vs.select(i)
            d = np.linalg.norm(position - point["position"])
            if d < distance:
                index = i
                distance = d
        return index

    def merge(self, other: 'Topology') -> List[int]:
        """Combines two topologies and returns the indices of the new points."""
        return [self.insert(vertex["position"]) for vertex in other._graph.vs]

    def remove(self, index: int) -> None:
        """Removes a point from the topology by index.
        
        Note:
            This operation invalidates existing indices.
            Indices are always continuous: deleting a point 
            will cause all points added chronologically 
            afterwards to be renumbered.
        """
        self._graph.delete_vertices(index)
