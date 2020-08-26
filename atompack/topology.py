from typing import List, Optional, Tuple

import numpy as np
from igraph import Graph

from atompack._cell import cell_contains, cell_enforce
from atompack.atom import Atom
from atompack.bond import Bond


class Topology(object):
    """An undirected graph of atoms."""

    def __init__(self) -> None:
        self._graph = Graph()

    def insert(self, atom: Atom) -> int:
        """Inserts an atom into the topology and returns its index."""
        self._graph.add_vertices(1)
        self._graph.vs[-1]["atom"] = atom
        return len(self._graph.vs) - 1

    def connect(self, a: int, b: int, bond: Optional[Bond] = None) -> None:
        """Creates an edge between indices `a` and `b`."""
        if bond is None:
            bond = Bond()
        self._graph.add_edges([(a, b)])
        self._graph.es[-1]["bond"] = bond

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
                vertex["atom"].position += translation
        elif shape == (n_atoms, 3):
            for i, vertex in enumerate(self._graph.vs):
                vertex["atom"].position += translation[i]
        else:
            raise ValueError

    def merge(self, other: 'Topology') -> List[int]:
        """Combines two topologies and returns the indices of the merged atoms.
        
        Args:
            other: The topology to merge in.
        """
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

    @property
    def bonds(self) -> List[Tuple[Tuple[int, int], Bond]]:
        """Returns a list of tuples of bond edge ids and bond objects."""
        return [((edge.source, edge.target), edge["bond"]) for edge in self._graph.es]


class BoundedTopology(Topology):
    """An undirected graph of atoms within a bounding paralellpiped cell.
    
    Args:
        cell: 3x3 matrix which defines the bounding area. 
    """

    def __init__(self, cell: np.ndarray) -> None:
        self.cell = cell
        super().__init__()

    def contains(self, position: np.ndarray, tolerance: float = 1.0e-6) -> bool:
        """Returns True if position is within the bounds of the topology."""
        return cell_contains(self.cell, position, tolerance)

    def check(self, tolerance: float = 1.0e-6) -> List[int]:
        """Returns the indices of atoms which are out of bounds."""
        indices = []
        for vertex in self._graph.vs:
            position = vertex["atom"].position
            if not self.contains(position, tolerance):
                indices.append(vertex.index)
        return indices

    def enforce(self, tolerance: float = 1.0e-6) -> None:
        """Enforces that all atoms in the topology are within the bounds.
            Any atoms found to be out of bounds will have their positions mutated to bring them back in.
        """
        for vertex in self._graph.vs:
            atom = vertex["atom"]
            cell_enforce(self.cell, atom.position, tolerance)
