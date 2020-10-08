from typing import Tuple

from _atompack import topology_free, topology_new

from atompack.atom import Atom
from atompack.bond import Bond


class Topology(object):
    """Internal abstraction for an undirected graph of atoms and bonds."""

    def __init__(self) -> None:
        self._obj = topology_new()

    def __del__(self) -> None:
        topology_free(self._obj)

    ####################
    #    Properties    #
    ####################

    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all atoms in the topology."""
        pass

    @property
    def bonds(self) -> List[Bond]:
        """Returns a list of all bonds in the topology."""
        pass

    ########################
    #    Public Methods    #
    ########################

    def insert_atom(self, atom: Atom) -> None:
        """Inserts an atom into the topology."""
        pass

    def remove_atom(self, index: int) -> Atom:
        """Removes an atom from the topology."""
        pass

    def select_atom(self, index: int) -> Atom:
        """Returns a reference to an atom in the topology."""
        pass

    def insert_bond(self, bond: Bond) -> None:
        """Inserts a bond into the topology."""
        pass

    def remove_bond(self, indices: Tuple[int, int]) -> Bond:
        """Removes a bond from the topology."""
        pass

    def select_bond(self, indices: Tuple[int, int]) -> Bond:
        """Returns a reference to a bond in the topology."""
        pass
