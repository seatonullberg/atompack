from typing import Tuple

from atompack.atom import Atom
from atompack.bond import Bond


class Topology(object):
    """Internal abstraction for a collection of atoms and bonds."""

    def __init__(self) -> None:
        self._atoms: List[Atom] = []
        self._bonds: List[Bond] = []

    ####################
    #    Properties    #
    ####################

    @property
    def atoms(self) -> List[Atom]:
        """Returns a list of all atoms in the topology."""
        return self._atoms

    @property
    def bonds(self) -> List[Bond]:
        """Returns a list of all bonds in the topology."""
        return self._bonds
