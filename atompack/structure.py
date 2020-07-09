from typing import List, Optional, Tuple

import numpy as np

from atompack.atom import Atom


class Structure(object):
    """Representation of an atomic structure with boundaries which may be periodic.
    
    Args:
        atoms: List of atoms in the structure.
            Mutating `atoms` is a logical error if the change results in atoms overlapping with eachother or atoms existing out of bounds.
        basis: 3x3 matrix defining the boundaries of the structure.
            Mutating `basis` is a logical error if the change forces atoms out of bounds.
        pbc: Boolean array that indicates which boundaries are considered periodic.
        tolerance: Radius of tolerance for operations on the structure.
            Mutating `tolerance` is a logical error if the change results in atoms overlapping with eachother or atoms existing out of bounds.
    """

    def __init__(self,
                 atoms: List[Atom],
                 basis: np.ndarray,
                 pbc: Optional[Tuple[bool, bool, bool]] = None,
                 tolerance: Optional[float] = None) -> None:
        self.atoms = atoms
        self.basis = basis
        if pbc is None:
            pbc = (False, False, False)
        self.pbc = pbc
        if tolerance is None:
            tolerance = 1.0e-6
        self.tolerance = tolerance

    def __iter__(self) -> 'Structure':
        self._iter = 0
        return self

    def __next__(self) -> Atom:
        if self._iter == len(self):
            raise StopIteration
        res = self.atoms[self._iter]
        self._iter += 1
        return res

    def __len__(self) -> int:
        return len(self.atoms)

    def insert(self, atom: Atom) -> None:
        """Inserts an atom into the structure safely with bounds checking.
        
        Args:
            atom: Atom inserted into the structure.

        Raises:
            `atompack.errors.PositionOccupiedError`: If the position is already occupied.
            `atompack.errors.PositionOutsideError`: If the position is out of bounds.
        """
        raise NotImplementedError

    def remove(self, position: np.ndarray) -> Atom:
        """Removes and returns an atom from the structure.
        
        Args:
            position: Position to remove the atom from.

        Raises:
            `atompack.errors.PositionUnoccupiedError`: If the position is not occupied.
        """
        raise NotImplementedError

    def select(self, position: np.ndarray) -> int:
        """Returns the index of an atom.
        
        Args:
            position: Position to select the atom from.

        Raises:
            `atompack.errors.PositionUnoccupiedError`: If the position is not occupied.
        """
        raise NotImplementedError
