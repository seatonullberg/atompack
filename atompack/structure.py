import copy
from typing import List, Optional

import numpy as np

from atompack import Atom
from atompack.bindings import cell_contains, load_libatompack, nearest_neighbor
from atompack.errors import (PositionOccupiedError, PositionOutsideError, PositionUnoccupiedError)


class Structure(object):

    def __init__(self,
                 atoms: Optional[List[Atom]] = None,
                 basis: Optional[np.ndarray] = None,
                 periodicity: Optional[np.ndarray] = None,
                 tolerance: float = 1.0e-6) -> None:
        if atoms is None:
            atoms = []
        self.atoms = atoms
        if basis is None:
            basis = np.identity(3)
        self.basis = basis
        if periodicity is None:
            periodicity = np.array([0, 0, 0])
        self.periodicity = periodicity
        self.tolerance = tolerance
        self._lib = load_libatompack()

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
        if not cell_contains(self._lib, self.basis, atom["position"], self.tolerance):
            raise PositionOutsideError(atom["position"])
        _, distance = nearest_neighbor(self._lib, atom["position"], self._positions(), self.basis, self.periodicity,
                                       self.tolerance)
        if distance < self.tolerance:
            raise PositionOccupiedError(atom["position"])
        self.atoms.append(atom)

    def remove(self, position: np.ndarray) -> Atom:
        index, distance = nearest_neighbor(self._lib, position, self._positions(), self.basis, self.periodicity,
                                           self.tolerance)
        if distance > self.tolerance:
            raise PositionUnoccupiedError(position)
        res = copy.deepcopy(self.atoms[index])
        del self.atoms[index]
        return res

    def select(self, position: np.ndarray) -> int:
        index, distance = nearest_neighbor(self._lib, position, self._positions(), self.basis, self.periodicity,
                                           self.tolerance)
        if distance > self.tolerance:
            raise PositionUnoccupiedError(position)
        return index

    def _positions(self):
        return np.array([atom["position"] for atom in self.atoms])
