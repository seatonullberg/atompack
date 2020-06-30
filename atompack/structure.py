from typing import List, Optional

import numpy as np

from atompack import Atom
from atompack.bindings import cell_contains, load_libatompack, nearest_neighbor


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
        pass

    def remove(self, position: np.ndarray) -> Atom:
        pass

    def select(self, position: np.ndarray) -> int:
        pass

    def _positions(self):
        return np.array([atom["position"] for atom in self.atoms])
