from typing import List, Optional, Tuple

import numpy as np

from atompack import Atom
from atompack.structure import Structure


class Crystal(Structure):

    def __init__(self,
                 atoms: Optional[List[Atom]] = None,
                 basis: Optional[np.ndarray] = None,
                 periodicity: Optional[np.ndarray] = None,
                 tolerance: float = 1.0e-6) -> None:
        super().__init__(atoms, basis, periodicity, tolerance)

    @classmethod
    def from_lattice_sites(
            cls,
            lattice_sites: List[Tuple[Atom, np.ndarray]],
            a: float,
            b: float,
            c: float,
            alpha: float,
            beta: float,
            gamma: float,
            duplicates: Tuple[int, int, int],
            orientation: np.ndarray,
    ) -> 'Crystal':
        pass

    @classmethod
    def from_spacegroup(
            cls,
            spacegroup: str,
            a: float,
            b: float,
            c: float,
            alpha: float,
            beta: float,
            gamma: float,
            duplicates: Tuple[int, int, int],
            orientation: np.ndarray,
    ) -> 'Crystal':
        pass

    @staticmethod
    def _build(lattice_sites: List[Tuple[Atom, np.ndarray]], a: float, b: float, c: float, alpha: float, beta: float,
               gamma: float, duplicates: Tuple[int, int,
                                               int], orientation: np.ndarray) -> Tuple[List[Atom], np.ndarray]:
        pass
