from atompack.atom import Atom, AtomCollection

import numpy as np
from typing import List, Optional, Tuple


class Crystal(AtomCollection):
    """A crystalline lattice."""

    def __init__(self,
                 a: float = 0,
                 b: float = 0,
                 c: float = 0,
                 alpha: float = 0,
                 beta: float = 0,
                 gamma: float = 0,
                 unit_cell: Optional[List[Tuple[Atom, np.ndarray]]] = None,
                 orientation: Optional[np.ndarray] = None,
                 rotation: Optional[np.ndarray] = None,
                 size: Optional[Tuple[int, int, int]] = None) -> None:
        self._a, self._b, self._c = a, b, c
        self._alpha, self._beta, self._gamma = alpha, beta, gamma
        if unit_cell is None:
            unit_cell = []
        self._unit_cell = unit_cell
        if orientation is None:
            orientation = np.array([])
        self._orientation = orientation
        if rotation is None:
            rotation = np.array([])
        self._rotation = rotation
        if size is None:
            size = (0, 0, 0)
        self._size = size

    @classmethod
    def triclinic(cls) -> 'Crystal':
        pass

    @classmethod
    def monoclinic(cls) -> 'Crystal':
        pass

    @classmethod
    def orthorhombic(cls) -> 'Crystal':
        pass

    @classmethod
    def tetragonal(cls) -> 'Crystal':
        pass

    @classmethod
    def rhombohedral(cls) -> 'Crystal':
        pass

    @classmethod
    def hexagonal(cls) -> 'Crystal':
        pass

    @classmethod
    def cubic(cls) -> 'Crystal':
        pass

    def _build(self) -> None:
        pass
