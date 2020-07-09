from typing import List, Optional, Tuple, Dict, Any

import numpy as np

from atompack.atom import Atom
from atompack.structure import Structure


class Crystal(Structure):
    """Representation of a crystalline lattice.
    
    Note:
        It is a logical error to mutate any of these attributes because the result of any changes will not be reflected in the shape of the basis or the content of the atoms list.

    Args:
        lattice_data: Atomic data used to initialize atoms at each lattice site.
        lattice_sites: Fractional coordinates of atoms in the lattice.
        a: Length of the x direction.
        b: Length of the y direction.
        c: Length of the z direction.
        alpha: Angle between y and z directions in radians.
        beta: Angle between x and z directions in radians.
        gamma: Angle between x and y directions in radians.
        duplicates: Number of duplications to apply to the finished structure along each direction.
        orientation: 3x3 matrix indicating the alignment of the lattice vectors.
    """

    def __init__(
            self,
            lattice_data: List[Dict[str, Any]],
            lattice_sites: np.ndarray,
            a: float,
            b: float,
            c: float,
            alpha: float,
            beta: float,
            gamma: float,
            duplicates: Optional[Tuple[int, int, int]] = None,
            orientation: Optional[np.ndarray] = None,
            pbc: Optional[Tuple[bool, bool, bool]] = None,
            tolerance: Optional[float] = None,
    ) -> None:
        self.lattice_data = lattice_data
        self.lattice_sites = lattice_sites
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self.duplicates = duplicates
        self.orientation = orientation
        atoms, basis = self._build()
        super().__init__(atoms, basis, pbc, tolerance)

    def _build(self) -> Tuple[List[Atom], np.ndarray]:
        raise NotImplementedError