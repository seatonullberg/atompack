from typing import List, Optional, Tuple

import numpy as np

from atompack.atom import Atom
from atompack.structure import Structure


class Molecule(Structure):
    """Representation of a molecule.
    
    Note:
        It is a logical error to mutate any of these attributes because the result 
        of any changes will not be reflected in the shape of the basis or the 
        content of the atoms list.
    """

    def __init__(
        self,
        atoms: List[Atom],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> None:
        pass

    @classmethod
    def linear2(
        cls,
        atoms: List[Atom],
        bond_distance: float,
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def linear3(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def trigonal_planar(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def bent(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        theta: float,
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def tetrahedral(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def trigonal_pyramidal(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        theta: float,
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def trigonal_bipyramidal(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass

    @classmethod
    def octahedral(
        cls,
        atoms: List[Atom],
        bond_distances: List[Tuple[str, str]],
        basis: Optional[np.ndarray],
        pbc: Optional[Tuple[bool, bool, bool]],
        tolerance: float = 1.0e-6,
    ) -> Molecule:
        pass