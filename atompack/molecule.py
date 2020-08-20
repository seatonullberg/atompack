from typing import Dict, List, Optional, Tuple

import numpy as np
from _cell import cell_contains

from atompack.atom import Atom
from atompack.errors import PositionOccupiedError, PositionOutsideError
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
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> None:
        super().__init__(atoms, basis, pbc, tolerance)

    @classmethod
    def linear2(
        cls,
        atoms: List[Atom],
        bond_distance: float,
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a pair/linear configuration."""
        # generate positions
        if bond_distance < tolerance:
            raise PositionOccupiedError(np.zeros(3))
        atoms[0].position = np.zeros(3)
        atoms[1].position = np.array([bond_distance, 0.0, 0.0])
        # validate basis
        if basis is None:
            # generate a minimally sized cubic basis
            basis = np.array([bond_distance, bond_distance, bond_distance])
        else:
            # ensure all positions fit in the basis
            for atom in atoms:
                if not cell_contains(basis, atom.position, tolerance):
                    raise PositionOutsideError(atom.position, basis)
        return cls(atoms, basis, pbc, tolerance)

    @classmethod
    def linear3(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a triplet/linear configuration."""
        pass

    @classmethod
    def trigonal_planar(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a trigonal planar configuration."""
        pass

    @classmethod
    def bent(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        theta: float,
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a bent configuration."""
        pass

    @classmethod
    def tetrahedral(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a tetrahedral configuration."""
        pass

    @classmethod
    def trigonal_pyramidal(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        theta: float,
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a trigonal pyrimidal configuration."""
        pass

    @classmethod
    def trigonal_bipyramidal(
        cls,
        core_atom: Atom,
        axial_atoms: List[Atom],
        lateral_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with a trigonal bipyramidal configuration."""
        pass

    @classmethod
    def octahedral(
        cls,
        core_atom: Atom,
        axial_atoms: List[Atom],
        lateral_atoms: List[Atom],
        bond_distances: Dict[Tuple[str, str], float],
        basis: Optional[np.ndarray] = None,
        pbc: Optional[Tuple[bool, bool, bool]] = None,
        tolerance: float = 1.0e-6,
    ) -> 'Molecule':
        """Initializes a molecule with an octahedral configuration."""
        pass
