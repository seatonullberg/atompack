from typing import List, Tuple

import numpy as np

from atompack.atom import Atom
from atompack.topology import Topology


class Molecule(Topology):
    """Representation of a molecule."""

    def __init__(self, atoms: List[Atom], bonds: List[Tuple[int, int]]) -> None:
        super().__init__()
        for atom in atoms:
            self.insert(atom)
        for a, b in bonds:
            self.connect(a, b)

    @classmethod
    def linear2(cls, atoms: List[Atom], bond_length: float) -> 'Molecule':
        """Initializes a molecule with an idealized pair/linear configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        # atoms[0].position = np.zeros(3)
        # atoms[1].position = np.array([bond_length, 0, 0])
        # bonds = [(0, 1)]
        # return cls(atoms, bonds)

    @classmethod
    def linear3(cls, core_atom: Atom, bonded_atoms: List[Atom], bond_lengths: List[float]) -> 'Molecule':
        """Initializes a molecule with an idealized triplet/linear configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        # bonded_atoms[0]._position = np.zeros(3)
        # core_atom._position = np.array([bond_lengths[0], 0, 0])
        # bonded_atoms[1]._position = np.array([bond_lengths[0] + bond_lengths[1], 0, 0])
        # atoms = bonded_atoms + core_atom
        # bonds = [(0, 2), (1, 2)]
        # return cls(atoms, bonds)

    @classmethod
    def trigonal_planar(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
    ) -> 'Molecule':
        """Initializes a molecule with an idealized trigonal planar configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass

    @classmethod
    def bent(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
        theta: float,
    ) -> 'Molecule':
        """Initializes a molecule with an idealized bent/angular configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass

    @classmethod
    def tetrahedral(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
    ) -> 'Molecule':
        """Initializes a molecule with an idealized tetrahedral configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass

    @classmethod
    def trigonal_pyramidal(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
        theta: float,
    ) -> 'Molecule':
        """Initializes a molecule with an idealized trigonal pyramidal configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass

    @classmethod
    def trigonal_bipyramidal(
        cls,
        core_atom: Atom,
        axial_atoms: List[Atom],
        lateral_atoms: List[Atom],
        axial_bond_lengths: List[float],
        lateral_bond_lengths: List[float],
    ) -> 'Molecule':
        """Initializes a molecule with an idealized trigonal bipyramidal configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass

    @classmethod
    def octahedral(
        cls,
        core_atom: Atom,
        axial_atoms: List[Atom],
        lateral_atoms: List[Atom],
        axial_bond_lengths: List[float],
        lateral_bond_lengths: List[float],
    ) -> 'Molecule':
        """Initializes a molecule with an idealized octahedral configuration.
        
        Note:
            The `position` attribute of each atom is ignored and overwritten.
        """
        pass
