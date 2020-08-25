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
        """Initializes a molecule with a pair/linear configuration."""
        pass

    @classmethod
    def linear3(cls, core_atom: Atom, bonded_atoms: List[Atom], bond_lengths: List[float]) -> 'Molecule':
        """Initializes a molecule with a triplet/linear configuration."""
        pass

    @classmethod
    def trigonal_planar(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
    ) -> 'Molecule':
        """Initializes a molecule with a trigonal planar configuration."""
        pass

    @classmethod
    def bent(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
        theta: float,
    ) -> 'Molecule':
        """Initializes a molecule with a bent/angular configuration."""
        pass

    @classmethod
    def tetrahedral(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[Atom],
    ) -> 'Molecule':
        """Initializes a molecule with a tetrahedral configuration."""
        pass

    @classmethod
    def trigonal_pyramidal(
        cls,
        core_atom: Atom,
        bonded_atoms: List[Atom],
        bond_lengths: List[float],
        theta: float,
    ) -> 'Molecule':
        """Initializes a molecule with a trigonal pyramidal configuration."""
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
        """Initializes a molecule with a trigonal bipyramidal configuration."""
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
        """Initializes a molecule with an octahedral configuration."""
        pass
