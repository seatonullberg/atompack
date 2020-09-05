"""
Module `crystal` provides abstractions to construct crystals.
At the lowest level, a generic lattice and atomic basis with no crystallographic constraints can be composed in a `UnitCell` object.
At the mid-level, bravais lattices with specific symmetry constraints can be constructed with classes such as `Cubic`, `Orthorhombic`, or `Monoclinic`.
At the highest level, users can create a number of common configurations such as `Fcc`, `Rutile`, or `Diamond`.
For transformations such as rotations, reorientations, and supercell duplication, the `Crystal` class wraps a `UnitCell` and modifies it with the desired transforms.
"""
from atompack.crystal.crystal import Crystal
from atompack.crystal.defect import (octahedral_defect, substitution_defect, tetrahedral_defect, vacancy_defect)
from atompack.crystal.unit_cell import (CsCl, Cubic, Diamond, Fcc, Hcp, Hexagonal, Monoclinic, NaCl, Orthorhombic,
                                        Rhombohedral, Rutile, Sc, Tetragonal, Triclinic, UnitCell, Wurtzite, ZincBlend)
from atompack.crystal.util import metric_tensor

__all__ = [
    "Crystal",
    "vacancy_defect",
    "substitution_defect",
    "octahedral_defect",
    "tetrahedral_defect",
    "CsCl",
    "Cubic",
    "Diamond",
    "Fcc",
    "Hcp",
    "Hexagonal",
    "Monoclinic",
    "NaCl",
    "Orthorhombic",
    "Rhombohedral",
    "Rutile",
    "Sc",
    "Tetragonal",
    "Triclinic",
    "UnitCell",
    "Wurtzite",
    "ZincBlend",
    "metric_tensor",
]
