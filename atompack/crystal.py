from collections.abc import MutableSequence
from typing import List, Tuple, Union

import numpy as np

from atompack.topology import Topology

DEG90 = np.pi / 2       # 90 degrees in radians
DEG120 = 2 * np.pi / 3  # 120 degrees in radians


class Basis(MutableSequence):
    """Crystalline basis."""

    def __init__(self, sites: List[Tuple[str, np.ndarray]]):
        self._sites = sites

    ########################################
    #    MutableSequence Implementation    #
    ########################################

    def __getitem__(self, index: int) -> Tuple[str, np.ndarray]:
        return self._sites[index]

    def __setitem__(self, index: int, value: Tuple[str, np.ndarray]) -> None:
        self._sites[index] = value

    def __delitem__(self, index: int) -> None:
        del self._sites[index]

    def __len__(self) -> int:
        return len(self._sites)

    def insert(self, index: int, value) -> None:
        self._sites.insert(index, value)

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def primitive(cls, specie: str) -> 'Basis':
        return cls([(specie, np.zeros(3))])


class Crystal(Topology):
    """Atomic structure with long range order.
    
    Args:
        basis: Asymmetric site occupancy.
        lattice_parameters: Lattice parameters object.
        spacegroup: Hermann Mauguin spacegroup symbol or international spacegroup number.
    """

    def __init__(self, basis: List[Tuple[str, np.ndarray]], lattice_parameters: LatticeParameters,
                 spacegroup: Union[int, str]) -> None:
        self._unit_cell = UnitCell(basis, lattice_parameters, spacegroup)
        self._lattice_vectors = self._unit_cell._lattice_parameters.metric_tensor
        self._extent: Optional[Tuple[int, int, int]] = None
        self._orientation: Optional[Orientation] = None
        self._projection_plane: Optional[Plane] = None
        self._orthogonalize: Optional[bool] = None
        self._cut_plane: Optional[Plane] = None
        super().__init__()
        self._build()

    ####################
    #    Properties    #
    ####################

    @property
    def lattice_vectors(self):
        """Returns the lattice vectors."""
        return self._lattice_vectors

    ########################
    #    Public Methods    #
    ########################

    def duplicate(self, extent: Tuple[int, int, int]) -> 'Crystal':
        """Duplicate a unit cell in 3 dimensions."""
        pass

    def orient(self, orientation: Orientation) -> 'Crystal':
        """Rotate the crystal."""
        pass

    def project(self, plane: Plane) -> 'Crystal':
        """Project the crystal onto a plane."""
        pass

    def orthogonalize(self) -> 'Crystal':
        """Transforms the crystal into an orthogonal representation.
        
        Note:
            This may produce very large structures for acute lattices.
        """
        pass

    def cut(self, plane: Plane) -> 'Crystal':
        """Cut the crystal in the direction normal to a plane."""
        pass

    #########################
    #    Private Methods    #
    #########################

    def _build(self) -> None:
        pass


class LatticeParameters(object):
    """Edge lengths and angles which define a lattice.
    
    Args:
        a: Length of the x lattice vector.
        b: Length of the y lattice vector.
        c: Length of the z lattice vector.
        alpha: Angle between the y and z directions (radians).
        beta: Angle between the x and z directions (radians).
        gamma: Angle between the x and y directions (radians).
    """

    def __init__(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def triclinic(cls, a, b, c, alpha, beta, gamma) -> 'LatticeParameters':
        return cls(a, b, c, alpha, beta, gamma)

    @classmethod
    def monoclinic(cls, a, b, c, beta) -> 'LatticeParameters':
        return cls(a, b, c, DEG90, beta, DEG90)

    @classmethod
    def orthorhombic(cls, a, b, c) -> 'LatticeParameters':
        return cls(a, b, c, DEG90, DEG90, DEG90)

    @classmethod
    def tetragonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG90)

    @classmethod
    def trigonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG120)

    @classmethod
    def rhombohedral(cls, a, alpha) -> 'LatticeParameters':
        return cls(a, a, a, alpha, alpha, alpha)

    @classmethod
    def hexagonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG120)

    @classmethod
    def cubic(cls, a) -> 'LatticeParameters':
        return cls(a, a, a, DEG90, DEG90, DEG90)

    ####################
    #    Properties    #
    ####################

    @property
    def metric_tensor(self) -> np.ndarray:
        """Returns the metric tensor defined by the lattice parameters."""
        return np.array([[self.a * self.a, self.a * self.b * np.cos(self.gamma), self.a * self.c * np.cos(self.beta)],
                         [self.a * self.b * np.cos(self.gamma), self.b * self.b, self.b * self.c * np.cos(self.alpha)],
                         [self.a * self.c * np.cos(self.beta), self.b * self.c * np.cos(self.alpha), self.c * self.c]])


class Orientation(object):
    pass


class Plane(object):
    pass


class UnitCell(Topology):
    """Minimal representation of a crystalline structure.

    Note:
        End users should not construct unit cell objects directly.
    
    Args:
        basis: Asymmetric site occupancy.
        lattice_parameters: Lattice parameters object.
        spacegroup: Hermann Mauguin spacegroup symbol or international spacegroup number.
    """

    def __init__(self, basis: Basis, lattice_parameters: LatticeParameters,
                 spacegroup: Union[int, str]) -> None:
        self._basis = basis
        self._lattice_parameters = lattice_parameters
        self._spacegroup = spacegroup
        super().__init__()
        self._build()

    ####################
    #    Properties    #
    ####################

    @property
    def basis(self) -> Basis:
        return self._basis

    @property
    def lattice_parameters(self) -> LatticeParameters:
        return self._lattice_parameters

    @property
    def spacegroup_number(self) -> int:
        pass

    @property
    def spacegroup_symbol(self) -> str:
        pass


    #########################
    #    Private Methods    #
    #########################

    def _build(self) -> None:
        pass
