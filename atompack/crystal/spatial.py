from typing import Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack.crystal.components import LatticeVectors


class MillerIndex(object):
    """Representation of a Miller index for describing crystallographic planes and directions."""

    def __init__(self, hkl: Tuple[int, int, int]) -> None:
        self.hkl = hkl

    @property
    def reciprocal(self) -> np.ndarray:
        """Returns the reciprocal index."""
        _max = max(self.hkl)
        return np.array([index / _max if index > 0 else np.inf for index in self.hkl])


class Orientation(Rotation):
    """Representation of a crystallographic orientation.
    This class inherits from `scipy.spatial.transform.Rotation`.
    """

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_miller_indices(cls, hkl: MillerIndex, uvw: MillerIndex) -> 'Orientation':
        """Initialize from Miller Indices.
        
        Args:
            hkl: Indices of the plane.
            uvw: Indices of the direction.
        """
        raise NotImplementedError

    ########################
    #    Public Methods    #
    ########################

    def as_miller_indices(self) -> Tuple[MillerIndex, MillerIndex]:
        """Represent as Miller Indices."""
        raise NotImplementedError


class Plane(object):
    """Representation of a crystallographic plane."""

    def __init__(self, coords: np.ndarray) -> None:
        self._coords = coords

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_miller_index(cls, miller_index: MillerIndex) -> 'Plane':
        """Initializes from a Miller index."""
        raise NotImplementedError

    @classmethod
    def from_axis_normal(cls, axis: int, length_fraction: float) -> 'Plane':
        """Initializes from an axis and a fractional position along that axis."""
        raise NotImplementedError

    ####################
    #    Properties    #
    ####################

    def centroid(self) -> np.ndarray:
        """Returns the center of the plane."""
        raise NotImplementedError

    ########################
    #    Public Methods    #
    ########################

    def coefficients(self, lattice_vectors: LatticeVectors) -> np.ndarray:
        """Returns the 'A', 'B', 'C', and 'D' coefficients of the equation of the plane."""
        raise NotImplementedError

    def coplanar_points(self, lattice_vectors: LatticeVectors) -> np.ndarray:
        """Returns 3 points which define the plane."""
        raise NotImplementedError

    def intercepts(self, lattice_vectors: LatticeVectors) -> np.ndarray:
        """Returns the 3 points at which the plane intercepts each axis.
        Each point is represented as a scalar since the off-axis components are by definition 0.
        """
        raise NotImplementedError
