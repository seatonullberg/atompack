from typing import Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack.crystal.components import LatticeVectors


class MillerIndex(object):
    """Representation of a Miller index for describing crystallographic planes and directions."""

    def __init__(self, hkl: Tuple[int, int, int]) -> None:
        self.hkl = hkl

    ####################
    #    Properties    #
    ####################

    @property
    def reciprocal(self) -> np.ndarray:
        """Returns the reciprocal index."""
        _max = max(self.hkl)
        return np.array([index / _max if index > 0 else np.inf for index in self.hkl])

    #########################
    #    Special Methods    #
    #########################

    def __eq__(self, other) -> bool:
        if not isinstance(other, MillerIndex):
            return NotImplemented
        return self.hkl == other.hkl


class Orientation(Rotation):
    """Representation of a crystallographic orientation.
    This class inherits from scipy's Rotation class for 
    efficient conversion between possible representations.
    """

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_miller_indices(cls, plane: MillerIndex, direction: MillerIndex) -> 'Orientation':
        """Initialize from Miller Indices.
        
        Args:
            plane: Indices of the plane.
            direction: Indices of the direction.
        """
        hkl = np.array(plane.hkl)
        uvw = np.array(direction.hkl)
        b_hat = uvw / np.linalg.norm(uvw)
        n_hat = hkl / np.linalg.norm(hkl)
        n_cross_b = np.cross(n_hat, b_hat)
        t_hat = n_cross_b / np.linalg.norm(n_cross_b)
        matrix = np.vstack((b_hat, t_hat))
        matrix = np.vstack((matrix, n_hat))
        return super().from_rotvec(matrix.T)

    ########################
    #    Public Methods    #
    ########################

    def as_miller_indices(self, tol: float = 1E-6) -> Tuple[MillerIndex, MillerIndex]:
        """Represent as Miller Indices."""
        matrix = self.as_rotvec()
        hkl = matrix[:,2]
        uvw = matrix[:,0]
        min_nonzero = lambda arr: np.min(arr[np.abs(arr) > tol])
        normalize = lambda arr: np.array([x / min_nonzero(arr) for x in arr])
        hkl = tuple(np.round(normalize(hkl)).astype(int))
        uvw = tuple(np.round(normalize(uvw)).astype(int))
        return MillerIndex(hkl), MillerIndex(uvw)


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
