from collections.abc import MutableSequence
from typing import List, Tuple, Union

import numpy as np
from scipy.spatial.transform import Rotation

from atompack.constants import DEG90, DEG120
from atompack.spacegroup import Spacegroup
from atompack.topology import Topology


class Basis(MutableSequence):
    """Crystalline basis.
    
    Example:
        >>> from atompack.crystal import Basis
        >>> import numpy as np
        >>>
        >>> # primitive basis of iron atoms
        >>> basis = Basis.primitive("Fe")
        >>> assert len(basis) == 1
        >>>
        >>> specie, site = basis[0]
        >>> assert specie == "Fe"
        >>> assert np.array_equal(site, np.zeros(3))
    """

    def __init__(self, basis: List[Tuple[str, np.ndarray]]) -> None:
        for _, site in basis:
            if max(site) > 1 or min(site) < 0:
                raise ValueError("basis sites must be represented in fractional coordinates")
        self._basis = basis

    ########################################
    #    MutableSequence Implementation    #
    ########################################

    def __getitem__(self, index: int) -> Tuple[str, np.ndarray]:
        return self._basis[index]

    def __setitem__(self, index: int, value: Tuple[str, np.ndarray]) -> None:
        self._basis[index] = value

    def __delitem__(self, index: int) -> None:
        del self._basis[index]

    def __len__(self) -> int:
        return len(self._basis)

    def insert(self, index: int, value) -> None:
        self._basis.insert(index, value)

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def primitive(cls, specie: str) -> 'Basis':
        """Returns a primitive basis."""
        return cls([(specie, np.zeros(3))])


class LatticeParameters(object):
    """Edge lengths and angles which define a lattice.
    
    Args:
        a: Length of the x lattice vector.
        b: Length of the y lattice vector.
        c: Length of the z lattice vector.
        alpha: Angle between the y and z directions (radians).
        beta: Angle between the x and z directions (radians).
        gamma: Angle between the x and y directions (radians).

    Example:
        >>> from atompack.crystal import LatticeParameters
        >>> import numpy as np
        >>>
        >>> # cubic lattice parameters
        >>> params = LatticeParameters.cubic(10)
        >>> assert params.a == params.b == params.c == 10
        >>> assert params.alpha == params.beta == params.gamma == np.pi / 2
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


class Orientation(Rotation):
    """Representation of a crystallographic orientation.
    This class inherits from `scipy.spatial.transform.Rotation`.
    """

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_miller_indices(cls, hkl: Tuple[int, int, int], uvw: Tuple[int, int, int]) -> 'Orientation':
        """Initialize from Miller Indices.
        
        Args:
            hkl: Indices of the plane.
            uvw: Indices of the direction.
        """
        pass

    ########################
    #    Public Methods    #
    ########################

    def as_miller_indices(self) -> Tuple[Tuple[int, int, int], Tuple[int, int, int]]:
        """Represent as Miller Indices."""
        pass


class Plane(object):

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_miller_indices(cls, hkl: Tuple[int, int, int]) -> 'Plane':
        pass


class UnitCell(Topology):
    """Minimal representation of a crystalline structure.

    Note:
        End users should not construct unit cell objects directly.
    
    Args:
        basis: Asymmetric site occupancy.
        lattice_parameters: Lattice parameters object.
        spacegroup: Spacegroup object.
    """

    def __init__(self, basis: Basis, lattice_parameters: LatticeParameters, spacegroup: Spacegroup) -> None:
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


class Crystal(Topology):
    """Atomic structure with long range order.
    
    Args:
        basis: Atomic basis set.
        lattice_parameters: Lattice parameters object.
        spacegroup: Hermann Mauguin spacegroup symbol or international spacegroup number.
    """

    def __init__(self, basis: List[Tuple[str, np.ndarray]], lattice_parameters: LatticeParameters,
                 spacegroup: Union[int, str]) -> None:
        # initialize attributes
        self._cut_plane: Optional[Plane] = None
        self._extent: Optional[Tuple[int, int, int]] = None
        self._orientation: Optional[Orientation] = None
        self._orthogonalize: Optional[bool] = None
        self._projection_plane: Optional[Plane] = None
        self._transformation_matrix: Optional[np.ndarray] = None
        # initialize superclass
        super().__init__()
        # copy attributes from unit cell
        self._unit_cell = UnitCell(basis, lattice_parameters, spacegroup)
        self._lattice_vectors = self._unit_cell._lattice_parameters.metric_tensor.copy()
        self.atoms = self._unit_cell.atoms.copy()

    ####################
    #    Properties    #
    ####################

    @property
    def lattice_vectors(self):
        """Returns the lattice vectors."""
        return self._lattice_vectors

    @property
    def unit_cell(self):
        """Returns the crystal's minimal representation."""
        return self._unit_cell

    ########################
    #    Public Methods    #
    ########################

    def general_transform(self, transformation: np.ndarray) -> 'Crystal':
        """Applies a general transformation matrix to the crystal.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            transformation: Transformation matrix.

        Note:
            This method is idempotent.
            The transform is not applied until the `finish` method is called.
        """
        self._transformation_matrix = transformation
        return self

    def supercell(self, extent: Tuple[int, int, int]) -> 'Crystal':
        """Creates a supercell by duplicating the crystal's unit cell in 3 dimensions.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            extent: Number of repeat units in each direction.

        Note:
            This method is idempotent.
            The transform is not applied until the `finish` method is called.
        """
        self._extent = extent
        return self

    def orient(self, orientation: Orientation) -> 'Crystal':
        """Change the crystal's orientation.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            orientation: Crystallographic orientation.

        Note:
            This method is idempotent.
            The transform is not applied until the `finish` method is called.
        """
        self._orientation = orientation
        return self

    def project(self, plane: Plane, orthogonalize: bool = False) -> 'Crystal':
        """Project the crystal onto a plane.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            plane: Projection plane.
            orthogonalize: Determines whether or not the projection is represented as an orthogonal lattice.

        Note:
            This method is idempotent.
            The transform is not applied until the `finish` method is called.
            Setting `orthogonalize` to True may result in very large structures for acute projections.
        """
        self._projection_plane = plane
        self._orthogonalize = orthogonalize
        return self

    def cut(self, plane: Plane) -> 'Crystal':
        """Cut the crystal along a plane.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            plane: Plane to cut along.

        Note:
            This method is idempotent.
            The transform is not applied until the `finish` method is called.
        """
        self._cut_plane = plane
        return self

    def finish(self) -> None:
        """Applies all active transforms to the crystal."""
        # TODO: call underlying implementations

        # reset attributes
        self._cut_plane = None
        self._extent = None
        self._orientation = None
        self._orthogonalize = None
        self._projection_plane = None
        self._transformation_matrix = None

    #########################
    #    Private Methods    #
    #########################

    def _general_transform(self) -> None:
        pass

    def _supercell(self) -> None:
        pass

    def _orient(self) -> None:
        pass

    def _project(self) -> None:
        pass

    def _cut(self) -> None:
        pass
