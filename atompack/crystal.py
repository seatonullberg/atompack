from typing import List, Tuple, Union

import numpy as np

from atompack.lattice import LatticeParameters
from atompack.topology import Topology
from atompack.unit_cell import UnitCell


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

    @property
    def lattice_vectors(self):
        """Returns the lattice vectors."""
        return self._lattice_vectors

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

    def _build(self) -> None:
        pass


class Orientation(object):
    pass


class Plane(object):
    pass
