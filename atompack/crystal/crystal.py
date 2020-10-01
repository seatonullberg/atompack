import copy
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack.atom import Atom
from atompack.crystal.unit_cell import UnitCell
from atompack.crystal.util import enforce_bounds
from atompack.topology import Topology


class Crystal(Topology):
    """Representation of a generic crystal."""

    def __init__(self,
                 unit_cell: UnitCell,
                 scale: Optional[Tuple[int, int, int]] = None,
                 orientation: Optional[np.ndarray] = None,
                 rotation: Optional[np.ndarray] = None,
                 tol: float = 1.0e-6) -> None:
        self._unit_cell = unit_cell
        if scale is None:
            scale = (1, 1, 1)
        self._scale = scale
        if orientation is None:
            orientation = np.identity(3)
        self._orientation = orientation
        if rotation is None:
            rotation = np.identity(3)
        self._rotation = rotation
        self._tol = tol
        atoms, self._lattice_vectors = self._build()
        super().__init__()
        for atom in atoms:
            self.insert(atom)

    @property
    def unit_cell(self) -> UnitCell:
        """Returns a copy of the crystal's unit cell."""
        return copy.deepcopy(self._unit_cell)

    @property
    def scale(self) -> Tuple[int, int, int]:
        """Returns a copy of the crystal's 3D scale factor."""
        return copy.deepcopy(self._scale)

    @property
    def orientation(self) -> np.ndarray:
        """Returns a copy of the crystal's orientation matrix."""
        return copy.deepcopy(self._orientation)

    @property
    def rotation(self) -> np.ndarray:
        """Returns a copy of the crystal's rotation matrix."""
        return copy.deepcopy(self._rotation)

    @property
    def lattice_vectors(self) -> np.ndarray:
        """Returns a copy of the crystal's lattice vectors."""
        return copy.deepcopy(self._lattice_vectors)

    def _build(self) -> Tuple[List[Atom], np.ndarray]:
        # transforms are applied in the following order:
        # - orientation
        # - rotation
        # - scale

        # calculate the magnitude of the lattice vectors
        lattice_vector_mags = np.linalg.norm(self._unit_cell.lattice_vectors, axis=0)

        # calculate the unit vector of each lattice vector
        lattice_unit_vectors = self._unit_cell.lattice_vectors / lattice_vector_mags

        # calculate the rotation matrix between the unoriented and oriented lattice vectors
        rotation = Rotation.align_vectors(lattice_unit_vectors, self._orientation)[0]

        # align the lattice vectors with the orientation
        oriented_lattice_vectors = np.matmul(self._orientation, self._unit_cell.lattice_vectors)

        # use QR decomposition to calculate an orthogonal representation
        # TODO: this should be optional
        _, r = np.linalg.qr(oriented_lattice_vectors.T)
        oriented_lattice_vectors = np.abs(r)  # removed multiplication by scale

        # calculate the magnitude of the oriented lattice vectors
        oriented_lattice_vector_mags = np.linalg.norm(oriented_lattice_vectors, axis=0)

        # determine smallest orthogonal size
        min_ortho_size = np.ceil(oriented_lattice_vector_mags / lattice_vector_mags)
        min_ortho_size = min_ortho_size.astype(int)

        # place the atoms
        atoms: List[Atom] = []
        for xsize in range(min_ortho_size[0]):
            for ysize in range(min_ortho_size[1]):
                for zsize in range(min_ortho_size[2]):
                    offset = np.matmul(np.array([xsize, ysize, zsize]), self._unit_cell.lattice_vectors)
                    for atom in self._unit_cell.atoms:

                        # calculate the cartesian position
                        position = atom.position + offset
                        position = rotation.apply(position)

                        # transform the position into the lattice
                        enforce_bounds(oriented_lattice_vector_mags, position, self._tol)

                        # accept the atom if the position is not yet occupied
                        positions = np.array([atom.position for atom in atoms])
                        is_occupied = False
                        for _position in positions:
                            if np.linalg.norm(position - _position) < self._tol:
                                is_occupied = True
                                break
                        if not is_occupied:
                            atom = copy.deepcopy(atom)
                            atom._position = position
                            atoms.append(atom)

        # TODO: apply a rotation matrix

        # tile the crystal in all directions
        current_atoms = copy.deepcopy(atoms)
        for xsize in range(self._scale[0]):
            for ysize in range(self._scale[1]):
                for zsize in range(self._scale[2]):
                    offset = np.matmul(np.array([xsize, ysize, zsize]), oriented_lattice_vectors)
                    if np.linalg.norm(offset) < self._tol:
                        continue
                    for atom in current_atoms:
                        atom = copy.deepcopy(atom)
                        atom._position += offset
                        atoms.append(atom)

        # multiply the oriented lattice vectors by the scale
        scaled_oriented_lattice_vectors = oriented_lattice_vectors * np.array(self._scale)

        return atoms, scaled_oriented_lattice_vectors
