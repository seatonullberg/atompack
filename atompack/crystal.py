from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from _pbc import pbc_nearest_neighbor
from atompack.atom import Atom
from atompack.structure import Structure


def metric_tensor(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    return np.array([[a * a, a * b * np.cos(gamma), a * c * np.cos(beta)],
                     [a * b * np.cos(gamma), b * b, b * c * np.cos(alpha)],
                     [a * c * np.cos(beta), b * c * np.cos(alpha), c * c]])


class Crystal(Structure):
    """Representation of a crystalline lattice.
    
    Note:
        It is a logical error to mutate any of these attributes because the result of any changes will not be reflected in the shape of the basis or the content of the atoms list.

    Args:
        lattice_data: Atomic data used to initialize atoms at each lattice site.
        lattice_sites: Fractional coordinates of atoms in the lattice.
        a: Length of the x direction.
        b: Length of the y direction.
        c: Length of the z direction.
        alpha: Angle between y and z directions in radians.
        beta: Angle between x and z directions in radians.
        gamma: Angle between x and y directions in radians.
        duplicates: Number of duplications to apply to the finished structure along each direction.
        orientation: 3x3 matrix indicating the alignment of the lattice vectors.
    """

    def __init__(
            self,
            lattice_data: List[Dict[str, Any]],
            lattice_sites: np.ndarray,
            a: float,
            b: float,
            c: float,
            alpha: float,
            beta: float,
            gamma: float,
            duplicates: Optional[Tuple[int, int, int]] = None,
            orientation: Optional[np.ndarray] = None,
            pbc: Optional[Tuple[bool, bool, bool]] = None,
            tolerance: float = 1.0e-6,
    ) -> None:
        if duplicates is None:
            duplicates = (1, 1, 1)
        if orientation is None:
            orientation = np.identity(3)
        if pbc is None:
            pbc = (False, False, False)
        atoms, basis = self._build(lattice_data, lattice_sites, a, b, c, alpha, beta, gamma, duplicates, orientation,
                                   pbc, tolerance)
        super().__init__(atoms, basis, pbc, tolerance)

    # TODO: Make a C extension to do this.
    @staticmethod
    def _build(
            lattice_data: List[Dict[str, Any]],
            lattice_sites: np.ndarray,
            a: float,
            b: float,
            c: float,
            alpha: float,
            beta: float,
            gamma: float,
            duplicates: Tuple[int, int, int],
            orientation: np.ndarray,
            pbc: Tuple[bool, bool, bool],
            tolerance: float,
    ) -> Tuple[List[Atom], np.ndarray]:
        """Construct a crystalline lattice and populate it with atoms."""
        # calculate lattice vectors from metric tensor
        lattice_vectors = np.sqrt(metric_tensor(a, b, c, alpha, beta, gamma))

        # calculate the magnitude of each lattice vector
        lattice_vectors_magnitude = np.linalg.norm(lattice_vectors, axis=0)

        # calculate the unit vector of each lattice vector
        lattice_vectors_hat = lattice_vectors / lattice_vectors_magnitude

        # calculate rotation matrix between natural orientation and desired orientation
        rotation = Rotation.align_vectors(lattice_vectors_hat, orientation)[0]

        # align the lattice vectors with the desired orientation by applying the
        lattice_vectors_oriented = np.matmul(orientation, lattice_vectors)

        # use QR decomposition to calculate an orthogonal representation of the oriented lattice vectors
        _, r = np.linalg.qr(lattice_vectors_oriented.T)
        lattice_vectors_oriented = np.abs(r) * np.array(duplicates)
        lattice_vectors_oriented_magnitude = np.linalg.norm(lattice_vectors_oriented, axis=0)

        # determine smallest possible size
        minimum_orthogonal_size = np.ceil(
            lattice_vectors_oriented_magnitude / lattice_vectors_magnitude) * np.array(duplicates)
        minimum_orthogonal_size = minimum_orthogonal_size.astype(int)

        atoms: List[Atom] = []
        for xsize in range(minimum_orthogonal_size[0]):
            for ysize in range(minimum_orthogonal_size[1]):
                for zsize in range(minimum_orthogonal_size[2]):
                    offset = np.matmul(np.array([xsize, ysize, zsize]), lattice_vectors)
                    for data, site in zip(lattice_data, lattice_sites):

                        # assign the cartesian position
                        position = np.matmul(site, lattice_vectors) + offset
                        position = rotation.apply(position)

                        # transform the position back into the bounding box
                        for i in range(3):
                            if position[i] < -tolerance:
                                new_position = position[i]
                                while new_position < -tolerance:
                                    new_position += lattice_vectors_oriented_magnitude[i]
                                position[i] = new_position
                            elif position[i] >= lattice_vectors_oriented_magnitude[i] - tolerance:
                                new_position = position[i]
                                while new_position >= lattice_vectors_oriented_magnitude[i] - tolerance:
                                    new_position -= lattice_vectors_oriented_magnitude[i]
                                position[i] = new_position

                        # accept the atom if the position is not occupied
                        positions = np.array([atom.position for atom in atoms])
                        distance, _ = pbc_nearest_neighbor(position, positions, lattice_vectors_oriented, pbc)
                        if distance > tolerance:
                            atom = Atom(position=position, **data)
                            atoms.append(atom)

        return atoms, lattice_vectors_oriented
