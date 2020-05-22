from atompack.atom import Atom, AtomCollection
from atompack.internal import metric_tensor, is_point_in_polyhedron

import numpy as np
from typing import List, Optional, Tuple


class Crystal(AtomCollection):
    """A crystalline lattice.
    
    Args:
        a: The \\(a\\) distance lattice parameter.
        b: The \\(b\\) distance lattice parameter.
        c: The \\(c\\) distance lattice parameter.
        alpha: The \\(\\alpha\\) angle lattice parameter.
        beta: The \\(\\beta\\) angle lattice parameter.
        gamma: The \\(\\gamma\\) angle lattice parameter.
    """

    def __init__(self,
                 a: float = 0,
                 b: float = 0,
                 c: float = 0,
                 alpha: float = 0,
                 beta: float = 0,
                 gamma: float = 0,
                 unit_cell: Optional[List[Tuple[Atom, np.ndarray]]] = None,
                 orientation: Optional[np.ndarray] = None,
                 rotation: Optional[np.ndarray] = None,
                 size: Optional[Tuple[int, int, int]] = None) -> None:
        """Initializes a new `Crystal`."""
        self._a, self._b, self._c = a, b, c
        self._alpha, self._beta, self._gamma = alpha, beta, gamma
        self._unit_cell = unit_cell
        self._orientation = orientation
        self._rotation = rotation
        self._size = size
        atoms, basis = self._build()
        super().__init__(atoms=atoms, basis=basis)

    @classmethod
    def triclinic(cls) -> 'Crystal':
        pass

    @classmethod
    def monoclinic(cls) -> 'Crystal':
        pass

    @classmethod
    def orthorhombic(cls) -> 'Crystal':
        pass

    @classmethod
    def tetragonal(cls) -> 'Crystal':
        pass

    @classmethod
    def rhombohedral(cls) -> 'Crystal':
        pass

    @classmethod
    def hexagonal(cls) -> 'Crystal':
        pass

    @classmethod
    def cubic(cls) -> 'Crystal':
        pass

    # TODO: Fix the `rotation_matrix`.
    def _build(self) -> None:
        # generate the basis set in its natural orientation from the metric tensor
        lattice_vectors = np.sqrt(
            metric_tensor(self._a, self._b, self._c, self._alpha, self._beta,
                          self._gamma))

        # process the size argument
        if self._size is None:
            self._size = (0, 0, 0)

        # use QR decomposition to find the appropriate basis if another orientation is provided
        oriented_vectors = np.zeros((3, 3))
        if self._orientation is None:
            oriented_vectors = lattice_vectors
        else:
            r = np.abs(np.linalg.qr(self._orientation, mode="r"))
            for i in range(3):
                for j in range(3):
                    r[i][j] *= np.linalg.norm(lattice_vectors[i])
            oriented_vectors = r

        # generate a rotation matrix to align the natural basis with the new orientation
        # https://stackoverflow.com/questions/50474886/finding-rotation-matrix-to-transform-one-3-vector-basis-to-another-in-3d
        rotation_matrix = np.linalg.solve(lattice_vectors, oriented_vectors).T

        # scale the basis such that it is (just) larger than necessary to accommodate the new orientation
        oriented_size = [0, 0, 0]
        for i in range(3):
            oriented_size[i] = np.ceil(
                np.linalg.norm(oriented_vectors[i]) /
                np.linalg.norm(lattice_vectors[i])) * self._size[i]
        oriented_size = [int(x) for x in oriented_size]

        # scale the oriented vectors by the size
        oriented_vectors *= np.array(self._size)

        # place the atoms
        atoms = []
        for x_size in range(oriented_size[0]):
            for y_size in range(oriented_size[1]):
                for z_size in range(oriented_size[2]):
                    offset = np.array([
                        np.linalg.norm(lattice_vectors[0]) * x_size,
                        np.linalg.norm(lattice_vectors[1]) * y_size,
                        np.linalg.norm(lattice_vectors[2]) * z_size,
                    ])

                    # replicate the unit cell
                    for atom, relative_position in self._unit_cell:
                        position = np.zeros((3,))
                        for i in range(3):
                            position[i] = (np.linalg.norm(lattice_vectors[i]) *
                                           relative_position[i]) + offset[i]
                        print("position: {}".format(position))
                        position = np.dot(rotation_matrix, position)
                        print("position rotated: {}\n".format(position))

                        # only add the atom to the collection if it falls within the oriented basis
                        if is_point_in_polyhedron(position, oriented_vectors):
                            atom.position = position
                            atoms.append(atom)

        # TODO: apply a rotation to the total collection

        return atoms, oriented_vectors
