import copy
from typing import List, Optional, Tuple

import numpy as np

from atompack.atom import Atom, AtomCollection
from atompack.error import (InvalidMonoclinicError, InvalidOrthorhombicError, InvalidTriclinicError)
from atompack.internal import (is_point_in_polyhedron, metric_tensor, rotation_matrix_from_vectors, search_for_atom)


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
    def triclinic(cls,
                  a: float,
                  b: float,
                  c: float,
                  alpha: float,
                  beta: float,
                  gamma: float,
                  unit_cell: List[Tuple[Atom, np.ndarray]],
                  orientation: Optional[np.ndarray] = None,
                  rotation: Optional[np.ndarray] = None,
                  size: Optional[Tuple[int, int, int]] = None) -> 'Crystal':
        """Initializes a `Crystal` with triclinic constraints."""
        if a == b or a == c or b == c:
            raise InvalidTriclinicError(a, b, c, alpha, beta, gamma)
        if alpha == beta or alpha == gamma or beta == gamma:
            raise InvalidTriclinicError(a, b, c, alpha, beta, gamma)
        return cls(a, b, c, alpha, beta, gamma, unit_cell, orientation, rotation, size)

    @classmethod
    def monoclinic(cls,
                   a: float,
                   b: float,
                   c: float,
                   beta: float,
                   unit_cell: List[Tuple[Atom, np.ndarray]],
                   orientation: Optional[np.ndarray] = None,
                   rotation: Optional[np.ndarray] = None,
                   size: Optional[Tuple[int, int, int]] = None) -> 'Crystal':
        """Initializes a `Crystal` with monoclinic constraints."""
        alpha, gamma = np.pi / 2, np.pi / 2
        if beta == np.pi / 2:
            raise InvalidMonoclinicError(a, b, c, alpha, beta, gamma)
        return cls(a, b, c, alpha, beta, gamma, unit_cell, orientation, rotation, size)

    @classmethod
    def orthorhombic(cls,
                     a: float,
                     b: float,
                     c: float,
                     unit_cell: List[Tuple[Atom, np.ndarray]],
                     orientation: Optional[np.ndarray] = None,
                     rotation: Optional[np.ndarray] = None,
                     size: Optional[Tuple[int, int, int]] = None) -> 'Crystal':
        """Initializes a `Crystal` with orthorhombic constraints."""
        alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
        if a == b or a == c or b == c:
            raise InvalidOrthorhombicError(a, b, c, alpha, beta, gamma)
        return cls(a, b, c, alpha, beta, gamma, unit_cell, orientation, rotation, size)

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

    # TODO: implement behavior for rotation
    def _build(self) -> Tuple[List[Atom], np.ndarray]:
        # process attributes
        if self._unit_cell is None:
            self._unit_cell = [(Atom(), np.zeros((3,)))]
        if self._orientation is None:
            self._orientation = np.identity(3)
        if self._rotation is None:
            self._rotation = np.identity(3)
        else:
            raise NotImplementedError()
        if self._size is None:
            self._size = np.array([1, 1, 1])

        # generate lattice vectors from metric tensor
        lattice_vectors = np.sqrt(metric_tensor(self._a, self._b, self._c, self._alpha, self._beta, self._gamma))

        # align lattice vectors with orientation
        oriented_lattice_vectors = np.matmul(self._orientation, lattice_vectors)

        # use QR decomposition to determine an orthogonal representation of the lattice vectors
        _, r = np.linalg.qr(oriented_lattice_vectors.T)
        q, _ = np.linalg.qr(oriented_lattice_vectors)
        oriented_lattice_vectors = np.abs(r) * self._size
        oriented_lattice_magnitudes = np.linalg.norm(oriented_lattice_vectors, axis=0)

        # determine the smallest required orthogonal cell
        minimum_orthogonal_size = (np.ceil(np.linalg.norm(r, axis=0)) * self._size).astype(int)

        # place atoms on the oriented lattice vectors
        atoms = []
        for xsize in range(minimum_orthogonal_size[0]):
            for ysize in range(minimum_orthogonal_size[1]):
                for zsize in range(minimum_orthogonal_size[2]):
                    offset = np.matmul(q, np.array([xsize, ysize, zsize]))
                    for atom, relative_position in self._unit_cell:
                        # assign cartesian position
                        position = np.abs(np.matmul(q, relative_position) + offset)
                        # reduce periodic images
                        reduced_position = copy.deepcopy(position)
                        for i in range(3):
                            diff = position[i] - oriented_lattice_magnitudes[i]
                            if diff >= -1e-6:
                                reduced_position[i] = diff
                        res = search_for_atom(atoms, reduced_position, 1e-6)
                        if res is not None:
                            continue
                        # ensure point is in the lattice
                        if is_point_in_polyhedron(reduced_position, oriented_lattice_vectors):
                            new_atom = copy.deepcopy(atom)
                            new_atom.position = reduced_position
                            atoms.append(new_atom)
        return atoms, oriented_lattice_vectors
