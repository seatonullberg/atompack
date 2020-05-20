from atompack.atom import Atom, AtomCollection
from atompack.internal import metric_tensor, is_point_in_polyhedron

import numpy as np
from typing import List, Optional, Tuple


class Crystal(AtomCollection):
    """A crystalline lattice."""

    # TODO: Call _build to place atoms and set basis
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
        self._a, self._b, self._c = a, b, c
        self._alpha, self._beta, self._gamma = alpha, beta, gamma
        self._unit_cell = unit_cell
        self._orientation = orientation
        self._rotation = rotation
        self._size = size
        # if unit_cell is None:
        #     unit_cell = []
        # self._unit_cell = unit_cell
        # if orientation is None:
        #     orientation = np.array([])
        # self._orientation = orientation
        # if rotation is None:
        #     rotation = np.array([])
        # self._rotation = rotation
        # if size is None:
        #     size = (0, 0, 0)
        # self._size = size
        self._build()
        super().__init__()

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

    # def _build(self) -> None:
    #     lattice_vectors = np.sqrt(metric_tensor(self._a, self._b, self._c, self._alpha, self._beta, self._gamma))
    #     if self._size is None:
    #         self._size = (0, 0, 0)
    #     oriented_vectors = np.zeros((3, 3))
    #     if self._orientation is None:
    #         oriented_vectors = lattice_vectors
    #     else:
    #         _, r = np.linalg.qr(self._orientation)
    #         for i in range(3):
    #             for j in range(3):
    #                 r[i][j] *= np.linalg.norm(lattice_vectors[i])
    #         oriented_vectors = r
    #     oriented_size = (0, 0, 0)
    #     for i in range(3):
    #         oriented_size[i] = np.linalg.norm(oriented_vectors[i]) / np.ceil(np.linalg.norm(lattice_vectors[i]))
    #         oriented_size[i] *= self._size[i]
    #     atoms = []
    #     for x_size in range(oriented_size[0]):
    #         for y_size in range(oriented_size[1]):
    #             for z_size in range(oriented_size[2]):
    #                 offset = np.array([
    #                     np.linalg.norm(oriented_vectors[0]) * x_size,
    #                     np.linalg.norm(oriented_vectors[1]) * y_size,
    #                     np.linalg.norm(oriented_vectors[2]) * z_size,
    #                 ])
    #                 for atom, relative_position in self._unit_cell:
    #                     position = np.zeros((3,))
    #                     for i in range(3):
    #                         position[i] = (np.linalg.norm(oriented_vectors[i]) * relative_position[i]) + offset[i]
    #                     if is_point_in_polyhedron(position, oriented_vectors):
    #                         atom.position = position
    #                         atoms.append(atom)
    #     # TODO: apply rotation

    def _build(self) -> None:
        lattice_vectors = np.sqrt(
            metric_tensor(self._a, self._b, self._c, self._alpha, self._beta,
                          self._gamma))
        if self._size is None:
            self._size = (0, 0, 0)
        oriented_vectors = np.zeros((3, 3))
        if self._orientation is None:
            oriented_vectors = lattice_vectors
        else:
            _, r = np.linalg.qr(self._orientation)
            for i in range(3):
                for j in range(3):
                    r[i][j] *= np.linalg.norm(lattice_vectors[i])
            oriented_vectors = r
        oriented_size = [0, 0, 0]
        for i in range(3):
            oriented_size[i] = np.linalg.norm(oriented_vectors[i]) / np.ceil(
                np.linalg.norm(lattice_vectors[i]))
            oriented_size[i] *= self._size[i]
        oriented_size = [int(x) for x in oriented_size]
        atoms = []
        for x_size in range(oriented_size[0]):
            for y_size in range(oriented_size[1]):
                for z_size in range(oriented_size[2]):
                    offset = np.array([
                        np.linalg.norm(lattice_vectors[0]) * x_size,
                        np.linalg.norm(lattice_vectors[1]) * y_size,
                        np.linalg.norm(lattice_vectors[2]) * z_size,
                    ])
                    for atom, relative_position in self._unit_cell:
                        position = np.zeros((3,))
                        for i in range(3):
                            position[i] = (np.linalg.norm(lattice_vectors[i]) *
                                           relative_position[i]) + offset[i]
                        atom.position = position
                        atoms.append(atom)
        for a in atoms:
            print(a.position)

        # TODO: rotate all atoms to given orientation
