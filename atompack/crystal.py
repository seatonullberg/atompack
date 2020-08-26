import copy
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack._cell import cell_enforce
from atompack.atom import Atom
from atompack.topology import Topology


def metric_tensor(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Returns the metric tensor given a set of lattice parameters."""
    return np.array([[a * a, a * b * np.cos(gamma), a * c * np.cos(beta)],
                     [a * b * np.cos(gamma), b * b, b * c * np.cos(alpha)],
                     [a * c * np.cos(beta), b * c * np.cos(alpha), c * c]])


class UnitCell(Topology):
    """Representation of a paralellpiped tileable unit cell.
    
    Args:
        atoms: The atoms in the cell.
        a: Length of the x lattice vector.
        b: Length of the y lattice vector.
        c: Length of the z lattice vector.
        alpha: Angle between y and z directions in radians.
        beta: Angle between x and z directions in radians.
        gamma: Angle between x and y directions in radians.
        fractional: (default True) Flag indicating the `position` attribute of incoming atoms 
            is in fractional coordinates. Fractional coordinates will be overwritten with 
            cartesian coordinates during initialization.
    """

    def __init__(self,
                 atoms: List[Atom],
                 a: float,
                 b: float,
                 c: float,
                 alpha: float,
                 beta: float,
                 gamma: float,
                 fractional: bool = True) -> None:
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self._lattice_vectors = np.sqrt(metric_tensor(self.a, self.b, self.c, self.alpha, self.beta, self.gamma))
        super().__init__()
        if fractional:
            for atom in atoms:
                atom.position = np.matmul(atom.position, self._lattice_vectors)
                self.insert(atom)
        else:
            for atom in atoms:
                self.insert(atom)

    @property
    def lattice_vectors(self):
        """Returns the lattice vectors of the unit cell."""
        return self._lattice_vectors

    @classmethod
    def triclinic(
        cls,
        atoms: List[Atom],
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
    ) -> 'UnitCell':
        """Initializes a unit cell with triclinic constraints.
        
        \\[a \\ne b \\ne c\\]

        \\[\\alpha \\ne \\beta \\ne \\gamma\\]
        """
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def monoclinic(
        cls,
        atoms: List[Atom],
        a: float,
        b: float,
        c: float,
        beta: float,
    ) -> 'UnitCell':
        """Initializes a crystal with monoclinic constraints.
        
        \\[a \\ne b \\ne c\\]
        
        \\[\\alpha \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
        
        \\[\\beta \\ne \\frac{\\pi}{2}\\]
        """
        alpha = gamma = np.pi / 2
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def orthorhombic(
        cls,
        atoms: List[Atom],
        a: float,
        b: float,
        c: float,
    ) -> 'UnitCell':
        """Initializes a crystal with orthorhombic constraints.
        
        \\[a \\ne b \\ne c\\]
        
        \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
        """
        alpha = beta = gamma = np.pi / 2
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def tetragonal(
        cls,
        atoms: List[Atom],
        a: float,
        c: float,
    ) -> 'UnitCell':
        """Initializes a crystal with tetragonal constraints.
        
        \\[a \\equiv b \\ne c\\]
        
        \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
        """
        b = a
        alpha = beta = gamma = np.pi / 2
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def rhombohedral(
        cls,
        atoms: List[Atom],
        a: float,
        alpha: float,
    ) -> 'UnitCell':
        """Initializes a crystal with rhombohedral constraints.
        
        \\[a \\equiv b \\equiv c\\]
        
        \\[\\alpha \\equiv \\beta \\equiv \\gamma \\ne \\frac{\\pi}{2}\\]
        """
        b = c = a
        beta = gamma = alpha
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def hexagonal(
        cls,
        atoms: List[Atom],
        a: float,
        c: float,
    ) -> 'UnitCell':
        """Initializes a crystal with hexagonal constraints.
        
        \\[a \\equiv b \\ne c\\]
        
        \\[\\alpha \\equiv \\beta \\equiv \\frac{\\pi}{2}\\]
        
        \\[\\gamma \\equiv \\frac{2\\pi}{3}\\]
        """
        b = a
        alpha = beta = np.pi / 2
        gamma = 2 * np.pi / 3
        return cls(atoms, a, b, c, alpha, beta, gamma)

    @classmethod
    def cubic(
        cls,
        atoms: List[Atom],
        a: float,
    ) -> 'UnitCell':
        """Initializes a crystal with cubic constraints.
        
        \\[a \\equiv b \\equiv c\\]
        
        \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
        """
        b = c = a
        alpha = beta = gamma = np.pi / 2
        return cls(atoms, a, b, c, alpha, beta, gamma)


class Crystal(Topology):
    """Representation of a generic crystal."""

    def __init__(self,
                 unit_cell: UnitCell,
                 scale: Optional[Tuple[int, int, int]] = None,
                 orientation: Optional[np.ndarray] = None,
                 rotation: Optional[np.ndarray] = None,
                 tolerance: float = 1.0e-6) -> None:
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
        self._tolerance = tolerance
        atoms, self._lattice_vectors = self._build()
        super().__init__()
        for atom in atoms:
            self.insert(atom)

    @property
    def lattice_vectors(self) -> np.ndarray:
        """Returns the lattice vectors of the crystal."""
        return self._lattice_vectors

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
                        cell_enforce(oriented_lattice_vectors, position, self._tolerance)

                        # accept the atom if the position is not yet occupied
                        positions = np.array([atom.position for atom in atoms])
                        is_occupied = False
                        for _position in positions:
                            if np.linalg.norm(position - _position) < self._tolerance:
                                is_occupied = True
                                break
                        if not is_occupied:
                            atom = copy.deepcopy(atom)
                            atom.position = position
                            atoms.append(atom)

        # TODO: apply a rotation matrix

        # tile the crystal in all directions
        current_atoms = copy.deepcopy(atoms)
        for xsize in range(self._scale[0]):
            for ysize in range(self._scale[1]):
                for zsize in range(self._scale[2]):
                    offset = np.matmul(np.array([xsize, ysize, zsize]), oriented_lattice_vectors)
                    if np.linalg.norm(offset) < self._tolerance:
                        continue
                    for atom in current_atoms:
                        atom = copy.deepcopy(atom)
                        atom.position += offset
                        atoms.append(atom)

        # multiply the oriented lattice vectors by the scale
        scaled_oriented_lattice_vectors = oriented_lattice_vectors * np.array(self._scale)

        return atoms, scaled_oriented_lattice_vectors
