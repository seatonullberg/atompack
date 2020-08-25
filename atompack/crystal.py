import copy
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack.atom import Atom
from atompack.topology import Topology


def metric_tensor(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Returns the metric tensor given a set of lattice parameters."""
    return np.array([[a * a, a * b * np.cos(gamma), a * c * np.cos(beta)],
                     [a * b * np.cos(gamma), b * b, b * c * np.cos(alpha)],
                     [a * c * np.cos(beta), b * c * np.cos(alpha), c * c]])


class UnitCell(Topology):
    """Representation of a paralellpiped tileable unit cell."""

    def __init__(self, atoms: List[Atom], a: float, b: float, c: float, alpha: float, beta: float,
                 gamma: float) -> None:
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self._vectors = np.sqrt(metric_tensor(self.a, self.b, self.c, self.alpha, self.beta, self.gamma))
        super().__init__()
        for atom in atoms:
            self.insert(atom)

    @property
    def vectors(self):
        return self._vectors

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

        Example:
            >>> from atompack import Atom, Crystal
            >>> import numpy as np
            >>> 
            >>> # Ag2F5: https://materialsproject.org/materials/mp-542298/#
            >>> # angstroms
            >>> a, b, c = 5.197, 7.589, 11.446
            >>> # radians
            >>> alpha, beta, gamma = 1.536, 1.557, 1.271
            >>> # fractional coordinates
            >>> atoms = [
            ...     Atom(np.array([0.0428, 0.2779, 0.7534]), symbol="Ag"),
            ...     Atom(np.array([0.046, 0.7848, 0.7656]),  symbol="Ag"),
            ...     Atom(np.array([0.5, 0.5, 0.0]),          symbol="Ag"),
            ...     Atom(np.array([0.5, 0.0, 0.5]),          symbol="Ag"),
            ...     Atom(np.array([0.5, 0.0, 0.0]),          symbol="Ag"),
            ...     Atom(np.array([0.5, 0.5, 0.5]),          symbol="Ag"),
            ...     Atom(np.array([0.954, 0.2152, 0.2344]),  symbol="Ag"),
            ...     Atom(np.array([0.9572, 0.7221, 0.2466]), symbol="Ag"),
            ...     Atom(np.array([0.0487, 0.0776, 0.3825]), symbol="F"),
            ...     Atom(np.array([0.056, 0.9765, 0.1535]),  symbol="F"),
            ...     Atom(np.array([0.1279, 0.6441, 0.9187]), symbol="F"),
            ...     Atom(np.array([0.1539, 0.5451, 0.6861]), symbol="F"),
            ...     Atom(np.array([0.2065, 0.6609, 0.399]),  symbol="F"),
            ...     Atom(np.array([0.2783, 0.5616, 0.1509]), symbol="F"),
            ...     Atom(np.array([0.3273, 0.2291, 0.9045]), symbol="F"),
            ...     Atom(np.array([0.3611, 0.122, 0.6541]),  symbol="F"),
            ...     Atom(np.array([0.4498, 0.2667, 0.4277]), symbol="F"),
            ...     Atom(np.array([0.4799, 0.1616, 0.1324]), symbol="F"),
            ...     Atom(np.array([0.5201, 0.8384, 0.8676]), symbol="F"),
            ...     Atom(np.array([0.5502, 0.7333, 0.5723]), symbol="F"),
            ...     Atom(np.array([0.6389, 0.878, 0.3459]),  symbol="F"),
            ...     Atom(np.array([0.6727, 0.7709, 0.0955]), symbol="F"),
            ...     Atom(np.array([0.7217, 0.4384, 0.8491]), symbol="F"),
            ...     Atom(np.array([0.7935, 0.3391, 0.601]),  symbol="F"),
            ...     Atom(np.array([0.8461, 0.4549, 0.3139]), symbol="F"),
            ...     Atom(np.array([0.8721, 0.3559, 0.0813]), symbol="F"),
            ...     Atom(np.array([0.944, 0.0235, 0.8465]),  symbol="F"),
            ...     Atom(np.array([0.9513, 0.9224, 0.6175]), symbol="F"),
            ... ]
            >>> unit_cell = UnitCell.triclinic(atoms, a, b, c, alpha, beta, gamma)
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
        atoms, self._vectors = self._build()
        super().__init__()
        for atom in atoms:
            self.insert(atom)

    @property
    def vectors(self) -> np.ndarray:
        """Returns the lattice vectors of the crystal."""
        return self._vectors

    def _build(self) -> Tuple[List[Atom], np.ndarray]:
        # transforms are applied in the following order:
        # - orientation
        # - rotation
        # - scale

        # calculate the magnitude of the lattice vectors
        lattice_vector_mags = np.linalg.norm(self._unit_cell.vectors, axis=0)

        # calculate the unit vector of each lattice vector
        lattice_unit_vectors = self._unit_cell.vectors / lattice_vector_mags

        # calculate the rotation matrix between the unoriented and oriented lattice vectors
        rotation = Rotation.align_vectors(lattice_unit_vectors, self._orientation)[0]

        # align the lattice vectors with the orientation
        oriented_lattice_vectors = np.matmul(self._orientation, self._unit_cell.vectors)

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
                    offset = np.matmul(np.array([xsize, ysize, zsize]), self._unit_cell.vectors)
                    for atom in self._unit_cell.atoms:

                        # calculate the cartesian position
                        position = np.matmul(atom.position, self._unit_cell.vectors) + offset
                        position = rotation.apply(position)

                        # TODO: transform the position into the lattice
                        for i, mag in enumerate(oriented_lattice_vector_mags):
                            if position[i] >= mag - self._tolerance:
                                position[i] -= mag
                            if position[i] <= -self._tolerance:
                                position[i] += mag

                        # TODO: accept the atom if the position is not yet occupied
                        positions = np.array([atom.position for atom in atoms])
                        is_occupied = False
                        for _position in positions:
                            if np.linalg.norm(position - _position) < self._tolerance:
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
                    if np.linalg.norm(offset) < self._tolerance:
                        continue
                    for atom in current_atoms:
                        atom = copy.deepcopy(atom)
                        atom._position += offset
                        atoms.append(atom)

        # multiply the oriented lattice vectors by the scale
        scaled_oriented_lattice_vectors = oriented_lattice_vectors * np.array(self._scale)

        return atoms, scaled_oriented_lattice_vectors
