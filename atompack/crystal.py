import copy
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from atompack._cell import cell_enforce
from atompack.atom import Atom
from atompack.elements import Element
from atompack.errors import CrystallographyError
from atompack.topology import Topology


def metric_tensor(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Returns the metric tensor for a set of lattice parameters."""
    return np.array([[a * a, a * b * np.cos(gamma), a * c * np.cos(beta)],
                     [a * b * np.cos(gamma), b * b, b * c * np.cos(alpha)],
                     [a * c * np.cos(beta), b * c * np.cos(alpha), c * c]])


def is_triclinic(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy triclinic constraints."
    if abs(a - b) < tol or abs(b - c) < tol or abs(a - c) < tol:
        return False
    if abs(alpha - beta) < tol or abs(beta - gamma) < tol or abs(alpha - gamma) < tol:
        return False
    return True


def is_monoclinic(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy monoclinic constraints."
    if abs(a - b) < tol or abs(b - c) < tol or abs(a - c) < tol:
        return False
    if abs(alpha - np.pi / 2) > tol or abs(gamma - np.pi / 2) > tol:
        return False
    if abs(beta - np.pi / 2) < tol:
        return False
    return True


def is_orthorhombic(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy orthorhombic constraints."
    if abs(a - b) < tol or abs(b - c) < tol or abs(a - c) < tol:
        return False
    if abs(alpha - np.pi / 2) > tol or abs(beta - np.pi / 2) > tol or abs(gamma - np.pi / 2) > tol:
        return False
    return True


def is_tetragonal(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy tetragonal constraints."
    if abs(a - b) > tol or abs(a - c) < tol or abs(b - c) < tol:
        return False
    if abs(alpha - np.pi / 2) > tol or abs(beta - np.pi / 2) > tol or abs(gamma - np.pi / 2) > tol:
        return False
    return True


def is_rhombohedral(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy rhombohedral constraints."
    if abs(a - b) > tol or abs(b - c) > tol or abs(a - c) > tol:
        return False
    if abs(alpha - beta) > tol or abs(beta - gamma) > tol or abs(alpha - gamma) > tol or abs(alpha - np.pi / 2) < tol:
        return False
    return True


def is_hexagonal(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy hexagonal constraints."
    if abs(a - b) > tol or abs(a - c) < tol or abs(b - c) < tol:
        return False
    if abs(alpha - beta) > tol or abs(alpha - np.pi / 2) > tol:
        return False
    if abs(gamma - 2 * np.pi / 3) > tol:
        return False
    return True


def is_cubic(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1.0e-6) -> bool:
    "Returns True if the lattice parameters satisfy cubic constraints."
    if abs(a - b) > tol or abs(b - c) > tol or abs(a - c) > tol:
        return False
    if abs(alpha - beta) > tol or abs(beta - gamma) > tol or abs(alpha - gamma) > tol or abs(alpha - np.pi / 2) > tol:
        return False
    return True


class UnitCell(Topology):
    """Representation of a paralellpiped tileable unit cell.
    
    Args:
        a: Length of the x lattice vector.
        b: Length of the y lattice vector.
        c: Length of the z lattice vector.
        alpha: Angle between the y and z directions in radians.
        beta: Angle between the x and z directions in radians.
        gamma: Angle between the x and y directions in radians.
        sites: List of fractional lattice sites.
        elements: Elemental data to associate with each site.
    """

    def __init__(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float, sites: np.ndarray,
                 elements: List[Optional[Element]]) -> None:
        self._a, self._b, self._c = a, b, c
        self._alpha, self._beta, self._gamma = alpha, beta, gamma
        self._sites = sites
        super().__init__()
        for site, element in zip(self.sites, elements):
            position = np.matmul(site, self.lattice_vectors)
            if element is None:
                atom = Atom(position)
            else:
                atom = Atom(position, **vars(element))
            self.insert(atom)

    @property
    def a(self) -> float:
        """Returns the length of the x lattice vector."""
        return self._a

    @property
    def b(self) -> float:
        """Returns the langth of the y lattice vector."""
        return self._b

    @property
    def c(self) -> float:
        """Returns the length of the z lattice vector."""
        return self._c

    @property
    def alpha(self) -> float:
        """Returns the angle between the y and z directions in radians."""
        return self._alpha

    @property
    def beta(self) -> float:
        """Returns the angle between the x and z directions in radians."""
        return self._beta

    @property
    def gamma(self) -> float:
        """Returns the angle between the x and y directions in radians."""
        return self._gamma

    @property
    def sites(self) -> np.ndarray:
        """Returns a copy of the fractional coordinates of each lattice site."""
        return copy.deepcopy(self._sites)

    @property
    def lattice_vectors(self) -> np.ndarray:
        """Returns the lattice vectors of the unit cell."""
        return np.sqrt(metric_tensor(self.a, self.b, self.c, self.alpha, self.beta, self.gamma))


class Triclinic(UnitCell):
    """Unit cell with triclinic constraints.
    
    \\[a \\ne b \\ne c\\]

    \\[\\alpha \\ne \\beta \\ne \\gamma\\]
    """

    def __init__(self,
                 a: float,
                 b: float,
                 c: float,
                 alpha: float,
                 beta: float,
                 gamma: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        if not is_triclinic(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a triclinic unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Monoclinic(UnitCell):
    """Unit cell with monoclinic constraints.
    
    \\[a \\ne b \\ne c\\]
        
    \\[\\alpha \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
        
    \\[\\beta \\ne \\frac{\\pi}{2}\\]
    """

    def __init__(self,
                 a: float,
                 b: float,
                 c: float,
                 beta: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        alpha = gamma = np.pi / 2
        if not is_monoclinic(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a monoclinic unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Orthorhombic(UnitCell):
    """Unit cell with orthorhombic constraints.
    
    \\[a \\ne b \\ne c\\]
        
    \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
    """

    def __init__(self,
                 a: float,
                 b: float,
                 c: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        alpha = beta = gamma = np.pi / 2
        if not is_orthorhombic(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for an orthorhombic unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Tetragonal(UnitCell):
    """Unit cell with tetragonal constraints.
    
    \\[a \\equiv b \\ne c\\]
        
    \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
    """

    def __init__(self,
                 a: float,
                 c: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        b = a
        alpha = beta = gamma = np.pi / 2
        if not is_tetragonal(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a Tetragonal unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Rutile(Tetragonal):
    """Rutile (TiO2) unit cell."""

    def __init__(self, a: float, c: float, ti: Optional[Element] = None, o: Optional[Element] = None) -> None:
        sites = np.array([
            [0.5, 0.5, 0.5],
            [0.0, 0.0, 0.0],
            [0.1954, 0.8046, 0.5],
            [0.8046, 0.1954, 0.50],
            [0.3046, 0.3046, 0.0],
            [0.6954, 0.6954, 0.0], 
        ])
        super(a, c, sites, [ti, ti, o, o, o, o])


class Rhombohedral(UnitCell):
    """Unit cell with rhombohedral constraints.
    
    \\[a \\equiv b \\equiv c\\]
        
    \\[\\alpha \\equiv \\beta \\equiv \\gamma \\ne \\frac{\\pi}{2}\\]
    """

    def __init__(self,
                 a: float,
                 alpha: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        b = c = a
        beta = gamma = alpha
        if not is_rhombohedral(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a rhombohedral unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Hexagonal(UnitCell):
    """Unit cell with hexagonal constraints.
    
    \\[a \\equiv b \\ne c\\]
        
    \\[\\alpha \\equiv \\beta \\equiv \\frac{\\pi}{2}\\]
        
    \\[\\gamma \\equiv \\frac{2\\pi}{3}\\]
    """

    def __init__(self,
                 a: float,
                 c: float,
                 sites: np.ndarray,
                 elements: List[Optional[Element]],
                 tol: float = 1.0e-6) -> None:
        b = a
        alpha = beta = np.pi / 2
        gamma = 2 * np.pi / 3
        if not is_hexagonal(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a hexagonal unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Hcp(Hexagonal):
    """Hexagonal close-packed unit cell."""

    def __init__(self, a, c, element: Optional[Element] = None) -> None:
        sites = np.ndarray([
            [2.0/3.0, 1.0/3.0, 0.75],
            [1.0/3.0, 2.0/3.0, 0.25],
        ])
        super().__init__(a, c, sites, [element, element])


class Wurtzite(Hexagonal):
    """Wurtzite (ZnS) unit cell."""

    def __init__(self, a, c, zn: Optional[Element] = None, s: Optional[Element] = None) -> None:
        sites = np.ndarray([
            [2.0/3.0, 1.0/3.0, 0.0],
            [1.0/3.0, 2.0/3.0, 0.5],
            [2.0/3.0, 1.0/3.0, 0.6259],
            [1.0/3.0, 2.0/3.0, 0.1259],
        ])
        super().__init__(a, c, sites, [zn, zn, s, s])

class Cubic(UnitCell):
    """Unit cell with cubic constraints.
    
    \\[a \\equiv b \\equiv c\\]
        
    \\[\\alpha \\equiv \\beta \\equiv \\gamma \\equiv \\frac{\\pi}{2}\\]
    """

    def __init__(self, a: float, sites: np.ndarray, elements: List[Optional[Element]], tol: float = 1.0e-6) -> None:
        b = c = a
        alpha = beta = gamma = np.pi / 2
        if not is_cubic(a, b, c, alpha, beta, gamma, tol):
            raise CrystallographyError("Invalid lattice parameters for a cubic unit cell.")
        super().__init__(a, b, c, alpha, beta, gamma, sites, elements)


class Sc(Cubic):
    """Simple cubic unit cell."""

    def __init__(self, a: float, element: Optional[Element] = None) -> None:
        sites = np.zeros(3)
        super().__init__(a, sites, [element])


class Bcc(Cubic):
    """Body-centered cubic unit cell."""

    def __init__(self, a: float, element: Optional[Element] = None) -> None:
        sites = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ])
        super().__init__(a, sites, [element, element])


class CsCl(Cubic):
    """Cesium Chloride unit cell."""

    def __init__(self, a: float, cs: Optional[Element] = None, cl: Optional[Element] = None):
        sites = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ])
        super().__init__(a, sites, [cs, cl])


class Fcc(Cubic):
    """Face-centered cubic unit cell."""

    def __init__(self, a: float, element: Optional[Element] = None) -> None:
        sites = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ])
        super().__init__(a, sites, [element, element, element, element])


class NaCl(Cubic):
    """Rock salt unit cell."""

    def __init__(self, a: float, na: Optional[Element] = None, cl: Optional[Element] = None) -> None:
        sites = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 0.5, 0.0],
        ])
        super().__init__(a, sites, [na, na, na, na, cl, cl, cl, cl])


class Diamond(Cubic):
    """Diamond unit cell."""

    def __init__(self, a: float, element: Optional[Element] = None) -> None:
        sites = np.array([
            [0.25, 0.75, 0.25],
            [0.0, 0.0, 0.5],
            [0.25, 0.25, 0.75],
            [0.0, 0.5, 0.0],
            [0.75, 0.75, 0.75],
            [0.5, 0.0, 0.0],
            [0.75, 0.25, 0.25],
            [0.5, 0.5, 0.5],
        ])
        super().__init__(a, sites, [element]*len(sites))


class ZincBlend(Cubic):
    """Zinc blend (ZnS) unit cell."""

    def __init__(self, a, zn: Optional[Element] = None, s: Optional[Element] = None) -> None:
        sites = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.25, 0.25, 0.75],
            [0.25, 0.75, 0.25],
            [0.75, 0.25, 0.25],
            [0.75, 0.75, 0.75],
        ])
        super().__init__(a, sites, [zn, zn, zn, zn, s, s, s, s])


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
                        cell_enforce(oriented_lattice_vectors, position, self._tol)

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
