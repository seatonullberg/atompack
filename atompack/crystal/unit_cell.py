import copy
from typing import List, Optional

import numpy as np

from atompack.atom import Atom
from atompack.crystal.util import (is_cubic, is_hexagonal, is_monoclinic, is_orthorhombic, is_rhombohedral,
                                   is_tetragonal, is_triclinic, metric_tensor)
from atompack.elements import Element
from atompack.errors import CrystallographyError
from atompack.topology import Topology


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
                atom = Atom(position, **element.as_dict())
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
        super().__init__(a, c, sites, [ti, ti, o, o, o, o])


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
        sites = np.array([
            [2.0 / 3.0, 1.0 / 3.0, 0.75],
            [1.0 / 3.0, 2.0 / 3.0, 0.25],
        ])
        super().__init__(a, c, sites, [element, element])


class Wurtzite(Hexagonal):
    """Wurtzite (ZnS) unit cell."""

    def __init__(self, a, c, zn: Optional[Element] = None, s: Optional[Element] = None) -> None:
        sites = np.array([
            [2.0 / 3.0, 1.0 / 3.0, 0.0],
            [1.0 / 3.0, 2.0 / 3.0, 0.5],
            [2.0 / 3.0, 1.0 / 3.0, 0.6259],
            [1.0 / 3.0, 2.0 / 3.0, 0.1259],
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
        super().__init__(a, sites, [element] * len(sites))


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
