import numpy as np

DEG90 = np.pi / 2  # 90 degrees in radians
DEG120 = 2 * np.pi / 3  # 120 degrees in radians


class LatticeParameters(object):
    """Edge lengths and angles which define a lattice.
    
    Args:
        a: Length of the x lattice vector.
        b: Length of the y lattice vector.
        c: Length of the z lattice vector.
        alpha: Angle between the y and z directions (radians).
        beta: Angle between the x and z directions (radians).
        gamma: Angle between the x and y directions (radians).
    """

    def __init__(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    @classmethod
    def triclinic(cls, a, b, c, alpha, beta, gamma) -> 'LatticeParameters':
        return cls(a, b, c, alpha, beta, gamma)

    @classmethod
    def monoclinic(cls, a, b, c, beta) -> 'LatticeParameters':
        return cls(a, b, c, DEG90, beta, DEG90)

    @classmethod
    def orthorhombic(cls, a, b, c) -> 'LatticeParameters':
        return cls(a, b, c, DEG90, DEG90, DEG90)

    @classmethod
    def tetragonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG90)

    @classmethod
    def trigonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG120)

    @classmethod
    def rhombohedral(cls, a, alpha) -> 'LatticeParameters':
        return cls(a, a, a, alpha, alpha, alpha)

    @classmethod
    def hexagonal(cls, a, c) -> 'LatticeParameters':
        return cls(a, a, c, DEG90, DEG90, DEG120)

    @classmethod
    def cubic(cls, a) -> 'LatticeParameters':
        return cls(a, a, a, DEG90, DEG90, DEG90)

    @property
    def metric_tensor(self) -> np.ndarray:
        """Returns the metric tensor defined by the lattice parameters."""
        return np.array([[self.a * self.a, self.a * self.b * np.cos(self.gamma), self.a * self.c * np.cos(self.beta)],
                         [self.a * self.b * np.cos(self.gamma), self.b * self.b, self.b * self.c * np.cos(self.alpha)],
                         [self.a * self.c * np.cos(self.beta), self.b * self.c * np.cos(self.alpha), self.c * self.c]])
