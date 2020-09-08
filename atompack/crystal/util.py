import numpy as np


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


def enforce_bounds(bounds: np.ndarray, position: np.ndarray, tol: float) -> None:
    """Mutates `position` such that it lies within the parallelpiped defined by `bounds`.

    Args:
        bounds: The magnitude of each vector of the parallelpiped boundary cell.
        position: The position to translate back into the cell.
        tol: The tolerance for float comparison.
    """
    for i in range(3):
        bound = bounds[i]
        tmpval = position[i]
        if tmpval > bound - tol:
            while tmpval > bound - tol:
                tmpval -= bound
        if tmpval < -tol:
            while tmpval < -tol:
                tmpval += bound
        position[i] = tmpval
