import numpy as np


def is_point_in_cell(cell: np.ndarray, point: np.ndarray, tol: float = 1E-6) -> bool:
    """Returns True if `point` is within the parallelepiped cell 
    defined by the 3x3 matrix `cell`."""
    bounds = np.linalg.norm(cell, axis=0)
    for i in range(3):
        if point[i] > (bounds[i] + tol):
            return False
        if point[i] < -tol:
            return False
    return True


def wrap_point_into_cell(cell: np.ndarray, point: np.ndarray, tol: float = 1E-6) -> np.ndarray:
    """Mutates `point` to reside within the parallelepiped cell defined by 
    the 3x3 matrix `cell`.
    Returns the mutated `point`.
    """
    bounds = np.linalg.norm(cell, axis=0)
    for i in range(3):
        bound = bounds[i]
        tmpval = point[i]
        if tmpval > (bound + tol):
            tmpval -= bound * (tmpval // bound)
        if tmpval < -tol:
            tmpval += bound * (1 + (-tmpval // bound))
        point[i] = tmpval
    return point
