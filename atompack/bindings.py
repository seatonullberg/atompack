from ctypes import POINTER, byref, c_double, c_int, c_size_t, cdll
from typing import Tuple

import numpy as np


def load_libatompack():
    # TODO: get names for other platforms
    name = "libatompack.so"
    lib = cdll.LoadLibrary(name)

    # cell_contains
    lib.cell_contains.restype = c_int
    lib.cell_contains.argtypes = [
        POINTER(POINTER(c_double)),
        POINTER(c_double),
        c_double,
    ]

    # metric_tensor
    lib.metric_tensor.restype = None
    lib.metric_tensor.argtypes = [
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        c_double,
        POINTER(POINTER(c_double)),
    ]

    # nearest_neighbor
    lib.nearest_neighbor.restype = c_size_t
    lib.nearest_neighbor.argtypes = [
        POINTER(c_double),
        POINTER(POINTER(c_double)), c_size_t,
        POINTER(POINTER(c_double)),
        POINTER(c_int), c_double,
        POINTER(c_double)
    ]
    return lib


def cell_contains(lib, cell: np.ndarray, position: np.ndarray, tolerance: float) -> bool:
    # process cell
    c_cell = cell.astype(np.float64)
    c_cell = c_cell.ctypes.data_as(POINTER(POINTER(c_double)))
    # process position
    c_position = position.astype(np.float64)
    c_position = c_position.ctypes.data_as(POINTER(c_double))
    # process tolerance
    c_tolerance = c_double(tolerance)
    res = lib.cell_contains(c_cell, c_position, c_tolerance)
    return bool(res)


def metric_tensor(lib, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    c_a = c_double(a)
    c_b = c_double(b)
    c_c = c_double(c)
    c_alpha = c_double(alpha)
    c_beta = c_double(beta)
    c_gamma = c_double(gamma)
    res = np.zeros((3, 3), dtype=np.float64)
    c_res = res.ctypes.data_as(POINTER(POINTER(c_double)))
    lib.metric_tensor(c_a, c_b, c_c, c_alpha, c_beta, c_gamma, c_res)
    return res


def nearest_neighbor(lib, position: np.ndarray, positions: np.ndarray, cell: np.ndarray, periodicity: np.ndarray,
                     tolerance: float) -> Tuple[int, float]:
    # process position
    c_position = position.astype(np.float64)
    c_position = c_position.ctypes.data_as(POINTER(c_double))
    # process positions
    c_positions = positions.astype(np.float64)
    c_positions = c_positions.ctypes.data_as(POINTER(POINTER(c_double)))
    # process length
    c_length = c_size_t(positions.shape[0])
    # process cell
    c_cell = cell.astype(np.float64)
    c_cell = c_cell.ctypes.data_as(POINTER(POINTER(c_double)))
    # process periodicity
    c_periodicity = periodicity.astype(np.int32)
    c_periodicity = c_periodicity.ctypes.data_as(POINTER(c_int))
    # process tolerance
    c_tolerance = c_double(tolerance)
    # process out
    out = c_double(0.0)
    c_out = byref(out)
    res = lib.nearest_neighbor(c_position, c_positions, c_length, c_cell, c_periodicity, c_tolerance, c_out)
    return int(res), float(out.value)
