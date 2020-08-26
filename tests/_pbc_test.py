import numpy as np

from atompack._pbc import pbc_nearest_neighbor


def test_pbc_nearest_neighbor_periodic():
    position = np.array([0.75, 0.75, 0.75])
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ])
    cell = np.identity(3)
    pbc = (True, True, True)

    (distance, index) = pbc_nearest_neighbor(position, positions, cell, pbc)

    assert distance == np.sqrt(3 * 0.25**2)
    assert index == 0


def test_pbc_nearest_neighbor_nonperiodic():
    position = np.array([0.75, 0.75, 0.75])
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ])
    cell = np.identity(3)
    pbc = (False, False, False)

    (distance, index) = pbc_nearest_neighbor(position, positions, cell, pbc)

    assert distance == np.sqrt(3 * 0.25**2)
    assert index == 1