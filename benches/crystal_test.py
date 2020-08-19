import numpy as np

from atompack.atom import Atom
from atompack.crystal import Crystal


def bench_crystal(lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, duplicates, orientation, pbc, tolerance):
    return Crystal(lattice_atoms,
                   lattice_sites,
                   a,
                   b,
                   c,
                   alpha,
                   beta,
                   gamma,
                   duplicates=duplicates,
                   orientation=orientation,
                   pbc=pbc,
                   tolerance=tolerance)


def test_bench_crystal_Fe_BCC_100_1x1x1(benchmark):
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    a, b, c, = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    duplicates = (1, 1, 1)
    orientation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    pbc = (False, False, False)
    tolerance = 1e-6
    args = (lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, duplicates, orientation, pbc, tolerance)
    res = benchmark.pedantic(bench_crystal, args=args, iterations=10, rounds=25)
    assert len(res) == 2


def test_bench_crystal_Fe_BCC_110_1x1x1(benchmark):
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    a, b, c, = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    duplicates = (1, 1, 1)
    orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
    pbc = (False, False, False)
    tolerance = 1e-6
    args = (lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, duplicates, orientation, pbc, tolerance)
    res = benchmark.pedantic(bench_crystal, args=args, iterations=10, rounds=25)
    assert len(res) == 4


def test_bench_crystal_Fe_BCC_111_1x1x1(benchmark):
    lattice_atoms = [Atom(symbol="Fe"), Atom(symbol="Fe")]
    lattice_sites = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    a, b, c, = 2.85, 2.85, 2.85
    alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
    duplicates = (1, 1, 1)
    orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
    pbc = (False, False, False)
    tolerance = 1e-6
    args = (lattice_atoms, lattice_sites, a, b, c, alpha, beta, gamma, duplicates, orientation, pbc, tolerance)
    res = benchmark.pedantic(bench_crystal, args=args, iterations=10, rounds=25)
    assert len(res) == 12
