# import numpy as np
# from numpy.testing import assert_array_almost_equal

# from atompack.atom import Atom
# from atompack.crystal import Crystal, UnitCell, metric_tensor

# def get_cubic_unit_cell():
#     a, b, c = 2.85, 2.85, 2.85
#     alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
#     sites = np.array([
#         [0.0, 0.0, 0.0],
#         [0.5, 0.5, 0.5],
#     ])
#     return UnitCell(a, b, c, alpha, beta, gamma, sites, [None, None])

# def test_crystal_metric_tensor():
#     a, b, c = 2.0, 2.0, 2.0
#     alpha, beta, gamma = np.pi / 2, np.pi / 2, np.pi / 2
#     res = metric_tensor(a, b, c, alpha, beta, gamma)
#     target = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
#     assert_array_almost_equal(res, target)

# def test_crystal_cubic_100_orientation():
#     unit_cell = get_cubic_unit_cell()
#     crystal = Crystal(unit_cell)
#     res_vectors = crystal.lattice_vectors
#     target_vectors = np.array([[2.85, 0.0, 0.0], [0.0, 2.85, 0.0], [0.0, 0.0, 2.85]])
#     assert_array_almost_equal(res_vectors, target_vectors)
#     res_positions = np.array([atom.position for atom in crystal.atoms])
#     target_positions = np.array([[0.0, 0.0, 0.0], [1.425, 1.425, 1.425]])
#     assert_array_almost_equal(res_positions, target_positions)

# def test_crystal_cubic_110_orientation():
#     unit_cell = get_cubic_unit_cell()
#     orientation = np.array([[-1, 1, 0], [0, 0, 1], [1, 1, 0]])
#     crystal = Crystal(unit_cell, orientation=orientation)
#     res_vectors = crystal.lattice_vectors
#     target_vectors = np.array([[np.sqrt(2) * 2.85, 0, 0], [0, 2.85, 0], [0, 0, np.sqrt(2) * 2.85]])
#     assert_array_almost_equal(res_vectors, target_vectors)
#     res_positions = np.array([atom.position for atom in crystal.atoms])
#     target_positions = np.array([
#         [0.0, 0.0, 0.0],
#         [0.0, 1.425, 2.015254],
#         [2.015254, 0.0, 2.015254],
#         [2.015254, 1.425, 0.0],
#     ])
#     assert_array_almost_equal(res_positions, target_positions)

# def test_crystal_cubic_111_orientation():
#     unit_cell = get_cubic_unit_cell()
#     orientation = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]])
#     crystal = Crystal(unit_cell, orientation=orientation)
#     res_vectors = crystal.lattice_vectors
#     target_vectors = np.array([[np.sqrt(2) * 2.85, 0, 0], [0, np.sqrt(6) * 2.85, 0], [0, 0, np.sqrt(3) * 2.85]])
#     assert_array_almost_equal(res_vectors, target_vectors)
#     res_positions = np.array([atom.position for atom in crystal.atoms])
#     target_positions = np.array([
#         [0.00000000, 0.00000000, 0.00000000],
#         [0.00000000, 0.00000000, 2.46817240],
#         [0.00000000, 4.65403051, 1.64544827],
#         [0.00000000, 4.65403051, 4.11362067],
#         [2.01525433, 1.16350763, 1.64544827],
#         [2.01525433, 1.16350763, 4.11362067],
#         [2.01525433, 5.81753814, 3.29089653],
#         [2.01525433, 5.81753814, 0.82272413],
#         [0.00000000, 2.32701526, 3.29089653],
#         [0.00000000, 2.32701526, 0.82272413],
#         [2.01525433, 3.49052288, 0.00000000],
#         [2.01525433, 3.49052288, 2.46817240],
#     ])
#     assert_array_almost_equal(res_positions, target_positions)

# def test_crystal_cubic_2x2x2_super_cell():
#     unit_cell = get_cubic_unit_cell()
#     scale = (2, 2, 2)
#     crystal = Crystal(unit_cell, scale=scale)
#     res_vectors = crystal.lattice_vectors
#     target_vectors = np.array([[5.7, 0, 0], [0, 5.7, 0], [0, 0, 5.7]])
#     assert_array_almost_equal(res_vectors, target_vectors)
#     res_positions = np.array([atom.position for atom in crystal.atoms])
#     target_positions = np.array([
#         [0.00000000, 0.00000000, 0.00000000],
#         [1.42500000, 1.42500000, 1.42500000],
#         [0.00000000, 0.00000000, 2.85000000],
#         [1.42500000, 1.42500000, 4.27500000],
#         [0.00000000, 2.85000000, 0.00000000],
#         [1.42500000, 4.27500000, 1.42500000],
#         [0.00000000, 2.85000000, 2.85000000],
#         [1.42500000, 4.27500000, 4.27500000],
#         [2.85000000, 0.00000000, 0.00000000],
#         [4.27500000, 1.42500000, 1.42500000],
#         [2.85000000, 0.00000000, 2.85000000],
#         [4.27500000, 1.42500000, 4.27500000],
#         [2.85000000, 2.85000000, 0.00000000],
#         [4.27500000, 4.27500000, 1.42500000],
#         [2.85000000, 2.85000000, 2.85000000],
#         [4.27500000, 4.27500000, 4.27500000],
#     ])
#     assert_array_almost_equal(res_positions, target_positions)
