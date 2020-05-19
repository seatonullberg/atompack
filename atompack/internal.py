import numpy as np


def search_for_atom(atoms, position, tolerance):
    index = None
    for i, atom in enumerate(atoms):
        if np.linalg.norm(atom.position - position) < tolerance:
            index = i
            break
    return index
