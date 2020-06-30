import numpy as np
import pytest

from atompack.errors import (PositionOccupiedError, PositionOutsideError, PositionUnoccupiedError)
from atompack.structure import Structure


def test_structure_iter():
    atoms = [{"symbol": "Fe"}, {"symbol": "Fe"}]
    structure = Structure(atoms)
    assert len(structure) == 2
    for atom in structure:
        assert atom["symbol"] == "Fe"


def test_structure_len():
    structure = Structure()
    assert len(structure) == 0
    structure.atoms.append({"symbol": "Fe"})
    assert len(structure) == 1


def test_structure_insert_occupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    with pytest.raises(PositionOccupiedError):
        structure.insert({"position": np.array([0.0, 0.0, 0.0])})


def test_structure_insert_unoccupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    structure.insert({"position": np.array([0.5, 0.5, 0.5])})
    assert len(structure) == 2


def test_structure_insert_outside():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    with pytest.raises(PositionOutsideError):
        structure.insert({"position": np.array([1.5, 1.5, 1.5])})


def test_structure_remove_occupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0]), "symbol": "Fe"}]
    structure = Structure(atoms)
    atom = structure.remove(np.array([0.0, 0.0, 0.0]))
    assert atom["symbol"] == "Fe"
    assert len(structure) == 0


def test_structure_remove_unoccupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    with pytest.raises(PositionUnoccupiedError):
        _ = structure.remove(np.array([0.5, 0.5, 0.5]))


def test_structure_select_occupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    index = structure.select(np.array([0.0, 0.0, 0.0]))
    assert index == 0


def test_structure_select_unoccupied():
    atoms = [{"position": np.array([0.0, 0.0, 0.0])}]
    structure = Structure(atoms)
    with pytest.raises(PositionUnoccupiedError):
        _ = structure.select(np.array([0.5, 0.5, 0.5]))
