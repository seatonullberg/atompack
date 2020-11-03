import numpy as np
import pytest

from atompack.atom import Atom
from atompack.bond import Bond
from atompack.topology import Topology

N_ATOMS = 10
N_BONDS = 5

#######################
#    Test Fixtures    #
#######################


@pytest.fixture
def topology():
    """Returns a pre-populated Topology of `N_ATOMS` and `N_BONDS`.
    The first atom is bonded to the next `N_BONDS` atoms.
    The remaining atoms are not bonded.
    """
    atoms = [Atom("TEST", np.zeros(3)) for _ in range(N_ATOMS)]
    bonds = [Bond((0, i + 1)) for i in range(N_BONDS)]
    res = Topology()
    res.insert_atoms(*atoms)
    # TODO: fix this in later release of retworkx
    for bond in bonds:
        res.insert_bond(bond)
    return res


########################
#    Topology Tests    #
########################


# TODO: is this test really necessary / helpful ?
def test_topology_insert_atoms(topology):
    # insert multiple atoms simultaneously
    indices = topology.insert_atoms(Atom("X", np.zeros(3)), Atom("Y", np.zeros(3)))
    assert len(topology.atoms) == N_ATOMS + len(indices)
    # insert just one atom
    indices += topology.insert_atoms(Atom("Z", np.zeros(3)))
    assert len(topology.atoms) == N_ATOMS + len(indices)


def test_topology_remove_atoms_invalid_index(topology):
    # try removing an invalid atom
    with pytest.raises(IndexError):
        _ = topology.remove_atoms(N_ATOMS + 1)


def test_topology_remove_bonded_atoms(topology):
    # remove the centrally bonded atom
    topology.remove_atoms(0)
    assert len(topology.atoms) == N_ATOMS - 1
    assert len(topology.bonds) == 0


def test_topology_remove_nonbonded_atoms(topology):
    # remove a nonbounded atom
    topology.remove_atoms(N_ATOMS - 1)
    assert len(topology.atoms) == N_ATOMS - 1
    assert len(topology.bonds) == N_BONDS


def test_topology_select_atoms_invalid_index(topology):
    # try selecting an invalid atom
    with pytest.raises(IndexError):
        _ = topology.select_atoms(N_ATOMS + 1)


def test_topology_select_atoms_mutability(topology):
    test_value = "TEST"
    # select all atoms
    atoms = topology.select_atoms(*range(N_ATOMS))
    # mutate each atom
    for atom in atoms:
        atom["test_value"] = test_value
    # select the atoms again
    atoms = topology.select_atoms(*range(N_ATOMS))
    # check each atom
    for atom in atoms:
        assert atom["test_value"] == test_value


def test_topology_to_from_json(topology):
    json_data = topology.to_json()
    new_topology = Topology.from_json(json_data)
    for i in range(N_ATOMS):
        assert new_topology.atoms[i].specie == topology.atoms[i].specie
        assert np.allclose(new_topology.atoms[i].position, topology.atoms[i].position)
    for i in range(N_BONDS):
        assert new_topology.bonds[i].indices == topology.bonds[i].indices


# TODO: tests for bond operations will be added after the retworkx update
