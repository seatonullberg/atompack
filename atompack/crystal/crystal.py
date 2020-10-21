import copy
import json
from typing import Optional, Tuple

import numpy as np

from atompack.atom import Atom
from atompack.bond import Bond
from atompack.crystal.components import (Basis, LatticeParameters, LatticeVectors)
from atompack.crystal.spatial import MillerIndex, Orientation, Plane
from atompack.symmetry import Spacegroup
from atompack.topology import Topology


class UnitCell(Topology):
    """Minimal representation of a crystalline structure.
    
    Args:
        basis: Atomic basis.
        lattice_parameters: Lattice parameters object.
        spacegroup: Spacegroup object.

    Example:
        >>> # primitive basis of iron
        >>> basis = Basis.primitive("Fe")
        >>> 
        >>> # cubic lattice parameters
        >>> lattparams = LatticeParameters.cubic(2.85)
        >>> 
        >>> # BCC spacegroup
        >>> spg = Spacegroup("I m -3 m")
        >>> 
        >>> # build the unit cell
        >>> unit_cell = UnitCell(basis, lattparams, spg)
        >>> assert len(unit_cell.atoms) == 2
    """

    def __init__(
        self,
        basis: Basis,
        lattice_parameters: LatticeParameters,
        spacegroup: Spacegroup,
        _build: bool = True,
    ) -> None:
        # initialize superclass
        super().__init__()
        # set attributes
        self._basis = basis
        self._lattice_parameters = lattice_parameters
        self._spacegroup = spacegroup
        # build the unit cell
        if _build:
            self._build()

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_json(cls, s: str) -> 'UnitCell':
        """Initializes from a JSON string."""
        topology = Topology.from_json(s)
        data = json.loads(s)
        basis = Basis.from_json(json.dumps(data["basis"]))
        params = LatticeParameters.from_json(json.dumps(data["lattice_parameters"]))
        spg = Spacegroup(int(data["spacegroup"]))
        # initialize unit cell without placing atoms
        res = cls(basis, params, spg, _build=False)
        res._graph = topology._graph
        return res

    ####################
    #    Properties    #
    ####################

    @property
    def basis(self) -> Basis:
        """Returns the basis."""
        return self._basis

    @property
    def lattice_parameters(self) -> LatticeParameters:
        """Returns the lattice parameters."""
        return self._lattice_parameters

    @property
    def spacegroup(self) -> Spacegroup:
        """Returns the spacegroup."""
        return self._spacegroup

    ########################
    #    Public Methods    #
    ########################

    def to_json(self) -> str:
        """Returns the JSON serialized representation."""
        topology_data = json.loads(super().to_json())
        return json.dumps({
            "basis": json.loads(self.basis.to_json()),
            "lattice_parameters": json.loads(self.lattice_parameters.to_json()),
            "spacegroup": self.spacegroup.international_number,
            "atoms": topology_data["atoms"],
            "bonds": topology_data["bonds"],
        })

    #########################
    #    Private Methods    #
    #########################

    def _build(self) -> None:
        vectors = LatticeVectors.from_lattice_parameters(self.lattice_parameters).vectors
        for specie, site in self.basis.apply_spacegroup(self.spacegroup):
            position = site * np.linalg.norm(vectors, axis=0)
            self.insert_atoms(Atom(specie, position))


class Crystal(Topology):
    """Atomic structure with long range order.
    
    Args:
        unit_cell: Minimal representation of a crystalline structure.
    """

    def __init__(
        self,
        unit_cell: UnitCell,
        _build: bool = True,
    ) -> None:
        # initialize superclass
        super().__init__()
        # initialize private attributes
        self._cut_plane: Optional[Plane] = None
        self._extent: Optional[Tuple[int, int, int]] = None
        self._orientation: Optional[Orientation] = None
        self._orthogonalize: Optional[bool] = None
        self._projection_plane: Optional[Plane] = None
        self._transformation_matrix: Optional[np.ndarray] = None
        # build the crystal
        self._unit_cell = unit_cell
        if _build:
            self._lattice_vectors = LatticeVectors.from_lattice_parameters(self._unit_cell.lattice_parameters)
            for atom in self._unit_cell.atoms:
                self.insert_atoms(copy.deepcopy(atom))
            for bond in self._unit_cell.bonds:
                self.insert_bond(copy.deepcopy(bond))

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_json(cls, s: str) -> 'Crystal':
        """Initializes from a JSON string."""
        topology = Topology.from_json(s)
        data = json.loads(s)
        unit_cell = UnitCell.from_json(json.dumps(data["unit_cell"]))
        res = cls(unit_cell, _build=False)
        res._lattice_vectors = LatticeVectors.from_json(json.dumps(data["lattice_vectors"]))
        res._graph = topology._graph
        return res

    ####################
    #    Properties    #
    ####################

    @property
    def lattice_vectors(self):
        """Returns the lattice vectors."""
        return self._lattice_vectors

    @property
    def unit_cell(self):
        """Returns the unit cell."""
        return self._unit_cell

    ########################
    #    Public Methods    #
    ########################

    def general_transform(self, transformation: np.ndarray) -> 'Crystal':
        """Applies a general transformation matrix to the crystal.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            transformation: Transformation matrix.

        Note:
            The transform is not applied until the `finish` method is called.
        """
        raise NotImplementedError
        self._transformation_matrix = transformation
        return self

    def supercell(self, extent: Tuple[int, int, int]) -> 'Crystal':
        """Creates a supercell by duplicating the crystal in 3 dimensions.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            extent: Number of repeat units in each direction.

        Note:
            The transform is not applied until the `finish` method is called.
        """
        self._extent = extent
        return self

    def orient(self, orientation: Orientation) -> 'Crystal':
        """Changes the crystal's orientation.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            orientation: Crystallographic orientation.

        Note:
            The transform is not applied until the `finish` method is called.
        """
        raise NotImplementedError
        self._orientation = orientation
        return self

    def project(self, plane: Plane, orthogonalize: bool = False) -> 'Crystal':
        """Projects the crystal onto a plane.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            plane: Projection plane.
            orthogonalize: Determines whether or not the projection is represented as an orthogonal lattice.

        Note:
            The transform is not applied until the `finish` method is called.
            Setting `orthogonalize` to True may result in very large structures for acute projections.
        """
        raise NotImplementedError
        self._projection_plane = plane
        self._orthogonalize = orthogonalize
        return self

    def cut(self, plane: Plane) -> 'Crystal':
        """Cuts the crystal along a plane.
        Mutates the crystal and returns a reference to itself to enable method chaining.

        Args:
            plane: Crystallographic plane to cut along.

        Note:
            The transform is not applied until the `finish` method is called.
        """
        raise NotImplementedError
        self._cut_plane = plane
        return self

    def finish(self) -> None:
        """Applies all active transformations to the crystal."""
        # call underlying implementations
        # TODO: determine best order
        self._project()
        self._general_transform()
        self._orient()
        self._supercell()
        self._cut()
        # reset attributes
        self._reset_attributes()

    def reset(self) -> None:
        """Undoes all transformations."""
        super().__init__()
        self._reset_attributes()
        self._lattice_vectors = LatticeVectors.from_lattice_parameters(self.unit_cell.lattice_parameters)
        for atom in self.unit_cell.atoms:
            self.insert_atoms(copy.deepcopy(atom))
        for bond in self.unit_cell.bonds:
            self.insert_bond(copy.deepcopy(bond))

    def to_json(self) -> str:
        """Returns the JSON serialized representation."""
        topology_data = json.loads(super().to_json())
        return json.dumps({
            "unit_cell": json.loads(self.unit_cell.to_json()),
            "lattice_vectors": json.loads(self.lattice_vectors.to_json()),
            "atoms": topology_data["atoms"],
            "bonds": topology_data["bonds"],
        })

    #########################
    #    Private Methods    #
    #########################

    def _reset_attributes(self) -> None:
        self._cut_plane = None
        self._extent = None
        self._orientation = None
        self._orthogonalize = None
        self._projection_plane = None
        self._transformation_matrix = None

    def _general_transform(self) -> None:
        if self._transformation_matrix is None:
            return

    def _supercell(self) -> None:
        if self._extent is None:
            return
        existing_atoms = self.atoms.copy()
        for x in range(self._extent[0]):
            for y in range(self._extent[1]):
                for z in range(self._extent[2]):
                    if x == y == z == 0:
                        continue
                    offset = np.matmul(np.array([x, y, z]), self.lattice_vectors.vectors)
                    for atom in existing_atoms:
                        _atom = copy.deepcopy(atom)
                        _atom.position += offset
                        self.insert_atoms(_atom)
        self.lattice_vectors.vectors *= self._extent

    def _orient(self) -> None:
        if self._orientation is None:
            return

    def _project(self) -> None:
        if self._projection_plane is None:
            return

    def _cut(self) -> None:
        if self._cut_plane is None:
            return
