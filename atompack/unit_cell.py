from typing import List, Tuple, Union

import numpy as np

from atompack.lattice import LatticeParameters
from atompack.topology import Topology


class UnitCell(Topology):
    """Minimal representation of a crystalline structure.

    Note:
        End users should not construct unit cell objects directly.
    
    Args:
        basis: Asymmetric site occupancy.
        lattice_parameters: Lattice parameters object.
        spacegroup: Hermann Mauguin spacegroup symbol or international spacegroup number.
    """

    def __init__(self, basis: List[Tuple[str, np.ndarray]], lattice_parameters: LatticeParameters,
                 spacegroup: Union[int, str]) -> None:
        self._basis = basis
        self._lattice_parameters = lattice_parameters
        self._spacegroup = spacegroup
        super().__init__()
        self._build()

    def _build(self) -> None:
        pass

