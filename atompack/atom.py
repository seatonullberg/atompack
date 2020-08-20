from typing import Optional

import numpy as np


class Atom(object):
    """Representation of an atom with dynamic attributes.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.
        This enables arbitrary atomic data to be associated with each instance in an
        elegant and pythonic way. `position` and `symbol` are left as optional arguments
        with default values so that each instance is guaranteed to have them as 
        attributes for use in internal operations.

    Args:
        position: Location vector in 3D cartesian space.
            Mutating `position` is a logical error if the change results in atoms 
            overlapping or existing out of bounds.
        symbol: IUPAC chemical symbol.

    Example:
        >>> from atompack import Atom
        >>> import numpy as np
        >>> 
        >>> atom = Atom(charge=-2)
        >>>
        >>> assert np.array_equal(atom.position, np.zeros(3))
        >>> assert atom.symbol == "Undefined"
        >>> assert atom.charge == -2
    """

    def __init__(self, position: Optional[np.ndarray] = None, symbol: Optional[str] = None, **kwargs) -> None:
        if position is None:
            position = np.zeros(3)
        self.position = position
        if symbol is None:
            symbol = "Undefined"
        self.symbol = symbol
        for k, v in kwargs.items():
            setattr(self, k, v)
