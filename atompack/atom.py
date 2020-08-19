from typing import Optional

import numpy as np


class Atom(object):
    """Representation of an atom with dynamic attributes.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.
        This enables arbitrary atomic data to be associated with each instance in an
        elegant and pythonic way. `position` is left as an optional positional 
        argument so that each instance is guaranteed to have it as an attribute.

    Args:
        position: Location vector in 3D cartesian space.
            Mutating `position` is a logical error if the change results in atoms 
            overlapping or existing out of bounds.

    Example:
        >>> from atompack import Atom
        >>> import numpy as np
        >>> 
        >>> position = np.array([1, 2, 3])
        >>> atom = Atom(position, symbol="O", charge=-2)
        >>>
        >>> assert np.array_equal(position, atom.position)
        >>> assert atom.symbol == "O"
        >>> assert atom.charge == -2
    """

    def __init__(self, position: Optional[np.ndarray] = None, **kwargs) -> None:
        if position is None:
            position = np.zeros(3)
        self.position = position
        for k, v in kwargs.items():
            setattr(self, k, v)
