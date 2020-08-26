from typing import Optional

import numpy as np


class Atom(object):
    """Container to store metadata about a single atom.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.

    Args:
        position: Position vector in 3D cartesian space.

    Example:
        >>> from atompack import Atom
        >>> import numpy as np
        >>> 
        >>> atom = Atom(charge=-2)
        >>>
        >>> assert np.array_equal(atom.position, np.zeros(3))
        >>> assert atom.charge == -2
    """

    def __init__(self, position: Optional[np.ndarray] = None, **kwargs) -> None:
        # set attributes dynamically
        for k, v in kwargs.items():
            setattr(self, k, v)
        # process position
        if position is None:
            position = np.zeros(3)
        self.position = position
