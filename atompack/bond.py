from typing import Optional

import numpy as np

from atompack.atom import Atom
from atompack.util import AttributeMap


class Bond(AttributeMap):
    """Container to store metadata about a bond between atoms.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.
    
    Example:
        >>> from atompack import Bond
        >>> import numpy as np
        >>> 
        >>> a = Atom(np.array([0.5, 0.5, 0.5]))
        >>> b = Atom(np.array([1.0, 2.0, 3.0]))
        >>> bond = Bond(a, b, kind="sigma")
        >>> assert bond.kind == "sigma"
        >>> assert bond["kind"] == "sigma"
        >>> assert np.array_equal(bond.vector, np.array([0.5, 1.5, 2.5]))
    """

    def __init__(self, a: Optional[Atom] = None, b: Optional[Atom] = None, **kwargs) -> None:
        super().__init__(**kwargs)
        if a is None:
            a = Atom()
        if b is None:
            b = Atom()
        self._a = a
        self._b = b

    @property
    def vector(self) -> np.ndarray:
        """Returns the vector between the bond's endpoints."""
        return self._b.position - self._a.position
