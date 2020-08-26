from atompack.util import AttributeMap


class Bond(AttributeMap):
    """Container to store metadata about a bond between atoms.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.
    
    Example:
        >>> from atompack import Bond
        >>> 
        >>> bond = Bond(kind="sigma")
        >>> assert bond.kind == "sigma"
        >>> assert bond["kind"] == "sigma"
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
