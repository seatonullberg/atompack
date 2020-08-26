class Bond(object):
    """Container to store metadata about a bond between atoms.
    
    Notes:
        Any `kwargs` passed to `__init__()` are dynamically set as instance variables.
    
    Example:
        >>> from atompack import Bond
        >>> 
        >>> bond = Bond(kind="sigma")
        >>> assert bond.kind == "sigma"
    """

    def __init__(self, **kwargs) -> None:
        # set attributes dynamically
        for k, v in kwargs.items():
            setattr(self, k, v)
