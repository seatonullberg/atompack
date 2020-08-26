from collections.abc import MutableMapping


class AttributeMap(MutableMapping):
    """A dict like container which converts all keys into attributes.
    
    Example:
        >>> from atompack.util import AttributeMap
        >>>
        >>> attmap = AttributeMap(name="test")
        >>> assert attmap.name == "test"
        >>> assert attmap["name"] == "test"
    """

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __getitem__(self, key):
        return vars(self)[key]

    def __delitem__(self, key):
        delattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __iter__(self):
        return iter(vars(self))

    def __len__(self):
        return len(vars(self))

    def __repr__(self):
        return repr(vars(self))
