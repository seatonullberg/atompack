import json
from collections.abc import Iterable, MutableMapping
from typing import Any, Tuple

from atompack.atom import Atom


class Bond(MutableMapping):
    """Dict like object containing arbitrary bond properties.
    
    Note:
        End users should not construct Bond objects directly.

    Args:
        indices: Index of each atom in the bond.
    """

    def __init__(self, indices: Tuple[int, int], **kwargs) -> None:
        self._attrs = {k: v for k, v in kwargs.items()}
        self._attrs["indices"] = indices

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_json(cls, s: str) -> 'Bond':
        """Initializes from a JSON string."""
        data = json.loads(s)
        # process indices
        indices = data.pop("indices")
        if indices is None:
            raise ValueError("`indices` is a required attribute")
        indices = tuple(indices)
        return cls(indices, **data)

    #######################################
    #    MutableMapping Implementation    #
    #######################################

    def __getitem__(self, key: str) -> Any:
        return self._attrs[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self._attrs[key] = value

    def __delitem__(self, key: str) -> None:
        del self._attrs[key]

    def __iter__(self) -> Iterable:
        return iter(self._attrs)

    def __len__(self) -> int:
        return len(self._attrs)

    ####################
    #    Properties    #
    ####################

    @property
    def indices(self) -> Tuple[int, int]:
        """Returns the index of each atom in the bond."""
        return self._attrs["indices"]

    ########################
    #    Public Methods    #
    ########################

    def to_json(self) -> str:
        """Returns the JSON serialized representation."""
        _attrs = self._attrs.copy()
        _attrs["indices"] = list(self.indices)
        return json.dumps(_attrs)
