from collections.abc import Iterable, MutableMapping
from typing import Any, Tuple

from atompack.atom import Atom


class Bond(MutableMapping):
    """Dict like object containing arbitrary bond properties.
    
    Note:
        End users should not construct bond objects directly.

    Args:
        endpoints: Atoms at each end of the bond.
    """

    #######################################
    #    MutableMapping Implementation    #
    #######################################

    def __init__(self, endpoints: Tuple[Atom, Atom], **kwargs) -> None:
        self._attrs = {k: v for k, v in kwargs.items()}
        self._attrs["endpoints"] = endpoints

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
    def endpoints(self) -> Tuple[Atom, Atom]:
        return self._attrs["endpoints"]

    @endpoints.setter
    def endpoints(self, value: Tuple[Atom, Atom]) -> None:
        self._attrs["endpoints"] = value
