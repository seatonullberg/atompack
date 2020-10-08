from collections.abc import Iterable, MutableMapping
from typing import Any

import numpy as np


class Atom(MutableMapping):
    """Dict like object containing arbitrary atomic properties.

    Note:
        End users should not construct atom objects directly.
    
    Args:
        position: 3D position in cartesian space.
        specie: Atomic specie.
    """

    def __init__(self, position: np.ndarray, specie: str, **kwargs) -> None:
        self._attrs = {k: v for k, v in kwargs.items()}
        self._attrs["position"] = position
        self._attrs["specie"] = specie

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
    def position(self) -> np.ndarray:
        return self._attrs["position"]

    @position.setter
    def position(self, value: np.ndarray) -> None:
        self._attrs["position"] = value

    @property
    def specie(self) -> str:
        return self._attrs["specie"]

    @specie.setter
    def specie(self, value: str) -> None:
        self._attrs["specie"] = value
