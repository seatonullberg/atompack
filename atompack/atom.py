import json
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

    ######################
    #    Constructors    #
    ######################

    @classmethod
    def from_json(cls, s: str) -> 'Atom':
        """Initializes from a JSON string."""
        data = json.loads(s)
        # process position
        position = data.pop("position")
        if position is None:
            raise ValueError("`position` is a required attribute")
        position = np.array(position)
        # process specie
        specie = data.pop("specie")
        if specie is None:
            raise ValueError("`specie` is a required attribute")
        return cls(position, specie, **data)

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

    ########################
    #    Public Methods    #
    ########################

    def to_json(self) -> str:
        """Returns the JSON serialized representation."""
        _attrs = self._attrs.copy()
        _attrs["position"] = self.position.tolist()
        return json.dumps(_attrs)
