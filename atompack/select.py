from typing import Callable, List, Optional

import numpy as np

from atompack.topology import Topology


class Selection(object):
    """Base class for complex atomic queries."""

    def __call__(self, topology) -> List[int]:
        raise NotImplementedError()


class Random(Selection):
    """Representation of a random query of atoms."""

    def __init__(self, number: Optional[int], percentage: Optional[float], func: Optional[Callable]) -> None:
        self._number = number
        self._percentage = percentage
        self._func = func

    def __call__(self, topology: Topology) -> List[int]:
        # get a reference to atoms
        atoms = topology.atoms
        # determine desired number of indices
        n = None
        if self._number is not None:
            n = self._number
        if self._percentage is not None:
            n = int(round(self._percentage * len(atoms)))
        if n is None:
            raise RuntimeError("random selection requires on of `number` of `percentage`")
        # process the callable filter
        func = lambda x: True
        if self._func is not None:
            func = self._func
        # determine the appropriate indices
        all_indices = np.arange(0, len(atoms))
        np.random.shuffle(all_indices)  # randomize
        res: List[int] = []
        for index in all_indices:
            if len(res) == n:
                break
            if func(atoms[index]):
                res.append(index)
        if len(res) < n:
            raise RuntimeError("not enough valid atoms to create the desired selection")
        return res

    @classmethod
    def from_number(cls, number: int, func: Optional[Callable] = None) -> 'Random':
        return cls(number, None, func)

    @classmethod
    def from_percentage(cls, percentage: float, func: Optional[Callable] = None) -> 'Random':
        return cls(None, percentage, func)
