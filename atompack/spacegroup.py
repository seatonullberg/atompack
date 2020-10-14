from typing import Callable, List, Union

import numpy as np


class Spacegroup(object):
    """Representation of a spacegroup.

    Args:
        spg: Hermann Mauguin symbol or International spacegroup number.
    """

    def __init__(self, spg: Union[int, str]) -> None:
        group = None
        if type(spg) is int:
            if spg > 230 or spg < 1:
                raise ValueError("`spg` must be in range 1..230")
            group = SPACEGROUPS[spg - 1]
        elif type(spg) is str:
            for i, _group in enumerate(SPACEGROUPS):
                if _group["hermann_mauguin"] == spg:
                    group = SPACEGROUPS[i]
                    break
            if group is None:
                raise ValueError("`spg` is not a valid Hermann Mauguin spacegroup symbol")
        else:
            raise TypeError("`spg` must be of type int or str")

        self._bravais_lattice = group["bravais_lattice"]
        self._international_number = group["international_number"]
        self._hermann_mauguin = group["hermann_mauguin"]
        self._genpos = group["genpos"]

    ####################
    #    Properties    #
    ####################

    @property
    def bravais_lattice(self) -> str:
        """Returns the bravais lattice name."""
        return self._bravais_lattice

    @property
    def international_number(self) -> int:
        """Returns the international spacegroup number."""
        return self._international_number

    @property
    def hermann_mauguin(self) -> str:
        """Returns the Hermann Mauguin spacegroup symbol."""
        return self._hermann_mauguin

    @property
    def genpos(self) -> List[str]:
        """Returns the general position expressions."""
        return self._genpos


SPACEGROUPS = [{
    "bravais_lattice": "triclinic",
    "international_number": 1,
    "hermann_mauguin": "P 1",
    "genpos": ["x,y,z"]
}, {
    "bravais_lattice": "triclinic",
    "international_number": 2,
    "hermann_mauguin": "P -1",
    "genpos": ["x,y,z", "-x,-y,-z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 3,
    "hermann_mauguin": "P 1 2 1",
    "genpos": ["x,y,z", "-x,y,-z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 4,
    "hermann_mauguin": "P 1 21 1",
    "genpos": ["x,y,z", "-x,y+1/2,-z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 5,
    "hermann_mauguin": "C 1 2 1",
    "genpos": ["x,y,z", "-x,y,-z", "x+1/2,y+1/2,z", "-x+1/2,y+1/2,-z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 6,
    "hermann_mauguin": "P 1 m 1",
    "genpos": ["x,y,z", "x,-y,z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 7,
    "hermann_mauguin": "P 1 c 1",
    "genpos": ["x,y,z", "x,-y,z+1/2"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 8,
    "hermann_mauguin": "C 1 m 1",
    "genpos": ["x,y,z", "x,-y,z", "x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 9,
    "hermann_mauguin": "C 1 c 1",
    "genpos": ["x,y,z", "x,-y,z+1/2", "x+1/2,y+1/2,z", "x+1/2,-y+1/2,z+1/2"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 10,
    "hermann_mauguin": "P 1 2/m 1",
    "genpos": ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 11,
    "hermann_mauguin": "P 1 21/m 1",
    "genpos": ["x,y,z", "-x,y+1/2,-z", "-x,-y,-z", "x,-y-1/2,z"]
}, {
    "bravais_lattice":
        "monoclinic",
    "international_number":
        12,
    "hermann_mauguin":
        "C 1 2/m 1",
    "genpos": [
        "x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "x+1/2,y+1/2,z", "-x+1/2,y+1/2,-z", "-x+1/2,-y+1/2,-z",
        "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 13,
    "hermann_mauguin": "P 1 2/c 1",
    "genpos": ["x,y,z", "-x,y,-z+1/2", "-x,-y,-z", "x,-y,z-1/2"]
}, {
    "bravais_lattice": "monoclinic",
    "international_number": 14,
    "hermann_mauguin": "P 1 21/c 1",
    "genpos": ["x,y,z", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x,-y-1/2,z-1/2"]
}, {
    "bravais_lattice":
        "monoclinic",
    "international_number":
        15,
    "hermann_mauguin":
        "C 1 2/c 1",
    "genpos": [
        "x,y,z", "-x,y,-z+1/2", "-x,-y,-z", "x,-y,z-1/2", "x+1/2,y+1/2,z", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,-y+1/2,-z",
        "x+1/2,-y+1/2,z-1/2"
    ]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 16,
    "hermann_mauguin": "P 2 2 2",
    "genpos": ["x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 17,
    "hermann_mauguin": "P 2 2 21",
    "genpos": ["x,y,z", "-x,-y,z+1/2", "x,-y,-z", "-x,y,-z+1/2"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 18,
    "hermann_mauguin": "P 21 21 2",
    "genpos": ["x,y,z", "-x,-y,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 19,
    "hermann_mauguin": "P 21 21 21",
    "genpos": ["x,y,z", "-x+1/2,-y,z+1/2", "x+1/2,-y+1/2,-z", "-x,y+1/2,-z+1/2"]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        20,
    "hermann_mauguin":
        "C 2 2 21",
    "genpos": [
        "x,y,z", "-x,-y,z+1/2", "x,-y,-z", "-x,y,-z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        21,
    "hermann_mauguin":
        "C 2 2 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,-z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        22,
    "hermann_mauguin":
        "F 2 2 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2",
        "-x,y+1/2,-z+1/2", "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2", "-x+1/2,y,-z+1/2", "x+1/2,y+1/2,z",
        "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        23,
    "hermann_mauguin":
        "I 2 2 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2",
        "-x+1/2,y+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        24,
    "hermann_mauguin":
        "I 21 21 21",
    "genpos": [
        "x,y,z", "-x,-y+1/2,z", "x,-y,-z+1/2", "-x,y+1/2,-z+1/2", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1,z+1/2",
        "x+1/2,-y+1/2,-z+1", "-x+1/2,y+1,-z+1"
    ]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 25,
    "hermann_mauguin": "P m m 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 26,
    "hermann_mauguin": "P m c 21",
    "genpos": ["x,y,z", "-x,-y,z+1/2", "-x,y,z", "x,-y,z+1/2"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 27,
    "hermann_mauguin": "P c c 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x,y,z+1/2", "x,-y,z+1/2"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 28,
    "hermann_mauguin": "P m a 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x+1/2,y,z", "x+1/2,-y,z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 29,
    "hermann_mauguin": "P c a 21",
    "genpos": ["x,y,z", "-x,-y,z+1/2", "-x+1/2,y,z+1/2", "x+1/2,-y,z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 30,
    "hermann_mauguin": "P n c 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x,y+1/2,z+1/2", "x,-y+1/2,z+1/2"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 31,
    "hermann_mauguin": "P m n 21",
    "genpos": ["x,y,z", "-x+1/2,-y,z+1/2", "-x,y,z", "x+1/2,-y,z+1/2"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 32,
    "hermann_mauguin": "P b a 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 33,
    "hermann_mauguin": "P n a 21",
    "genpos": ["x,y,z", "-x,-y,z+1/2", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y+1/2,z"]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 34,
    "hermann_mauguin": "P n n 2",
    "genpos": ["x,y,z", "-x,-y,z", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2"]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        35,
    "hermann_mauguin":
        "C m m 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        36,
    "hermann_mauguin":
        "C m c 21",
    "genpos": [
        "x,y,z", "-x,-y,z+1/2", "-x,y,z", "x,-y,z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,z",
        "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        37,
    "hermann_mauguin":
        "C c c 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z+1/2", "x,-y,z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "-x+1/2,y+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        38,
    "hermann_mauguin":
        "A m m 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "-x,y+1/2,z+1/2", "x,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        39,
    "hermann_mauguin":
        "A e m 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y+1/2,z", "x,-y+1/2,z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "-x,y+1,z+1/2",
        "x,-y+1,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        40,
    "hermann_mauguin":
        "A m a 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x+1/2,y,z", "x+1/2,-y,z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        41,
    "hermann_mauguin":
        "A e a 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "-x+1/2,y+1,z+1/2",
        "x+1/2,-y+1,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        42,
    "hermann_mauguin":
        "F m m 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "-x,y+1/2,z+1/2", "x,-y+1/2,z+1/2",
        "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "-x+1/2,y,z+1/2", "x+1/2,-y,z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z",
        "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        43,
    "hermann_mauguin":
        "F d d 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x+1/4,y+1/4,z+1/4", "x+3/4,-y+3/4,z+1/4", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2",
        "-x+1/4,y+3/4,z+3/4", "x+3/4,-y+5/4,z+3/4", "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "-x+3/4,y+1/4,z+3/4",
        "x+5/4,-y+3/4,z+3/4", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "-x+3/4,y+3/4,z+1/4", "x+5/4,-y+5/4,z+1/4"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        44,
    "hermann_mauguin":
        "I m m 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        45,
    "hermann_mauguin":
        "I b a 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x,y,z+1/2", "x,-y,z+1/2", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "-x+1/2,y+1/2,z+1",
        "x+1/2,-y+1/2,z+1"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        46,
    "hermann_mauguin":
        "I m a 2",
    "genpos": [
        "x,y,z", "-x,-y,z", "-x+1/2,y,z", "x+1/2,-y,z", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "-x+1,y+1/2,z+1/2",
        "x+1,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 47,
    "hermann_mauguin": "P m m m",
    "genpos": ["x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z"]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        48,
    "hermann_mauguin":
        "P n n n",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x+1/2,-y+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 49,
    "hermann_mauguin": "P c c m",
    "genpos": ["x,y,z", "-x,-y,z", "x,-y,-z+1/2", "-x,y,-z+1/2", "-x,-y,-z", "x,y,-z", "-x,y,z-1/2", "x,-y,z-1/2"]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        50,
    "hermann_mauguin":
        "P b a n",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x+1/2,-y+1/2,-z", "x+1/2,y+1/2,-z", "-x+1/2,y+1/2,z",
        "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice": "orthorhombic",
    "international_number": 51,
    "hermann_mauguin": "P m m a",
    "genpos": ["x,y,z", "-x+1/2,-y,z", "x+1/2,-y,-z", "-x,y,-z", "-x,-y,-z", "x-1/2,y,-z", "-x-1/2,y,z", "x,-y,z"]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        52,
    "hermann_mauguin":
        "P n n a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z", "x,-y+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x-1/2,y,-z", "-x,y-1/2,z-1/2",
        "x-1/2,-y-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        53,
    "hermann_mauguin":
        "P m n a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x,-y,-z", "-x+1/2,y,-z+1/2", "-x,-y,-z", "x-1/2,y,-z-1/2", "-x,y,z",
        "x-1/2,-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        54,
    "hermann_mauguin":
        "P c c a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z", "x+1/2,-y,-z+1/2", "-x,y,-z+1/2", "-x,-y,-z", "x-1/2,y,-z", "-x-1/2,y,z-1/2",
        "x,-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        55,
    "hermann_mauguin":
        "P b a m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z", "-x,-y,-z", "x,y,-z", "-x-1/2,y-1/2,z",
        "x-1/2,-y-1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        56,
    "hermann_mauguin":
        "P c c n",
    "genpos": [
        "x,y,z", "-x+1/2,-y+1/2,z", "x+1/2,-y,-z+1/2", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x-1/2,y-1/2,-z",
        "-x-1/2,y,z-1/2", "x,-y-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        57,
    "hermann_mauguin":
        "P b c m",
    "genpos": [
        "x,y,z", "-x,-y,z+1/2", "x,-y+1/2,-z", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x,y,-z-1/2", "-x,y-1/2,z",
        "x,-y-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        58,
    "hermann_mauguin":
        "P n n m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x+1/2,-y+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x,y,-z", "-x-1/2,y-1/2,z-1/2",
        "x-1/2,-y-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        59,
    "hermann_mauguin":
        "P m m n",
    "genpos": [
        "x,y,z", "-x,-y,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z", "-x+1/2,-y+1/2,-z", "x+1/2,y+1/2,-z", "-x,y,z",
        "x,-y,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        60,
    "hermann_mauguin":
        "P b c n",
    "genpos": [
        "x,y,z", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z", "-x,y,-z+1/2", "-x,-y,-z", "x-1/2,y-1/2,-z-1/2",
        "-x-1/2,y-1/2,z", "x,-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        61,
    "hermann_mauguin":
        "P b c a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x+1/2,-y+1/2,-z", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x-1/2,y,-z-1/2",
        "-x-1/2,y-1/2,z", "x,-y-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        62,
    "hermann_mauguin":
        "P n m a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x+1/2,-y+1/2,-z+1/2", "-x,y+1/2,-z", "-x,-y,-z", "x-1/2,y,-z-1/2",
        "-x-1/2,y-1/2,z-1/2", "x,-y-1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        63,
    "hermann_mauguin":
        "C m c m",
    "genpos": [
        "x,y,z", "-x,-y,z+1/2", "x,-y,-z", "-x,y,-z+1/2", "-x,-y,-z", "x,y,-z-1/2", "-x,y,z", "x,-y,z-1/2",
        "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,-y+1/2,-z",
        "x+1/2,y+1/2,-z-1/2", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        64,
    "hermann_mauguin":
        "C m c e",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x,-y,-z", "-x+1/2,y,-z+1/2", "-x,-y,-z", "x-1/2,y,-z-1/2", "-x,y,z",
        "x-1/2,-y,z-1/2", "x+1/2,y+1/2,z", "-x+1,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z", "-x+1,y+1/2,-z+1/2",
        "-x+1/2,-y+1/2,-z", "x,y+1/2,-z-1/2", "-x+1/2,y+1/2,z", "x,-y+1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        65,
    "hermann_mauguin":
        "C m m m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "x+1/2,y+1/2,z",
        "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z", "-x+1/2,-y+1/2,-z", "x+1/2,y+1/2,-z", "-x+1/2,y+1/2,z",
        "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        66,
    "hermann_mauguin":
        "C c c m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z+1/2", "-x,y,-z+1/2", "-x,-y,-z", "x,y,-z", "-x,y,z-1/2", "x,-y,z-1/2",
        "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,-y+1/2,-z",
        "x+1/2,y+1/2,-z", "-x+1/2,y+1/2,z-1/2", "x+1/2,-y+1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        67,
    "hermann_mauguin":
        "C m m a",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z", "x,-y,-z", "-x+1/2,y,-z", "-x,-y,-z", "x-1/2,y,-z", "-x,y,z", "x-1/2,-y,z",
        "x+1/2,y+1/2,z", "-x+1,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1,y+1/2,-z", "-x+1/2,-y+1/2,-z", "x,y+1/2,-z",
        "-x+1/2,y+1/2,z", "x,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        68,
    "hermann_mauguin":
        "C c c e",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x+1/2,-y,-z+1/2", "x+1/2,y,-z+1/2", "-x+1/2,y,z+1/2",
        "x+1/2,-y,z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z",
        "-x+1,-y+1/2,-z+1/2", "x+1,y+1/2,-z+1/2", "-x+1,y+1/2,z+1/2", "x+1,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        69,
    "hermann_mauguin":
        "F m m m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "x,y+1/2,z+1/2",
        "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2", "-x,y+1/2,-z+1/2", "-x,-y+1/2,-z+1/2", "x,y+1/2,-z+1/2", "-x,y+1/2,z+1/2",
        "x,-y+1/2,z+1/2", "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2", "-x+1/2,y,-z+1/2", "-x+1/2,-y,-z+1/2",
        "x+1/2,y,-z+1/2", "-x+1/2,y,z+1/2", "x+1/2,-y,z+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,-z", "-x+1/2,-y+1/2,-z", "x+1/2,y+1/2,-z", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        70,
    "hermann_mauguin":
        "F d d d",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x+1/4,-y+1/4,-z+1/4", "x+1/4,y+1/4,-z+1/4", "-x+1/4,y+1/4,z+1/4",
        "x+1/4,-y+1/4,z+1/4", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2", "-x,y+1/2,-z+1/2",
        "-x+1/4,-y+3/4,-z+3/4", "x+1/4,y+3/4,-z+3/4", "-x+1/4,y+3/4,z+3/4", "x+1/4,-y+3/4,z+3/4", "x+1/2,y,z+1/2",
        "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2", "-x+1/2,y,-z+1/2", "-x+3/4,-y+1/4,-z+3/4", "x+3/4,y+1/4,-z+3/4",
        "-x+3/4,y+1/4,z+3/4", "x+3/4,-y+1/4,z+3/4", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,-z", "-x+3/4,-y+3/4,-z+1/4", "x+3/4,y+3/4,-z+1/4", "-x+3/4,y+3/4,z+1/4", "x+3/4,-y+3/4,z+1/4"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        71,
    "hermann_mauguin":
        "I m m m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "x+1/2,y+1/2,z+1/2",
        "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2", "-x+1/2,-y+1/2,-z+1/2",
        "x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        72,
    "hermann_mauguin":
        "I b a m",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z+1/2", "-x,y,-z+1/2", "-x,-y,-z", "x,y,-z", "-x,y,z-1/2", "x,-y,z-1/2",
        "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z+1", "-x+1/2,y+1/2,-z+1", "-x+1/2,-y+1/2,-z+1/2",
        "x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z", "x+1/2,-y+1/2,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        73,
    "hermann_mauguin":
        "I b c a",
    "genpos": [
        "x,y,z", "-x,-y+1/2,z", "x,-y,-z+1/2", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x,y-1/2,-z", "-x,y,z-1/2",
        "x,-y-1/2,z-1/2", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1,z+1/2", "x+1/2,-y+1/2,-z+1", "-x+1/2,y+1,-z+1",
        "-x+1/2,-y+1/2,-z+1/2", "x+1/2,y,-z+1/2", "-x+1/2,y+1/2,z", "x+1/2,-y,z"
    ]
}, {
    "bravais_lattice":
        "orthorhombic",
    "international_number":
        74,
    "hermann_mauguin":
        "I m m a",
    "genpos": [
        "x,y,z", "-x,-y+1/2,z", "x,-y,-z", "-x,y+1/2,-z", "-x,-y,-z", "x,y-1/2,-z", "-x,y,z", "x,-y-1/2,z",
        "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1,z+1/2", "x+1/2,-y+1/2,-z+1/2", "-x+1/2,y+1,-z+1/2", "-x+1/2,-y+1/2,-z+1/2",
        "x+1/2,y,-z+1/2", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y,z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 75,
    "hermann_mauguin": "P 4",
    "genpos": ["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 76,
    "hermann_mauguin": "P 41",
    "genpos": ["x,y,z", "-y,x,z+1/4", "-x,-y,z+1/2", "y,-x,z+3/4"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 77,
    "hermann_mauguin": "P 42",
    "genpos": ["x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 78,
    "hermann_mauguin": "P 43",
    "genpos": ["x,y,z", "-y,x,z+3/4", "-x,-y,z+1/2", "y,-x,z+1/4"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        79,
    "hermann_mauguin":
        "I 4",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2",
        "y+1/2,-x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        80,
    "hermann_mauguin":
        "I 41",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4",
        "-x+1,-y+1,z+1", "y+1,-x+1/2,z+5/4"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 81,
    "hermann_mauguin": "P -4",
    "genpos": ["x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        82,
    "hermann_mauguin":
        "I -4",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x+1/2,y+1/2,z+1/2", "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2",
        "-y+1/2,x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 83,
    "hermann_mauguin": "P 4/m",
    "genpos": ["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 84,
    "hermann_mauguin": "P 42/m",
    "genpos": ["x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "-x,-y,-z", "y,-x,-z-1/2", "x,y,-z", "-y,x,-z-1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        85,
    "hermann_mauguin":
        "P 4/n",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z", "-x,-y,z", "y+1/2,-x+1/2,z", "-x+1/2,-y+1/2,-z", "y,-x,-z", "x+1/2,y+1/2,-z",
        "-y,x,-z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        86,
    "hermann_mauguin":
        "P 42/n",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "-x+1/2,-y+1/2,-z+1/2", "y,-x,-z",
        "x+1/2,y+1/2,-z+1/2", "-y,x,-z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        87,
    "hermann_mauguin":
        "I 4/m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "x+1/2,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "-x+1/2,-y+1/2,-z+1/2",
        "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-y+1/2,x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        88,
    "hermann_mauguin":
        "I 41/a",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "-x,-y+1/2,-z+1/4", "y,-x,-z",
        "x-1/2,y,-z-1/4", "-y-1/2,x+1/2,-z-1/2", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4", "-x+1,-y+1,z+1",
        "y+1,-x+1/2,z+5/4", "-x+1/2,-y+1,-z+3/4", "y+1/2,-x+1/2,-z+1/2", "x,y+1/2,-z+1/4", "-y,x+1,-z"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 89,
    "hermann_mauguin": "P 4 2 2",
    "genpos": ["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        90,
    "hermann_mauguin":
        "P 4 21 2",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z", "-x,-y,z", "y+1/2,-x+1/2,z", "x+1/2,-y+1/2,-z", "y,x,-z", "-x+1/2,y+1/2,-z",
        "-y,-x,-z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        91,
    "hermann_mauguin":
        "P 41 2 2",
    "genpos": [
        "x,y,z", "-y,x,z+1/4", "-x,-y,z+1/2", "y,-x,z+3/4", "x,-y,-z+1/2", "y,x,-z+3/4", "-x,y,-z", "-y,-x,-z+1/4"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        92,
    "hermann_mauguin":
        "P 41 21 2",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/4", "-x,-y,z+1/2", "y+1/2,-x+1/2,z+3/4", "x+1/2,-y+1/2,-z+3/4", "y,x,-z",
        "-x+1/2,y+1/2,-z+1/4", "-y,-x,-z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 93,
    "hermann_mauguin": "P 42 2 2",
    "genpos": ["x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "x,-y,-z", "y,x,-z+1/2", "-x,y,-z", "-y,-x,-z+1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        94,
    "hermann_mauguin":
        "P 42 21 2",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y,x,-z",
        "-x+1/2,y+1/2,-z+1/2", "-y,-x,-z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        95,
    "hermann_mauguin":
        "P 43 2 2",
    "genpos": [
        "x,y,z", "-y,x,z+3/4", "-x,-y,z+1/2", "y,-x,z+1/4", "x,-y,-z+1/2", "y,x,-z+1/4", "-x,y,-z", "-y,-x,-z+3/4"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        96,
    "hermann_mauguin":
        "P 43 21 2",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+3/4", "-x,-y,z+1/2", "y+1/2,-x+1/2,z+1/4", "x+1/2,-y+1/2,-z+1/4", "y,x,-z",
        "-x+1/2,y+1/2,-z+3/4", "-y,-x,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        97,
    "hermann_mauguin":
        "I 4 2 2",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "x+1/2,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y+1/2,-z+1/2", "-y+1/2,-x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        98,
    "hermann_mauguin":
        "I 41 2 2",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "x,-y+1/2,-z+1/4", "y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y,-z+3/4", "-y,-x,-z", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4", "-x+1,-y+1,z+1", "y+1,-x+1/2,z+5/4",
        "x+1/2,-y+1,-z+3/4", "y+1,x+1,-z+1", "-x+1,y+1/2,-z+5/4", "-y+1/2,-x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 99,
    "hermann_mauguin": "P 4 m m",
    "genpos": ["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,y,z", "-y,-x,z", "x,-y,z", "y,x,z"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        100,
    "hermann_mauguin":
        "P 4 b m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x+1/2,y+1/2,z", "-y+1/2,-x+1/2,z", "x+1/2,-y+1/2,z", "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 101,
    "hermann_mauguin": "P 42 c m",
    "genpos": ["x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "-x,y,z+1/2", "-y,-x,z", "x,-y,z+1/2", "y,x,z"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        102,
    "hermann_mauguin":
        "P 42 n m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2", "-y,-x,z",
        "x+1/2,-y+1/2,z+1/2", "y,x,z"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 103,
    "hermann_mauguin": "P 4 c c",
    "genpos": ["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,y,z+1/2", "-y,-x,z+1/2", "x,-y,z+1/2", "y,x,z+1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        104,
    "hermann_mauguin":
        "P 4 n c",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2",
        "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 105,
    "hermann_mauguin": "P 42 m c",
    "genpos": ["x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "-x,y,z", "-y,-x,z+1/2", "x,-y,z", "y,x,z+1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        106,
    "hermann_mauguin":
        "P 42 b c",
    "genpos": [
        "x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "-x+1/2,y+1/2,z", "-y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,z",
        "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        107,
    "hermann_mauguin":
        "I 4 m m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,y,z", "-y,-x,z", "x,-y,z", "y,x,z", "x+1/2,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2", "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        108,
    "hermann_mauguin":
        "I 4 c m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "-x,y,z+1/2", "-y,-x,z+1/2", "x,-y,z+1/2", "y,x,z+1/2",
        "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,z+1",
        "-y+1/2,-x+1/2,z+1", "x+1/2,-y+1/2,z+1", "y+1/2,x+1/2,z+1"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        109,
    "hermann_mauguin":
        "I 41 m d",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "-x,y,z", "-y,-x+1/2,z+1/4",
        "x+1/2,-y+1/2,z+1/2", "y+1/2,x,z+3/4", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4", "-x+1,-y+1,z+1",
        "y+1,-x+1/2,z+5/4", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1,z+3/4", "x+1,-y+1,z+1", "y+1,x+1/2,z+5/4"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        110,
    "hermann_mauguin":
        "I 41 c d",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "-x,y,z+1/2", "-y,-x+1/2,z+3/4",
        "x+1/2,-y+1/2,z", "y+1/2,x,z+1/4", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4", "-x+1,-y+1,z+1", "y+1,-x+1/2,z+5/4",
        "-x+1/2,y+1/2,z+1", "-y+1/2,-x+1,z+5/4", "x+1,-y+1,z+1/2", "y+1,x+1/2,z+3/4"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 111,
    "hermann_mauguin": "P -4 2 m",
    "genpos": ["x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z", "-y,-x,z", "-x,y,-z", "y,x,z"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 112,
    "hermann_mauguin": "P -4 2 c",
    "genpos": ["x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z+1/2", "-y,-x,z+1/2", "-x,y,-z+1/2", "y,x,z+1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        113,
    "hermann_mauguin":
        "P -4 21 m",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x+1/2,-y+1/2,-z", "-y+1/2,-x+1/2,z", "-x+1/2,y+1/2,-z",
        "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        114,
    "hermann_mauguin":
        "P -4 21 c",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x+1/2,-y+1/2,-z+1/2", "-y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 115,
    "hermann_mauguin": "P -4 m 2",
    "genpos": ["x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x,y,z", "y,x,-z", "x,-y,z", "-y,-x,-z"]
}, {
    "bravais_lattice": "tetragonal",
    "international_number": 116,
    "hermann_mauguin": "P -4 c 2",
    "genpos": ["x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x,y,z+1/2", "y,x,-z+1/2", "x,-y,z+1/2", "-y,-x,-z+1/2"]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        117,
    "hermann_mauguin":
        "P -4 b 2",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x+1/2,y+1/2,z", "y+1/2,x+1/2,-z", "x+1/2,-y+1/2,z",
        "-y+1/2,-x+1/2,-z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        118,
    "hermann_mauguin":
        "P -4 n 2",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x+1/2,y+1/2,z+1/2", "y+1/2,x+1/2,-z+1/2", "x+1/2,-y+1/2,z+1/2",
        "-y+1/2,-x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        119,
    "hermann_mauguin":
        "I -4 m 2",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x,y,z", "y,x,-z", "x,-y,z", "-y,-x,-z", "x+1/2,y+1/2,z+1/2",
        "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2", "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "y+1/2,x+1/2,-z+1/2",
        "x+1/2,-y+1/2,z+1/2", "-y+1/2,-x+1/2,-z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        120,
    "hermann_mauguin":
        "I -4 c 2",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "-x,y,z+1/2", "y,x,-z+1/2", "x,-y,z+1/2", "-y,-x,-z+1/2",
        "x+1/2,y+1/2,z+1/2", "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2", "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z+1",
        "y+1/2,x+1/2,-z+1", "x+1/2,-y+1/2,z+1", "-y+1/2,-x+1/2,-z+1"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        121,
    "hermann_mauguin":
        "I -4 2 m",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z", "-y,-x,z", "-x,y,-z", "y,x,z", "x+1/2,y+1/2,z+1/2",
        "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2", "-y+1/2,x+1/2,-z+1/2", "x+1/2,-y+1/2,-z+1/2",
        "-y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2", "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        122,
    "hermann_mauguin":
        "I -4 2 d",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y+1/2,-z+1/4", "-y+1/2,-x,z+3/4", "-x,y+1/2,-z+1/4",
        "y+1/2,x,z+3/4", "x+1/2,y+1/2,z+1/2", "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2", "-y+1/2,x+1/2,-z+1/2",
        "x+1/2,-y+1,-z+3/4", "-y+1,-x+1/2,z+5/4", "-x+1/2,y+1,-z+3/4", "y+1,x+1/2,z+5/4"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        123,
    "hermann_mauguin":
        "P 4/m m m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "-x,-y,-z", "y,-x,-z",
        "x,y,-z", "-y,x,-z", "-x,y,z", "-y,-x,z", "x,-y,z", "y,x,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        124,
    "hermann_mauguin":
        "P 4/m c c",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z+1/2", "y,x,-z+1/2", "-x,y,-z+1/2", "-y,-x,-z+1/2", "-x,-y,-z",
        "y,-x,-z", "x,y,-z", "-y,x,-z", "-x,y,z-1/2", "-y,-x,z-1/2", "x,-y,z-1/2", "y,x,z-1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        125,
    "hermann_mauguin":
        "P 4/n b m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "-x+1/2,-y+1/2,-z",
        "y+1/2,-x+1/2,-z", "x+1/2,y+1/2,-z", "-y+1/2,x+1/2,-z", "-x+1/2,y+1/2,z", "-y+1/2,-x+1/2,z", "x+1/2,-y+1/2,z",
        "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        126,
    "hermann_mauguin":
        "P 4/n n c",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "-x+1/2,-y+1/2,-z+1/2",
        "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2",
        "x+1/2,-y+1/2,z+1/2", "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        127,
    "hermann_mauguin":
        "P 4/m b m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x+1/2,-y+1/2,-z", "y+1/2,x+1/2,-z", "-x+1/2,y+1/2,-z",
        "-y+1/2,-x+1/2,-z", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "-x-1/2,y-1/2,z", "-y-1/2,-x-1/2,z",
        "x-1/2,-y-1/2,z", "y-1/2,x-1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        128,
    "hermann_mauguin":
        "P 4/m n c",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x+1/2,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "-y+1/2,-x+1/2,-z+1/2", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "-x-1/2,y-1/2,z-1/2", "-y-1/2,-x-1/2,z-1/2",
        "x-1/2,-y-1/2,z-1/2", "y-1/2,x-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        129,
    "hermann_mauguin":
        "P 4/n m m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z", "-x,-y,z", "y+1/2,-x+1/2,z", "x+1/2,-y+1/2,-z", "y,x,-z", "-x+1/2,y+1/2,-z",
        "-y,-x,-z", "-x+1/2,-y+1/2,-z", "y,-x,-z", "x+1/2,y+1/2,-z", "-y,x,-z", "-x,y,z", "-y+1/2,-x+1/2,z", "x,-y,z",
        "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        130,
    "hermann_mauguin":
        "P 4/n c c",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z", "-x,-y,z", "y+1/2,-x+1/2,z", "x+1/2,-y+1/2,-z+1/2", "y,x,-z+1/2",
        "-x+1/2,y+1/2,-z+1/2", "-y,-x,-z+1/2", "-x+1/2,-y+1/2,-z", "y,-x,-z", "x+1/2,y+1/2,-z", "-y,x,-z", "-x,y,z-1/2",
        "-y+1/2,-x+1/2,z-1/2", "x,-y,z-1/2", "y+1/2,x+1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        131,
    "hermann_mauguin":
        "P 42/m m c",
    "genpos": [
        "x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "x,-y,-z", "y,x,-z+1/2", "-x,y,-z", "-y,-x,-z+1/2", "-x,-y,-z",
        "y,-x,-z-1/2", "x,y,-z", "-y,x,-z-1/2", "-x,y,z", "-y,-x,z-1/2", "x,-y,z", "y,x,z-1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        132,
    "hermann_mauguin":
        "P 42/m c m",
    "genpos": [
        "x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "x,-y,-z+1/2", "y,x,-z", "-x,y,-z+1/2", "-y,-x,-z", "-x,-y,-z",
        "y,-x,-z-1/2", "x,y,-z", "-y,x,-z-1/2", "-x,y,z-1/2", "-y,-x,z", "x,-y,z-1/2", "y,x,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        133,
    "hermann_mauguin":
        "P 42/n b c",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x,-y,-z+1/2", "y+1/2,x+1/2,-z", "-x,y,-z+1/2",
        "-y+1/2,-x+1/2,-z", "-x+1/2,-y+1/2,-z+1/2", "y,-x,-z", "x+1/2,y+1/2,-z+1/2", "-y,x,-z", "-x+1/2,y+1/2,z",
        "-y,-x,z+1/2", "x+1/2,-y+1/2,z", "y,x,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        134,
    "hermann_mauguin":
        "P 42/n n m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x,-y,-z", "y+1/2,x+1/2,-z+1/2", "-x,y,-z",
        "-y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,-z+1/2", "y,-x,-z", "x+1/2,y+1/2,-z+1/2", "-y,x,-z",
        "-x+1/2,y+1/2,z+1/2", "-y,-x,z", "x+1/2,-y+1/2,z+1/2", "y,x,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        135,
    "hermann_mauguin":
        "P 42/m b c",
    "genpos": [
        "x,y,z", "-y,x,z+1/2", "-x,-y,z", "y,-x,z+1/2", "x+1/2,-y+1/2,-z", "y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,-z",
        "-y+1/2,-x+1/2,-z+1/2", "-x,-y,-z", "y,-x,-z-1/2", "x,y,-z", "-y,x,-z-1/2", "-x-1/2,y-1/2,z",
        "-y-1/2,-x-1/2,z-1/2", "x-1/2,-y-1/2,z", "y-1/2,x-1/2,z-1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        136,
    "hermann_mauguin":
        "P 42/m n m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y,x,-z",
        "-x+1/2,y+1/2,-z+1/2", "-y,-x,-z", "-x,-y,-z", "y-1/2,-x-1/2,-z-1/2", "x,y,-z", "-y-1/2,x-1/2,-z-1/2",
        "-x-1/2,y-1/2,z-1/2", "-y,-x,z", "x-1/2,-y-1/2,z-1/2", "y,x,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        137,
    "hermann_mauguin":
        "P 42/n m c",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y,x,-z",
        "-x+1/2,y+1/2,-z+1/2", "-y,-x,-z", "-x+1/2,-y+1/2,-z+1/2", "y,-x,-z", "x+1/2,y+1/2,-z+1/2", "-y,x,-z", "-x,y,z",
        "-y+1/2,-x+1/2,z+1/2", "x,-y,z", "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        138,
    "hermann_mauguin":
        "P 42/n c m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z", "y,x,-z+1/2",
        "-x+1/2,y+1/2,-z", "-y,-x,-z+1/2", "-x+1/2,-y+1/2,-z+1/2", "y,-x,-z", "x+1/2,y+1/2,-z+1/2", "-y,x,-z",
        "-x,y,z+1/2", "-y+1/2,-x+1/2,z", "x,-y,z+1/2", "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        139,
    "hermann_mauguin":
        "I 4/m m m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "-x,-y,-z", "y,-x,-z",
        "x,y,-z", "-y,x,-z", "-x,y,z", "-y,-x,z", "x,-y,z", "y,x,z", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1/2,z+1/2",
        "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "-y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,-z+1/2", "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2",
        "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2", "y+1/2,x+1/2,z+1/2"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        140,
    "hermann_mauguin":
        "I 4/m c m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z+1/2", "y,x,-z+1/2", "-x,y,-z+1/2", "-y,-x,-z+1/2", "-x,-y,-z",
        "y,-x,-z", "x,y,-z", "-y,x,-z", "-x,y,z-1/2", "-y,-x,z-1/2", "x,-y,z-1/2", "y,x,z-1/2", "x+1/2,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1", "y+1/2,x+1/2,-z+1",
        "-x+1/2,y+1/2,-z+1", "-y+1/2,-x+1/2,-z+1", "-x+1/2,-y+1/2,-z+1/2", "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2",
        "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z", "-y+1/2,-x+1/2,z", "x+1/2,-y+1/2,z", "y+1/2,x+1/2,z"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        141,
    "hermann_mauguin":
        "I 41/a m d",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "x,-y+1/2,-z+1/4", "y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y,-z+3/4", "-y,-x,-z", "-x,-y+1/2,-z+1/4", "y,-x,-z", "x-1/2,y,-z-1/4", "-y-1/2,x+1/2,-z-1/2", "-x,y,z",
        "-y-1/2,-x,z-1/4", "x-1/2,-y+1/2,z-1/2", "y,x+1/2,z+1/4", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4",
        "-x+1,-y+1,z+1", "y+1,-x+1/2,z+5/4", "x+1/2,-y+1,-z+3/4", "y+1,x+1,-z+1", "-x+1,y+1/2,-z+5/4",
        "-y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1,-z+3/4", "y+1/2,-x+1/2,-z+1/2", "x,y+1/2,-z+1/4", "-y,x+1,-z",
        "-x+1/2,y+1/2,z+1/2", "-y,-x+1/2,z+1/4", "x,-y+1,z", "y+1/2,x+1,z+3/4"
    ]
}, {
    "bravais_lattice":
        "tetragonal",
    "international_number":
        142,
    "hermann_mauguin":
        "I 41/a c d",
    "genpos": [
        "x,y,z", "-y,x+1/2,z+1/4", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x,z+3/4", "x+1/2,-y,-z+1/4", "y,x,-z+1/2",
        "-x,y+1/2,-z+3/4", "-y+1/2,-x+1/2,-z", "-x,-y+1/2,-z+1/4", "y,-x,-z", "x-1/2,y,-z-1/4", "-y-1/2,x+1/2,-z-1/2",
        "-x-1/2,y+1/2,z", "-y,-x+1/2,z-1/4", "x,-y,z-1/2", "y-1/2,x,z+1/4", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1,z+3/4",
        "-x+1,-y+1,z+1", "y+1,-x+1/2,z+5/4", "x+1,-y+1/2,-z+3/4", "y+1/2,x+1/2,-z+1", "-x+1/2,y+1,-z+5/4",
        "-y+1,-x+1,-z+1/2", "-x+1/2,-y+1,-z+3/4", "y+1/2,-x+1/2,-z+1/2", "x,y+1/2,-z+1/4", "-y,x+1,-z", "-x,y+1,z+1/2",
        "-y+1/2,-x+1,z+1/4", "x+1/2,-y+1/2,z", "y,x+1/2,z+3/4"
    ]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 143,
    "hermann_mauguin": "P 3",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 144,
    "hermann_mauguin": "P 31",
    "genpos": ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 145,
    "hermann_mauguin": "P 32",
    "genpos": ["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3"]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        146,
    "hermann_mauguin":
        "R 3",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "x+2/3,y+1/3,z+1/3", "-y+2/3,x-y+1/3,z+1/3", "-x+y+2/3,-x+1/3,z+1/3",
        "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3"
    ]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 147,
    "hermann_mauguin": "P -3",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z"]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        148,
    "hermann_mauguin":
        "R -3",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z", "x+2/3,y+1/3,z+1/3",
        "-y+2/3,x-y+1/3,z+1/3", "-x+y+2/3,-x+1/3,z+1/3", "-x+2/3,-y+1/3,-z+1/3", "y+2/3,-x+y+1/3,-z+1/3",
        "x-y+2/3,x+1/3,-z+1/3", "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3",
        "-x+1/3,-y+2/3,-z+2/3", "y+1/3,-x+y+2/3,-z+2/3", "x-y+1/3,x+2/3,-z+2/3"
    ]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 149,
    "hermann_mauguin": "P 3 1 2",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,-z", "x,x-y,-z", "-x+y,y,-z"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 150,
    "hermann_mauguin": "P 3 2 1",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z", "-x,-x+y,-z", "x-y,-y,-z"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 151,
    "hermann_mauguin": "P 31 1 2",
    "genpos": ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3", "-y,-x,-z+2/3", "x,x-y,-z", "-x+y,y,-z+1/3"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 152,
    "hermann_mauguin": "P 31 2 1",
    "genpos": ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3", "y,x,-z", "-x,-x+y,-z+1/3", "x-y,-y,-z+2/3"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 153,
    "hermann_mauguin": "P 32 1 2",
    "genpos": ["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3", "-y,-x,-z+1/3", "x,x-y,-z", "-x+y,y,-z+2/3"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 154,
    "hermann_mauguin": "P 32 2 1",
    "genpos": ["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3", "y,x,-z", "-x,-x+y,-z+2/3", "x-y,-y,-z+1/3"]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        155,
    "hermann_mauguin":
        "R 3 2",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z", "-x,-x+y,-z", "x-y,-y,-z", "x+2/3,y+1/3,z+1/3",
        "-y+2/3,x-y+1/3,z+1/3", "-x+y+2/3,-x+1/3,z+1/3", "y+2/3,x+1/3,-z+1/3", "-x+2/3,-x+y+1/3,-z+1/3",
        "x-y+2/3,-y+1/3,-z+1/3", "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3",
        "y+1/3,x+2/3,-z+2/3", "-x+1/3,-x+y+2/3,-z+2/3", "x-y+1/3,-y+2/3,-z+2/3"
    ]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 156,
    "hermann_mauguin": "P 3 m 1",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,z", "x,x-y,z", "-x+y,y,z"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 157,
    "hermann_mauguin": "P 3 1 m",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z", "-x,-x+y,z", "x-y,-y,z"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 158,
    "hermann_mauguin": "P 3 c 1",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,z+1/2", "x,x-y,z+1/2", "-x+y,y,z+1/2"]
}, {
    "bravais_lattice": "trigonal",
    "international_number": 159,
    "hermann_mauguin": "P 3 1 c",
    "genpos": ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z+1/2", "-x,-x+y,z+1/2", "x-y,-y,z+1/2"]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        160,
    "hermann_mauguin":
        "R 3 m",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,z", "x,x-y,z", "-x+y,y,z", "x+2/3,y+1/3,z+1/3", "-y+2/3,x-y+1/3,z+1/3",
        "-x+y+2/3,-x+1/3,z+1/3", "-y+2/3,-x+1/3,z+1/3", "x+2/3,x-y+1/3,z+1/3", "-x+y+2/3,y+1/3,z+1/3",
        "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3", "-y+1/3,-x+2/3,z+2/3",
        "x+1/3,x-y+2/3,z+2/3", "-x+y+1/3,y+2/3,z+2/3"
    ]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        161,
    "hermann_mauguin":
        "R 3 c",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,z+1/2", "x,x-y,z+1/2", "-x+y,y,z+1/2", "x+2/3,y+1/3,z+1/3",
        "-y+2/3,x-y+1/3,z+1/3", "-x+y+2/3,-x+1/3,z+1/3", "-y+2/3,-x+1/3,z+5/6", "x+2/3,x-y+1/3,z+5/6",
        "-x+y+2/3,y+1/3,z+5/6", "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3",
        "-y+1/3,-x+2/3,z+7/6", "x+1/3,x-y+2/3,z+7/6", "-x+y+1/3,y+2/3,z+7/6"
    ]
}, {
    "bravais_lattice":
        "trigonal",
    "international_number":
        162,
    "hermann_mauguin":
        "P -3 1 m",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,-z", "x,x-y,-z", "-x+y,y,-z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z",
        "y,x,z", "-x,-x+y,z", "x-y,-y,z"
    ]
}, {
    "bravais_lattice":
        "trigonal",
    "international_number":
        163,
    "hermann_mauguin":
        "P -3 1 c",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "-y,-x,-z+1/2", "x,x-y,-z+1/2", "-x+y,y,-z+1/2", "-x,-y,-z", "y,-x+y,-z",
        "x-y,x,-z", "y,x,z-1/2", "-x,-x+y,z-1/2", "x-y,-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "trigonal",
    "international_number":
        164,
    "hermann_mauguin":
        "P -3 m 1",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z", "-x,-x+y,-z", "x-y,-y,-z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z",
        "-y,-x,z", "x,x-y,z", "-x+y,y,z"
    ]
}, {
    "bravais_lattice":
        "trigonal",
    "international_number":
        165,
    "hermann_mauguin":
        "P -3 c 1",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z+1/2", "-x,-x+y,-z+1/2", "x-y,-y,-z+1/2", "-x,-y,-z", "y,-x+y,-z",
        "x-y,x,-z", "-y,-x,z-1/2", "x,x-y,z-1/2", "-x+y,y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        166,
    "hermann_mauguin":
        "R -3 m",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z", "-x,-x+y,-z", "x-y,-y,-z", "-x,-y,-z", "y,-x+y,-z", "x-y,x,-z",
        "-y,-x,z", "x,x-y,z", "-x+y,y,z", "x+2/3,y+1/3,z+1/3", "-y+2/3,x-y+1/3,z+1/3", "-x+y+2/3,-x+1/3,z+1/3",
        "y+2/3,x+1/3,-z+1/3", "-x+2/3,-x+y+1/3,-z+1/3", "x-y+2/3,-y+1/3,-z+1/3", "-x+2/3,-y+1/3,-z+1/3",
        "y+2/3,-x+y+1/3,-z+1/3", "x-y+2/3,x+1/3,-z+1/3", "-y+2/3,-x+1/3,z+1/3", "x+2/3,x-y+1/3,z+1/3",
        "-x+y+2/3,y+1/3,z+1/3", "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3", "-x+y+1/3,-x+2/3,z+2/3",
        "y+1/3,x+2/3,-z+2/3", "-x+1/3,-x+y+2/3,-z+2/3", "x-y+1/3,-y+2/3,-z+2/3", "-x+1/3,-y+2/3,-z+2/3",
        "y+1/3,-x+y+2/3,-z+2/3", "x-y+1/3,x+2/3,-z+2/3", "-y+1/3,-x+2/3,z+2/3", "x+1/3,x-y+2/3,z+2/3",
        "-x+y+1/3,y+2/3,z+2/3"
    ]
}, {
    "bravais_lattice":
        "rhombohedral",
    "international_number":
        167,
    "hermann_mauguin":
        "R -3 c",
    "genpos": [
        "x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z+1/2", "-x,-x+y,-z+1/2", "x-y,-y,-z+1/2", "-x,-y,-z", "y,-x+y,-z",
        "x-y,x,-z", "-y,-x,z-1/2", "x,x-y,z-1/2", "-x+y,y,z-1/2", "x+2/3,y+1/3,z+1/3", "-y+2/3,x-y+1/3,z+1/3",
        "-x+y+2/3,-x+1/3,z+1/3", "y+2/3,x+1/3,-z+5/6", "-x+2/3,-x+y+1/3,-z+5/6", "x-y+2/3,-y+1/3,-z+5/6",
        "-x+2/3,-y+1/3,-z+1/3", "y+2/3,-x+y+1/3,-z+1/3", "x-y+2/3,x+1/3,-z+1/3", "-y+2/3,-x+1/3,z-1/6",
        "x+2/3,x-y+1/3,z-1/6", "-x+y+2/3,y+1/3,z-1/6", "x+1/3,y+2/3,z+2/3", "-y+1/3,x-y+2/3,z+2/3",
        "-x+y+1/3,-x+2/3,z+2/3", "y+1/3,x+2/3,-z+7/6", "-x+1/3,-x+y+2/3,-z+7/6", "x-y+1/3,-y+2/3,-z+7/6",
        "-x+1/3,-y+2/3,-z+2/3", "y+1/3,-x+y+2/3,-z+2/3", "x-y+1/3,x+2/3,-z+2/3", "-y+1/3,-x+2/3,z+1/6",
        "x+1/3,x-y+2/3,z+1/6", "-x+y+1/3,y+2/3,z+1/6"
    ]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 168,
    "hermann_mauguin": "P 6",
    "genpos": ["x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 169,
    "hermann_mauguin": "P 61",
    "genpos": ["x,y,z", "x-y,x,z+1/6", "-y,x-y,z+1/3", "-x,-y,z+1/2", "-x+y,-x,z+2/3", "y,-x+y,z+5/6"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 170,
    "hermann_mauguin": "P 65",
    "genpos": ["x,y,z", "x-y,x,z+5/6", "-y,x-y,z+2/3", "-x,-y,z+1/2", "-x+y,-x,z+1/3", "y,-x+y,z+1/6"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 171,
    "hermann_mauguin": "P 62",
    "genpos": ["x,y,z", "x-y,x,z+1/3", "-y,x-y,z+2/3", "-x,-y,z", "-x+y,-x,z+1/3", "y,-x+y,z+2/3"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 172,
    "hermann_mauguin": "P 64",
    "genpos": ["x,y,z", "x-y,x,z+2/3", "-y,x-y,z+1/3", "-x,-y,z", "-x+y,-x,z+2/3", "y,-x+y,z+1/3"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 173,
    "hermann_mauguin": "P 63",
    "genpos": ["x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2"]
}, {
    "bravais_lattice": "hexagonal",
    "international_number": 174,
    "hermann_mauguin": "P -6",
    "genpos": ["x,y,z", "-x+y,-x,-z", "-y,x-y,z", "x,y,-z", "-x+y,-x,z", "-y,x-y,-z"]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        175,
    "hermann_mauguin":
        "P 6/m",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "-x,-y,-z", "-x+y,-x,-z", "y,-x+y,-z",
        "x,y,-z", "x-y,x,-z", "-y,x-y,-z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        176,
    "hermann_mauguin":
        "P 63/m",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "-x,-y,-z", "-x+y,-x,-z-1/2",
        "y,-x+y,-z", "x,y,-z-1/2", "x-y,x,-z", "-y,x-y,-z-1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        177,
    "hermann_mauguin":
        "P 6 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "-y,-x,-z", "x-y,-y,-z", "x,x-y,-z",
        "y,x,-z", "-x+y,y,-z", "-x,-x+y,-z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        178,
    "hermann_mauguin":
        "P 61 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z+1/6", "-y,x-y,z+1/3", "-x,-y,z+1/2", "-x+y,-x,z+2/3", "y,-x+y,z+5/6", "-y,-x,-z+5/6",
        "x-y,-y,-z", "x,x-y,-z+1/6", "y,x,-z+1/3", "-x+y,y,-z+1/2", "-x,-x+y,-z+2/3"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        179,
    "hermann_mauguin":
        "P 65 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z+5/6", "-y,x-y,z+2/3", "-x,-y,z+1/2", "-x+y,-x,z+1/3", "y,-x+y,z+1/6", "-y,-x,-z+1/6",
        "x-y,-y,-z", "x,x-y,-z+5/6", "y,x,-z+2/3", "-x+y,y,-z+1/2", "-x,-x+y,-z+1/3"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        180,
    "hermann_mauguin":
        "P 62 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z+1/3", "-y,x-y,z+2/3", "-x,-y,z", "-x+y,-x,z+1/3", "y,-x+y,z+2/3", "-y,-x,-z+2/3", "x-y,-y,-z",
        "x,x-y,-z+1/3", "y,x,-z+2/3", "-x+y,y,-z", "-x,-x+y,-z+1/3"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        181,
    "hermann_mauguin":
        "P 64 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z+2/3", "-y,x-y,z+1/3", "-x,-y,z", "-x+y,-x,z+2/3", "y,-x+y,z+1/3", "-y,-x,-z+1/3", "x-y,-y,-z",
        "x,x-y,-z+2/3", "y,x,-z+1/3", "-x+y,y,-z", "-x,-x+y,-z+2/3"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        182,
    "hermann_mauguin":
        "P 63 2 2",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "-y,-x,-z+1/2", "x-y,-y,-z",
        "x,x-y,-z+1/2", "y,x,-z", "-x+y,y,-z+1/2", "-x,-x+y,-z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        183,
    "hermann_mauguin":
        "P 6 m m",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "y,x,z", "-x+y,y,z", "-x,-x+y,z", "-y,-x,z",
        "x-y,-y,z", "x,x-y,z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        184,
    "hermann_mauguin":
        "P 6 c c",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "y,x,z+1/2", "-x+y,y,z+1/2",
        "-x,-x+y,z+1/2", "-y,-x,z+1/2", "x-y,-y,z+1/2", "x,x-y,z+1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        185,
    "hermann_mauguin":
        "P 63 c m",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "y,x,z", "-x+y,y,z+1/2",
        "-x,-x+y,z", "-y,-x,z+1/2", "x-y,-y,z", "x,x-y,z+1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        186,
    "hermann_mauguin":
        "P 63 m c",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "y,x,z+1/2", "-x+y,y,z",
        "-x,-x+y,z+1/2", "-y,-x,z", "x-y,-y,z+1/2", "x,x-y,z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        187,
    "hermann_mauguin":
        "P -6 m 2",
    "genpos": [
        "x,y,z", "-x+y,-x,-z", "-y,x-y,z", "x,y,-z", "-x+y,-x,z", "-y,x-y,-z", "-y,-x,-z", "-x+y,y,z", "x,x-y,-z",
        "-y,-x,z", "-x+y,y,-z", "x,x-y,z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        188,
    "hermann_mauguin":
        "P -6 c 2",
    "genpos": [
        "x,y,z", "-x+y,-x,-z+1/2", "-y,x-y,z", "x,y,-z+1/2", "-x+y,-x,z", "-y,x-y,-z+1/2", "-y,-x,-z", "-x+y,y,z+1/2",
        "x,x-y,-z", "-y,-x,z+1/2", "-x+y,y,-z", "x,x-y,z+1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        189,
    "hermann_mauguin":
        "P -6 2 m",
    "genpos": [
        "x,y,z", "-x+y,-x,-z", "-y,x-y,z", "x,y,-z", "-x+y,-x,z", "-y,x-y,-z", "y,x,z", "x-y,-y,-z", "-x,-x+y,z",
        "y,x,-z", "x-y,-y,z", "-x,-x+y,-z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        190,
    "hermann_mauguin":
        "P -6 2 c",
    "genpos": [
        "x,y,z", "-x+y,-x,-z+1/2", "-y,x-y,z", "x,y,-z+1/2", "-x+y,-x,z", "-y,x-y,-z+1/2", "y,x,z+1/2", "x-y,-y,-z",
        "-x,-x+y,z+1/2", "y,x,-z", "x-y,-y,z+1/2", "-x,-x+y,-z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        191,
    "hermann_mauguin":
        "P 6/m m m",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "-y,-x,-z", "x-y,-y,-z", "x,x-y,-z",
        "y,x,-z", "-x+y,y,-z", "-x,-x+y,-z", "-x,-y,-z", "-x+y,-x,-z", "y,-x+y,-z", "x,y,-z", "x-y,x,-z", "-y,x-y,-z",
        "y,x,z", "-x+y,y,z", "-x,-x+y,z", "-y,-x,z", "x-y,-y,z", "x,x-y,z"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        192,
    "hermann_mauguin":
        "P 6/m c c",
    "genpos": [
        "x,y,z", "x-y,x,z", "-y,x-y,z", "-x,-y,z", "-x+y,-x,z", "y,-x+y,z", "-y,-x,-z+1/2", "x-y,-y,-z+1/2",
        "x,x-y,-z+1/2", "y,x,-z+1/2", "-x+y,y,-z+1/2", "-x,-x+y,-z+1/2", "-x,-y,-z", "-x+y,-x,-z", "y,-x+y,-z",
        "x,y,-z", "x-y,x,-z", "-y,x-y,-z", "y,x,z-1/2", "-x+y,y,z-1/2", "-x,-x+y,z-1/2", "-y,-x,z-1/2", "x-y,-y,z-1/2",
        "x,x-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        193,
    "hermann_mauguin":
        "P 63/m c m",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "-y,-x,-z", "x-y,-y,-z+1/2",
        "x,x-y,-z", "y,x,-z+1/2", "-x+y,y,-z", "-x,-x+y,-z+1/2", "-x,-y,-z", "-x+y,-x,-z-1/2", "y,-x+y,-z",
        "x,y,-z-1/2", "x-y,x,-z", "-y,x-y,-z-1/2", "y,x,z", "-x+y,y,z-1/2", "-x,-x+y,z", "-y,-x,z-1/2", "x-y,-y,z",
        "x,x-y,z-1/2"
    ]
}, {
    "bravais_lattice":
        "hexagonal",
    "international_number":
        194,
    "hermann_mauguin":
        "P 63/m m c",
    "genpos": [
        "x,y,z", "x-y,x,z+1/2", "-y,x-y,z", "-x,-y,z+1/2", "-x+y,-x,z", "y,-x+y,z+1/2", "-y,-x,-z+1/2", "x-y,-y,-z",
        "x,x-y,-z+1/2", "y,x,-z", "-x+y,y,-z+1/2", "-x,-x+y,-z", "-x,-y,-z", "-x+y,-x,-z-1/2", "y,-x+y,-z",
        "x,y,-z-1/2", "x-y,x,-z", "-y,x-y,-z-1/2", "y,x,z-1/2", "-x+y,y,z", "-x,-x+y,z-1/2", "-y,-x,z", "x-y,-y,z-1/2",
        "x,x-y,z"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        195,
    "hermann_mauguin":
        "P 2 3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        196,
    "hermann_mauguin":
        "F 2 3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2", "-x,y+1/2,-z+1/2", "z,x+1/2,y+1/2",
        "-z,-x+1/2,y+1/2", "z,-x+1/2,-y+1/2", "-z,x+1/2,-y+1/2", "y,z+1/2,x+1/2", "y,-z+1/2,-x+1/2", "-y,z+1/2,-x+1/2",
        "-y,-z+1/2,x+1/2", "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2", "-x+1/2,y,-z+1/2", "z+1/2,x,y+1/2",
        "-z+1/2,-x,y+1/2", "z+1/2,-x,-y+1/2", "-z+1/2,x,-y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z,-x+1/2", "-y+1/2,z,-x+1/2",
        "-y+1/2,-z,x+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z", "z+1/2,x+1/2,y",
        "-z+1/2,-x+1/2,y", "z+1/2,-x+1/2,-y", "-z+1/2,x+1/2,-y", "y+1/2,z+1/2,x", "y+1/2,-z+1/2,-x", "-y+1/2,z+1/2,-x",
        "-y+1/2,-z+1/2,x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        197,
    "hermann_mauguin":
        "I 2 3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "z+1/2,x+1/2,y+1/2", "-z+1/2,-x+1/2,y+1/2", "z+1/2,-x+1/2,-y+1/2", "-z+1/2,x+1/2,-y+1/2", "y+1/2,z+1/2,x+1/2",
        "y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,-x+1/2", "-y+1/2,-z+1/2,x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        198,
    "hermann_mauguin":
        "P 21 3",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x+1/2,-y+1/2,-z", "-x,y+1/2,-z+1/2", "z,x,y", "-z+1/2,-x,y+1/2", "z+1/2,-x+1/2,-y",
        "-z,x+1/2,-y+1/2", "y,z,x", "y+1/2,-z+1/2,-x", "-y,z+1/2,-x+1/2", "-y+1/2,-z,x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        199,
    "hermann_mauguin":
        "I 21 3",
    "genpos": [
        "x,y,z", "-x,-y+1/2,z", "x,-y,-z+1/2", "-x,y+1/2,-z+1/2", "z,x,y", "-z,-x+1/2,y", "z,-x,-y+1/2",
        "-z,x+1/2,-y+1/2", "y,z,x", "y,-z,-x+1/2", "-y,z+1/2,-x+1/2", "-y+1/2,-z,x+1/2", "x+1/2,y+1/2,z+1/2",
        "-x+1/2,-y+1,z+1/2", "x+1/2,-y+1/2,-z+1", "-x+1/2,y+1,-z+1", "z+1/2,x+1/2,y+1/2", "-z+1/2,-x+1,y+1/2",
        "z+1/2,-x+1/2,-y+1", "-z+1/2,x+1,-y+1", "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1", "-y+1/2,z+1,-x+1",
        "-y+1,-z+1/2,x+1"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        200,
    "hermann_mauguin":
        "P m -3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "-z,-x,-y", "z,x,-y", "-z,x,y", "z,-x,y",
        "-y,-z,-x", "-y,z,x", "y,-z,x", "y,z,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        201,
    "hermann_mauguin":
        "P n -3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "-x+1/2,-y+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2",
        "-z+1/2,-x+1/2,-y+1/2", "z+1/2,x+1/2,-y+1/2", "-z+1/2,x+1/2,y+1/2", "z+1/2,-x+1/2,y+1/2",
        "-y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,x+1/2", "y+1/2,z+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        202,
    "hermann_mauguin":
        "F m -3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "-z,-x,-y", "z,x,-y", "-z,x,y", "z,-x,y",
        "-y,-z,-x", "-y,z,x", "y,-z,x", "y,z,-x", "x,y+1/2,z+1/2", "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2",
        "-x,y+1/2,-z+1/2", "z,x+1/2,y+1/2", "-z,-x+1/2,y+1/2", "z,-x+1/2,-y+1/2", "-z,x+1/2,-y+1/2", "y,z+1/2,x+1/2",
        "y,-z+1/2,-x+1/2", "-y,z+1/2,-x+1/2", "-y,-z+1/2,x+1/2", "-x,-y+1/2,-z+1/2", "x,y+1/2,-z+1/2", "-x,y+1/2,z+1/2",
        "x,-y+1/2,z+1/2", "-z,-x+1/2,-y+1/2", "z,x+1/2,-y+1/2", "-z,x+1/2,y+1/2", "z,-x+1/2,y+1/2", "-y,-z+1/2,-x+1/2",
        "-y,z+1/2,x+1/2", "y,-z+1/2,x+1/2", "y,z+1/2,-x+1/2", "x+1/2,y,z+1/2", "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2",
        "-x+1/2,y,-z+1/2", "z+1/2,x,y+1/2", "-z+1/2,-x,y+1/2", "z+1/2,-x,-y+1/2", "-z+1/2,x,-y+1/2", "y+1/2,z,x+1/2",
        "y+1/2,-z,-x+1/2", "-y+1/2,z,-x+1/2", "-y+1/2,-z,x+1/2", "-x+1/2,-y,-z+1/2", "x+1/2,y,-z+1/2", "-x+1/2,y,z+1/2",
        "x+1/2,-y,z+1/2", "-z+1/2,-x,-y+1/2", "z+1/2,x,-y+1/2", "-z+1/2,x,y+1/2", "z+1/2,-x,y+1/2", "-y+1/2,-z,-x+1/2",
        "-y+1/2,z,x+1/2", "y+1/2,-z,x+1/2", "y+1/2,z,-x+1/2", "x+1/2,y+1/2,z", "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z",
        "-x+1/2,y+1/2,-z", "z+1/2,x+1/2,y", "-z+1/2,-x+1/2,y", "z+1/2,-x+1/2,-y", "-z+1/2,x+1/2,-y", "y+1/2,z+1/2,x",
        "y+1/2,-z+1/2,-x", "-y+1/2,z+1/2,-x", "-y+1/2,-z+1/2,x", "-x+1/2,-y+1/2,-z", "x+1/2,y+1/2,-z", "-x+1/2,y+1/2,z",
        "x+1/2,-y+1/2,z", "-z+1/2,-x+1/2,-y", "z+1/2,x+1/2,-y", "-z+1/2,x+1/2,y", "z+1/2,-x+1/2,y", "-y+1/2,-z+1/2,-x",
        "-y+1/2,z+1/2,x", "y+1/2,-z+1/2,x", "y+1/2,z+1/2,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        203,
    "hermann_mauguin":
        "F d -3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "-x+1/4,-y+1/4,-z+1/4", "x+1/4,y+1/4,-z+1/4", "-x+1/4,y+1/4,z+1/4", "x+1/4,-y+1/4,z+1/4",
        "-z+1/4,-x+1/4,-y+1/4", "z+1/4,x+1/4,-y+1/4", "-z+1/4,x+1/4,y+1/4", "z+1/4,-x+1/4,y+1/4",
        "-y+1/4,-z+1/4,-x+1/4", "-y+1/4,z+1/4,x+1/4", "y+1/4,-z+1/4,x+1/4", "y+1/4,z+1/4,-x+1/4", "x,y+1/2,z+1/2",
        "-x,-y+1/2,z+1/2", "x,-y+1/2,-z+1/2", "-x,y+1/2,-z+1/2", "z,x+1/2,y+1/2", "-z,-x+1/2,y+1/2", "z,-x+1/2,-y+1/2",
        "-z,x+1/2,-y+1/2", "y,z+1/2,x+1/2", "y,-z+1/2,-x+1/2", "-y,z+1/2,-x+1/2", "-y,-z+1/2,x+1/2",
        "-x+1/4,-y+3/4,-z+3/4", "x+1/4,y+3/4,-z+3/4", "-x+1/4,y+3/4,z+3/4", "x+1/4,-y+3/4,z+3/4",
        "-z+1/4,-x+3/4,-y+3/4", "z+1/4,x+3/4,-y+3/4", "-z+1/4,x+3/4,y+3/4", "z+1/4,-x+3/4,y+3/4",
        "-y+1/4,-z+3/4,-x+3/4", "-y+1/4,z+3/4,x+3/4", "y+1/4,-z+3/4,x+3/4", "y+1/4,z+3/4,-x+3/4", "x+1/2,y,z+1/2",
        "-x+1/2,-y,z+1/2", "x+1/2,-y,-z+1/2", "-x+1/2,y,-z+1/2", "z+1/2,x,y+1/2", "-z+1/2,-x,y+1/2", "z+1/2,-x,-y+1/2",
        "-z+1/2,x,-y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z,-x+1/2", "-y+1/2,z,-x+1/2", "-y+1/2,-z,x+1/2",
        "-x+3/4,-y+1/4,-z+3/4", "x+3/4,y+1/4,-z+3/4", "-x+3/4,y+1/4,z+3/4", "x+3/4,-y+1/4,z+3/4",
        "-z+3/4,-x+1/4,-y+3/4", "z+3/4,x+1/4,-y+3/4", "-z+3/4,x+1/4,y+3/4", "z+3/4,-x+1/4,y+3/4",
        "-y+3/4,-z+1/4,-x+3/4", "-y+3/4,z+1/4,x+3/4", "y+3/4,-z+1/4,x+3/4", "y+3/4,z+1/4,-x+3/4", "x+1/2,y+1/2,z",
        "-x+1/2,-y+1/2,z", "x+1/2,-y+1/2,-z", "-x+1/2,y+1/2,-z", "z+1/2,x+1/2,y", "-z+1/2,-x+1/2,y", "z+1/2,-x+1/2,-y",
        "-z+1/2,x+1/2,-y", "y+1/2,z+1/2,x", "y+1/2,-z+1/2,-x", "-y+1/2,z+1/2,-x", "-y+1/2,-z+1/2,x",
        "-x+3/4,-y+3/4,-z+1/4", "x+3/4,y+3/4,-z+1/4", "-x+3/4,y+3/4,z+1/4", "x+3/4,-y+3/4,z+1/4",
        "-z+3/4,-x+3/4,-y+1/4", "z+3/4,x+3/4,-y+1/4", "-z+3/4,x+3/4,y+1/4", "z+3/4,-x+3/4,y+1/4",
        "-y+3/4,-z+3/4,-x+1/4", "-y+3/4,z+3/4,x+1/4", "y+3/4,-z+3/4,x+1/4", "y+3/4,z+3/4,-x+1/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        204,
    "hermann_mauguin":
        "I m -3",
    "genpos": [
        "x,y,z", "-x,-y,z", "x,-y,-z", "-x,y,-z", "z,x,y", "-z,-x,y", "z,-x,-y", "-z,x,-y", "y,z,x", "y,-z,-x",
        "-y,z,-x", "-y,-z,x", "-x,-y,-z", "x,y,-z", "-x,y,z", "x,-y,z", "-z,-x,-y", "z,x,-y", "-z,x,y", "z,-x,y",
        "-y,-z,-x", "-y,z,x", "y,-z,x", "y,z,-x", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2",
        "-x+1/2,y+1/2,-z+1/2", "z+1/2,x+1/2,y+1/2", "-z+1/2,-x+1/2,y+1/2", "z+1/2,-x+1/2,-y+1/2", "-z+1/2,x+1/2,-y+1/2",
        "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,-x+1/2", "-y+1/2,-z+1/2,x+1/2",
        "-x+1/2,-y+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2",
        "-z+1/2,-x+1/2,-y+1/2", "z+1/2,x+1/2,-y+1/2", "-z+1/2,x+1/2,y+1/2", "z+1/2,-x+1/2,y+1/2",
        "-y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,x+1/2", "y+1/2,z+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        205,
    "hermann_mauguin":
        "P a -3",
    "genpos": [
        "x,y,z", "-x+1/2,-y,z+1/2", "x+1/2,-y+1/2,-z", "-x,y+1/2,-z+1/2", "z,x,y", "-z+1/2,-x,y+1/2", "z+1/2,-x+1/2,-y",
        "-z,x+1/2,-y+1/2", "y,z,x", "y+1/2,-z+1/2,-x", "-y,z+1/2,-x+1/2", "-y+1/2,-z,x+1/2", "-x,-y,-z",
        "x-1/2,y,-z-1/2", "-x-1/2,y-1/2,z", "x,-y-1/2,z-1/2", "-z,-x,-y", "z-1/2,x,-y-1/2", "-z-1/2,x-1/2,y",
        "z,-x-1/2,y-1/2", "-y,-z,-x", "-y-1/2,z-1/2,x", "y,-z-1/2,x-1/2", "y-1/2,z,-x-1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        206,
    "hermann_mauguin":
        "I a -3",
    "genpos": [
        "x,y,z", "-x,-y+1/2,z", "x,-y,-z+1/2", "-x,y+1/2,-z+1/2", "z,x,y", "-z,-x+1/2,y", "z,-x,-y+1/2",
        "-z,x+1/2,-y+1/2", "y,z,x", "y,-z,-x+1/2", "-y,z+1/2,-x+1/2", "-y+1/2,-z,x+1/2", "-x,-y,-z", "x,y-1/2,-z",
        "-x,y,z-1/2", "x,-y-1/2,z-1/2", "-z,-x,-y", "z,x-1/2,-y", "-z,x,y-1/2", "z,-x-1/2,y-1/2", "-y,-z,-x",
        "-y,z,x-1/2", "y,-z-1/2,x-1/2", "y-1/2,z,-x-1/2", "x+1/2,y+1/2,z+1/2", "-x+1/2,-y+1,z+1/2", "x+1/2,-y+1/2,-z+1",
        "-x+1/2,y+1,-z+1", "z+1/2,x+1/2,y+1/2", "-z+1/2,-x+1,y+1/2", "z+1/2,-x+1/2,-y+1", "-z+1/2,x+1,-y+1",
        "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1", "-y+1/2,z+1,-x+1", "-y+1,-z+1/2,x+1", "-x+1/2,-y+1/2,-z+1/2",
        "x+1/2,y,-z+1/2", "-x+1/2,y+1/2,z", "x+1/2,-y,z", "-z+1/2,-x+1/2,-y+1/2", "z+1/2,x,-y+1/2", "-z+1/2,x+1/2,y",
        "z+1/2,-x,y", "-y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,x", "y+1/2,-z,x", "y,z+1/2,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        207,
    "hermann_mauguin":
        "P 4 3 2",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        208,
    "hermann_mauguin":
        "P 42 3 2",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x,-y,-z", "y+1/2,x+1/2,-z+1/2", "-x,y,-z",
        "-y+1/2,-x+1/2,-z+1/2", "z,x,y", "-x+1/2,z+1/2,y+1/2", "-z,-x,y", "x+1/2,-z+1/2,y+1/2", "z,-x,-y",
        "x+1/2,z+1/2,-y+1/2", "-z,x,-y", "-x+1/2,-z+1/2,-y+1/2", "y,z,x", "y,-z,-x", "z+1/2,y+1/2,-x+1/2", "-y,z,-x",
        "-z+1/2,-y+1/2,-x+1/2", "-y,-z,x", "z+1/2,-y+1/2,x+1/2", "-z+1/2,y+1/2,x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        209,
    "hermann_mauguin":
        "F 4 3 2",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "x,y+1/2,z+1/2", "-y,x+1/2,z+1/2", "-x,-y+1/2,z+1/2",
        "y,-x+1/2,z+1/2", "x,-y+1/2,-z+1/2", "y,x+1/2,-z+1/2", "-x,y+1/2,-z+1/2", "-y,-x+1/2,-z+1/2", "z,x+1/2,y+1/2",
        "-x,z+1/2,y+1/2", "-z,-x+1/2,y+1/2", "x,-z+1/2,y+1/2", "z,-x+1/2,-y+1/2", "x,z+1/2,-y+1/2", "-z,x+1/2,-y+1/2",
        "-x,-z+1/2,-y+1/2", "y,z+1/2,x+1/2", "y,-z+1/2,-x+1/2", "z,y+1/2,-x+1/2", "-y,z+1/2,-x+1/2", "-z,-y+1/2,-x+1/2",
        "-y,-z+1/2,x+1/2", "z,-y+1/2,x+1/2", "-z,y+1/2,x+1/2", "x+1/2,y,z+1/2", "-y+1/2,x,z+1/2", "-x+1/2,-y,z+1/2",
        "y+1/2,-x,z+1/2", "x+1/2,-y,-z+1/2", "y+1/2,x,-z+1/2", "-x+1/2,y,-z+1/2", "-y+1/2,-x,-z+1/2", "z+1/2,x,y+1/2",
        "-x+1/2,z,y+1/2", "-z+1/2,-x,y+1/2", "x+1/2,-z,y+1/2", "z+1/2,-x,-y+1/2", "x+1/2,z,-y+1/2", "-z+1/2,x,-y+1/2",
        "-x+1/2,-z,-y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z,-x+1/2", "z+1/2,y,-x+1/2", "-y+1/2,z,-x+1/2", "-z+1/2,-y,-x+1/2",
        "-y+1/2,-z,x+1/2", "z+1/2,-y,x+1/2", "-z+1/2,y,x+1/2", "x+1/2,y+1/2,z", "-y+1/2,x+1/2,z", "-x+1/2,-y+1/2,z",
        "y+1/2,-x+1/2,z", "x+1/2,-y+1/2,-z", "y+1/2,x+1/2,-z", "-x+1/2,y+1/2,-z", "-y+1/2,-x+1/2,-z", "z+1/2,x+1/2,y",
        "-x+1/2,z+1/2,y", "-z+1/2,-x+1/2,y", "x+1/2,-z+1/2,y", "z+1/2,-x+1/2,-y", "x+1/2,z+1/2,-y", "-z+1/2,x+1/2,-y",
        "-x+1/2,-z+1/2,-y", "y+1/2,z+1/2,x", "y+1/2,-z+1/2,-x", "z+1/2,y+1/2,-x", "-y+1/2,z+1/2,-x", "-z+1/2,-y+1/2,-x",
        "-y+1/2,-z+1/2,x", "z+1/2,-y+1/2,x", "-z+1/2,y+1/2,x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        210,
    "hermann_mauguin":
        "F 41 3 2",
    "genpos": [
        "x,y,z", "-y+1/4,x+1/4,z+1/4", "-x,-y+1/2,z+1/2", "y+3/4,-x+1/4,z+3/4", "x,-y,-z", "y+1/4,x+1/4,-z+1/4",
        "-x,y+1/2,-z+1/2", "-y+3/4,-x+1/4,-z+3/4", "z,x,y", "-x+1/4,z+1/4,y+1/4", "-z,-x+1/2,y+1/2",
        "x+3/4,-z+1/4,y+3/4", "z,-x,-y", "x+1/4,z+1/4,-y+1/4", "-z,x+1/2,-y+1/2", "-x+3/4,-z+1/4,-y+3/4", "y,z,x",
        "y+1/2,-z,-x+1/2", "z+1/4,y+3/4,-x+3/4", "-y+1/2,z+1/2,-x", "-z+1/4,-y+1/4,-x+1/4", "-y,-z,x",
        "z+1/4,-y+3/4,x+3/4", "-z+3/4,y+3/4,x+1/4", "x,y+1/2,z+1/2", "-y+1/4,x+3/4,z+3/4", "-x,-y+1,z+1",
        "y+3/4,-x+3/4,z+5/4", "x,-y+1/2,-z+1/2", "y+1/4,x+3/4,-z+3/4", "-x,y+1,-z+1", "-y+3/4,-x+3/4,-z+5/4",
        "z,x+1/2,y+1/2", "-x+1/4,z+3/4,y+3/4", "-z,-x+1,y+1", "x+3/4,-z+3/4,y+5/4", "z,-x+1/2,-y+1/2",
        "x+1/4,z+3/4,-y+3/4", "-z,x+1,-y+1", "-x+3/4,-z+3/4,-y+5/4", "y,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1",
        "z+1/4,y+5/4,-x+5/4", "-y+1/2,z+1,-x+1/2", "-z+1/4,-y+3/4,-x+3/4", "-y,-z+1/2,x+1/2", "z+1/4,-y+5/4,x+5/4",
        "-z+3/4,y+5/4,x+3/4", "x+1/2,y,z+1/2", "-y+3/4,x+1/4,z+3/4", "-x+1/2,-y+1/2,z+1", "y+5/4,-x+1/4,z+5/4",
        "x+1/2,-y,-z+1/2", "y+3/4,x+1/4,-z+3/4", "-x+1/2,y+1/2,-z+1", "-y+5/4,-x+1/4,-z+5/4", "z+1/2,x,y+1/2",
        "-x+3/4,z+1/4,y+3/4", "-z+1/2,-x+1/2,y+1", "x+5/4,-z+1/4,y+5/4", "z+1/2,-x,-y+1/2", "x+3/4,z+1/4,-y+3/4",
        "-z+1/2,x+1/2,-y+1", "-x+5/4,-z+1/4,-y+5/4", "y+1/2,z,x+1/2", "y+1,-z,-x+1", "z+3/4,y+3/4,-x+5/4",
        "-y+1,z+1/2,-x+1/2", "-z+3/4,-y+1/4,-x+3/4", "-y+1/2,-z,x+1/2", "z+3/4,-y+3/4,x+5/4", "-z+5/4,y+3/4,x+3/4",
        "x+1/2,y+1/2,z", "-y+3/4,x+3/4,z+1/4", "-x+1/2,-y+1,z+1/2", "y+5/4,-x+3/4,z+3/4", "x+1/2,-y+1/2,-z",
        "y+3/4,x+3/4,-z+1/4", "-x+1/2,y+1,-z+1/2", "-y+5/4,-x+3/4,-z+3/4", "z+1/2,x+1/2,y", "-x+3/4,z+3/4,y+1/4",
        "-z+1/2,-x+1,y+1/2", "x+5/4,-z+3/4,y+3/4", "z+1/2,-x+1/2,-y", "x+3/4,z+3/4,-y+1/4", "-z+1/2,x+1,-y+1/2",
        "-x+5/4,-z+3/4,-y+3/4", "y+1/2,z+1/2,x", "y+1,-z+1/2,-x+1/2", "z+3/4,y+5/4,-x+3/4", "-y+1,z+1,-x",
        "-z+3/4,-y+3/4,-x+1/4", "-y+1/2,-z+1/2,x", "z+3/4,-y+5/4,x+3/4", "-z+5/4,y+5/4,x+1/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        211,
    "hermann_mauguin":
        "I 4 3 2",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "x+1/2,y+1/2,z+1/2", "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2",
        "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,-z+1/2",
        "-y+1/2,-x+1/2,-z+1/2", "z+1/2,x+1/2,y+1/2", "-x+1/2,z+1/2,y+1/2", "-z+1/2,-x+1/2,y+1/2", "x+1/2,-z+1/2,y+1/2",
        "z+1/2,-x+1/2,-y+1/2", "x+1/2,z+1/2,-y+1/2", "-z+1/2,x+1/2,-y+1/2", "-x+1/2,-z+1/2,-y+1/2", "y+1/2,z+1/2,x+1/2",
        "y+1/2,-z+1/2,-x+1/2", "z+1/2,y+1/2,-x+1/2", "-y+1/2,z+1/2,-x+1/2", "-z+1/2,-y+1/2,-x+1/2",
        "-y+1/2,-z+1/2,x+1/2", "z+1/2,-y+1/2,x+1/2", "-z+1/2,y+1/2,x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        212,
    "hermann_mauguin":
        "P 43 3 2",
    "genpos": [
        "x,y,z", "-y+3/4,x+1/4,z+3/4", "-x+1/2,-y,z+1/2", "y+3/4,-x+3/4,z+1/4", "x+1/2,-y+1/2,-z", "y+1/4,x+3/4,-z+3/4",
        "-x,y+1/2,-z+1/2", "-y+1/4,-x+1/4,-z+1/4", "z,x,y", "-x+3/4,z+1/4,y+3/4", "-z+1/2,-x,y+1/2",
        "x+3/4,-z+3/4,y+1/4", "z+1/2,-x+1/2,-y", "x+1/4,z+3/4,-y+3/4", "-z,x+1/2,-y+1/2", "-x+1/4,-z+1/4,-y+1/4",
        "y,z,x", "y+1/2,-z+1/2,-x", "z+1/4,y+3/4,-x+3/4", "-y,z+1/2,-x+1/2", "-z+1/4,-y+1/4,-x+1/4", "-y+1/2,-z,x+1/2",
        "z+3/4,-y+3/4,x+1/4", "-z+3/4,y+1/4,x+3/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        213,
    "hermann_mauguin":
        "P 41 3 2",
    "genpos": [
        "x,y,z", "-y+1/4,x+3/4,z+1/4", "-x+1/2,-y,z+1/2", "y+1/4,-x+1/4,z+3/4", "x+1/2,-y+1/2,-z", "y+3/4,x+1/4,-z+1/4",
        "-x,y+1/2,-z+1/2", "-y+3/4,-x+3/4,-z+3/4", "z,x,y", "-x+1/4,z+3/4,y+1/4", "-z+1/2,-x,y+1/2",
        "x+1/4,-z+1/4,y+3/4", "z+1/2,-x+1/2,-y", "x+3/4,z+1/4,-y+1/4", "-z,x+1/2,-y+1/2", "-x+3/4,-z+3/4,-y+3/4",
        "y,z,x", "y+1/2,-z+1/2,-x", "z+3/4,y+1/4,-x+1/4", "-y,z+1/2,-x+1/2", "-z+3/4,-y+3/4,-x+3/4", "-y+1/2,-z,x+1/2",
        "z+1/4,-y+1/4,x+3/4", "-z+1/4,y+3/4,x+1/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        214,
    "hermann_mauguin":
        "I 41 3 2",
    "genpos": [
        "x,y,z", "-y+1/4,x+3/4,z+1/4", "-x+1/2,-y,z+1/2", "y+1/4,-x+1/4,z+3/4", "x,-y,-z+1/2", "y+1/4,x+3/4,-z+3/4",
        "-x+1/2,y,-z", "-y+1/4,-x+1/4,-z+1/4", "z,x,y", "-x+1/4,z+3/4,y+1/4", "-z+1/2,-x,y+1/2", "x+1/4,-z+1/4,y+3/4",
        "z,-x,-y+1/2", "x+1/4,z+3/4,-y+3/4", "-z+1/2,x,-y", "-x+1/4,-z+1/4,-y+1/4", "y,z,x", "y+1/2,-z+1/2,-x",
        "z+3/4,y+1/4,-x+1/4", "-y,z+1/2,-x+1/2", "-z+1/4,-y+1/4,-x+1/4", "-y+1/2,-z,x+1/2", "z+3/4,-y+3/4,x+1/4",
        "-z+3/4,y+1/4,x+3/4", "x+1/2,y+1/2,z+1/2", "-y+3/4,x+5/4,z+3/4", "-x+1,-y+1/2,z+1", "y+3/4,-x+3/4,z+5/4",
        "x+1/2,-y+1/2,-z+1", "y+3/4,x+5/4,-z+5/4", "-x+1,y+1/2,-z+1/2", "-y+3/4,-x+3/4,-z+3/4", "z+1/2,x+1/2,y+1/2",
        "-x+3/4,z+5/4,y+3/4", "-z+1,-x+1/2,y+1", "x+3/4,-z+3/4,y+5/4", "z+1/2,-x+1/2,-y+1", "x+3/4,z+5/4,-y+5/4",
        "-z+1,x+1/2,-y+1/2", "-x+3/4,-z+3/4,-y+3/4", "y+1/2,z+1/2,x+1/2", "y+1,-z+1,-x+1/2", "z+5/4,y+3/4,-x+3/4",
        "-y+1/2,z+1,-x+1", "-z+3/4,-y+3/4,-x+3/4", "-y+1,-z+1/2,x+1", "z+5/4,-y+5/4,x+3/4", "-z+5/4,y+3/4,x+5/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        215,
    "hermann_mauguin":
        "P -4 3 m",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z", "-y,-x,z", "-x,y,-z", "y,x,z", "z,x,y", "x,-z,-y",
        "-z,-x,y", "-x,z,-y", "z,-x,-y", "-x,-z,y", "-z,x,-y", "x,z,y", "y,z,x", "y,-z,-x", "-z,-y,x", "-y,z,-x",
        "z,y,x", "-y,-z,x", "-z,y,-x", "z,-y,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        216,
    "hermann_mauguin":
        "F -4 3 m",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z", "-y,-x,z", "-x,y,-z", "y,x,z", "z,x,y", "x,-z,-y",
        "-z,-x,y", "-x,z,-y", "z,-x,-y", "-x,-z,y", "-z,x,-y", "x,z,y", "y,z,x", "y,-z,-x", "-z,-y,x", "-y,z,-x",
        "z,y,x", "-y,-z,x", "-z,y,-x", "z,-y,-x", "x,y+1/2,z+1/2", "y,-x+1/2,-z+1/2", "-x,-y+1/2,z+1/2",
        "-y,x+1/2,-z+1/2", "x,-y+1/2,-z+1/2", "-y,-x+1/2,z+1/2", "-x,y+1/2,-z+1/2", "y,x+1/2,z+1/2", "z,x+1/2,y+1/2",
        "x,-z+1/2,-y+1/2", "-z,-x+1/2,y+1/2", "-x,z+1/2,-y+1/2", "z,-x+1/2,-y+1/2", "-x,-z+1/2,y+1/2",
        "-z,x+1/2,-y+1/2", "x,z+1/2,y+1/2", "y,z+1/2,x+1/2", "y,-z+1/2,-x+1/2", "-z,-y+1/2,x+1/2", "-y,z+1/2,-x+1/2",
        "z,y+1/2,x+1/2", "-y,-z+1/2,x+1/2", "-z,y+1/2,-x+1/2", "z,-y+1/2,-x+1/2", "x+1/2,y,z+1/2", "y+1/2,-x,-z+1/2",
        "-x+1/2,-y,z+1/2", "-y+1/2,x,-z+1/2", "x+1/2,-y,-z+1/2", "-y+1/2,-x,z+1/2", "-x+1/2,y,-z+1/2", "y+1/2,x,z+1/2",
        "z+1/2,x,y+1/2", "x+1/2,-z,-y+1/2", "-z+1/2,-x,y+1/2", "-x+1/2,z,-y+1/2", "z+1/2,-x,-y+1/2", "-x+1/2,-z,y+1/2",
        "-z+1/2,x,-y+1/2", "x+1/2,z,y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z,-x+1/2", "-z+1/2,-y,x+1/2", "-y+1/2,z,-x+1/2",
        "z+1/2,y,x+1/2", "-y+1/2,-z,x+1/2", "-z+1/2,y,-x+1/2", "z+1/2,-y,-x+1/2", "x+1/2,y+1/2,z", "y+1/2,-x+1/2,-z",
        "-x+1/2,-y+1/2,z", "-y+1/2,x+1/2,-z", "x+1/2,-y+1/2,-z", "-y+1/2,-x+1/2,z", "-x+1/2,y+1/2,-z", "y+1/2,x+1/2,z",
        "z+1/2,x+1/2,y", "x+1/2,-z+1/2,-y", "-z+1/2,-x+1/2,y", "-x+1/2,z+1/2,-y", "z+1/2,-x+1/2,-y", "-x+1/2,-z+1/2,y",
        "-z+1/2,x+1/2,-y", "x+1/2,z+1/2,y", "y+1/2,z+1/2,x", "y+1/2,-z+1/2,-x", "-z+1/2,-y+1/2,x", "-y+1/2,z+1/2,-x",
        "z+1/2,y+1/2,x", "-y+1/2,-z+1/2,x", "-z+1/2,y+1/2,-x", "z+1/2,-y+1/2,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        217,
    "hermann_mauguin":
        "I -4 3 m",
    "genpos": [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z", "-y,-x,z", "-x,y,-z", "y,x,z", "z,x,y", "x,-z,-y",
        "-z,-x,y", "-x,z,-y", "z,-x,-y", "-x,-z,y", "-z,x,-y", "x,z,y", "y,z,x", "y,-z,-x", "-z,-y,x", "-y,z,-x",
        "z,y,x", "-y,-z,x", "-z,y,-x", "z,-y,-x", "x+1/2,y+1/2,z+1/2", "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1/2,z+1/2",
        "-y+1/2,x+1/2,-z+1/2", "x+1/2,-y+1/2,-z+1/2", "-y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1/2,-z+1/2", "y+1/2,x+1/2,z+1/2",
        "z+1/2,x+1/2,y+1/2", "x+1/2,-z+1/2,-y+1/2", "-z+1/2,-x+1/2,y+1/2", "-x+1/2,z+1/2,-y+1/2", "z+1/2,-x+1/2,-y+1/2",
        "-x+1/2,-z+1/2,y+1/2", "-z+1/2,x+1/2,-y+1/2", "x+1/2,z+1/2,y+1/2", "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1/2",
        "-z+1/2,-y+1/2,x+1/2", "-y+1/2,z+1/2,-x+1/2", "z+1/2,y+1/2,x+1/2", "-y+1/2,-z+1/2,x+1/2", "-z+1/2,y+1/2,-x+1/2",
        "z+1/2,-y+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        218,
    "hermann_mauguin":
        "P -4 3 n",
    "genpos": [
        "x,y,z", "y+1/2,-x+1/2,-z+1/2", "-x,-y,z", "-y+1/2,x+1/2,-z+1/2", "x,-y,-z", "-y+1/2,-x+1/2,z+1/2", "-x,y,-z",
        "y+1/2,x+1/2,z+1/2", "z,x,y", "x+1/2,-z+1/2,-y+1/2", "-z,-x,y", "-x+1/2,z+1/2,-y+1/2", "z,-x,-y",
        "-x+1/2,-z+1/2,y+1/2", "-z,x,-y", "x+1/2,z+1/2,y+1/2", "y,z,x", "y,-z,-x", "-z+1/2,-y+1/2,x+1/2", "-y,z,-x",
        "z+1/2,y+1/2,x+1/2", "-y,-z,x", "-z+1/2,y+1/2,-x+1/2", "z+1/2,-y+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        219,
    "hermann_mauguin":
        "F -4 3 c",
    "genpos": [
        "x,y,z", "y+1/2,-x,-z", "-x+1/2,-y+1/2,z", "-y,x+1/2,-z", "x,-y,-z", "-y+1/2,-x,z", "-x+1/2,y+1/2,-z",
        "y,x+1/2,z", "z,x,y", "x+1/2,-z,-y", "-z+1/2,-x+1/2,y", "-x,z+1/2,-y", "z,-x,-y", "-x+1/2,-z,y",
        "-z+1/2,x+1/2,-y", "x,z+1/2,y", "y,z,x", "y,-z+1/2,-x+1/2", "-z,-y,x+1/2", "-y+1/2,z,-x+1/2", "z+1/2,y,x",
        "-y,-z,x", "-z,y,-x+1/2", "z+1/2,-y+1/2,-x+1/2", "x,y+1/2,z+1/2", "y+1/2,-x+1/2,-z+1/2", "-x+1/2,-y+1,z+1/2",
        "-y,x+1,-z+1/2", "x,-y+1/2,-z+1/2", "-y+1/2,-x+1/2,z+1/2", "-x+1/2,y+1,-z+1/2", "y,x+1,z+1/2", "z,x+1/2,y+1/2",
        "x+1/2,-z+1/2,-y+1/2", "-z+1/2,-x+1,y+1/2", "-x,z+1,-y+1/2", "z,-x+1/2,-y+1/2", "-x+1/2,-z+1/2,y+1/2",
        "-z+1/2,x+1,-y+1/2", "x,z+1,y+1/2", "y,z+1/2,x+1/2", "y,-z+1,-x+1", "-z,-y+1/2,x+1", "-y+1/2,z+1/2,-x+1",
        "z+1/2,y+1/2,x+1/2", "-y,-z+1/2,x+1/2", "-z,y+1/2,-x+1", "z+1/2,-y+1,-x+1", "x+1/2,y,z+1/2", "y+1,-x,-z+1/2",
        "-x+1,-y+1/2,z+1/2", "-y+1/2,x+1/2,-z+1/2", "x+1/2,-y,-z+1/2", "-y+1,-x,z+1/2", "-x+1,y+1/2,-z+1/2",
        "y+1/2,x+1/2,z+1/2", "z+1/2,x,y+1/2", "x+1,-z,-y+1/2", "-z+1,-x+1/2,y+1/2", "-x+1/2,z+1/2,-y+1/2",
        "z+1/2,-x,-y+1/2", "-x+1,-z,y+1/2", "-z+1,x+1/2,-y+1/2", "x+1/2,z+1/2,y+1/2", "y+1/2,z,x+1/2",
        "y+1/2,-z+1/2,-x+1", "-z+1/2,-y,x+1", "-y+1,z,-x+1", "z+1,y,x+1/2", "-y+1/2,-z,x+1/2", "-z+1/2,y,-x+1",
        "z+1,-y+1/2,-x+1", "x+1/2,y+1/2,z", "y+1,-x+1/2,-z", "-x+1,-y+1,z", "-y+1/2,x+1,-z", "x+1/2,-y+1/2,-z",
        "-y+1,-x+1/2,z", "-x+1,y+1,-z", "y+1/2,x+1,z", "z+1/2,x+1/2,y", "x+1,-z+1/2,-y", "-z+1,-x+1,y", "-x+1/2,z+1,-y",
        "z+1/2,-x+1/2,-y", "-x+1,-z+1/2,y", "-z+1,x+1,-y", "x+1/2,z+1,y", "y+1/2,z+1/2,x", "y+1/2,-z+1,-x+1/2",
        "-z+1/2,-y+1/2,x+1/2", "-y+1,z+1/2,-x+1/2", "z+1,y+1/2,x", "-y+1/2,-z+1/2,x", "-z+1/2,y+1/2,-x+1/2",
        "z+1,-y+1,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        220,
    "hermann_mauguin":
        "I -4 3 d",
    "genpos": [
        "x,y,z", "y+1/4,-x+3/4,-z+1/4", "-x,-y+1/2,z", "-y+3/4,x+3/4,-z+1/4", "x,-y,-z+1/2", "-y+1/4,-x+3/4,z+3/4",
        "-x,y+1/2,-z+1/2", "y+3/4,x+3/4,z+3/4", "z,x,y", "x+1/4,-z+3/4,-y+1/4", "-z,-x+1/2,y", "-x+3/4,z+3/4,-y+1/4",
        "z,-x,-y+1/2", "-x+1/4,-z+3/4,y+3/4", "-z,x+1/2,-y+1/2", "x+3/4,z+3/4,y+3/4", "y,z,x", "y,-z,-x+1/2",
        "-z+1/4,-y+3/4,x+3/4", "-y,z+1/2,-x+1/2", "z+1/4,y+1/4,x+1/4", "-y+1/2,-z,x+1/2", "-z+1/4,y+1/4,-x+3/4",
        "z+3/4,-y+1/4,-x+3/4", "x+1/2,y+1/2,z+1/2", "y+3/4,-x+5/4,-z+3/4", "-x+1/2,-y+1,z+1/2", "-y+5/4,x+5/4,-z+3/4",
        "x+1/2,-y+1/2,-z+1", "-y+3/4,-x+5/4,z+5/4", "-x+1/2,y+1,-z+1", "y+5/4,x+5/4,z+5/4", "z+1/2,x+1/2,y+1/2",
        "x+3/4,-z+5/4,-y+3/4", "-z+1/2,-x+1,y+1/2", "-x+5/4,z+5/4,-y+3/4", "z+1/2,-x+1/2,-y+1", "-x+3/4,-z+5/4,y+5/4",
        "-z+1/2,x+1,-y+1", "x+5/4,z+5/4,y+5/4", "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1", "-z+3/4,-y+5/4,x+5/4",
        "-y+1/2,z+1,-x+1", "z+3/4,y+3/4,x+3/4", "-y+1,-z+1/2,x+1", "-z+3/4,y+3/4,-x+5/4", "z+5/4,-y+3/4,-x+5/4"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        221,
    "hermann_mauguin":
        "P m -3 m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "-x,y,z", "-y,-x,z",
        "x,-y,z", "y,x,z", "-z,-x,-y", "x,-z,-y", "z,x,-y", "-x,z,-y", "-z,x,y", "-x,-z,y", "z,-x,y", "x,z,y",
        "-y,-z,-x", "-y,z,x", "-z,-y,x", "y,-z,x", "z,y,x", "y,z,-x", "-z,y,-x", "z,-y,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        222,
    "hermann_mauguin":
        "P n -3 n",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "-x+1/2,-y+1/2,-z+1/2", "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2",
        "-y+1/2,x+1/2,-z+1/2", "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2", "y+1/2,x+1/2,z+1/2",
        "-z+1/2,-x+1/2,-y+1/2", "x+1/2,-z+1/2,-y+1/2", "z+1/2,x+1/2,-y+1/2", "-x+1/2,z+1/2,-y+1/2",
        "-z+1/2,x+1/2,y+1/2", "-x+1/2,-z+1/2,y+1/2", "z+1/2,-x+1/2,y+1/2", "x+1/2,z+1/2,y+1/2", "-y+1/2,-z+1/2,-x+1/2",
        "-y+1/2,z+1/2,x+1/2", "-z+1/2,-y+1/2,x+1/2", "y+1/2,-z+1/2,x+1/2", "z+1/2,y+1/2,x+1/2", "y+1/2,z+1/2,-x+1/2",
        "-z+1/2,y+1/2,-x+1/2", "z+1/2,-y+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        223,
    "hermann_mauguin":
        "P m -3 n",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x,-y,-z", "y+1/2,x+1/2,-z+1/2", "-x,y,-z",
        "-y+1/2,-x+1/2,-z+1/2", "z,x,y", "-x+1/2,z+1/2,y+1/2", "-z,-x,y", "x+1/2,-z+1/2,y+1/2", "z,-x,-y",
        "x+1/2,z+1/2,-y+1/2", "-z,x,-y", "-x+1/2,-z+1/2,-y+1/2", "y,z,x", "y,-z,-x", "z+1/2,y+1/2,-x+1/2", "-y,z,-x",
        "-z+1/2,-y+1/2,-x+1/2", "-y,-z,x", "z+1/2,-y+1/2,x+1/2", "-z+1/2,y+1/2,x+1/2", "-x,-y,-z",
        "y-1/2,-x-1/2,-z-1/2", "x,y,-z", "-y-1/2,x-1/2,-z-1/2", "-x,y,z", "-y-1/2,-x-1/2,z-1/2", "x,-y,z",
        "y-1/2,x-1/2,z-1/2", "-z,-x,-y", "x-1/2,-z-1/2,-y-1/2", "z,x,-y", "-x-1/2,z-1/2,-y-1/2", "-z,x,y",
        "-x-1/2,-z-1/2,y-1/2", "z,-x,y", "x-1/2,z-1/2,y-1/2", "-y,-z,-x", "-y,z,x", "-z-1/2,-y-1/2,x-1/2", "y,-z,x",
        "z-1/2,y-1/2,x-1/2", "y,z,-x", "-z-1/2,y-1/2,-x-1/2", "z-1/2,-y-1/2,-x-1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        224,
    "hermann_mauguin":
        "P n -3 m",
    "genpos": [
        "x,y,z", "-y+1/2,x+1/2,z+1/2", "-x,-y,z", "y+1/2,-x+1/2,z+1/2", "x,-y,-z", "y+1/2,x+1/2,-z+1/2", "-x,y,-z",
        "-y+1/2,-x+1/2,-z+1/2", "z,x,y", "-x+1/2,z+1/2,y+1/2", "-z,-x,y", "x+1/2,-z+1/2,y+1/2", "z,-x,-y",
        "x+1/2,z+1/2,-y+1/2", "-z,x,-y", "-x+1/2,-z+1/2,-y+1/2", "y,z,x", "y,-z,-x", "z+1/2,y+1/2,-x+1/2", "-y,z,-x",
        "-z+1/2,-y+1/2,-x+1/2", "-y,-z,x", "z+1/2,-y+1/2,x+1/2", "-z+1/2,y+1/2,x+1/2", "-x+1/2,-y+1/2,-z+1/2",
        "y,-x,-z", "x+1/2,y+1/2,-z+1/2", "-y,x,-z", "-x+1/2,y+1/2,z+1/2", "-y,-x,z", "x+1/2,-y+1/2,z+1/2", "y,x,z",
        "-z+1/2,-x+1/2,-y+1/2", "x,-z,-y", "z+1/2,x+1/2,-y+1/2", "-x,z,-y", "-z+1/2,x+1/2,y+1/2", "-x,-z,y",
        "z+1/2,-x+1/2,y+1/2", "x,z,y", "-y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,x+1/2", "-z,-y,x", "y+1/2,-z+1/2,x+1/2",
        "z,y,x", "y+1/2,z+1/2,-x+1/2", "-z,y,-x", "z,-y,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        225,
    "hermann_mauguin":
        "F m -3 m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "-x,y,z", "-y,-x,z",
        "x,-y,z", "y,x,z", "-z,-x,-y", "x,-z,-y", "z,x,-y", "-x,z,-y", "-z,x,y", "-x,-z,y", "z,-x,y", "x,z,y",
        "-y,-z,-x", "-y,z,x", "-z,-y,x", "y,-z,x", "z,y,x", "y,z,-x", "-z,y,-x", "z,-y,-x", "x,y+1/2,z+1/2",
        "-y,x+1/2,z+1/2", "-x,-y+1/2,z+1/2", "y,-x+1/2,z+1/2", "x,-y+1/2,-z+1/2", "y,x+1/2,-z+1/2", "-x,y+1/2,-z+1/2",
        "-y,-x+1/2,-z+1/2", "z,x+1/2,y+1/2", "-x,z+1/2,y+1/2", "-z,-x+1/2,y+1/2", "x,-z+1/2,y+1/2", "z,-x+1/2,-y+1/2",
        "x,z+1/2,-y+1/2", "-z,x+1/2,-y+1/2", "-x,-z+1/2,-y+1/2", "y,z+1/2,x+1/2", "y,-z+1/2,-x+1/2", "z,y+1/2,-x+1/2",
        "-y,z+1/2,-x+1/2", "-z,-y+1/2,-x+1/2", "-y,-z+1/2,x+1/2", "z,-y+1/2,x+1/2", "-z,y+1/2,x+1/2",
        "-x,-y+1/2,-z+1/2", "y,-x+1/2,-z+1/2", "x,y+1/2,-z+1/2", "-y,x+1/2,-z+1/2", "-x,y+1/2,z+1/2", "-y,-x+1/2,z+1/2",
        "x,-y+1/2,z+1/2", "y,x+1/2,z+1/2", "-z,-x+1/2,-y+1/2", "x,-z+1/2,-y+1/2", "z,x+1/2,-y+1/2", "-x,z+1/2,-y+1/2",
        "-z,x+1/2,y+1/2", "-x,-z+1/2,y+1/2", "z,-x+1/2,y+1/2", "x,z+1/2,y+1/2", "-y,-z+1/2,-x+1/2", "-y,z+1/2,x+1/2",
        "-z,-y+1/2,x+1/2", "y,-z+1/2,x+1/2", "z,y+1/2,x+1/2", "y,z+1/2,-x+1/2", "-z,y+1/2,-x+1/2", "z,-y+1/2,-x+1/2",
        "x+1/2,y,z+1/2", "-y+1/2,x,z+1/2", "-x+1/2,-y,z+1/2", "y+1/2,-x,z+1/2", "x+1/2,-y,-z+1/2", "y+1/2,x,-z+1/2",
        "-x+1/2,y,-z+1/2", "-y+1/2,-x,-z+1/2", "z+1/2,x,y+1/2", "-x+1/2,z,y+1/2", "-z+1/2,-x,y+1/2", "x+1/2,-z,y+1/2",
        "z+1/2,-x,-y+1/2", "x+1/2,z,-y+1/2", "-z+1/2,x,-y+1/2", "-x+1/2,-z,-y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z,-x+1/2",
        "z+1/2,y,-x+1/2", "-y+1/2,z,-x+1/2", "-z+1/2,-y,-x+1/2", "-y+1/2,-z,x+1/2", "z+1/2,-y,x+1/2", "-z+1/2,y,x+1/2",
        "-x+1/2,-y,-z+1/2", "y+1/2,-x,-z+1/2", "x+1/2,y,-z+1/2", "-y+1/2,x,-z+1/2", "-x+1/2,y,z+1/2", "-y+1/2,-x,z+1/2",
        "x+1/2,-y,z+1/2", "y+1/2,x,z+1/2", "-z+1/2,-x,-y+1/2", "x+1/2,-z,-y+1/2", "z+1/2,x,-y+1/2", "-x+1/2,z,-y+1/2",
        "-z+1/2,x,y+1/2", "-x+1/2,-z,y+1/2", "z+1/2,-x,y+1/2", "x+1/2,z,y+1/2", "-y+1/2,-z,-x+1/2", "-y+1/2,z,x+1/2",
        "-z+1/2,-y,x+1/2", "y+1/2,-z,x+1/2", "z+1/2,y,x+1/2", "y+1/2,z,-x+1/2", "-z+1/2,y,-x+1/2", "z+1/2,-y,-x+1/2",
        "x+1/2,y+1/2,z", "-y+1/2,x+1/2,z", "-x+1/2,-y+1/2,z", "y+1/2,-x+1/2,z", "x+1/2,-y+1/2,-z", "y+1/2,x+1/2,-z",
        "-x+1/2,y+1/2,-z", "-y+1/2,-x+1/2,-z", "z+1/2,x+1/2,y", "-x+1/2,z+1/2,y", "-z+1/2,-x+1/2,y", "x+1/2,-z+1/2,y",
        "z+1/2,-x+1/2,-y", "x+1/2,z+1/2,-y", "-z+1/2,x+1/2,-y", "-x+1/2,-z+1/2,-y", "y+1/2,z+1/2,x", "y+1/2,-z+1/2,-x",
        "z+1/2,y+1/2,-x", "-y+1/2,z+1/2,-x", "-z+1/2,-y+1/2,-x", "-y+1/2,-z+1/2,x", "z+1/2,-y+1/2,x", "-z+1/2,y+1/2,x",
        "-x+1/2,-y+1/2,-z", "y+1/2,-x+1/2,-z", "x+1/2,y+1/2,-z", "-y+1/2,x+1/2,-z", "-x+1/2,y+1/2,z", "-y+1/2,-x+1/2,z",
        "x+1/2,-y+1/2,z", "y+1/2,x+1/2,z", "-z+1/2,-x+1/2,-y", "x+1/2,-z+1/2,-y", "z+1/2,x+1/2,-y", "-x+1/2,z+1/2,-y",
        "-z+1/2,x+1/2,y", "-x+1/2,-z+1/2,y", "z+1/2,-x+1/2,y", "x+1/2,z+1/2,y", "-y+1/2,-z+1/2,-x", "-y+1/2,z+1/2,x",
        "-z+1/2,-y+1/2,x", "y+1/2,-z+1/2,x", "z+1/2,y+1/2,x", "y+1/2,z+1/2,-x", "-z+1/2,y+1/2,-x", "z+1/2,-y+1/2,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        226,
    "hermann_mauguin":
        "F m -3 c",
    "genpos": [
        "x,y,z", "-y+1/2,x,z", "-x+1/2,-y+1/2,z", "y,-x+1/2,z", "x,-y,-z", "y+1/2,x,-z", "-x+1/2,y+1/2,-z",
        "-y,-x+1/2,-z", "z,x,y", "-x+1/2,z,y", "-z+1/2,-x+1/2,y", "x,-z+1/2,y", "z,-x,-y", "x+1/2,z,-y",
        "-z+1/2,x+1/2,-y", "-x,-z+1/2,-y", "y,z,x", "y,-z+1/2,-x+1/2", "z,y,-x+1/2", "-y+1/2,z,-x+1/2", "-z+1/2,-y,-x",
        "-y,-z,x", "z,-y,x+1/2", "-z+1/2,y+1/2,x+1/2", "-x,-y,-z", "y-1/2,-x,-z", "x-1/2,y-1/2,-z", "-y,x-1/2,-z",
        "-x,y,z", "-y-1/2,-x,z", "x-1/2,-y-1/2,z", "y,x-1/2,z", "-z,-x,-y", "x-1/2,-z,-y", "z-1/2,x-1/2,-y",
        "-x,z-1/2,-y", "-z,x,y", "-x-1/2,-z,y", "z-1/2,-x-1/2,y", "x,z-1/2,y", "-y,-z,-x", "-y,z-1/2,x-1/2",
        "-z,-y,x-1/2", "y-1/2,-z,x-1/2", "z-1/2,y,x", "y,z,-x", "-z,y,-x-1/2", "z-1/2,-y-1/2,-x-1/2", "x,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1,z+1/2", "y,-x+1,z+1/2", "x,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y+1,-z+1/2", "-y,-x+1,-z+1/2", "z,x+1/2,y+1/2", "-x+1/2,z+1/2,y+1/2", "-z+1/2,-x+1,y+1/2",
        "x,-z+1,y+1/2", "z,-x+1/2,-y+1/2", "x+1/2,z+1/2,-y+1/2", "-z+1/2,x+1,-y+1/2", "-x,-z+1,-y+1/2", "y,z+1/2,x+1/2",
        "y,-z+1,-x+1", "z,y+1/2,-x+1", "-y+1/2,z+1/2,-x+1", "-z+1/2,-y+1/2,-x+1/2", "-y,-z+1/2,x+1/2", "z,-y+1/2,x+1",
        "-z+1/2,y+1,x+1", "-x,-y+1/2,-z+1/2", "y-1/2,-x+1/2,-z+1/2", "x-1/2,y,-z+1/2", "-y,x,-z+1/2", "-x,y+1/2,z+1/2",
        "-y-1/2,-x+1/2,z+1/2", "x-1/2,-y,z+1/2", "y,x,z+1/2", "-z,-x+1/2,-y+1/2", "x-1/2,-z+1/2,-y+1/2",
        "z-1/2,x,-y+1/2", "-x,z,-y+1/2", "-z,x+1/2,y+1/2", "-x-1/2,-z+1/2,y+1/2", "z-1/2,-x,y+1/2", "x,z,y+1/2",
        "-y,-z+1/2,-x+1/2", "-y,z,x", "-z,-y+1/2,x", "y-1/2,-z+1/2,x", "z-1/2,y+1/2,x+1/2", "y,z+1/2,-x+1/2",
        "-z,y+1/2,-x", "z-1/2,-y,-x", "x+1/2,y,z+1/2", "-y+1,x,z+1/2", "-x+1,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2",
        "x+1/2,-y,-z+1/2", "y+1,x,-z+1/2", "-x+1,y+1/2,-z+1/2", "-y+1/2,-x+1/2,-z+1/2", "z+1/2,x,y+1/2", "-x+1,z,y+1/2",
        "-z+1,-x+1/2,y+1/2", "x+1/2,-z+1/2,y+1/2", "z+1/2,-x,-y+1/2", "x+1,z,-y+1/2", "-z+1,x+1/2,-y+1/2",
        "-x+1/2,-z+1/2,-y+1/2", "y+1/2,z,x+1/2", "y+1/2,-z+1/2,-x+1", "z+1/2,y,-x+1", "-y+1,z,-x+1", "-z+1,-y,-x+1/2",
        "-y+1/2,-z,x+1/2", "z+1/2,-y,x+1", "-z+1,y+1/2,x+1", "-x+1/2,-y,-z+1/2", "y,-x,-z+1/2", "x,y-1/2,-z+1/2",
        "-y+1/2,x-1/2,-z+1/2", "-x+1/2,y,z+1/2", "-y,-x,z+1/2", "x,-y-1/2,z+1/2", "y+1/2,x-1/2,z+1/2",
        "-z+1/2,-x,-y+1/2", "x,-z,-y+1/2", "z,x-1/2,-y+1/2", "-x+1/2,z-1/2,-y+1/2", "-z+1/2,x,y+1/2", "-x,-z,y+1/2",
        "z,-x-1/2,y+1/2", "x+1/2,z-1/2,y+1/2", "-y+1/2,-z,-x+1/2", "-y+1/2,z-1/2,x", "-z+1/2,-y,x", "y,-z,x",
        "z,y,x+1/2", "y+1/2,z,-x+1/2", "-z+1/2,y,-x", "z,-y-1/2,-x", "x+1/2,y+1/2,z", "-y+1,x+1/2,z", "-x+1,-y+1,z",
        "y+1/2,-x+1,z", "x+1/2,-y+1/2,-z", "y+1,x+1/2,-z", "-x+1,y+1,-z", "-y+1/2,-x+1,-z", "z+1/2,x+1/2,y",
        "-x+1,z+1/2,y", "-z+1,-x+1,y", "x+1/2,-z+1,y", "z+1/2,-x+1/2,-y", "x+1,z+1/2,-y", "-z+1,x+1,-y",
        "-x+1/2,-z+1,-y", "y+1/2,z+1/2,x", "y+1/2,-z+1,-x+1/2", "z+1/2,y+1/2,-x+1/2", "-y+1,z+1/2,-x+1/2",
        "-z+1,-y+1/2,-x", "-y+1/2,-z+1/2,x", "z+1/2,-y+1/2,x+1/2", "-z+1,y+1,x+1/2", "-x+1/2,-y+1/2,-z", "y,-x+1/2,-z",
        "x,y,-z", "-y+1/2,x,-z", "-x+1/2,y+1/2,z", "-y,-x+1/2,z", "x,-y,z", "y+1/2,x,z", "-z+1/2,-x+1/2,-y",
        "x,-z+1/2,-y", "z,x,-y", "-x+1/2,z,-y", "-z+1/2,x+1/2,y", "-x,-z+1/2,y", "z,-x,y", "x+1/2,z,y",
        "-y+1/2,-z+1/2,-x", "-y+1/2,z,x-1/2", "-z+1/2,-y+1/2,x-1/2", "y,-z+1/2,x-1/2", "z,y+1/2,x", "y+1/2,z+1/2,-x",
        "-z+1/2,y+1/2,-x-1/2", "z,-y,-x-1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        227,
    "hermann_mauguin":
        "F d -3 m",
    "genpos": [
        "x,y,z", "-y+1/4,x+1/4,z+1/4", "-x,-y+1/2,z+1/2", "y+3/4,-x+1/4,z+3/4", "x,-y,-z", "y+1/4,x+1/4,-z+1/4",
        "-x,y+1/2,-z+1/2", "-y+3/4,-x+1/4,-z+3/4", "z,x,y", "-x+1/4,z+1/4,y+1/4", "-z,-x+1/2,y+1/2",
        "x+3/4,-z+1/4,y+3/4", "z,-x,-y", "x+1/4,z+1/4,-y+1/4", "-z,x+1/2,-y+1/2", "-x+3/4,-z+1/4,-y+3/4", "y,z,x",
        "y+1/2,-z,-x+1/2", "z+1/4,y+3/4,-x+3/4", "-y+1/2,z+1/2,-x", "-z+1/4,-y+1/4,-x+1/4", "-y,-z,x",
        "z+1/4,-y+3/4,x+3/4", "-z+3/4,y+3/4,x+1/4", "-x+1/4,-y+1/4,-z+1/4", "y,-x,-z", "x+1/4,y-1/4,-z-1/4",
        "-y-1/2,x,-z-1/2", "-x+1/4,y+1/4,z+1/4", "-y,-x,z", "x+1/4,-y-1/4,z-1/4", "y-1/2,x,z-1/2",
        "-z+1/4,-x+1/4,-y+1/4", "x,-z,-y", "z+1/4,x-1/4,-y-1/4", "-x-1/2,z,-y-1/2", "-z+1/4,x+1/4,y+1/4", "-x,-z,y",
        "z+1/4,-x-1/4,y-1/4", "x-1/2,z,y-1/2", "-y+1/4,-z+1/4,-x+1/4", "-y-1/4,z+1/4,x-1/4", "-z,-y-1/2,x-1/2",
        "y-1/4,-z-1/4,x+1/4", "z,y,x", "y+1/4,z+1/4,-x+1/4", "-z,y-1/2,-x-1/2", "z-1/2,-y-1/2,-x", "x,y+1/2,z+1/2",
        "-y+1/4,x+3/4,z+3/4", "-x,-y+1,z+1", "y+3/4,-x+3/4,z+5/4", "x,-y+1/2,-z+1/2", "y+1/4,x+3/4,-z+3/4",
        "-x,y+1,-z+1", "-y+3/4,-x+3/4,-z+5/4", "z,x+1/2,y+1/2", "-x+1/4,z+3/4,y+3/4", "-z,-x+1,y+1",
        "x+3/4,-z+3/4,y+5/4", "z,-x+1/2,-y+1/2", "x+1/4,z+3/4,-y+3/4", "-z,x+1,-y+1", "-x+3/4,-z+3/4,-y+5/4",
        "y,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1", "z+1/4,y+5/4,-x+5/4", "-y+1/2,z+1,-x+1/2", "-z+1/4,-y+3/4,-x+3/4",
        "-y,-z+1/2,x+1/2", "z+1/4,-y+5/4,x+5/4", "-z+3/4,y+5/4,x+3/4", "-x+1/4,-y+3/4,-z+3/4", "y,-x+1/2,-z+1/2",
        "x+1/4,y+1/4,-z+1/4", "-y-1/2,x+1/2,-z", "-x+1/4,y+3/4,z+3/4", "-y,-x+1/2,z+1/2", "x+1/4,-y+1/4,z+1/4",
        "y-1/2,x+1/2,z", "-z+1/4,-x+3/4,-y+3/4", "x,-z+1/2,-y+1/2", "z+1/4,x+1/4,-y+1/4", "-x-1/2,z+1/2,-y",
        "-z+1/4,x+3/4,y+3/4", "-x,-z+1/2,y+1/2", "z+1/4,-x+1/4,y+1/4", "x-1/2,z+1/2,y", "-y+1/4,-z+3/4,-x+3/4",
        "-y-1/4,z+3/4,x+1/4", "-z,-y,x", "y-1/4,-z+1/4,x+3/4", "z,y+1/2,x+1/2", "y+1/4,z+3/4,-x+3/4", "-z,y,-x",
        "z-1/2,-y,-x+1/2", "x+1/2,y,z+1/2", "-y+3/4,x+1/4,z+3/4", "-x+1/2,-y+1/2,z+1", "y+5/4,-x+1/4,z+5/4",
        "x+1/2,-y,-z+1/2", "y+3/4,x+1/4,-z+3/4", "-x+1/2,y+1/2,-z+1", "-y+5/4,-x+1/4,-z+5/4", "z+1/2,x,y+1/2",
        "-x+3/4,z+1/4,y+3/4", "-z+1/2,-x+1/2,y+1", "x+5/4,-z+1/4,y+5/4", "z+1/2,-x,-y+1/2", "x+3/4,z+1/4,-y+3/4",
        "-z+1/2,x+1/2,-y+1", "-x+5/4,-z+1/4,-y+5/4", "y+1/2,z,x+1/2", "y+1,-z,-x+1", "z+3/4,y+3/4,-x+5/4",
        "-y+1,z+1/2,-x+1/2", "-z+3/4,-y+1/4,-x+3/4", "-y+1/2,-z,x+1/2", "z+3/4,-y+3/4,x+5/4", "-z+5/4,y+3/4,x+3/4",
        "-x+3/4,-y+1/4,-z+3/4", "y+1/2,-x,-z+1/2", "x+3/4,y-1/4,-z+1/4", "-y,x,-z", "-x+3/4,y+1/4,z+3/4",
        "-y+1/2,-x,z+1/2", "x+3/4,-y-1/4,z+1/4", "y,x,z", "-z+3/4,-x+1/4,-y+3/4", "x+1/2,-z,-y+1/2",
        "z+3/4,x-1/4,-y+1/4", "-x,z,-y", "-z+3/4,x+1/4,y+3/4", "-x+1/2,-z,y+1/2", "z+3/4,-x-1/4,y+1/4", "x,z,y",
        "-y+3/4,-z+1/4,-x+3/4", "-y+1/4,z+1/4,x+1/4", "-z+1/2,-y-1/2,x", "y+1/4,-z-1/4,x+3/4", "z+1/2,y,x+1/2",
        "y+3/4,z+1/4,-x+3/4", "-z+1/2,y-1/2,-x", "z,-y-1/2,-x+1/2", "x+1/2,y+1/2,z", "-y+3/4,x+3/4,z+1/4",
        "-x+1/2,-y+1,z+1/2", "y+5/4,-x+3/4,z+3/4", "x+1/2,-y+1/2,-z", "y+3/4,x+3/4,-z+1/4", "-x+1/2,y+1,-z+1/2",
        "-y+5/4,-x+3/4,-z+3/4", "z+1/2,x+1/2,y", "-x+3/4,z+3/4,y+1/4", "-z+1/2,-x+1,y+1/2", "x+5/4,-z+3/4,y+3/4",
        "z+1/2,-x+1/2,-y", "x+3/4,z+3/4,-y+1/4", "-z+1/2,x+1,-y+1/2", "-x+5/4,-z+3/4,-y+3/4", "y+1/2,z+1/2,x",
        "y+1,-z+1/2,-x+1/2", "z+3/4,y+5/4,-x+3/4", "-y+1,z+1,-x", "-z+3/4,-y+3/4,-x+1/4", "-y+1/2,-z+1/2,x",
        "z+3/4,-y+5/4,x+3/4", "-z+5/4,y+5/4,x+1/4", "-x+3/4,-y+3/4,-z+1/4", "y+1/2,-x+1/2,-z", "x+3/4,y+1/4,-z-1/4",
        "-y,x+1/2,-z-1/2", "-x+3/4,y+3/4,z+1/4", "-y+1/2,-x+1/2,z", "x+3/4,-y+1/4,z-1/4", "y,x+1/2,z-1/2",
        "-z+3/4,-x+3/4,-y+1/4", "x+1/2,-z+1/2,-y", "z+3/4,x+1/4,-y-1/4", "-x,z+1/2,-y-1/2", "-z+3/4,x+3/4,y+1/4",
        "-x+1/2,-z+1/2,y", "z+3/4,-x+1/4,y-1/4", "x,z+1/2,y-1/2", "-y+3/4,-z+3/4,-x+1/4", "-y+1/4,z+3/4,x-1/4",
        "-z+1/2,-y,x-1/2", "y+1/4,-z+1/4,x+1/4", "z+1/2,y+1/2,x", "y+3/4,z+3/4,-x+1/4", "-z+1/2,y,-x-1/2", "z,-y,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        228,
    "hermann_mauguin":
        "F d -3 c",
    "genpos": [
        "x,y,z", "-y+1/4,x+1/4,z+1/4", "-x,-y+1/2,z+1/2", "y+3/4,-x+1/4,z+3/4", "x,-y,-z", "y+1/4,x+1/4,-z+1/4",
        "-x,y+1/2,-z+1/2", "-y+3/4,-x+1/4,-z+3/4", "z,x,y", "-x+1/4,z+1/4,y+1/4", "-z,-x+1/2,y+1/2",
        "x+3/4,-z+1/4,y+3/4", "z,-x,-y", "x+1/4,z+1/4,-y+1/4", "-z,x+1/2,-y+1/2", "-x+3/4,-z+1/4,-y+3/4", "y,z,x",
        "y+1/2,-z,-x+1/2", "z+1/4,y+3/4,-x+3/4", "-y+1/2,z+1/2,-x", "-z+1/4,-y+1/4,-x+1/4", "-y,-z,x",
        "z+1/4,-y+3/4,x+3/4", "-z+3/4,y+3/4,x+1/4", "-x+3/4,-y+1/4,-z+1/4", "y+1/2,-x,-z", "x+3/4,y-1/4,-z-1/4",
        "-y,x,-z-1/2", "-x+3/4,y+1/4,z+1/4", "-y+1/2,-x,z", "x+3/4,-y-1/4,z-1/4", "y,x,z-1/2", "-z+3/4,-x+1/4,-y+1/4",
        "x+1/2,-z,-y", "z+3/4,x-1/4,-y-1/4", "-x,z,-y-1/2", "-z+3/4,x+1/4,y+1/4", "-x+1/2,-z,y", "z+3/4,-x-1/4,y-1/4",
        "x,z,y-1/2", "-y+3/4,-z+1/4,-x+1/4", "-y+1/4,z+1/4,x-1/4", "-z+1/2,-y-1/2,x-1/2", "y+1/4,-z-1/4,x+1/4",
        "z+1/2,y,x", "y+3/4,z+1/4,-x+1/4", "-z+1/2,y-1/2,-x-1/2", "z,-y-1/2,-x", "x,y+1/2,z+1/2", "-y+1/4,x+3/4,z+3/4",
        "-x,-y+1,z+1", "y+3/4,-x+3/4,z+5/4", "x,-y+1/2,-z+1/2", "y+1/4,x+3/4,-z+3/4", "-x,y+1,-z+1",
        "-y+3/4,-x+3/4,-z+5/4", "z,x+1/2,y+1/2", "-x+1/4,z+3/4,y+3/4", "-z,-x+1,y+1", "x+3/4,-z+3/4,y+5/4",
        "z,-x+1/2,-y+1/2", "x+1/4,z+3/4,-y+3/4", "-z,x+1,-y+1", "-x+3/4,-z+3/4,-y+5/4", "y,z+1/2,x+1/2",
        "y+1/2,-z+1/2,-x+1", "z+1/4,y+5/4,-x+5/4", "-y+1/2,z+1,-x+1/2", "-z+1/4,-y+3/4,-x+3/4", "-y,-z+1/2,x+1/2",
        "z+1/4,-y+5/4,x+5/4", "-z+3/4,y+5/4,x+3/4", "-x+3/4,-y+3/4,-z+3/4", "y+1/2,-x+1/2,-z+1/2", "x+3/4,y+1/4,-z+1/4",
        "-y,x+1/2,-z", "-x+3/4,y+3/4,z+3/4", "-y+1/2,-x+1/2,z+1/2", "x+3/4,-y+1/4,z+1/4", "y,x+1/2,z",
        "-z+3/4,-x+3/4,-y+3/4", "x+1/2,-z+1/2,-y+1/2", "z+3/4,x+1/4,-y+1/4", "-x,z+1/2,-y", "-z+3/4,x+3/4,y+3/4",
        "-x+1/2,-z+1/2,y+1/2", "z+3/4,-x+1/4,y+1/4", "x,z+1/2,y", "-y+3/4,-z+3/4,-x+3/4", "-y+1/4,z+3/4,x+1/4",
        "-z+1/2,-y,x", "y+1/4,-z+1/4,x+3/4", "z+1/2,y+1/2,x+1/2", "y+3/4,z+3/4,-x+3/4", "-z+1/2,y,-x", "z,-y,-x+1/2",
        "x+1/2,y,z+1/2", "-y+3/4,x+1/4,z+3/4", "-x+1/2,-y+1/2,z+1", "y+5/4,-x+1/4,z+5/4", "x+1/2,-y,-z+1/2",
        "y+3/4,x+1/4,-z+3/4", "-x+1/2,y+1/2,-z+1", "-y+5/4,-x+1/4,-z+5/4", "z+1/2,x,y+1/2", "-x+3/4,z+1/4,y+3/4",
        "-z+1/2,-x+1/2,y+1", "x+5/4,-z+1/4,y+5/4", "z+1/2,-x,-y+1/2", "x+3/4,z+1/4,-y+3/4", "-z+1/2,x+1/2,-y+1",
        "-x+5/4,-z+1/4,-y+5/4", "y+1/2,z,x+1/2", "y+1,-z,-x+1", "z+3/4,y+3/4,-x+5/4", "-y+1,z+1/2,-x+1/2",
        "-z+3/4,-y+1/4,-x+3/4", "-y+1/2,-z,x+1/2", "z+3/4,-y+3/4,x+5/4", "-z+5/4,y+3/4,x+3/4", "-x+5/4,-y+1/4,-z+3/4",
        "y+1,-x,-z+1/2", "x+5/4,y-1/4,-z+1/4", "-y+1/2,x,-z", "-x+5/4,y+1/4,z+3/4", "-y+1,-x,z+1/2",
        "x+5/4,-y-1/4,z+1/4", "y+1/2,x,z", "-z+5/4,-x+1/4,-y+3/4", "x+1,-z,-y+1/2", "z+5/4,x-1/4,-y+1/4", "-x+1/2,z,-y",
        "-z+5/4,x+1/4,y+3/4", "-x+1,-z,y+1/2", "z+5/4,-x-1/4,y+1/4", "x+1/2,z,y", "-y+5/4,-z+1/4,-x+3/4",
        "-y+3/4,z+1/4,x+1/4", "-z+1,-y-1/2,x", "y+3/4,-z-1/4,x+3/4", "z+1,y,x+1/2", "y+5/4,z+1/4,-x+3/4",
        "-z+1,y-1/2,-x", "z+1/2,-y-1/2,-x+1/2", "x+1/2,y+1/2,z", "-y+3/4,x+3/4,z+1/4", "-x+1/2,-y+1,z+1/2",
        "y+5/4,-x+3/4,z+3/4", "x+1/2,-y+1/2,-z", "y+3/4,x+3/4,-z+1/4", "-x+1/2,y+1,-z+1/2", "-y+5/4,-x+3/4,-z+3/4",
        "z+1/2,x+1/2,y", "-x+3/4,z+3/4,y+1/4", "-z+1/2,-x+1,y+1/2", "x+5/4,-z+3/4,y+3/4", "z+1/2,-x+1/2,-y",
        "x+3/4,z+3/4,-y+1/4", "-z+1/2,x+1,-y+1/2", "-x+5/4,-z+3/4,-y+3/4", "y+1/2,z+1/2,x", "y+1,-z+1/2,-x+1/2",
        "z+3/4,y+5/4,-x+3/4", "-y+1,z+1,-x", "-z+3/4,-y+3/4,-x+1/4", "-y+1/2,-z+1/2,x", "z+3/4,-y+5/4,x+3/4",
        "-z+5/4,y+5/4,x+1/4", "-x+5/4,-y+3/4,-z+1/4", "y+1,-x+1/2,-z", "x+5/4,y+1/4,-z-1/4", "-y+1/2,x+1/2,-z-1/2",
        "-x+5/4,y+3/4,z+1/4", "-y+1,-x+1/2,z", "x+5/4,-y+1/4,z-1/4", "y+1/2,x+1/2,z-1/2", "-z+5/4,-x+3/4,-y+1/4",
        "x+1,-z+1/2,-y", "z+5/4,x+1/4,-y-1/4", "-x+1/2,z+1/2,-y-1/2", "-z+5/4,x+3/4,y+1/4", "-x+1,-z+1/2,y",
        "z+5/4,-x+1/4,y-1/4", "x+1/2,z+1/2,y-1/2", "-y+5/4,-z+3/4,-x+1/4", "-y+3/4,z+3/4,x-1/4", "-z+1,-y,x-1/2",
        "y+3/4,-z+1/4,x+1/4", "z+1,y+1/2,x", "y+5/4,z+3/4,-x+1/4", "-z+1,y,-x-1/2", "z+1/2,-y,-x"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        229,
    "hermann_mauguin":
        "I m -3 m",
    "genpos": [
        "x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z", "x,-y,-z", "y,x,-z", "-x,y,-z", "-y,-x,-z", "z,x,y", "-x,z,y",
        "-z,-x,y", "x,-z,y", "z,-x,-y", "x,z,-y", "-z,x,-y", "-x,-z,-y", "y,z,x", "y,-z,-x", "z,y,-x", "-y,z,-x",
        "-z,-y,-x", "-y,-z,x", "z,-y,x", "-z,y,x", "-x,-y,-z", "y,-x,-z", "x,y,-z", "-y,x,-z", "-x,y,z", "-y,-x,z",
        "x,-y,z", "y,x,z", "-z,-x,-y", "x,-z,-y", "z,x,-y", "-x,z,-y", "-z,x,y", "-x,-z,y", "z,-x,y", "x,z,y",
        "-y,-z,-x", "-y,z,x", "-z,-y,x", "y,-z,x", "z,y,x", "y,z,-x", "-z,y,-x", "z,-y,-x", "x+1/2,y+1/2,z+1/2",
        "-y+1/2,x+1/2,z+1/2", "-x+1/2,-y+1/2,z+1/2", "y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,-z+1/2", "y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y+1/2,-z+1/2", "-y+1/2,-x+1/2,-z+1/2", "z+1/2,x+1/2,y+1/2", "-x+1/2,z+1/2,y+1/2", "-z+1/2,-x+1/2,y+1/2",
        "x+1/2,-z+1/2,y+1/2", "z+1/2,-x+1/2,-y+1/2", "x+1/2,z+1/2,-y+1/2", "-z+1/2,x+1/2,-y+1/2",
        "-x+1/2,-z+1/2,-y+1/2", "y+1/2,z+1/2,x+1/2", "y+1/2,-z+1/2,-x+1/2", "z+1/2,y+1/2,-x+1/2", "-y+1/2,z+1/2,-x+1/2",
        "-z+1/2,-y+1/2,-x+1/2", "-y+1/2,-z+1/2,x+1/2", "z+1/2,-y+1/2,x+1/2", "-z+1/2,y+1/2,x+1/2",
        "-x+1/2,-y+1/2,-z+1/2", "y+1/2,-x+1/2,-z+1/2", "x+1/2,y+1/2,-z+1/2", "-y+1/2,x+1/2,-z+1/2",
        "-x+1/2,y+1/2,z+1/2", "-y+1/2,-x+1/2,z+1/2", "x+1/2,-y+1/2,z+1/2", "y+1/2,x+1/2,z+1/2", "-z+1/2,-x+1/2,-y+1/2",
        "x+1/2,-z+1/2,-y+1/2", "z+1/2,x+1/2,-y+1/2", "-x+1/2,z+1/2,-y+1/2", "-z+1/2,x+1/2,y+1/2", "-x+1/2,-z+1/2,y+1/2",
        "z+1/2,-x+1/2,y+1/2", "x+1/2,z+1/2,y+1/2", "-y+1/2,-z+1/2,-x+1/2", "-y+1/2,z+1/2,x+1/2", "-z+1/2,-y+1/2,x+1/2",
        "y+1/2,-z+1/2,x+1/2", "z+1/2,y+1/2,x+1/2", "y+1/2,z+1/2,-x+1/2", "-z+1/2,y+1/2,-x+1/2", "z+1/2,-y+1/2,-x+1/2"
    ]
}, {
    "bravais_lattice":
        "cubic",
    "international_number":
        230,
    "hermann_mauguin":
        "I a -3 d",
    "genpos": [
        "x,y,z", "-y+1/4,x+3/4,z+1/4", "-x+1/2,-y,z+1/2", "y+1/4,-x+1/4,z+3/4", "x,-y,-z+1/2", "y+1/4,x+3/4,-z+3/4",
        "-x+1/2,y,-z", "-y+1/4,-x+1/4,-z+1/4", "z,x,y", "-x+1/4,z+3/4,y+1/4", "-z+1/2,-x,y+1/2", "x+1/4,-z+1/4,y+3/4",
        "z,-x,-y+1/2", "x+1/4,z+3/4,-y+3/4", "-z+1/2,x,-y", "-x+1/4,-z+1/4,-y+1/4", "y,z,x", "y+1/2,-z+1/2,-x",
        "z+3/4,y+1/4,-x+1/4", "-y,z+1/2,-x+1/2", "-z+1/4,-y+1/4,-x+1/4", "-y+1/2,-z,x+1/2", "z+3/4,-y+3/4,x+1/4",
        "-z+3/4,y+1/4,x+3/4", "-x,-y,-z", "y-1/4,-x-3/4,-z-1/4", "x-1/2,y,-z-1/2", "-y-1/4,x-1/4,-z-3/4", "-x,y,z-1/2",
        "-y-1/4,-x-3/4,z-3/4", "x-1/2,-y,z", "y-1/4,x-1/4,z-1/4", "-z,-x,-y", "x-1/4,-z-3/4,-y-1/4", "z-1/2,x,-y-1/2",
        "-x-1/4,z-1/4,-y-3/4", "-z,x,y-1/2", "-x-1/4,-z-3/4,y-3/4", "z-1/2,-x,y", "x-1/4,z-1/4,y-1/4", "-y,-z,-x",
        "-y-1/2,z-1/2,x", "-z-3/4,-y-1/4,x-1/4", "y,-z-1/2,x-1/2", "z-1/4,y-1/4,x-1/4", "y-1/2,z,-x-1/2",
        "-z-3/4,y-3/4,-x-1/4", "z-3/4,-y-1/4,-x-3/4", "x+1/2,y+1/2,z+1/2", "-y+3/4,x+5/4,z+3/4", "-x+1,-y+1/2,z+1",
        "y+3/4,-x+3/4,z+5/4", "x+1/2,-y+1/2,-z+1", "y+3/4,x+5/4,-z+5/4", "-x+1,y+1/2,-z+1/2", "-y+3/4,-x+3/4,-z+3/4",
        "z+1/2,x+1/2,y+1/2", "-x+3/4,z+5/4,y+3/4", "-z+1,-x+1/2,y+1", "x+3/4,-z+3/4,y+5/4", "z+1/2,-x+1/2,-y+1",
        "x+3/4,z+5/4,-y+5/4", "-z+1,x+1/2,-y+1/2", "-x+3/4,-z+3/4,-y+3/4", "y+1/2,z+1/2,x+1/2", "y+1,-z+1,-x+1/2",
        "z+5/4,y+3/4,-x+3/4", "-y+1/2,z+1,-x+1", "-z+3/4,-y+3/4,-x+3/4", "-y+1,-z+1/2,x+1", "z+5/4,-y+5/4,x+3/4",
        "-z+5/4,y+3/4,x+5/4", "-x+1/2,-y+1/2,-z+1/2", "y+1/4,-x-1/4,-z+1/4", "x,y+1/2,-z", "-y+1/4,x+1/4,-z-1/4",
        "-x+1/2,y+1/2,z", "-y+1/4,-x-1/4,z-1/4", "x,-y+1/2,z+1/2", "y+1/4,x+1/4,z+1/4", "-z+1/2,-x+1/2,-y+1/2",
        "x+1/4,-z-1/4,-y+1/4", "z,x+1/2,-y", "-x+1/4,z+1/4,-y-1/4", "-z+1/2,x+1/2,y", "-x+1/4,-z-1/4,y-1/4",
        "z,-x+1/2,y+1/2", "x+1/4,z+1/4,y+1/4", "-y+1/2,-z+1/2,-x+1/2", "-y,z,x+1/2", "-z-1/4,-y+1/4,x+1/4",
        "y+1/2,-z,x", "z+1/4,y+1/4,x+1/4", "y,z+1/2,-x", "-z-1/4,y-1/4,-x+1/4", "z-1/4,-y+1/4,-x-1/4"
    ]
}]