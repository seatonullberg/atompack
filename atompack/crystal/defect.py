import copy

from atompack.atom import Atom
from atompack.crystal.crystal import Crystal
from atompack.elements import Element
from atompack.select import Selection


def substitution_defect(crystal: Crystal, selection: Selection, element: Element) -> Crystal:
    res = copy.deepcopy(crystal)
    indices = selection(crystal)
    kwargs = {k: getattr(element, k) for k in dir(element) if not k.startswith("_")}
    for index in indices:
        atom = copy.deepcopy(res.atoms[index])
        new_atom = Atom(atom.position, **kwargs)
        res.remove(index)
        res.insert(new_atom)
    return res


def vacancy_defect(crystal: Crystal, selection: Selection) -> Crystal:
    res = copy.deepcopy(crystal)
    indices = selection(crystal)
    for index in indices:
        res.remove(index)
    return res


# TODO
def tetrahedral_defect(crystal: Crystal, selection: Selection) -> Crystal:
    pass


# TODO
def octahedral_defect(crystal: Crystal, selection: Selection) -> Crystal:
    pass
