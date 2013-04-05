from structures import *
import itertools as its


def widen(structure, ex, ey, ez) :
    """Returns an (ex X ey X ez) widenened AtomStruct
    (AtomStruct, int, int, int) -> AtomStruct
    """
    coord_x = ex*structure.coordx
    coord_y = ex*structure.coordy
    coord_z = ex*structure.coordz
    #ats = [Atom(at.species, xx*at.x, yy*at.y, zz*at.z) for at in \
    #       structure.atoms for xx in range(1, ex+1) for yy in range(1, ey+1) \
    #       for zz in range(1, ez+1)]
    ats = []
    for xx, yy, zz in its.product(range(1, ex+1), range(1, ey+1), 
                                  range(1,ez+1)):
        print xx, yy, zz
        for at in structure.atoms:
            at_ = Atom(at.species, at.x*xx, at.y*yy, at.z*zz)
            at_.Charge(at.charge)
            ats.append(at_)

    newstruct = AtomStruct(ats, (coord_x, coord_y, coord_z, structure.alpha,
                           structure.beta, structure.gamma))
    return newstruct
