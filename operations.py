from structures import *
import itertools as its


def repeat(structure, ex, ey, ez) :
    """Returns an (ex X ey X ez) widenened AtomStruct
    (AtomStruct, int, int, int) -> AtomStruct
    """
    coord_x = ex*structure.coordx
    coord_y = ey*structure.coordy
    coord_z = ez*structure.coordz
    ats = []
    for xx, yy, zz in its.product(range(ex), range(ey), 
                                  range(ez)):
        # print xx, yy, zz
        for at in structure.atoms:
            at_ = Atom(at.species, at.x + xx*structure.coordx, 
                       at.y + yy*structure.coordy, 
                       at.z + zz*structure.coordz)
            at_.Charge(at.charge)
            ats.append(at_)

    newstruct = AtomStruct(ats, (coord_x, coord_y, coord_z, structure.alpha,
                           structure.beta, structure.gamma))
    return newstruct
