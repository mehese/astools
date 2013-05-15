from structures import *
import math
import itertools as its
#from __future__ import division

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

def expand(structure, X = (0., 0.), Y = (0., 0.), 
           Z= (0., 0.)) :
    """Expands a structure by repeating its unit cell until its dimensions 
    reach the given parameters
    (AtomStruct, (float, float), (float, float), (float, float)) -> AtomStruct
    
    """

    (xn, xp) = map(int, X) 
    dnx = math.modf(X[0])[0]
    dpx = math.modf(X[1])[0]
    (yn, yp) = map(int, Y) 
    dny = math.modf(Y[0])[0]
    dpy = math.modf(Y[1])[0]
    (zn, zp) = map(int, Z) 
    dnz = math.modf(Z[0])[0]
    dpz = math.modf(Z[1])[0]
    #print xn, xp, range(xn, xp+1)
    #print yn, yp, range(yn, yp+1)
    #print zn, zp, range(zn, zp+1)
    #print dnx, dpx, dny, dpy, dnz, dpz
    # Initialising atom list
    ats  = []
    for xx, yy, zz in its.product(range(xn, xp+1), range(yn, yp+1), 
                                  range(zn, zp+1)):
        print xx, yy, zz
        for at in structure.atoms:
            at_ = Atom(at.species, at.x + xx*structure.coordx, 
                       at.y + yy*structure.coordy, 
                       at.z + zz*structure.coordz)
            
            try :
                at_.Charge(at.charge)
            except AttributeError:
                at_.Charge(0.)
            ats.append(at_)


    newstruct = AtomStruct(ats, (2., 2., 2., 90., 90., 
                           90.))

    return newstruct
