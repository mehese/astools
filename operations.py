#!/usr/bin/env python

from structures import *
from ReadWrite import *
import math
import itertools as its
#from __future__ import division

##############################################################################
#                                                                            #
#    Contains functions that modify structures defined in structures.py      #
#                                                                            #
##############################################################################

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
            at_.tags = list(at.tags)
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

    # find multiplication integers (xn, xp, etc)
    # also find boundaries we need to expand
    (xn, xp) = map(int, X) 
    (dnx, dpx) = (math.modf(X[0])[0]*structure.coordx,
                  math.modf(X[1])[0]*structure.coordx)
    (yn, yp) = map(int, Y) 
    (dny, dpy) = (math.modf(Y[0])[0]*structure.coordy, 
                  math.modf(Y[1])[0]*structure.coordy)
    (zn, zp) = map(int, Z) 
    (dnz, dpz) = (math.modf(Z[0])[0]*structure.coordz, 
                  math.modf(Z[1])[0]*structure.coordz)
    #print xn, xp, range(xn, xp+1)
    #print yn, yp, range(yn, yp+1)
    #print zn, zp, range(zn, zp+1)
    #print dnx, dpx, dny, dpy, dnz, dpz
    # Initialising atom list
    ats  = []
    # Expanding unit cell
    for xx, yy, zz in its.product(range(xn, xp+1), range(yn, yp+1), 
                                  range(zn, zp+1)):
        #print xx, yy, zz
        for at in structure.atoms:
            at_ = Atom(at.species, at.x + xx*structure.coordx, 
                       at.y + yy*structure.coordy, 
                       at.z + zz*structure.coordz)
            
            if (xx, yy, zz) == (0, 0, 0):
                at_.tags = at.tags
            else:
                at_.tags = [tag for tag in at.tags if tag != 'original']
            try :
                at_.Charge(at.charge)
            except AttributeError:
                at_.Charge(0.)
            ats.append(at_)
    # Expanding franctional boundaries
    newats = [] 
    for at in ats:
        # lower atoms
        if at.x < structure.coordx*xn + dpx:
            newats.append(Atom(at.species, 
                               at.x + structure.coordx*(-1*xn+xp+1.),
                               at.y, at.z))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
        # higher atoms
        if at.x > structure.coordx*(xp+1.) + dnx:
            newats.append(Atom(at.species, 
                               at.x - structure.coordx*(-1*xn+xp+1.),
                               at.y, at.z))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
    ats.extend(newats)
    newats = []
    for at in ats:
        # lower atoms
        if at.y < structure.coordy*yn + dpy:
            newats.append(Atom(at.species, at.x,
                               at.y + structure.coordy*(-1*yn+yp+1.), at.z))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
        # higher atoms
        if at.y > structure.coordy*(yp+1.) + dny:
            newats.append(Atom(at.species, at.x,
                               at.y - structure.coordy*(-1*yn+yp+1.), at.z))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
    ats.extend(newats)
    newats = []
    for at in ats:
        # lower atoms
        if at.z < structure.coordz*zn + dpz:
            newats.append(Atom(at.species, at.x, at.y,
                               at.z + structure.coordz*(-1*zn+zp+1.)))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
        # higher atoms
        if at.z > structure.coordz*(zp+1.) + dnz:
            newats.append(Atom(at.species, at.x, at.y,
                               at.z - structure.coordz*(-1*zn+zp+1.)))
            newats[-1].tags = [tag for tag in at.tags if tag != 'original']
    ats.extend(newats)
    for at in ats:
        at.x = at.x + structure.coordx*(-1*xn) + (-1)*dnx
        at.y = at.y + structure.coordy*(-1*yn) + (-1)*dny
        at.z = at.z + structure.coordz*(-1*zn) + (-1)*dnz
    coord_tuple = (structure.coordx*(-1*xn+xp+1.) + (-1)*dnx + dpx, 
                   structure.coordy*(-1*yn+yp+1.) + (-1)*dny + dpy,
                   structure.coordz*(-1*zn+zp+1.) + (-1)*dnz + dpz, 
                   structure.alpha, structure.beta, structure.gamma)
        
    newstruct = AtomStruct(ats, coord_tuple)
    return newstruct

def slab_hfo2(structure, X, Y, Z) :
    ats = []
    cx = structure.coordx*X
    cy = structure.coordy*Y
    # sin(9.419 deg) = 0.1636531
    # cos(9.419 deg) = 0.9865179
    h = structure.coordz*0.9865179
    dx = structure.coordz*0.1636531

    for i, j in its.product(range(X),
                            range(Y)):
        #print i, j
        for at in structure.atoms:
            x = at.x + i*structure.coordx
            y = at.y + j*structure.coordx
            z = at.z
            spec = at.species
            at_ = Atom(spec, x, y, z)
            ats.append(at_)
    structure.atoms = ats[:]
    for i in range(1,int(Z)) :
        for at in structure.atoms:
            x = at.x - dx*i  
            y = at.y
            z = at.z + h*i
            spec = at.species
            at_ = Atom(spec, x, y, z)
            ats.append(at_)

    cz = h*Z
    crystal = AtomStruct(ats, (cx, cy, cz, 90., 90., 90.))
    return crystal

def main():
    HfO2 = ReadStruct('HfO2_out', style='crystal_out')
    HfO2 = slab_hfo2(HfO2, 2, 2, 1.0)
    PrintStruct(HfO2, 'crystal_inp', name='INPUT_newcell_hfo2')
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
