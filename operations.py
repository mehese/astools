#!/usr/bin/env python

from structures import *
from analysis import *
from ReadWrite import *
import math
import numpy as np
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
        for at in structure.atoms:
            at_ = Atom(at.species, at.x + xx*structure.coordx, 
                       at.y + yy*structure.coordy, 
                       at.z + zz*structure.coordz)
            at_.tags = list(at.tags)
            # try to give the atoms the same charge
            try :
                at_.Charge(at.charge)
            except AttributeError:
                at_.Charge(0.0)
                at.Charge(0.0)
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

def getH(at1, at2, verbose=False):
    """Returns an Atom object mirroring the position of at2 with respect to at1
    (Atom, Atom) -> Atom
    """
    
    if verbose:
        print 'Computing coordinates for H atom'
    if at1.species == 'Si':
        dHr = 1.48
    elif at1.species == 'O':
        dHr = 1.09
    else :
        dHr = 1.3

    dx = at1.x - at2.x
    dy = at1.y - at2.y
    dz = at1.z - at2.z
    
    dr = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    dHx = dHr*dx/dr 
    dHy = dHr*dy/dr 
    dHz = dHr*dz/dr 

    H = Atom('H', at1.x + dHx, at1.y + dHy, at1.z + dHz)

    if verbose:
        print '--> DISTANCE from main atom:', distance(at1, H), ' Angstroem'
        print '--> Colinearity area of points (the closer to 0 the better) :',
        # Check (18) @ http://mathworld.wolfram.com/TriangleArea.html 
        m1 = np.array([[at1.y, at1.z,1],[at2.y,at2.z,1],[H.y,H.z,1]])
        m2 = np.array([[at1.z, at1.x,1],[at2.z,at2.x,1],[H.z,H.x,1]])
        m3 = np.array([[at1.x, at1.y,1],[at2.x,at2.y,1],[H.x,H.y,1]])
        delta = 0.5*np.sqrt(np.linalg.det(m1)**2 + np.linalg.det(m2)**2 + 
                            np.linalg.det(m3)**2) 
        print delta, ' A^3'
        print 
    return H

def get_tet_centre(at0, at1, at2, at3):
    """Given four atoms it will find the coordinates for the projection of at0
    on the plane determined by at1, at2 & at3
    (Atom, Atom, Atom, Atom) -> numpy.array
    """
    proj = np.zeros(3)
    c0 = np.array([at0.x, at0.y, at0.z])
    c1 = np.array([at1.x, at1.y, at1.z])
    c2 = np.array([at2.x, at2.y, at2.z])
    c3 = np.array([at3.x, at3.y, at3.z])
    # Get the 3 points through the plane
    # http://en.wikipedia.org/wiki/Plane_(geometry), let d = -10
    D = np.array([c1, c2, c3]).transpose()
    # -d/D -> fact
    fact = -10/np.linalg.det(D)
    # calculating the main factors a, b, c
    a = np.copy(D); a[:,0] = np.ones(3); a = fact*np.linalg.det(a) 
    b = np.copy(D); b[:,1] = np.ones(3); b = fact*np.linalg.det(b) 
    c = np.copy(D); c[:,2] = np.ones(3); c = fact*np.linalg.det(c) 
    n = np.array([a, b, c])
    n = n/np.linalg.norm(n)
    # another way of calculating the normal to a plane, this shit might be safer
    n2 = np.cross(c2-c1, c3-c1)
    n2 = n2/np.linalg.norm(n2)

    # pray this shit is correct
    # http://stackoverflow.com/questions/8942950/how-do-i-find-the-orthogonal-projection-of-a-point-onto-a-plane 
    proj = c0 - np.dot(c0 - c3, n2)*n2

    return proj 

def passivateSi(s):
    """Passivates all undercoordinated Si atoms with a hydrogen
    AtomStruct -> AtomStruct

    """
    import copy

    newstruct = copy.deepcopy(s) 

    bls = {'Si-O':1.62, 'Si-Si':2.39, 'Si-H':1.48}

    s.ClearTags()

    # Iterate though all the Si atoms in the structure
    for a, i in zip(s.atoms, range(1, len(s)+1)):
        if a.species == 'Si':
            nbs = get_neighbours(a, s, dmax=5., verbose=False)
            # Let's suppose our atoms is fully coordinated
            uncoord = 0
            # Iterate through all the neighbours of atom a
            for k in nbs :
                ideal = bls['Si-'+k.at.species] 
                err = abs(ideal - k.length)/ideal
                # if the bond length is overstretched
                if err > 0.3 :
                    # Our uncoordinations increases
                    uncoord += 1

            a.tags.append('u'+str(uncoord))

            if uncoord == 1:
                cds = tuple(get_tet_centre(a, nbs[0].at, nbs[1].at, nbs[2].at))
                #x = sum([nbs[i].at.x for i in range(3)])/3  
                #y = sum([nbs[i].at.y for i in range(3)])/3  
                #z = sum([nbs[i].at.z for i in range(3)])/3  
                #H = getH(a, Atom('X', x, y, z), verbose=False) 
                H = getH(a, Atom('X', *cds), verbose=False) 
                newstruct.atoms.append(H)
                #newstruct.atoms.append(Atom('X', *cds))

            if uncoord == 2:
                # Ok, let's refer everything to the origin of our axis, let
                # T_vec be the translation vector
                T_vec = np.array([a.x, a.y, a.z])
                v1 = np.array([nbs[0].at.x, nbs[0].at.y, nbs[0].at.z]) - T_vec
                v2 = np.array([nbs[1].at.x, nbs[1].at.y, nbs[1].at.z]) - T_vec
                # The cross product will give us something perpendicular to the
                # plane of the 3 atoms
                cp = np.cross(v1, v2)
                cp = (cp * 1.48)/np.linalg.norm(cp)
                # Don't forget to translate it back up
                H1 = Atom('H', *tuple(cp+T_vec))
                proj_coords = get_tet_centre(a, H1, nbs[0].at, nbs[1].at)
                X = Atom('X', *tuple(proj_coords))
                #newstruct.atoms.append(X)
                H2 = getH(a, X, verbose=False)
                newstruct.atoms.extend([H1, H2])

            if uncoord > 2:
                print 'WARNING!!!! '*4
                print
                print 'Silicon that has a coordination number of 1 (or lower) found!!!'
                print 'Atom id = ', i
                print 'Atom info: ', at, 'tags: ', at.tags
                print 'WARNING!!!! '*4
                print
                print 'WARNING!!!! '*4
                exit()
    
    return newstruct


def passivateO(s):
    """Passivates all undercoordinated O atoms with a hydrogen
    AtomStruct -> AtomStruct

    See passivateSi for more details

    """
    import copy
    bl_opt = 1.62

    newstruct = copy.deepcopy(s) 

    for a, i in zip(s.atoms, range(1, len(s)+1)):
        if a.species == 'O':
            nbs = get_neighbours(a, s, dmax=5., verbose=False)
            uncoord = 0
            for k in nbs :
                ideal = bl_opt 
                err = abs(ideal - k.length)/ideal
                # if the bond length is overstretched
                if err > 0.3 :
                    # Our uncoordinations increases
                    uncoord += 1

            a.tags.append('u'+str(uncoord))

            if uncoord == 1:
                H1 = getH(a, nbs[0].at, verbose=True)
                newstruct.atoms.append(H1)

    return newstruct



def main():
    import matplotlib.pylab as plt
    struct = ReadStruct('INPUT_nopass', style='crystal')
    ns = passivateSi(struct)
    ns = passivateO(ns)
    #x, y = rdf(struct, 300, dist=5.)
    #HfO2 = slab_hfo2(HfO2, 2, 2, 1.0)
    #plt.plot(x, y, 'r-')
    #plt.show()
    PrintStruct(ns, 'crystal_inp', name='INPUT_passivated')
    print 'Done!'

# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
