# -*- coding: utf-8 -*-p

import numpy as np
from astools.operations import *
from astools.ReadWrite import PrintStruct, ReadStruct
from operator import itemgetter, attrgetter
from itertools import product


global coordination_dict 
coordination_dict = {'H': 1, 'O':2, 'Si': 4, 'Hf':4}

# dictionary containing the atomic mass units of elemets, taken from Wolfram
# Alpha
mass_dict = {'H':1.00749, 'O':15.9994, 'Si':28.0855, 'Hf':178.49}

class dist :
    def __init__(self, species, r) :
        self.atom_type = species
        self.length = r

def distance(at1, at2):
    """Returns the distance between two atoms
    (Atom, Atom) -> float
    """
    return np.sqrt((at1.x - at2.x)*(at1.x - at2.x) + 
                   (at1.y - at2.y)*(at1.y - at2.y) +
                   (at1.z - at2.z)*(at1.z - at2.z))

class neighbour:
    def __init__(self, atom, r):
        self.at = atom
        self.length = r
    def __str__(self):
        return '{:2} d = {:10.5f} A'.format(self.at.species, self.length)


class at_neigbors:
    """Class that contains the species of an atoms, and a list of dist with
    its N neighbors, where N is found in coordination_dict
    """
    def __init__(self, species):
        self.atom_type = species
        self.neighbors = [dist('H', 9e+10)]*coordination_dict[species]

    def __str__(self):
        s = 'Main atom: {} {} \n'.format(self.atom_type, 
                                         species2Z[self.atom_type])
        for i, d in zip(range(1, len(self.neighbors)+1), self.neighbors) :
            s = s+'\t neighbor {:2} : {:2} {:3}  length: {:9.6f} Å \n'.format( \
            i, d.atom_type, species2Z[d.atom_type], d.length)

        return s

def rdf(structure, nbin, dist=25.):
    """Returns the pair distribution function for a system in nbin bins
    (AtomStruct, int, float) -> (list, list)
    """
    for atom in structure.atoms:
        atom.tags.append('original')
    dr = dist/float(2*nbin)
    r = [(i*dist/nbin + dist/(2*nbin)) for i in range(nbin)]
    g = [0]*len(r)
    i = structure.atoms[0]
    newstr = expand(structure,
                    X=(-dist/structure.coordx, dist/structure.coordx),
                    Y=(-dist/structure.coordy, dist/structure.coordy),
                    Z=(-dist/structure.coordz, dist/structure.coordz))
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]
    
    # counting the number of atoms at certain distance intervals

    for d in [distance(at1, at2) for at1 in newstr[:len(structure.atoms)] \
              for at2 in newstr if 1e-5 < distance(at1, at2) < dist]:
        for i in range(len(r)):
            if (r[i] - dr) <= d < (r[i] + dr):
                g[i] += 1.

    # normalising with respect to the radius

    for i in range(len(r)):
        g[i] = g[i]/(6*r[i]*r[i]*dr + 2*dr*dr*dr)/len(structure)

    return (r, g)

def rdf_triclinic(structure, nbin, dist=25., verbose=False):
    """Returns the pair distribution function for a system in nbin bins
    (AtomStruct, int, float) -> (list, list)
    Works for triclinic cells
    """

    at_set = set([x.species for x in structure.atoms])
    pair_set = list(set(['-'.join(sorted([x, y])) for x in at_set for y in \
               at_set]))

    if verbose:
        print 'Cell parameters:\na     = {:7.4f}    b = {:7.4f}     c =\
        {:7.4f}'.format(structure.coordx, structure.coordy, structure.coordz)
        print 'alpha = {:7.4f} beta = {:7.4f} gamma =\
        {:5.3f}'.format(structure.alpha, structure.beta, structure.gamma)
    for atom in structure.atoms:
        atom.tags.append('original')
    dr = dist/float(2*nbin)
    r = [(i*dist/nbin + dist/(2*nbin)) for i in range(nbin)]
    g = []
    for i in range(len(pair_set)):
        g.append([0]*len(r))
    i = structure.atoms[0]
    max_xyz = triclinic2xyz(Atom('X', 1., 1., 1.), structure)
    if verbose:
        print 'Direction vectors...'
        print xyz2triclinic(Atom('X', dist, 0., 0.), structure)
        print xyz2triclinic(Atom('X', 0., dist, 0.), structure)
        print xyz2triclinic(Atom('X', 0., 0., dist), structure)
    xx, yy, zz = tuple(map(lambda k: np.modf(dist/k)[1]+3,max_xyz))
    if verbose:
        print 'Expansion coef:\n{:10.6f} {:10.6f}\n{:10.6f} {:10.6f}\n{:10.6f} {:10.6f}'.format(-xx, xx, -yy, yy, -zz, zz)

    newstr = expand2(structure,
                    X=(-xx, xx),
                    Y=(-yy, yy),
                    Z=(-zz, zz))
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]
    
    # counting the number of atoms at certain distance intervals

    #for d in [distance(at1, at2) for at1 in newstr[:len(structure.atoms)] \
    #          for at2 in newstr if 1e-5 < distance(at1, at2) < dist]:
    #    for i in range(len(r)):
    #        if (r[i] - dr) <= d < (r[i] + dr):
    #            g[i] += 1.
    if verbose:
        print 'Starting to check pairs...'
    for i in range(len(pair_set)):
        o1, o2 = tuple(pair_set[i].split('-'))
        #l1 = [at for at in newstr[:len(structure.atoms)] if at.species == o1]
        #l2 = [at for at in newstr[] if at.species == o2]
        for d in [distance(at1, at2) \
                  for at1 in newstr[:len(structure.atoms)] \
                  if at1.species == o1 \
                  for at2 in newstr if (1e-5 < distance(at1, at2) < dist and \
                  at2.species == o2)]:
            for j in range(len(r)):
                if (r[j] - dr) <= d < (r[j] + dr):
                    g[i][j] += 1.
                    if o1 != o2 :
                        g[i][j] += 1
    if verbose:
        print 'Pair checking done, scaling final function'

    # normalising with respect to the radius

    # normalising with respect to the radius

    #for i in range(len(r)):
    #    g[i] = g[i]/(6*r[i]*r[i]*dr + 2*dr*dr*dr)/len(structure)
    for i in range(len(g)):
        for j in range(len(r)):
            g[i][j] = g[i][j]/(6*r[j]*r[j]*dr + 2*dr*dr*dr)/len(structure)
        g[i].append(pair_set[i])

    if verbose:
        print 'RDF done!'

    return (r, g)

def rdf2(structure, nbin, dist=25.):
    """Returns the pair distribution function for a system in nbin bins
    (AtomStruct, int, float) -> (list, list)
    """
    #print 'O Hai'

    at_set = set([x.species for x in structure.atoms])
    pair_set = list(set(['-'.join(sorted([x, y])) for x in at_set for y in \
               at_set]))

    for atom in structure.atoms:
        atom.tags.append('original')
    dr = dist/float(2*nbin)
    r = [(i*dist/nbin + dist/(2*nbin)) for i in range(nbin)]
    g = []
    for i in range(len(pair_set)):
        g.append([0]*len(r))
    #print len(r), len(g[0])
    i = structure.atoms[0]
    newstr = expand(structure,
                    X=(-dist/structure.coordx, dist/structure.coordx),
                    Y=(-dist/structure.coordy, dist/structure.coordy),
                    Z=(-dist/structure.coordz, dist/structure.coordz))
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]
    
    # counting the number of atoms at certain distance intervals


    for i in range(len(pair_set)):
        o1, o2 = tuple(pair_set[i].split('-'))
        #l1 = [at for at in newstr[:len(structure.atoms)] if at.species == o1]
        #l2 = [at for at in newstr[] if at.species == o2]
        for d in [distance(at1, at2) \
                  for at1 in newstr[:len(structure.atoms)] \
                  if at1.species == o1 \
                  for at2 in newstr if (1e-5 < distance(at1, at2) < dist and \
                  at2.species == o2)]:
            for j in range(len(r)):
                if (r[j] - dr) <= d < (r[j] + dr):
                    g[i][j] += 1.
                    if o1 != o2 :
                        g[i][j] += 1

    # normalising with respect to the radius

    for i in range(len(g)):
        for j in range(len(r)):
            g[i][j] = g[i][j]/(6*r[j]*r[j]*dr + 2*dr*dr*dr)/len(structure)
        g[i].append(pair_set[i])

    return (r, g)


def pairof4(structure, dmax = 10.):
    """Returns nearest 4 neighbors, and next nearest 4 neighbors
    (AtomStruct, int) -> ((list, list), (list, list))

    """

    for atom in structure.atoms:
        atom.tags.append('original')
    dlist = []
    d2list = []

    # expand the structure all the way to dmax
    newstr = expand(structure,
                    X=(-dmax/structure.coordx, dmax/structure.coordx),
                    Y=(-dmax/structure.coordx, dmax/structure.coordx),
                    Z=(-dmax/structure.coordx, dmax/structure.coordx))
    #PrintStruct(newstr, 'crystal_inp') 
    # redistribute atoms such as the ones in the original structure come first 
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]

    # this for only iterates though atoms in the original structure
    for atx in newstr[:len(structure.atoms)] :
        distances = [(99999., ''), (99999., ''), (99999., ''), (99999., ''),  
                    (99999., ''), (99999., ''), (99999., ''), (99999., '')]
        # this one iterates through all atoms that are not atx 
        for aty in [at for at in newstr if at != atx] : 
            max_ = max(distances, key = itemgetter(0))
            if distance(atx, aty) < max_[0] :
                if ((atx.species, aty.species) == ('Si', 'O') or \
                   (aty.species, atx.species == 'Si', 'O')):
                    tag = 'SiO'
                if (atx.species == aty.species == 'Si'):
                    tag = 'SiSi' 
                if (atx.species == aty.species == 'O'):
                    tag = 'OO'
                distances.remove(max_)
                distances.append((distance(atx, aty), tag))
        distances = sorted(distances, key = itemgetter(0))

        # add first 4 neighbors to first 4 neighbors list
        dlist.extend(distances[:4])
        # add next 4 neighbors to next 4 neighbors list
        d2list.extend(distances[4:])

    return dlist, d2list 

def get_neighbours(at_main, structure, dmax=10., verbose=False):
    """Returns a list of neighbours for the atom in the structure

    (Atom, AtomStruct, float) -> at_neighbours
    """
    for atom in structure.atoms:
        if 'original' not in atom.tags:
            atom.tags.append('original')
        if (abs(atom.x - at_main.x) < 1e-4 and
            abs(atom.y - at_main.y) < 1e-4 and
            abs(atom.z - at_main.z) < 1e-4):
            atom.tags.append('main')
        
    neighbor_lst = []

    if verbose:
        print 'Expanding structure...'

    # expand the structure all the way to dmax
    newstr = expand(structure,
                    X=(-dmax/structure.coordx, dmax/structure.coordx),
                    Y=(-dmax/structure.coordy, dmax/structure.coordy),
                    Z=(-dmax/structure.coordz, dmax/structure.coordz))

    # redistribute atoms such as the ones in the original structure come first 
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]

    if verbose:
        print 'Done!'
        print 'Identifying atom in structure...'
    # this for only iterates though atoms in the original structure
    for at in newstr[:len(structure.atoms)] :
        if 'main' in at.tags :
            # We identify our main atom in the expanded structure and
            # we also calculate dx, dy, dz to help us to trace back the 
            # coordinates of its neighbours
            atx = at
            if verbose:
                print 'Done identifying!'
            dx = atx.x - at_main.x
            dy = atx.y - at_main.y
            dz = atx.z - at_main.z
    new_elem = at_neigbors(atx.species)
   
    if verbose:
        print 'Creating initial neighbours list...'

    # i will iterate through the coordination number as to create an initial
    # list that has an appropriate length (i.e 2 for oxygen, 4 for Si and Hf)
    i, voisins = 0, []
    while len(voisins) < coordination_dict[atx.species]:
        if newstr[i] != atx:
            d = distance(atx, newstr[i])
            voisins.append(neighbour(newstr[i], d))
            # One less neighbour to worry about
        i = np.random.randint(len(newstr)) 

    if verbose:
        print 'Initial list:', [str(l) for l in voisins]
        print 'Finding closest atoms...'

    # this one iterates through all atoms that are not atx 
    for aty in [at for at in newstr if at != atx] : 
        # Element in voisins list that is the farthest away from atx
        max_ = max(voisins, key = lambda x: x.length)
        # If there's an atom at a smaller distance than than the maximum one
        # in the list of nearest neighbors for atom atx
        newD = distance(atx, aty)  
        if newD < max_.length :
            # replace that element in the list with the new neighbor
            voisins.remove(max_) 
            voisins.append(neighbour(aty, newD)) 

    # Make sure you give the coordinates of the atoms with respect to the
    # atom given initially
    for v in voisins:
        v.at.x -= dx
        v.at.y -= dy
        v.at.z -= dz

    voisins = sorted(voisins, key = lambda x: x.length)
        

    # remove that main tag, if you don't you will get problems when calling 
    # this function again

    for at in structure.atoms:
        if 'main' in at.tags:
            at.tags.remove('main')

    if verbose:
        print 'Done with finding the nearest neighbours of atom ', at_main

    return voisins



def nearest_neighbors(structure, dmax = 10., verbose=False):
    """Returns nearest neighbors
    (AtomStruct, float) -> list of at_neigbors

    """
    
    for atom in structure.atoms:
        if 'original' not in atom.tags:
            atom.tags.append('original')

    neighbor_lst = []

    # expand the structure all the way to dmax
    newstr = expand(structure,
                    X=(-dmax/structure.coordx, dmax/structure.coordx),
                    Y=(-dmax/structure.coordy, dmax/structure.coordy),
                    Z=(-dmax/structure.coordz, dmax/structure.coordz))

    # redistribute atoms such as the ones in the original structure come first 
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]

    # this for only iterates though atoms in the original structure
    for atx in newstr[:len(structure.atoms)] :

        new_elem = at_neigbors(atx.species)
        # this one iterates through all atoms that are not atx 
        for aty in [at for at in newstr if at != atx] : 
            max_ = max(new_elem.neighbors, key = lambda x: x.length)
            # If there's an atom at a smaller distance than than the maximum one
            # in the list of nearest neighbors for atom atx
            newD = distance(atx, aty)  
            if (newD < 1.5) and verbose :
                print 'WARNING2: This atom is too damn close to your main one'
                print at 
                print
            if newD < max_.length :
                # replace that element in the list with the new neighbor
                new_elem.neighbors.remove(max_) 
                new_elem.neighbors.append(dist(aty.species, newD)) 

        new_elem.neighbors = sorted(new_elem.neighbors, 
                                    key = lambda x: x.length)
        
        neighbor_lst.append(new_elem)

    return neighbor_lst 

def neighbor_statistics(structure, dmax = 6., limit=0.20, verbose=True):
    """Returns nearest neighbors
    (AtomStruct, float, float, bool) -> floar, dict

    """

    bond_lengths = np.zeros((100, 100))
    bond_lengths[8][14] = bond_lengths[14][8] = 1.543 # Si-O bonds in Si
    bond_lengths[14][14] = 2.34 # Si-Si in bulk Si
    bond_lengths[8][8] = 2.46 # O-O bonds in SiO2
    

    neighbor_lst = nearest_neighbors(structure)

    if verbose:
        print '  Total number of atoms in structure: {:4}'.format(
               len(structure.atoms))

    no_ats = len(structure.atoms)
    # dictionary species: (broken bonds, total expected)
    at_dict = {x:[0, coordination_dict[x]*len([y for y in structure.atoms if \
               y.species == x])] for x in set([at.species for at \
               in structure.atoms])}

    flawed_ats = []
    for at1 in neighbor_lst :
        #print at1
        for at2 in at1.neighbors:
            if at1.atom_type == at2.atom_type == 'O':
                at_dict[at1.atom_type][0] += 1 
                flawed_ats.append(at1)

            #if at2.length > \
            #bond_lengths[species2Z[at1.atom_type]][species2Z[at2.atom_type]]\
            #*(1+limit):
            if (at2.length > \
            bond_lengths[species2Z[at1.atom_type]][species2Z[at2.atom_type]]\
            *(1+limit)) and (not (at1.atom_type == at2.atom_type == 'O')):
                #print 'overcoordinated'  
                at_dict[at1.atom_type][0] += 1 
                flawed_ats.append(at1)

    fraction_flawed = len(set(flawed_ats))/float(len(structure.atoms))

    if verbose:
        print '  Fraction of atoms with defects: {:10.5f}%'.format(
               fraction_flawed*100)

    if verbose:
        for at, lst in at_dict.items():
            # calculate fraction of broken bonds
            f = float(lst[0])/lst[1]
            print '  For {:2} atoms the fraction of broken bonds is {:10.5f}%'\
                   .format(at, f*100)

    return fraction_flawed, at_dict 

def vertical_density_profile(struct, deltaz, no_points=200, full=False):
    """ Takes in a structure, returns a set of points along the z axis, and the
    corresponding densities
    (AtomStruct, float) -> (float, float)
    density is defined as (mass of atoms with z-deltaz < at.z < z+deltaz) /
    2*deltaz*struct.coordx*struct.coordy
    """
    # go half a structure up, and half a structure up so  
    z = np.linspace(-0.5*struct.coordz, 1.5*struct.coordz, 2*no_points)
    rho = np.zeros(2*no_points)

    work_struct = repeat(struct, 1, 1, 3)
    for at in work_struct.atoms:
        at.z = at.z - struct.coordz

    for i in range(len(z)):
        for at in work_struct.atoms:
            if z[i] - deltaz < at.z < z[i] + deltaz:
                rho[i] = rho[i] + mass_dict[at.species]

    # works for orthorhombic cells only, 1.66 -- conversion to g/cm^3 
    volume = 2*deltaz*struct.coordx*struct.coordy/1.66

    # full returns the whole array, full=False only the z=0, z=zmax one
    if full:
        return z, rho/volume
    else:
        return (z[int(0.5*no_points): int(1.5*no_points)], 
                rho[int(0.5*no_points): int(1.5*no_points)]/volume)

    #z, rho = np.linspace(0., struct.coordz, no_points), np.zeros(no_points)
    #
    #work_struct = repeat(struct, 1, 1, 3)
    #for at in work_struct.atoms:
    #    at.z = at.z - struct.coordz

    #for i in range(len(z)):
    #    for at in work_struct.atoms:
    #        if z[i] - deltaz < at.z < z[i] + deltaz:
    #            rho[i] = rho[i] + mass_dict[at.species]

    ## works for orthorhombic cells only 
    #volume = 2*deltaz*struct.coordx*struct.coordy 
    #return z, 1.66*rho/volume

def main():
    fmin = 0.9
    j = 1
    positions = range(-1, -2001, -50)
    #positions = range(-1, -2, -50)
    struct = ReadStruct('INPUT_Si', 
                        style='crystal')

    struct = repeat(struct, 3, 3 ,3)

    PrintStruct(struct, 'crystal_inp')

    f, _ = neighbor_statistics(struct, verbose=True, limit=0.2)

    print 'Done!'

if __name__ == "__main__":
    main()
