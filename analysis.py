#! /usr/bin/env python
# -*- coding: utf-8 -*-p

import numpy as np
import matplotlib.pylab as plt
from operations import *
from operator import itemgetter, attrgetter
from itertools import product


global coordination_dict 
coordination_dict = {'H': 1, 'O':2, 'Si': 4, 'Hf':4}

class dist :
    def __init__(self, species, r) :
        self.atom_type = species
        self.length = r

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

def distance(at1, at2):
    """Returns the distance between two atoms
    (Atom, Atom) -> float
    """
    return np.sqrt((at1.x - at2.x)*(at1.x - at2.x) + 
                   (at1.y - at2.y)*(at1.y - at2.y) +
                   (at1.z - at2.z)*(at1.z - at2.z))

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
                    Y=(-dist/structure.coordx, dist/structure.coordx),
                    Z=(-dist/structure.coordx, dist/structure.coordx))
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
        g[i] = g[i]/(6*r[i]*r[i]*dr + 2*dr*dr*dr)

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

def nearest_neighbors(structure, dmax = 10.):
    """Returns nearest neighbors
    (AtomStruct, float) -> list of at_neigbors

    """
    
    for atom in structure.atoms:
        atom.tags.append('original')

    neighbor_lst = []

    # expand the structure all the way to dmax
    newstr = expand(structure,
                    X=(-dmax/structure.coordx, dmax/structure.coordx),
                    Y=(-dmax/structure.coordx, dmax/structure.coordx),
                    Z=(-dmax/structure.coordx, dmax/structure.coordx))

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

def main():
    lO1, lO2, lSi1, lSi2 = [], [], [], []
    fmin = 0.9
    j = 1
    positions = range(-1, -2001, -50)
    #positions = range(-1, -2, -50)
    for i in positions : 
        temp_range =  range(3000, 6000, 500)
        #temp_range =  range(3000, 3001, 500)
        for k in temp_range:
            print '\nCalculating structure {:4} out of {:4} ...\n'.format(
                   j, len(positions)*len(temp_range)*2)
            struct = ReadStruct(
                     '/home/eric/lammps-6Dec12/Si/SiO2melt/cmp/dump.SiO2Simelt'+\
                     str(k)+'t200r5000', style='lmp_dump', pos=i)
        
            f, _ = neighbor_statistics(struct, verbose=True, limit=0.2)
            if f < fmin:
                fmin, structmin = f, struct
                PrintStruct(structmin, 'crystal_inp', name='INPUT_mindefects')
            j += 1
            print '\n'+60*'-'+'\n'

            print '\nCalculating structure {:4} out of {:4} ...\n'.format(
                   j, len(positions)*len(temp_range)*2)
            struct = ReadStruct(
                     '/home/eric/lammps-6Dec12/Si/SiO2melt/cmp/dump.SiO2Simelt'+\
                     str(k)+'t200r10000', style='lmp_dump', pos=i)
        
            f, _ = neighbor_statistics(struct, verbose=True, limit=0.2)
            if f < fmin:
                fmin, structmin = f, struct
                PrintStruct(structmin, 'crystal_inp', name='INPUT_mindefects')
            j += 1
            print '\n'+60*'-'+'\n'

    print 'Min percentage of defects {:10.5f} %'.format(fmin)

    print 'Done!'

if __name__ == "__main__":
    main()
