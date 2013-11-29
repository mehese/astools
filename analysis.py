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

    # print sorted([(at, distance(i, atom)) for at in structure.atoms], 
    #         key=itemgetter(1))
    return (r, g)

def pairof4(structure, dmax = 10.):
    """Returns nearest 4 neighbors, and next nearest 4 neighbors
    (AtomStruct, int) -> ((list, list), (list, list))

    """
    #class adistance:
    #    def __init__(self, r, label) :
    #        self.d = r
    #        self.label = label

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
    #PrintStruct(newstr, 'crystal_inp') 
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

def neighbor_statistics(structure, dmax = 6., limit=0.20):
    """Returns nearest neighbors
    (AtomStruct, float) -> list of at_neigbors

    """

    bond_lengths = np.zeros((100, 100))
    bond_lengths[8][14] = bond_lengths[14][8] = 1.543
    bond_lengths[14][14] = 2.34
    
    print bond_lengths

    neighbor_lst = nearest_neighbors(structure)

    print 'got the neighbors'

    for at1 in neighbor_lst :
        #print at1.atom_type, species2Z[at1.atom_type]
        print at1
        #for at2 in at1.neighbors:
        #    print '\t {} {}'.format(at2.atom_type, species2Z[at.atom_type])

    return neighbor_lst 




def main():
    lO1, lO2, lSi1, lSi2 = [], [], [], []
    for i in range(-1, -2, -1):
        print i
        struct = ReadStruct('dump.SiO2Simelt', style='lmp_dump', pos=i)
        
        stats = neighbor_statistics(struct)

        #n_lst = nearest_neighbors(struct, dmax = 6.0)
        #for at1 in n_lst:
        #    for at2 in at1.neighbors:
        #        if at1.atom_type == 'O' and at2.atom_type == 'O':
        #            lO1.append(at2.length)
        #        if at1.atom_type == 'O' and at2.atom_type == 'Si':
        #            lO2.append(at2.length)
        #        if at1.atom_type == 'Si' and at2.atom_type == 'Si':
        #            lSi1.append(at2.length)
        #        if at1.atom_type == 'Si' and at2.atom_type == 'O':
        #            lSi2.append(at2.length)
             
    #plt.title('Nearest 2 neighbours, O centred')
    #plt.hist(lO1, 100, label='O-O bonds')
    #plt.hist(lO2, 100, label='O-Si bonds')
    #plt.xlabel('pair distance $(\\AA)$')
    #plt.legend()
    #plt.figure()
    #plt.title('Nearest 4 neighbours, Si centred')
    #plt.hist(lSi1, 100, label='Si-Si bonds')
    #plt.hist(lSi2, 100, label='Si-O bonds')
    #plt.xlabel('pair distance $(\\AA)$')
    #plt.legend()
    #plt.show()


    #sio = plt.hist(p1, 100, label='Si-O bonds')
    #sisi = plt.hist(p2, 100, label='Si-Si bonds')
    #oo = plt.hist(p3, 100, label = 'O-O bonds')
    #plt.title('Nearest 4 neighbors histogram SiO$_2$ (40 snapshot)', fontsize=20)
    #plt.xlabel('Pair distance $(\\AA)$', fontsize=18)
    #plt.legend(prop={'size':25})
    #plt.yticks([])
    #plt.figure()
    #pairs1 = [dist[0] for dist in rlist2 if dist[1] == 'SiO']
    #pairs2 = [dist[0] for dist in rlist2 if dist[1] == 'SiSi']
    #pairs3 = [dist[0] for dist in rlist2 if dist[1] == 'OO']
    #sio = plt.hist(pairs1, 100, label='SiO bonds')
    #sisi = plt.hist(pairs2, 100, label='SiSi bonds')
    #oo = plt.hist(pairs3, 100, label = 'OO bonds')
    #plt.title('neighbors 5-8 histogram SiO2Si, run 1')
    #plt.xticks(fontsize=17)

    print 'Done!'

if __name__ == "__main__":
    main()
