#! /usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from operations import *
from operator import itemgetter, attrgetter
from itertools import product

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
    """Returns nearest 4 neighbours, and next nearest 4 neighbours
    (AtomStruct, int) -> ((list, list), (list, list))

    """
    class adistance:
        def __init__(self, r, label) :
            self.d = r
            self.label = label

    for atom in structure.atoms:
        atom.tags.append('original')
    dlist = []
    d2list = []

    # expand the structure all the way to dmax
    newstr = expand(structure,
                    X=(-dmax/structure.coordx, dmax/structure.coordx),
                    Y=(-dmax/structure.coordx, dmax/structure.coordx),
                    Z=(-dmax/structure.coordx, dmax/structure.coordx))
    
    # redistribute atoms such as the ones in the original structure come first 
    newstr = [tup[0] for tup in sorted([(at, 'original' in at.tags) \
              for at in newstr.atoms], key = itemgetter(1), 
              reverse = True)]
    for atx in newstr[:len(structure.atoms)] :
        distances = [(99999., ''), (99999., ''), (99999., ''), (99999., ''),  
                    (99999., ''), (99999., ''), (99999., ''), (99999., '')]
        for aty in [at for at in newstr if at != atx] : 
            #print '--', max(distances, key = itemgetter(0))[0]
            max_ = max(distances, key = itemgetter(0))
            if distance(atx, aty) < max_[0] :
                if ((atx.species, aty.species) == ('Si', 'O') or \
                   (aty.species, atx.species == 'Si', 'O')):
                    tag = 'SiO'
                    #pass
                if (atx.species == aty.species == 'Si'):
                    tag = 'SiSi' 
                if (atx.species == aty.species == 'O'):
                    tag = 'OO'
                distances.remove(max_)
                distances.append((distance(atx, aty), tag))
        # add first 4 neighbours to first 4 neighbours list
        dlist.extend(distances[:4])
        # add next 4 neighbours to next 4 neighbours list
        d2list.extend(distances[4:])

    return dlist, d2list 

def main():
    #atlist = [Atom('Si', 0.6779, 0.6779, 0.6779), 
    #          Atom('Si', 3.3879, 3.3879, 0.6779),
    #          Atom('Si', 2.0329, 4.7429, 2.0329),
    #          Atom('Si', 4.7429, 2.0329, 2.0329),
    #          Atom('Si', 0.6779, 3.3879, 3.3879),
    #          Atom('Si', 3.3879, 0.6779, 3.3879),
    #          Atom('Si', 2.0329, 2.0329, 4.7429),
    #          Atom('Si', 4.7429, 4.7429, 4.7429)]
    #SiBulk = AtomStruct(atlist, (5.42, 5.42, 5.42, 90.0, 90.0, 90.0))
    SiO2Si = ReadStruct('INPUT_castep', style='crystal')

    rlist, _ = pairof4(SiO2Si)
    pairs1 = [dist[0] for dist in rlist if dist[1] == 'SiO']
    pairs2 = [dist[0] for dist in rlist if dist[1] == 'SiSi']
    pairs3 = [dist[0] for dist in rlist if dist[1] == 'OO']
    sio = plt.hist(pairs1, 100, label='SiO bonds')
    sisi = plt.hist(pairs2, 100, label='SiSi bonds')
    oo = plt.hist(pairs3, 100, label = 'OO bonds')
    plt.title('Nearest 4 neighbours histogram SiO2Si, run 1')
    plt.xlabel('pair distance $(\\AA)$')
    plt.legend()
    plt.show()
    print 'Done!'

if __name__ == "__main__":
    main()
