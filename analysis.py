import numpy as np
import matplotlib.pylab as plt
from operations import *

def distance(at1, at2):
    """Returns the distance between two atoms
    (Atom, Atom) -> float
    """
    return np.sqrt((at1.x - at2.x)*(at1.x - at2.x) + 
                   (at1.y - at2.y)*(at1.y - at2.y) +
                   (at1.z - at2.z)*(at1.z - at2.z))

def rdf(structure, nbin, dist=15.):
    """Returns the pair distribution function for a system in nbin bins
    (AtomStruct, int) -> (list, list)
    """
    for atom in structure.atoms:
        atom.tags.append('original')
    r = [(i*dist/nbin + dist/(2*nbin)) for i in range(nbin)]
    print r
    g = []
    return (r, g)
