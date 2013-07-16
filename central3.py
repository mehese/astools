#!/usr/bin/env python

import numpy as np
import random as rn
import itertools as it

from structures import *
from ReadWrite import *
from operations import *
from analysis import *

def chargey(struct):
    for at in struct.atoms:
        if at.species == 'Si':
            at.Charge(2.0)
        elif at.species == 'O':
            at.Charge(-1.0)

def main() :
    SiO2Bulk = ReadStruct('CrystalCell', style='crystal')
    SiO2Bulk.normalise()

    # Creating bulk Si cell

    atlist = [Atom('Si', 0.6779, 0.6779, 0.6779), 
              Atom('Si', 3.3879, 3.3879, 0.6779),
              Atom('Si', 2.0329, 4.7429, 2.0329),
              Atom('Si', 4.7429, 2.0329, 2.0329),
              Atom('Si', 0.6779, 3.3879, 3.3879),
              Atom('Si', 3.3879, 0.6779, 3.3879),
              Atom('Si', 2.0329, 2.0329, 4.7429),
              Atom('Si', 4.7429, 4.7429, 4.7429)]
    SiBulk = AtomStruct(atlist, (5.42, 5.42, 5.42, 90.0, 90.0, 90.0))
    for atom in SiBulk.atoms:
        atom.Charge(0.0)

    # Flipping coordinates of cell
    for at in SiO2Bulk.atoms:
        tmp = at.z
        at.z = at.y
        at.y = tmp
    tmp = SiO2Bulk.coordz
    SiO2Bulk.coordz = SiO2Bulk.coordy
    SiO2Bulk.coordy = tmp
    #SiO2Bulk = expand(SiO2Bulk, Z = (0., 1.))
    SiO2Bulk.normalise()
    # Setting the cell charge
    #SiO2Bulk = expand(SiO2Bulk, X = (0., 1.), Y=(0., 1.), Z=(0., 0.5))
    chargey(SiO2Bulk)
    SiO2Bulk = repeat(SiO2Bulk, 2, 2, 1)
    SiO2Bulk = expand(SiO2Bulk, Z=(0., 0.5))
    chargey(SiO2Bulk)
    print SiO2Bulk.charge()

    PrintStruct(SiO2Bulk, 'crystal_inp', name='INPUT_oxide')
    PrintStruct(SiO2Bulk, 'lmp_data', name='data.SiO2tomelt')

    #SiBulk = expand(SiBulk, Z = (-1., 1.))


    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
