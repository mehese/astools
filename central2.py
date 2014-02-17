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
    SiO2Bulk = ReadStruct('INPUT_mindefects', style='crystal')
    SiO2Bulk.normalise()
    for at in SiO2Bulk.atoms:
        at.tags.append('oxide')
    chargey(SiO2Bulk)

    deltax = 0.5

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
    SiBulk = repeat(SiBulk, 2, 2, 1)
    SiBulk = expand(SiBulk, Z=(0., 0.5))
    for atom in SiBulk.atoms:
        atom.Charge(0.0)

    SiO2Bulk.normalise()
    cds = (SiBulk.coordx, SiBulk.coordy,
           SiBulk.coordz + SiO2Bulk.coordz + deltax,
           90., 90., 90.)
    newstruct = AtomStruct(SiBulk.atoms+SiO2Bulk.atoms, cds)
    for at in newstruct.atoms:
        if 'oxide' in at.tags:
            at.z = SiBulk.coordz + at.z + deltax - 0.2

    for at in newstruct.atoms:
        if 'oxide' in at.tags:
            #if at.species == 'Si':
            #    at.Charge(2.0)
            #elif at.species == 'O':
            #    at.Charge(-1.0)
            #else :
            #    at.Charge(0.0)
            at.x = (at.x * newstruct.coordx)/ SiO2Bulk.coordx
            at.y = (at.y * newstruct.coordy)/ SiO2Bulk.coordy
    PrintStruct(newstruct, 'castep_inp', name='SiO2Si.cell')
    PrintStruct(newstruct, 'crystal_inp', name='INPUT_SiO2Si')
    PrintStruct(newstruct, 'lmp_data', name='data.SiO2Sinewstruct')
    print 'Structure has', len(newstruct.atoms), 'atoms'
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
