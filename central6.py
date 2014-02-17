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
    SiO2Bulk = ReadStruct('dump.SiO2tomelt', style='lmp_dump')
    SiO2Bulk.normalise()
    for at in SiO2Bulk.atoms:
        at.tags.append('oxide')

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
    SiBulk = repeat(SiBulk, 2, 2, 2)

    SiO2Bulk.normalise()
    cds = (SiBulk.coordx, SiBulk.coordy,
           SiO2Bulk.coordz,
           90., 90., 90.)
    newstruct = AtomStruct(SiO2Bulk.atoms, cds)
    for at in newstruct.atoms:
        if 'oxide' in at.tags:
            at.z = SiBulk.coordz + at.z
    HfO2 = ReadStruct('HfO2_out', style='crystal_out')
    HfO2 = slab_hfo2(HfO2, 2, 2, 2.0)
    HfO2.normalise()
    for at in HfO2.atoms:
        at.tags.append('hfo2')
    newstruct.atoms.extend(HfO2.atoms) 
    newstruct.coordz += HfO2.coordz 
    for at in newstruct.atoms:
        if 'hfo2' in at.tags:
            at.z += newstruct.coordz - 5.
            at.x = (at.x * newstruct.coordx)/ HfO2.coordx
            at.y = (at.y * newstruct.coordy)/ HfO2.coordy
        elif 'oxide' in at.tags:
            at.x = (at.x * newstruct.coordx)/ SiO2Bulk.coordx
            at.y = (at.y * newstruct.coordy)/ SiO2Bulk.coordy
            at.z = at.z - 5.
        at.z = at.z - 5.
    PrintStruct(newstruct, 'castep_inp', name='HfSiO2.cell')
    print 'Structure has', len(newstruct.atoms), 'atoms'
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
