#!/usr/bin/env python

import numpy as np
import random as rn
import itertools as it

from structures import *
from ReadWrite import *
from operations import *
from analysis import *

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

    ###################################################################
    
    ##################### CREATING STRUCTURE ########################## 

    ###################################################################

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
    for atom in SiO2Bulk.atoms :
        atom.z = SiO2Bulk.coordz - atom.z
        atom.x = atom.x*(SiBulk.coordx/SiO2Bulk.coordx) 
        atom.y = atom.y*(SiBulk.coordy/SiO2Bulk.coordy) 
        atom.tags.append('oxide')
        if atom.species == 'Si':
            atom.Charge(2.0)
        elif atom.species == 'O':
            atom.Charge(-1.0)
    
    #PrintStruct(SiO2Bulk, 'crystal_inp', name='INPUT_oxide')
    
    #SiBulk = expand(SiBulk, Z = (-1., 1.))
    new_struct = AtomStruct([x for x in SiBulk.atoms+SiO2Bulk.atoms],
                            (SiBulk.coordx, SiBulk.coordy, 
                             SiBulk.coordz + SiO2Bulk.coordz + 1., 
                             SiBulk.alpha, SiBulk.beta,SiBulk.gamma), 
                             pb='bulk')
    for at in new_struct.atoms:
        if 'oxide' in at.tags:
            at.z = at.z + SiBulk.coordz + 0.1
            

    for val in [distance(at_a, at_b) < .5 for at_a in new_struct.atoms \
                for at_b in new_struct.atoms if at_a != at_b] :
        if val == True :
            print 'print WARNING', val
    #new_struct.atoms = set(new_struct.atoms)
    oldx = new_struct.coordx
    oldy = new_struct.coordy
    new_struct.coordz = 19.
    #new_struct = repeat(new_struct, 3, 3, 1)

    print 'Total number of atoms in structure', len(new_struct.atoms)
    new_str_ats_cp = new_struct.atoms[:]
    #PrintStruct(new_struct, 'crystal_inp', name='INPUT_newstruct')
    j = 1
    k = 0
    #for x, y in it.product(np.arange(0., 
    #                                 oldx/2., 0.05),
    #                       np.arange(0., 
    #                                 oldy/2., 0.05)): 
    for x, y in it.product(np.arange(0., 
                                     oldx/2., 0.3),
                           np.arange(0., 
                                     oldy/2., 0.3)):
        #j += 1
        for i in range(len(new_struct.atoms)):
            if 'oxide' in new_struct.atoms[i].tags :
                new_struct.atoms[i].x = new_str_ats_cp[i].x + x
        str_ = '{:.2f}_{:.2f}_SiO2'.format(x, y).replace('.',
                                                         'p')

        PrintStruct(new_struct, 'castep_inp', 
                    name=str_+'.cell')
        new_struct_lmp = repeat(new_struct, 3, 3, 1)
        #PrintStruct(new_struct_lmp, 'lmp_data', 
        #            name='data.SiO2_offset_{:.2f}_{:.2f}'.format(x, y))
        k += 1
        print k, x, y, str_
    #PrintStruct(new_struct, 'lmp_data', 
    #            name='data.SiO2_test',
    #            nocharge=False)
    PrintStruct(new_struct, 'crystal_inp', name='INPUT_interface')
    #PrintStruct(new_struct, 'castep_inp', name='SiO2Si_topt.cell')
    #new_struct = repeat(new_struct, 3, 3, 1)
    #PrintStruct(new_struct, 'lmp_data', name='data.SiO2Si_topt')
    #r, y = rdf(new_struct, 10)
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
