#!/usr/bin/env python

from structures import *
from ReadWrite import *
from operations import *

def main() :
    SiO2Bulk = ReadStruct('CrystalCell', style='crystal')
    SiO2Bulk.normalise()
    #SiO2Bulk = repeat(SiO2Bulk, 3, 3, 3) 
    #PrintStruct(SiO2Bulk, 'lmp_data', name='data.SiO2_veloc')
    PrintStruct(SiO2Bulk, 'crystal_inp', name='INPUT_astools')

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
    #PrintStruct(SiBulk, 'crystal_inp', name='INPUT_SiBulk')
    #newSi = expand(SiBulk, X=(-0.5, 0.5), Y=(-0.5, 0.5), Z=(-1.5, 1.3))
    #PrintStruct(newSi, 'crystal_inp', name='INPUT_newSi')

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
    SiO2Bulk = expand(SiO2Bulk, Z = (0., 1.))
    # Setting the cell charge
    for atom in SiO2Bulk.atoms :
        atom.tags.append('oxide')
        if atom.species == 'Si':
            atom.Charge(2.0)
        elif atom.species == 'O':
            atom.Charge(-1.0)
    
    
    SiBulk = expand(SiBulk, Z = (-1., 1.))
    new_struct = AtomStruct([x for x in SiBulk.atoms+SiO2Bulk.atoms],
                            (SiBulk.coordx, SiBulk.coordy, 
                             SiBulk.alpha), pb='slab')
    for at in new_struct.atoms:
        if 'oxide' in at.tags:
            at.z = at.z + SiBulk.coordz + 0.5
            

    for val in [distance(at_a, at_b) < .5 for at_a in new_struct.atoms \
                for at_b in new_struct.atoms if at_a != at_b] :
        if val == True :
            print 'print WARNING', val
    new_struct.atoms = set(new_struct.atoms)
    print 'Total number of atoms in structure', len(new_struct.atoms)
    #PrintStruct(new_struct, 'crystal_inp', name='INPUT_newstruct')
    PrintStruct(new_struct, 'lmp_data', name='data.SiO2_veloc')

    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__" :
    main()
