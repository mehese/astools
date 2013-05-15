#!/usr/bin/env python

from structures import *
from ReadWrite import *
from operations import *

def main() :
    Str1 = ReadStruct('CrystalCell', style='crystal')
    # Setting the cell charge
    for atom in Str1.atoms :
        if atom.species == 'Si':
            atom.Charge(2.0)
        elif atom.species == 'O':
            atom.Charge(-1.0)
    #print Str1.charge()
    Str1 = repeat(Str1, 3, 3, 3) 
    #PrintStruct(Str1, 'crystal_inp', name='INPUT_astools')
    PrintStruct(Str1, 'lmp_data', name='data.SiO2_veloc')

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
    expand(SiBulk, X=(-2.1, 2.1))
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__" :
    main()
