#!/usr/bin/env python

from structures import *
from ReadWrite import *
from operations import *

def main() :
    #struct1.elems = list(set([x for x in struct1.elems if x.z > 2.] + [x for x in struct1.elems if x.label == 'layer']))
    ## removing lowerst Si atoms => ability to create thicker layers
    #for elem in struct1.elems :
    #   if elem.label == 'bulk' and (elem.z < 4.) :
    #       struct1.elems.remove(elem)
    #
    #struct2 = gener_(struct1, 28.5)
    #struct2.separate(11.2)
    ## x.45 expansion of SiO2 yields good stoicheometry
    #struct2.expand(.8, 1., 5.4289, 8.5126)
    #struct2.widen(2, 2)
    #struct2.printlmp()
    #struct2.check_charge()struct1 = StructIn('interface5')

    Str1 = ReadStruct('CrystalCell', style='crystal')
    # Setting the cell charge
    for atom in Str1.atoms :
        if atom.species == 'Si':
            atom.Charge(2.0)
        elif atom.species == 'O':
            atom.Charge(-1.0)
    #print Str1.charge()
    #PrintStruct(Str1, 'lmp_data', name='data.astools')
    Str1 = widen(Str1, 3, 3, 3) 
    PrintStruct(Str1, 'crystal_inp', name='INPUT_astools')
    #PrintStruct(Str2, 'lmp_data', name='data.astools')
    print 'ze end!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__" :
    main()
