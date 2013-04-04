#!/usr/bin/env python

from structures import *

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
    #struct2.check_ratio()
    #struct2.printy()
    #struct2.printlmp()
    #struct2.check_charge()struct1 = StructIn('interface5')
    print AtomStruct.__doc__
    print 'Hello World!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__" :
    main()
