#!/usr/bin/env python

import numpy as np
import random as rn
import itertools as it
from copy import deepcopy
from math import sqrt, sin, cos, pi

from structures import *
from ReadWrite import *
from operations import *
from analysis import *

def main() :
    cx, cy, cz = 13., 13., 13.
   
    bl = 1.9770
    a = 4* bl / sqrt(6)

    neutSi = Atom('Si', cx/2., cy/2., 2*bl + bl/3.) 
    neutSi.Charge(0.0)
    chargedSi = Atom('Si', cx/2., cy/2., cz) 
    chargedSi.Charge(2.0)
    Simolec = Atom('Si', a/2., a*cos(pi/6.)/3, bl/3.) 
    Simolec.Charge(2.0)
    Simolec.tags = ['oxide']
    O1 = Atom('O', 0.0, 0.0, 0.0) 
    O1.Charge(-1.0)
    O1.tags = ['oxide']
    O2 = Atom('O', a, 0.0, 0.0) 
    O2.Charge(-1.0)
    O2.tags = ['oxide']
    O3 = Atom('O', a/2., a*cos(pi/6.),  0.0) 
    O3.Charge(-1.0)
    O3.tags = ['oxide']
    O4 = Atom('O', a/2., a*cos(pi/6.)/3, bl/3. + bl) 
    O4.Charge(-1.0)
    O4.tags = ['oxide']
   
    S = AtomStruct([neutSi, chargedSi, Simolec, O1, O2, O3, O4],
                   (cx, cy, cz, 90., 90., 90.))
    Sb = deepcopy(S)
    at_ref = deepcopy(Simolec)
    
    for at in S.atoms:
        if 'oxide' in at.tags:
            at.x = at.x - at_ref.x
            at.y = at.y - at_ref.y + cy/2.
            #at.z = at.z - at_ref.z + 1.5


    for x, y in it.product(np.arange(0., 
                                     cx, .1),
                           [0.00]): 
        for i in range(2, len(S.atoms)):
            S.atoms[i].x = Sb.atoms[i].x + x
        PrintStruct(S, 'lmp_data', 
                    name='data.SiO2_offset_{:.2f}_{:.2f}'.format(x,y),
                    nocharge=False)
        print x, y
    PrintStruct(S, 'crystal_inp', name='INPUT_test')
    print 'Done!'


    # make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
