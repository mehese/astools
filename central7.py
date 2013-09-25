#!/usr/bin/env python

import numpy as np
import random as rn
import itertools as it

from structures import *
from ReadWrite import *
from operations import *
from analysis import *

def main() :
    HfO2 = ReadStruct('HfO2_out', style='crystal_out')
    #HfO2 = slab_hfo2(HfO2, 2, 2, 2.0)
    #HfO2.normalise()
    for at in HfO2.atoms:
        at.tags.append('hfo2')
    PrintStruct(HfO2, 'castep_inp', name='HfO2.cell')
    print 'Structure has', len(HfO2.atoms), 'atoms'
    print 'Done!'


# make sure the codelines don't get executed if we insert this as a module
if __name__ == "__main__":
    main()
