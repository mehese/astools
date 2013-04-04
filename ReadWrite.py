from structures import *

def ReadStruct(filename, style='crystal'):
    """ Reads and input file and returns an AtomStruct
        (str, str) -> AtomStruct
    """
    if style == 'crystal' :
        f = [line.split() for line in open(filename, 'r').readlines()]
        print f[1][0]


    return 'Lalala'
