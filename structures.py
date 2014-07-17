from operator import attrgetter
        
class Atom:
    """
    basically is an object that defines an atom by it ABSOLUTE coordinates, species and (optional) label
    """
    def __init__(self, species, x, y, z) :
        self.species, self.x, self.y, self.z = species, x, y, z
        self.tags = []
    def Label(self, text) :
        self.label = text
    def Charge(self, val) :
        self.charge = val
    def Addtag(tag):
        self.tags.append(tag)
    def __str__(self):
        return '{:2} {:10.5f} {:10.5f} {:10.5f}'.format(self.species,
                self.x, self.y, self.z)

class AtomStruct:
    """ Defines a structure of atoms

    Keyword arguments:

    at_list -- a list of Atom objects
    coords -- a tuple (x_dimension, y_dimension, z_dimension,
                       alpha_angle, beta_angle, gamma_angle)
    coordstyle -- (default='Angles') a string that defines how to read 
        the coords tuple. If coordstyle='vecs' then coors tuple is read
        as (x0, x1, y0, y1, z0, z1)

    Contains:

    self.coord[xyz] = cell cordinates
    alpha, beta, gamma = cell angles
    self.atoms = atoms list

    """

    def __init__(self, at_list, coords, coordstyle='Angles', pb='bulk') :
        if coordstyle == 'Angles' and pb=='bulk':
            (self.coordx, self.coordy, self.coordz, 
             self.alpha, self.beta, self.gamma) = coords
        elif coordstyle == 'Angles' and pb=='slab':
            (self.coordx, self.coordy, 
             self.alpha) = coords
        else :
            # Print and error message for unimplemented stuff and exit
            ERR = 'Can\'t do ( Angles = '+ str(coordstyle)+' and bp = '+ \
                  str(pb) + ')'
            exit()
        self.atoms = list(at_list)
        self.periodicity = pb
    
    def charge(self):
        return sum([atom.charge for atom in self.atoms])
    
    def ClearTags(self):
        for at in self.atoms:
            at.tags = []

    def normalise(self):
        """ Removes negative coordinates, all coordinates must fit box
        (AtomStruct) -> None
        """
        for at in self.atoms:
            if at.x < 0. :
                at.x = self.coordx + at.x
            if at.y < 0. :
                at.y = self.coordy + at.y
            if at.z < 0. :
                at.z = self.coordz + at.z
    def __len__(self):
        return len(self.atoms)

