from operator import attrgetter
import numpy as np
        
def triclinic2xyz(at, struct):
    """
    http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_cartesian_coordinates
    """
    alpha, beta, gamma = np.radians(struct.alpha), np.radians(struct.beta), \
                         np.radians(struct.gamma) 
    abc = np.matrix([[at.x], [at.y], [at.z]])
    a, b, c = struct.coordx, struct.coordy, struct.coordz
    v = np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 +
                2*np.cos(gamma)*np.cos(beta)*np.cos(gamma))
    # transformation matrix as shown on the reference
    transf_mat = np.matrix([
    [a, b*np.cos(gamma), c*np.cos(beta)], 
    [0, b*np.sin(gamma), c*(np.cos(alpha) - (np.cos(beta)*np.cos(gamma)))/np.sin(gamma)], 
    [0, 0, c*(v/np.sin(gamma))]])
    xyz = np.dot(transf_mat, abc)
    return tuple(map(float, xyz))

def xyz2triclinic(at, struct):
    """
    http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_from_cartesian_coordinates
    """
    alpha, beta, gamma = np.radians(struct.alpha), np.radians(struct.beta), \
                         np.radians(struct.gamma) 
    xyz = np.matrix([[at.x], [at.y], [at.z]])
    a, b, c = struct.coordx, struct.coordy, struct.coordz
    v = np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 +
                2*np.cos(gamma)*np.cos(beta)*np.cos(gamma))
    # transformation matrix as shown on the reference
    transf_mat = np.matrix([
    [ 1./a,
     -1*(np.cos(gamma)/(a*np.sin(gamma))), 
     (np.cos(alpha)*np.cos(gamma) - np.cos(beta))/(a*v*np.sin(gamma))
    ], 
    [0,    
     1./(b*np.sin(gamma)), 
     (np.cos(beta)*np.cos(gamma)-np.cos(alpha))/(b*v*np.sin(gamma))
    ], 
    [0,    
     0,
     np.sin(gamma)/(c*v)     
    ]])
    abc = np.dot(transf_mat, xyz)
    return tuple(map(float, abc))

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
    def Addtag(self,tag):
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

    def __init__(self, at_list, coords, coordstyle='Angles', pb='bulk',
                 frac=False) :
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
        self.periodicity = pb
        self.atoms = list(at_list)
        
        
        # if coordinates are fractional compute compute the Cartesian values
        if frac:
            for at in self.atoms:
                if at.x < 0. :
                    at.x = 1. + at.x
                if at.y < 0. :
                    at.y = 1. + at.y
                if at.z < 0. :
                    at.z = 1. + at.z
                at.x, at.y, at.z = triclinic2xyz(at, self)
    
    def charge(self):
        return sum([atom.charge for atom in self.atoms])

    def volume(self):
        """Cell is considered triclic
        http://webmineral.com/help/CellDimensions.shtml#.VNpIhjasWBY
        """
        _alpha = np.radians(self.alpha)
        _beta  = np.radians(self.beta)
        _gamma = np.radians(self.gamma)
        return self.coordx*self.coordy*self.coordz*np.sqrt(
               2*np.cos(_alpha)*np.cos(_beta)*np.cos(_gamma) -\
               np.cos(_alpha)**2 - np.cos(_beta)**2 - np.cos(_gamma)**2 + 1)
    
    def ClearTags(self):
        for at in self.atoms:
            at.tags = []

    def normalise(self):
        """ Removes negative coordinates, all coordinates must fit box
        (AtomStruct) -> None
        NB: Not tested for triclinic cells
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

