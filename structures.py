from operator import attrgetter
        
class Atom:
    """
    basically is an object that defines an atom by it ABSOLUTE coordinates, species and (optional) label
    """
    def __init__(self, species, x, y, z) :
        self.species, self.x, self.y, self.z = species, x, y, z
    def Label(self, text) :
        self.label = text
    def Charge(self, val) :
        self.charge = val

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
        else :
            # Print and error message for unimplemented stuff and exit
            ERR = 'Can\'t do ( Angles = '+ str(coordstyle)+' and bp = '+ \
                  str(pb) + ')'
            exit()
        self.atoms = list(at_list)
    
    def charge(self):
        return sum([atom.charge for atom in self.atoms])

class StructIn:
    """
    reads in a slab structure, gets its x and y dimensions 
    also creates the atom list using the structure Atom
    """
    def __init__(self, path) :
        f = open(path,'r').readlines()
        first_lines = [x.split() for x in f[0:5]]
        self.coordx = float(first_lines[3][0])
        first_lines[3][0] = str(self.coordx)
        self.coordy = float(first_lines[3][1])
        first_lines[3][1] = str(self.coordy)
        no_ats = int(first_lines[-1][0])
        last_lines = f[no_ats+5:]
        self.elems = []
        for line in f[5:no_ats+5] :
            line = [line.split()[0]] + map(lambda x: float(x), line.replace('E', 'e').split()[1:4]) + line.split()[4:]
            for i in range(1,3) :
                if line[i] < 0. :
                    line[i] = 1. + line[i]
                    elem = Atom(line[0], line[1]*self.coordx, line[2]*self.coordy, line[3]) 
            if line[-1] == '$' :
                elem.Label('layer')
            else :
                elem.Label('bulk')
            self.elems.append(elem)

class gener_:
    """
    generates an interface from a structin structure
    """
    def __init__(self, mold, z) :
        self.coordx = mold.coordx
        self.coordy = mold.coordy
        max_ = max(mold.elems, key = lambda p: p.z).z
        self.coordz = z
        self.elems = mold.elems
        #self.delta1 = max([q for q in struct1.elems if q.label != 'bulk'], key = lambda p: p.z).z - min([q for q in struct1.elems if q.label != 'bulk'], key = lambda p: p.z).z 
        #self.delta2 = max([q for q in struct1.elems if q.label != 'layer'], key = lambda p: p.z).z - min([q for q in struct1.elems if q.label != 'layer'], key = lambda p: p.z).z 
    def separate(self, dz) :
        for elem in self.elems :
            if elem.label == 'layer' :
                elem.z = elem.z + dz
        #self.elems = [x for x in self.elems if x.label == 'bulk']
        #self.coordz = delta1 + delta2 + dz
        #self.coordz = 23.
    def expand(self, fracSi, fracSiO, uSi, uSiO ) :
        """ 
        uSiO2 = 8.5126
        uSi = 5.4289
        """ 
        highbulk = max([q for q in struct1.elems if q.label == 'bulk'], key = lambda p: p.z).z 
        lowlayer = min([q for q in struct1.elems if q.label == 'layer'], key = lambda p: p.z).z 
        newelems = []
        for elem in [q for q in self.elems if q.label == 'bulk'] : 
            if elem.z > highbulk - uSi*(fracSi % 1) :
                newelem = Atom(elem.species, elem.x, elem.y, elem.z - uSi*int(1 + fracSi))
                newelem.Label(elem.label)
                newelems.append(newelem)
            for i in range(1, int(fracSi) + 1) :
                newelem = Atom(elem.species, elem.x, elem.y, elem.z - uSi*i)
                newelem.Label(elem.label)
                newelems.append(newelem)
                
        for elem in [q for q in self.elems if q.label == 'layer'] : 
            if elem.z < lowlayer + uSiO*(fracSiO % 1 ):
                newelem = Atom(elem.species, elem.x, elem.y, elem.z + uSiO*int(1 + fracSiO))
                newelem.Label(elem.label)
                newelems.append(newelem)
            for i in range(1, int(fracSiO) + 1) :
                newelem = Atom(elem.species, elem.x, elem.y, elem.z + uSiO*i)
                newelem.Label(elem.label)
                newelems.append(newelem)
        self.elems.extend(newelems)
    def widen(self, wx, wy) :
        newelems = []
        for xx in range(0, wx) :
            for yy in range(0, wy) :
                if (xx != 0 or yy != 0) :
                        for elem in self.elems :
                            newelem = Atom(elem.species, elem.x + xx*self.coordx, elem.y + yy*self.coordy, elem.z)
                            newelem.Label(elem.label)
                            newelems.append(newelem)

        self.elems.extend(newelems)
        self.coordx = wx*self.coordx
        self.coordy = wy*self.coordy
    def check_ratio(self) :
        print float(len([p for p in self.elems if p.label == 'layer' and p.species == '108']))/len([p for p in self.elems if p.label == 'layer' and p.species == '14'])
    def printy(self) :
        print 'Created with surfdo.py\nCRYSTAL\n0  0  0\n1'
        print '%f %f %f 90.0   90.0   90.0'%(self.coordx, self.coordy, self.coordz)
        print '%i'%(len(self.elems))
        for elem in self.elems :
            print '%4s %15.8f %15.8f %15.8f %s'%(elem.species, elem.x/self.coordx, elem.y/self.coordy, elem.z/self.coordz, elem.label)      
            #print elem.species, elem.x/self.coordx, elem.y/self.coordy, elem.z/self.coordz     
        print 'EXTPRT\nCOORPRT\nSTOP'
    def printlmp(self) :
        f2 = open('data.cr2lmp', 'w')
        f2.write('Created with cr2lmp.py\n\n')
        atom_types = []
        for e in self.elems :
            if e.species not in atom_types :
                atom_types.append(e.species)
        f2.write('%17i atoms\n \t\t0 bonds\n\t\t0 angles\n\t\t0 dihedrals\n\t\t0 impropers\n\n'%(len(self.elems)))
        f2.write('%17i atom types\n \t\t0 bond types\n\t\t0 angle types\n\t\t0 dihedral types\n\t\t0 improper types\n\n'%(len(atom_types)))
        f2.write('\t%18.8f\t%18.8f\txlo xhi\n'%(0.0, self.coordx))
        f2.write('\t%18.8f\t%18.8f\tylo yhi\n'%(0.0, self.coordy))
        f2.write('\t%18.8f\t%18.8f\tzlo zhi\n'%(0.0, self.coordz))
        for elem in self.elems :
            if elem.label == 'bulk' :
                elem.Charge(0.0)
            elif elem.species == '14' :
                elem.Charge(2.0)
            else :
                elem.Charge(-1.0)
                
        f2.write('\nAtoms\n\n')
        for i in range(len(self.elems)) :
            #f2.write('%6i\t%6i\t %15.8f\t %15.8f\t %15.8f\t\n'%(i+1, atom_types.index(self.elems[i].species) + 1, self.elems[i].x, self.elems[i].y, self.elems[i].z))
            f2.write('%6i\t%6i\t%7.4f %15.8f\t %15.8f\t %15.8f\t\n'%(i+1, atom_types.index(self.elems[i].species) + 1, self.elems[i].charge, self.elems[i].x, self.elems[i].y, self.elems[i].z))
        """
        for elem in elems:
            f2.write('%6i\t%6i\t %15.8f\t %15.8f\t %15.8f\t %s\n'%(elem[0], elem[1], elem[2], elem[3], elem[4], elem[-1]))
        """
        f2.close()
    def check_charge(self) :
        print sum([p.charge for p in self.elems])


