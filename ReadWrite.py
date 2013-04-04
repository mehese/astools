from structures import *

# Dictionary that changes the Z to a string defining the atom species
# You should probably try completing it
Z2species = {1:'H', 2:'He', 8:'O', 14:'Si', 72:'Hf'}

def ReadStruct(filename, style='crystal'):
    """ Reads and input file and returns an AtomStruct
    (str, str) -> AtomStruct
    """
    if style == 'crystal':
        f = [line.split() for line in open(filename, 'r').readlines()]
        if f[1][0][:2] == 'CR' :
            # Periodic boundaries
            # Coordinates
            cds = tuple(map(lambda x: float(x), f[4]))
            # Initializing atom list
            atoms = []
            # Going through the atom list
            for line in f[6:int(f[5][0])+6] :
                spec = Z2species[int(line[0][-2:])] 
                x = float(line[1])*cds[0] # Atom object stores 
                y = float(line[2])*cds[1] # absolute distances, not
                z = float(line[3])*cds[2] # fractional coordinates
                at = Atom(spec, x, y, z)
                #appending the atom object to the atom list
                atoms.append(at)
            crystal = AtomStruct(atoms, cds, coordstyle='Angles', pb='bulk')
            return crystal
        else:
            print 'Bad File !!!'
            exit()
    else:
    
        print 'Bad File !!!'
        exit()

def PrintStruct(structure, filetype, name='PrintStruct.out'):
    """Prints an AtomStruct to an ASCII file of desired format
    (AtomStruct, str) -> None
    """
    if filetype == 'lmp_data':
        print 'file name :', name
        f2 = open(name, 'w')
        f2.write('Created with cr2lmp.py\n\n')
        atom_types = []
        for e in structure.atoms:
            if e.species not in atom_types:
                atom_types.append(e.species)
        f2.write('%17i atoms\n \t\t0 bonds\n\t\t0 angles\n\t\t0 dihedrals \
                 \n\t\t0 impropers\n\n'%(len(structure.atoms)))
        f2.write('%17i atom types\n \t\t0 bond types\n\t\t0 angle types\n \
                 \t\t0 dihedral types\n\t\t0 improper types\n\n' \
                 %(len(atom_types)))
        f2.write('\t%18.8f\t%18.8f\txlo xhi\n'%(0.0, structure.coordx))
        f2.write('\t%18.8f\t%18.8f\tylo yhi\n'%(0.0, structure.coordy))
        f2.write('\t%18.8f\t%18.8f\tzlo zhi\n'%(0.0, structure.coordz))
        f2.write('\nAtoms\n\n')
        for i in range(len(structure.atoms)) :
            f2.write('%6i\t%6i\t%7.4f %15.8f\t %15.8f\t %15.8f\t\n' \
                     %(i+1, atom_types.index(structure.atoms[i].species) + 1, 
                       structure.atoms[i].charge, structure.atoms[i].x, 
                       structure.atoms[i].y, structure.atoms[i].z))
        f2.close()

