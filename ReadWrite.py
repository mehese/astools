#!/usr/bin/env python

from structures import *

# Dictionary that changes the Z to a string defining the atom species
# You should probably try completing it
Z2species = {1:'H', 2:'He', 8:'O', 14:'Si', 72:'Hf'}
species2Z = {v:k for k, v in Z2species.items()}

def ReadStruct(filename, style='crystal', pos=-1):
    """ Reads and input file and returns an AtomStruct
    (str, str) -> AtomStruct

    Available styles:

     crystal     -- CRYSTAL input files
     crystal_out -- CRYSTAL output files
     lmp_dump    -- LAMMPS dump file

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
    if style == 'crystal_out':
        f = open(filename, 'r').read()
        print 'opened', filename
        c = f.split('CARTESIAN COORDINATES')[-2].split('PRIMITIVE CELL')[-1]
        cds = tuple(map(float, c.split('\n')[2].split()))
        c = f.split('CARTESIAN COORDINATES')[-1].split('\n')[4:50]
        i = 0
        atoms = []
        while c[i].split() != [] :
            l = c[i].split()
            spec = l[2][0] + l[2][1:].lower()
            x, y, z = tuple(map(float, l[3:]))
            at = Atom(spec, x, y, z)
            atoms.append(at)
            i += 1
        crystal = AtomStruct(atoms, cds, coordstyle='Angles', pb='bulk')
        return crystal
    if style == 'lmp_dump':
        f = open(filename, 'r').read()
        # Get the number of atoms at given position
        no_ats = int(f.split('ITEM: BOX BOUNDS')[pos-1].split('ATOMS'
                 )[-1].split()[0])
        dats=[a for a in f.split('ITEM: BOX BOUNDS')[pos].split('\n')[1:] if \
              a != '']
        cx = map(float, dats[0].split())
        cx = cx[1] - cx[0]
        cy = map(float, dats[1].split())
        cy = cy[1] - cy[0]
        cz = map(float, dats[2].split())
        cz = cz[1] - cz[0]
        #print cx, cy, cz
        atoms = []
        # Iterate through the lines with the atoms info
        for line in dats[4:no_ats+4] :
            a = line.split()
            spec = Z2species[int(float(a[1]))/2]
            x, y, z, q = tuple(map(float, a[2:])) 
            at = Atom(spec, x, y, z)
            at.Charge(q)
            atoms.append(at)
        #print len(atoms)
        cds = (cx, cy, cz, 90., 90., 90.)
        crystal = AtomStruct(atoms, cds, coordstyle='Angles', pb='bulk')
        return crystal
    if style == 'castep_inp' :
        f = open(filename, 'r').read()
        cx, cy, cz, ang1, ang2, ang3 = tuple(map(float, 
                                    f.split('LATTICE_ABC')[1].split()[:-1]))
        cds = (cx, cy, cz, ang1, ang2, ang3)
        atoms = []
        for line in f.split('positions_frac' \
                            )[1][1:-len('%endblock  ')].split('\n'):
            spec = line.split()[0]
            x, y, z = tuple(map(float, line.split()[1:]))
            x, y, z = x*cx, y*cy, z*cz
            at = Atom(spec, x, y, z)
            atoms.append(at)
        crystal = AtomStruct(atoms, cds, coordstyle='Angles', pb='bulk')
        return crystal

    else:
        print 'Filetype not recognized !!!'
        exit()

def PrintStruct(structure, filetype, name='PrintStruct.out', nocharge=False):
    """Prints an AtomStruct to an ASCII file of desired format
    (AtomStruct, str) -> None
    """
    if filetype == 'lmp_data':
        f2 = open(name, 'w')
        f2.write('Created with astools\n\n')
        atom_types = []
        for e in structure.atoms:
            if e.species not in atom_types:
                atom_types.append(e.species)
        f2.write('%17i atoms\n\t\t0 bonds\n\t\t0 angles\n\t\t0 dihedrals' \
                 '\n\t\t0 impropers\n\n'%(len(structure.atoms)))
        f2.write('%17i atom types\n\t\t0 bond types\n\t\t0 angle types' \
                 '\n\t\t0 dihedral types\n\t\t0 improper types\n\n' \
                 %(len(atom_types)))
        f2.write('\t%18.8f\t%18.8f\txlo xhi\n'%(0.0, structure.coordx))
        f2.write('\t%18.8f\t%18.8f\tylo yhi\n'%(0.0, structure.coordy))
        f2.write('\t%18.8f\t%18.8f\tzlo zhi\n'%(0.0, structure.coordz))
        f2.write('\nAtoms\n\n')
        if nocharge :
            for i in range(len(structure.atoms)) :
                f2.write('%6i\t%6i %15.8f\t %15.8f\t %15.8f\t\n' \
                         %(i+1, atom_types.index(structure.atoms[i].species) \
                           + 1, structure.atoms[i].x, structure.atoms[i].y, 
                           structure.atoms[i].z))
        else :
            for i in range(len(structure.atoms)) :
                f2.write('%6i\t%6i\t%7.4f %15.8f\t %15.8f\t %15.8f\t\n' \
                         %(i+1, atom_types.index(structure.atoms[i].species) \
                           + 1, structure.atoms[i].charge, 
                           structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z))
        f2.close()
    elif filetype == 'crystal_inp':
        if structure.periodicity == 'bulk':
            f2 = open(name, 'w')
            f2.write('Created with astools\nCRYSTAL\n0 0 0\n1\n')
            f2.write('{:-g} {:-g} {:-g} {:-3.5f} {:-3.5f} {:-3.5f}\n'.format(
                     structure.coordx, structure.coordy, structure.coordz,
                     structure.alpha, structure.beta, structure.gamma))
            f2.write('{}\n'.format(len(structure.atoms)))
            for atom in structure.atoms :
                x_ = atom.x / structure.coordx # CRYSTAL input files require 
                y_ = atom.y / structure.coordy # fractional coordinates
                z_ = atom.z / structure.coordz #
                f2.write('{:3} {: 15.9f} {: 15.9f} {: 15.9f}\n'.format(
                         species2Z[atom.species], x_, y_, z_))
            f2.write('ENDGEOM\nSTOP\n')
            f2.close()
        if structure.periodicity == 'slab':
            f2 = open(name, 'w')
            f2.write('Created with astools\nSLAB\n1\n')
            f2.write('{:-g} {:-g} {:-3.5f}\n'.format(
                     structure.coordx, structure.coordy, structure.alpha))
            f2.write('{}\n'.format(len(structure.atoms)))
            for atom in structure.atoms :
                x_ = atom.x / structure.coordx # CRYSTAL input files require 
                y_ = atom.y / structure.coordy # fractional coordinates
                z_ = atom.z                    # but not for z coordinates in
                                               # the case of slabs
                f2.write('{:3} {: 15.9f} {: 15.9f} {: 15.9f}\n'.format(
                         species2Z[atom.species], x_, y_, z_))
            f2.write('ENDGEOM\nSTOP\n')
            f2.close()
    elif filetype == 'castep_inp':
        if structure.periodicity == 'bulk':
            f2 = open(name, 'w')
            types = set([nm.species for nm in structure.atoms])
            f2.write('%block LATTICE_ABC\n')
            f2.write('{:f} {:f} {:f}\n'.format(structure.coordx,
                                               structure.coordy,
                                               structure.coordz))
            f2.write('{:f} {:f} {:f}\n'.format(structure.alpha,
                                               structure.beta,
                                               structure.gamma))
            f2.write('%endblock LATTICE_ABC\n\n')
            f2.write('%block positions_frac\n')
            for at in structure.atoms:
                f2.write('{}    {:f}    {:f}    {:f}\n'.format(
                at.species.rjust(2), at.x/structure.coordx, 
                at.y/structure.coordy, at.z/structure.coordz))
            f2.write('%endblock positions_frac\n')
            f2.write('\nfix_all_cell : true\n\n')
            f2.write('kpoints_mp_grid 4 4 2\n\n')
            f2.write('SYMMETRY_GENERATE\n\n')
            f2.write('%block species_pot\n')
            for s in types:
                f2.write('{}   {}_00PBE_OP.recpot\n'.format(s, s))
            f2.write('%endblock species_pot\n')
            f2.close()

def main():
    str1 = ReadStruct('dump.SiO2tomelt', style='lmp_dump')
    str2 = ReadStruct('dump.SiO2tomelt', style='lmp_dump', pos=-2)
    str3 = ReadStruct('dump.SiO2tomelt', style='lmp_dump', pos=-3)
    str4 = ReadStruct('dump.SiO2tomelt', style='lmp_dump', pos=-4)
    

    #sio2 = ReadStruct('SiO2Si.cell', style='castep_inp')
    PrintStruct(str1, 'crystal_inp', name='INPUT_1')
    PrintStruct(str2, 'crystal_inp', name='INPUT_2')
    PrintStruct(str3, 'crystal_inp', name='INPUT_3')
    PrintStruct(str4, 'crystal_inp', name='INPUT_4')
    print 'All done!'

if __name__ == "__main__" :
    main()
