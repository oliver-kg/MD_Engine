from vpython import vector
import math

atom_type = []
atom_pos = []
atom_index = []
atom_charge = []

#bond values
a_a = []
a_b = []
r_0 = []
k_engine = []

# bond angles:
a_i = []
a_j = []
a_k = []
b_angle = []
k_ang = []

# dihedrals
d_i = []
d_j = []
d_k = []
d_l = []
ph = []
k_dih = []
n = []


#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/1,6-HEXANEDIAMINE/1,6-HEXANEDIAMINE_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Butanol/butanol_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/C10H23NO6/C10H23NO6_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/BigMolecule/BigMolecule_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_PDB.txt
f = open(r"C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_PDB.txt", "r")             # opens the pos file of the molecule

for i in f:
    if "CONECT" in i:
        break

    if ("AUTHOR" in i) or ("TITLE" in i) or ("HEADER" in i):
        continue

    atom_index.append(int(i[6:11].strip())-1)                                                   # gets the index of the atom
    atom_type.append(i[76:78].strip())                                                          # removes all digits from string, giving the atom type
    x_pos = float(i[30:38].strip())                                                             # get x,y,z positions
    y_pos = float(i[39:46].strip())
    z_pos = float(i[47:54].strip())
    atom_pos.append(vector(x_pos/10, y_pos/10, z_pos/10))                                                # put position vector into array

f.close()

#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/1,6-HEXANEDIAMINE/1,6-HEXANEDIAMINE_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Butanol/butanol_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/C10H23NO6/C10H23NO6_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/BigMolecule/BigMolecule_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_ITP.txt
f = open(r"C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_ITP.txt")


reading_atoms = False
reading_bonds = False
reading_angles = False
reading_dihedrals = False

for i in f:
    if "[ atoms ]" in i:                                                                 # start reading atoms list
        reading_atoms = True
        continue
    
    if "[ bonds ]" in i:                                                                # start reading bonds from the bonds list       
        reading_bonds = True
        reading_atoms = False
        continue

    if "[ angles ]" in i:
        reading_bonds = False
        reading_angles = True
        continue
    
    if "[ dihedrals ]" in i:
        reading_angles = False
        reading_dihedrals = True
        continue
    
    if reading_atoms and "; total charge of the molecule:" in i:                        # skips the last line of the atoms list, and stops
        reading_atoms = False
        continue

    if reading_bonds and "[ pairs ]" in i:                                              # skips the last line of the bonds list, and stops
        reading_bonds = False
        continue

    if reading_angles and "[ dihedrals ]" in i:
        reading_angles = False
        continue
    
    if reading_dihedrals and "[ exclusions ]" in i:
        reading_dihedrals = False
        continue
    
    if (reading_atoms or reading_bonds or reading_angles or reading_dihedrals) and ";" in i:                 # skips the first line
        continue

    # ----------------------------------------
    # ATOMS
    # ----------------------------------------

    if reading_atoms:
       atom_charge.append(i[36:45].strip())                                             # gets atom charge

    # ----------------------------------------
    # BONDS
    # ----------------------------------------

    if reading_bonds:
        a_a.append(int(i[1:5].strip())-1)                                               # gets atom positions in a bond
        a_b.append(int(i[6:10].strip())-1)
        
        r = float(i[17:24].strip())                                                     # r value from file
        r_0.append(r)                                                      

        k_file = float(i[27:37].strip())                                                # k value from file

        funct = float(i[14:15].strip())                                                 # gets the function number for the bond

        if funct == 1:
            k_engine.append(k_file / 2.0)                                               # for the GROMACS harmonic bond

        elif funct == 2:
            k_engine.append(k_file * r**2)                                              # for the GROMOS-96 fourth-power bond

        else:
            raise ValueError("Unsupported bond function type ", funct)

    # ----------------------------------------
    # ANGLES
    # ----------------------------------------

    if reading_angles:
        a_i.append(int(i[1:5].strip())-1)                                               # get the three molecules involved in the angle
        a_j.append(int(i[6:10].strip())-1)
        a_k.append(int(i[11:15].strip())-1)

        theta0_deg = float(i[23:30].strip())
        b_angle.append(theta0_deg)                                         # get ideal angle
        theta0_rad = math.radians(theta0_deg)

        k_file = float(i[32:39].strip())                                                # get k value

        funct = float(i[19:20].strip())                                                 # gets the function number for the bond

        if funct == 1:
            k_ang.append(k_file)                                                        # for the normal GROMACS harmonic angle

        elif funct == 2:
            k_ang.append(k_file * math.sin(theta0_rad)**2)                                              # for the GROMOS-96 fourth-power bond

        else:
            raise ValueError("Unsupported bond function type ", funct)

    # ----------------------------------------
    # DIHEDRALS
    # ----------------------------------------
        
    if reading_dihedrals:
        d_i.append(int(i[1:5].strip())-1)                                               # get the four molecules involved in the dihedral
        d_j.append(int(i[6:10].strip())-1)
        d_k.append(int(i[11:15].strip())-1)
        d_l.append(int(i[16:20].strip())-1)
        ph.append(float(i[28:35].strip()))                                              # get phi angle
        k_dih.append(float(i[38:44].strip()))                                           # get k val
        n.append(int(i[47:49].strip()))                                               # get n
