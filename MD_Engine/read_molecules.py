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

# LJ force-field information
lj_c6 = []
lj_c12 = []

# ff_nonbonded information
ff_atom_types = []



#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/1,6-HEXANEDIAMINE/1,6-HEXANEDIAMINE_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Butanol/butanol_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/C10H23NO6/C10H23NO6_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/BigMolecule/BigMolecule_PDB.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_PDB.txt
PDB_file = open(r"C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_PDB.txt", "r")             # opens the pos file of the molecule

for i in PDB_file:
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

PDB_file.close()

#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/1,6-HEXANEDIAMINE/1,6-HEXANEDIAMINE_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Butanol/butanol_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/C10H23NO6/C10H23NO6_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/BigMolecule/BigMolecule_ITP.txt
#C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Water/Water_ITP.txt
ITP_file = open(r"C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/Cholestane/Cholestane_ITP.txt")


reading_atoms = False
reading_bonds = False
reading_angles = False
reading_dihedrals = False

for i in ITP_file:
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

       ff_atom_types.append(i[6:11].replace(" ", ""))                              # strip only the white spaces


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

ITP_file.close()

# -------------------------------------------------------------------------------------
# FORCE FIELD STUFF FOR INDIVIDUAL LJ_SIGMA and LJ_EPSILON
# -------------------------------------------------------------------------------------

ff_file = open(r"C:/Dev/projects/MD_Engine/MD_Engine/Molecule_Files/ffnonbonded_ITP.txt")

used_ff_types = set(ff_atom_types)

required_lj_pairs = set()

for type_i in used_ff_types:
    for type_j in used_ff_types:

        pair = tuple(sorted((type_i, type_j)))

        required_lj_pairs.add(pair)


ff_atomtypes = {}
ff_nonbond_params = {}

reading_atomtypes = False
reading_nonbond_params = False

with ff_file as f:

    for i, line in enumerate(f):

        line = line.strip()

        if not line or line.startswith(";"):
            continue

        # atom types
        if line == "[ atomtypes ]":
            reading_atomtypes = True
            continue

        if line.startswith("[") and line != "[ atomtypes ]":
            reading_atomtypes = False

        # nonbonded params
        if line == "[ nonbond_params ]":
            reading_nonbond_params = True
            continue

        if line.startswith("[") and line != "[ nonbond_params ]":
            reading_nonbond_params = False

        # atom types
        if reading_atomtypes:

            fields = line.split()

            if len(fields) < 7:
                continue
            if fields[0] not in used_ff_types:
                continue



            C6 = float(fields[5])
            C12 = float(fields[6])

            ff_atomtypes[fields[0]] = (C6, C12)

        # nonbonded params
        if reading_nonbond_params:
        
            fields = line.split()

            if len(fields) < 5:
                continue

            type_i = fields[0]
            type_j = fields[1]

            choice = tuple(sorted((type_i, type_j)))


            if choice not in required_lj_pairs:
                continue

            C6 = float(fields[3])
            C12 = float(fields[4])

            pair = tuple(sorted((type_i, type_j)))

            ff_nonbond_params[pair] = (C6, C12)

        
resolved_lj_params = {}

for pair in required_lj_pairs:

    if pair in ff_nonbond_params:
        resolved_lj_params[pair] = ff_nonbond_params[pair]
        source = "nonbond_params"

    else:
        type_i, type_j = pair

        if type_i == type_j:
            if type_i in ff_atomtypes:
                resolved_lj_params[pair] = ff_atomtypes[type_i]
                source = "atomtypes"
            else:
                raise ValueError(
                    f"No atomtype parameters for {type_i}"
                )

        else:
            source = "NOT RESOLVED"