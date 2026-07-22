from vpython import *                   # imports all I use
import tkinter as tk
import math
import read_molecules
import time

# -- CLASS STRUCTURES --

# dictionary of different elements
ELEMENTS = {
    "H": {"mass": 1.0080,  "radius": 0.03, "color": vector(1, 1, 1)},
    "C": {"mass": 12.0110, "radius": 0.05, "color": vector(0.2, 0.2, 0.2)},
    "O": {"mass": 15.9994, "radius": 0.04, "color": vector(1, 0, 0)},
    "N": {"mass": 14.0067, "radius": 0.047, "color": vector(0, 0, 1)},
    "S": {"mass": 32.0, "radius": 0.15, "color": vector(1, 1, 0)},
    "Cl": {"mass": 35.453, "radius": 0.05, "color": vector(0, 1, 0)}
}

# atom class
class Atom:
    def __init__(self, pos, vel, force, charge, element):
        props = ELEMENTS[element]
        self.mass = props["mass"]
        self.pos = pos
        self.vel = vel
        self.force = force
        self.charge = charge
        self.radius = props["radius"]
        self.element = element
        self.base_colour = props["color"]

# molecule class
class Molecule:
    def __init__(self, atom_indices):
        self.atom_indices = atom_indices

# bonds class
class Bond:
    def __init__(self, a1, a2, ideal_dist, K):
        self.a1 = a1
        self.a2 = a2
        self.ideal_dist = ideal_dist
        self.K = K

#bond angle class
class BondAngle:
    def __init__(self, ideal_ang_deg, theta, first, second, third, k):
        self.ideal_ang = math.radians(ideal_ang_deg)
        self.theta = theta
        self.first = first
        self.second = second
        self.third = third
        self.k = k
    
# torsion angle class
class TorsionAngle:
    def __init__(self, a1, a2, a3, a4, terms):
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.psi = 0.0
        self.terms = []
        
        # create structure for having more than one term for 3-term torsion equations
        for k, n, delta_deg in terms:
            self.terms.append((k, n, math.radians(delta_deg)))

# molecule properties arrays
atoms = []
balls = []
bonds = []
bond_angles = []
torsion_angles = []
bond_visuals = []

molecules = [
     Molecule([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]),
     Molecule([15])
    ]

# simulation pararmeters
pot_e = 0
k_e = 0
total_e = 0
TIME_STEP = 0.0004
LJ_SIGMA = 0.3
LJ_EPSILON = 0.02
COULOMB_CONSTANT = 138.935456
PBC_BOX_LENGTH = 10


time_taken_all_foces = 0
time_taken_exclusions = 0
time_taken_bonds = 0
time_taken_bond_angles = 0
time_taken_torsions = 0
time_taken_non_bonded = 0
time_taken_graphics = 0
time_taken_total = 0
time_other = 0


# -- SETUP --

# Get screen dimensions
root = tk.Tk()
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
root.destroy()

# Create a canvas that fills the screen 
scene = canvas(width=screen_width,
                height=screen_height)    

# draws a faint box of the simulation area
L = PBC_BOX_LENGTH
box_visual = box(
    pos=vector(0,0,0),
    size=vector(L, L, L),
    opacity=0.02,        # faint
    color=vector(1,1,1)
)

# Lighting
scene.lights = []      # Remove all default lights
scene.ambient = color.gray(0.22)   # Soft fill light
distant_light(
    direction=vector(-1,-1,-1),
    color=color.gray(0.35)
)

distant_light(
    direction=vector(1,0.5,1),
    color=color.gray(0.15)
)

# keys for camera movement
keys = {"w": False, "a": False, "s": False, "d": False, "q": False, "e": False}

def keydown(evt):
    key = evt.key.lower()
    if key in keys:
        keys[key] = True

def keyup(evt):
    key = evt.key.lower()
    if key in keys:
        keys[key] = False

def render_initial_model():
    # create balls
    for a in atoms:
        b = sphere(pos=a.pos, radius=a.radius, color=a.base_colour, make_trail=False)
        balls.append(b)

    # Draw lines between bonds
    for bond in bonds:
        c = cylinder(
            pos=atoms[bond.a1].pos,
            axis=atoms[bond.a2].pos - atoms[bond.a1].pos,
            radius=0.01,
            color=vector(1,1,1)
        )
        bond_visuals.append(c)

def build_molecule_from_file():
    # create atoms from the read_molecules python file
    for ind in read_molecules.atom_index:
        atoms.append(Atom(read_molecules.atom_pos[ind], vector(0,0,0), vector(0.1,0.1,0.1), read_molecules.atom_charge[ind], read_molecules.atom_type[ind]))

    # create bonds from the read_molecules python file    
    for ind in range(len(read_molecules.a_a)):    
        bonds.append(Bond(read_molecules.a_a[ind], read_molecules.a_b[ind], read_molecules.r_0[ind], read_molecules.k_engine[ind]))

    # create angles from the read_molecules python file
    for ind in range(len(read_molecules.a_i)):
        bond_angles.append(BondAngle(read_molecules.b_angle[ind], 0, read_molecules.a_i[ind], read_molecules.a_j[ind], read_molecules.a_k[ind], read_molecules.k_ang[ind]))

    # create dihedrals from the read_molecules python file
    for ind in range(len(read_molecules.d_i)):
        torsion_angles.append(TorsionAngle(read_molecules.d_i[ind], read_molecules.d_j[ind], read_molecules.d_k[ind], read_molecules.d_l[ind], [(read_molecules.k_dih[ind], read_molecules.n[ind], read_molecules.ph[ind])]))

#create test molecules

build_molecule_from_file()
atoms.append(Atom(vector(-0.51,0.51,0.81), vector(0,0,0), vector(0.1,0.1,0.1), 0.005, "H")),

# create number of atoms after all atoms have been made
no_balls = len(atoms)

# Lennard Jones Force Equation (VDW's)
def calc_LJ_force(dist, direction):
    return  -((24*LJ_EPSILON*((2*((LJ_SIGMA/dist)**12)) - ((LJ_SIGMA/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def calc_KE():
    k_e = 0
    for i in range(no_balls): 
        k_e = k_e + (0.5*atoms[i].mass*(dot(atoms[i].vel, atoms[i].vel))) # use dot product as v^2 to calc ke
    return k_e

# calc Coulombs
def calc_coulombs(q1, q2, r, direction, k):
    return -(k*((q1*q2)/(r**2))) * direction

# builds the neibour list structure
def build_neighbours():
    neighbours = {i: set() for i in range(no_balls)}        # set makes sure no duplicates and fast to check
    for bond in bonds:
        neighbours[bond.a1].add(bond.a2)                    # bonding works both ways so adds the bond to the specific atom you are looking at at the moment
        neighbours[bond.a2].add(bond.a1)

    return neighbours

# filters out all the bonds into 1-2, 1-3 or 1-4 bonds
def build_exclusion_sets():
    exclusions_start = time.perf_counter()

    neighbours = build_neighbours()

    bonded_12 = set()                                       # get ready to create set of bonds which are 1-2, 1-3 and 1-4 bonded
    bonded_13 = set()
    bonded_14 = set()

    # for 1-2 pairs
    for i in neighbours:                                    
        for j in neighbours[i]:                             # if j is in the neighbour list at pos i, then there is a bond between i and j
            bonded_12.add(tuple(sorted((i, j))))             # tuple and sorted means that each pair only gets added once, so no duplicates

    # for 1-3 pairs
    for i in neighbours:                                    # same thing as 1-2 bonds
        for j in neighbours[i]:
            for k in neighbours[j]:                         
                if k != i:                                  # cant loop back to the same atom, so skips it
                    bonded_13.add(tuple(sorted((i, k))))
    
    # for 1-4 pairs
    for i in neighbours:
        for j in neighbours[i]:
            for k in neighbours[j]:
                if k != i:
                    for l in neighbours[k]:
                        if l != j and l != i:
                            pair = tuple(sorted((i, l)))

                            if pair not in bonded_12 and pair not in bonded_13: # makes sure to skip over already bonded 1-2 and 1-3 pairs
                                bonded_14.add(pair)

    #calculate time taken
    exclusions_end = time.perf_counter()
    global time_taken_exclusions 
    time_taken_exclusions = exclusions_end - exclusions_start

    return bonded_12, bonded_13, bonded_14

bonded_12, bonded_13, bonded_14 = build_exclusion_sets()    # builds the sets to be used

# checks if a molecule has reached the box boundery and needs warping
def wrap_molecules():
    half = PBC_BOX_LENGTH / 2

    for mol in molecules:
        ref_pos = atoms[mol.atom_indices[0]].pos
        shift = vector(0, 0, 0)

        if ref_pos.x > half:
            shift.x -= PBC_BOX_LENGTH
        elif ref_pos.x < -half:
            shift.x += PBC_BOX_LENGTH

        if ref_pos.y > half:
            shift.y -= PBC_BOX_LENGTH
        elif ref_pos.y < -half:
            shift.y += PBC_BOX_LENGTH

        if ref_pos.z > half:
            shift.z -= PBC_BOX_LENGTH
        elif ref_pos.z < -half:
            shift.z += PBC_BOX_LENGTH

        if shift.x != 0 or shift.y != 0 or shift.z != 0:
            for i in mol.atom_indices:
                atoms[i].pos += shift

# minimum image periodic correction - turns the vector into the nearest image vector
def PBC_box_for_vectors(r_vector):                                  
    if r_vector.x > PBC_BOX_LENGTH/2:
        r_vector.x = r_vector.x - PBC_BOX_LENGTH
    if r_vector.x < -PBC_BOX_LENGTH/2:
        r_vector.x = r_vector.x + PBC_BOX_LENGTH

    if r_vector.y > PBC_BOX_LENGTH/2:
        r_vector.y = r_vector.y - PBC_BOX_LENGTH
    if r_vector.y < -PBC_BOX_LENGTH/2:
        r_vector.y = r_vector.y + PBC_BOX_LENGTH

    if r_vector.z > PBC_BOX_LENGTH/2:
        r_vector.z = r_vector.z - PBC_BOX_LENGTH
    if r_vector.z < -PBC_BOX_LENGTH/2:
        r_vector.z = r_vector.z + PBC_BOX_LENGTH

# simple camera movement via key press
def update_camera(dt):
    speed = 20 

    forward = norm(scene.forward)
    right = norm(cross(forward, scene.up))
    up = scene.up

    move = vector(0, 0, 0)

    if keys["w"]:
        move += forward
    if keys["s"]:
        move -= forward
    if keys["a"]:
        move -= right
    if keys["d"]:
        move += right
    if keys["e"]:
        move += up
    if keys["q"]:
        move -= up

    if mag(move) > 0:
        move = norm(move) * speed * dt
        scene.camera.pos += move

# caculate physics stuff - the heart <3
def calc_physics():
    pot_e_total = 0.0
    start_all_foces = time.perf_counter()

    
    # bond forces
    start_bonds = time.perf_counter()

    for bond in bonds:
        a = bond.a1
        b = bond.a2
        r0 = bond.ideal_dist
        k = bond.K

        # compute the bond vector and bond lengths
        r_vec = atoms[b].pos - atoms[a].pos
        PBC_box_for_vectors(r_vec)
        r = mag(r_vec)

        if r == 0:                              # avoid devide by zero just in case
            continue

        r_hat = norm(r_vec)                     # bond direction

        F = 2 * k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds

        atoms[a].force += F                          # apply Newtons third law - equal and opposite forces
        atoms[b].force -= F

        pot_e_total += k * (r - r0)**2                # calc the harmonic bond potential in function 2 type bond

    end_bonds = time.perf_counter()
    global time_taken_bonds
    time_taken_bonds = end_bonds - start_bonds


    # bond angles
    start_bond_angles = time.perf_counter()

    for angle in bond_angles:

        theta = angle.theta
        ideal_ang = angle.ideal_ang
        atomA = angle.first
        atomB = angle.second
        atomC = angle.third
        k = angle.k

        BA = atoms[atomA].pos - atoms[atomB].pos                                # find vector of B to A
        BC = atoms[atomC].pos - atoms[atomB].pos                                # find vector of B to C

        PBC_box_for_vectors(BA)                                                 # update the coppy vectors
        PBC_box_for_vectors(BC)
        
        r_BA = mag(BA)                                                          # get the lengths, as force depends on direction, and how long the arms are
        r_BC = mag(BC)
        
        if r_BA < 1e-12 or r_BC < 1e-12:
            continue
       
        theta_cos = (dot(BA, BC))/(r_BA*r_BC)                                           # calculate the dot product hrer
        theta_cos = max(-1.0, min(1.0, theta_cos))                                      # clamp for safety
        theta = math.acos(theta_cos)                                                    
        angle.theta = theta                                                             # update theta

        sin_theta = math.sin(theta)                                                     # to not divide by zero
        if abs(sin_theta) < 1e-8:
            continue

        dU_dtheta = k * (theta - ideal_ang)                                             # the error of the current angle vs ideal bond angle

        f_a = -(dU_dtheta / sin_theta) * ((theta_cos / (r_BA * r_BA)) * BA - (1.0 / (r_BA * r_BC)) * BC) # build the part of the force on A that changes the angle without streaching BA

        f_c = -(dU_dtheta / sin_theta) * ((theta_cos / (r_BC * r_BC)) * BC - (1.0 / (r_BA * r_BC)) * BA) # build the part of the force on C that changes the angle without streaching BC

        f_b = -(f_a + f_c)                                                      # force to be applied on b to conserve energy
    
        atoms[atomA].force += f_a                                               # update the forces
        atoms[atomB].force += f_b
        atoms[atomC].force += f_c

        pot_e_total += 0.5 * k * (theta - ideal_ang)**2                # update the potential energy change from the angle

    end_bond_angles = time.perf_counter()
    global time_taken_bond_angles
    time_taken_bond_angles = end_bond_angles - start_bond_angles


    # torsion angles
    start_torsions = time.perf_counter()
    
    for angle in torsion_angles:
        a1 = angle.a1                           # 4 atoms involved in the torsion bond
        a2 = angle.a2
        a3 = angle.a3
        a4 = angle.a4


        b1 = atoms[a2].pos - atoms[a1].pos      # directions of the three bonds
        b2 = atoms[a3].pos - atoms[a2].pos
        b3 = atoms[a4].pos - atoms[a3].pos

        PBC_box_for_vectors(b1)                 # update the vectors of coppies
        PBC_box_for_vectors(b2)
        PBC_box_for_vectors(b3)

        n1 = cross(b1, b2)                      # normal to plane ABC
        n2 = cross(b2, b3)                      # normal to plane BCD

        eps = 1e-12
        n1_sq = dot(n1, n1)
        n2_sq = dot(n2, n2)
        b2_sq = dot(b2, b2)
        b2_mag = mag(b2)

        if n1_sq < eps or n2_sq < eps or b2_sq < eps:    # makes sure small values dont blow the system up
            continue
        
        x = dot(n1, n2)
        y = dot(cross(n1, n2), norm(b2))

        psi = math.atan2(y, x)                  # caluclate current psi angle
        angle.psi = psi 

        torsion_e = 0
        dV_dpsi = 0
        
        for k, n, delta in angle.terms:                             # repeats for each term, updating the cos graph
            torsion_e += k * (1 + math.cos(n * psi - delta))   # calculate the potential energy
            dV_dpsi -= k * n * math.sin(n * psi - delta)            # how strongly the torsion wants to rotate
        
            
        fa_pref = dV_dpsi * (b2_mag / n1_sq)           # calculates the geometric scallings of the force
        fd_pref = -dV_dpsi * (b2_mag / n2_sq)

        f_a = fa_pref * n1                              # aligning the forces with the direction to the plane
        f_d = fd_pref * n2

        c1 = dot(b1, b2) / b2_sq                        # calculates how much b and c lean along the middle bond
        c2 = dot(b3, b2) / b2_sq

        f_b = -(1.0 + c1) * f_a + c2 * f_d                # calculate the final forces of b and c, by taking into account the forces of a and d
        f_c = f_c = -(f_a + f_b + f_d)

        atoms[a1].force += f_a                         # apply force to atom a
        atoms[a2].force += f_b                         # apply force to atom b
        atoms[a3].force += f_c                         # apply force to atom c
        atoms[a4].force += f_d                         # apply force to atom d
        
        pot_e_total += torsion_e  # calculate the potential energy

        end_bonds = time.perf_counter()

    end_torsions = time.perf_counter()
    global time_taken_torsions
    time_taken_torsions = end_torsions - start_torsions

    
    #LJ forces and Coulombs (non bonded loop)
    start_non_bonded = time.perf_counter()

    for i in range(no_balls):                   # loop over each unique atom pair once
        for j in range(i + 1, no_balls):

            scale = 1.0
            pair = tuple(sorted((i, j)))

            if pair in bonded_12 or pair in bonded_13:  # compleately skips 1-2 and 1-3 bonds
                continue
            
            if pair in bonded_14:                   # dapens 1-4 LJ forces
                scale = 0.5
            else:
                scale = 1

            r_vec = atoms[j].pos - atoms[i].pos     # find pair distances ect
            PBC_box_for_vectors(r_vec)              # update "copy vectors"
            r = mag(r_vec)
            
            if r == 0:             # skip invalid
                continue

            r_hat = norm(r_vec)

            F_LJ = calc_LJ_force(r, r_hat)                # calculate LJ force
            
            q1 = float(atoms[i].charge)                    # get atom charges
            q2 = float(atoms[j].charge)
            
            F_C = calc_coulombs(q1, q2, r, r_hat, COULOMB_CONSTANT)

            total_force = (scale*F_LJ) + (scale*F_C)
            atoms[i].force += total_force                      # Newtons third law - apply to both atoms
            atoms[j].force -= total_force

            pot_e_total += scale * 4 * LJ_EPSILON * ((LJ_SIGMA / r)**12 - (LJ_SIGMA / r)**6) # add the LJ potential energy on
            pot_e_total += scale * (COULOMB_CONSTANT * ((q1*q2)/ r))

    end_non_bonded = time.perf_counter()
    global time_taken_non_bonded
    time_taken_non_bonded = end_non_bonded - start_non_bonded

    end_all_foces = time.perf_counter()
    global time_taken_all_foces
    time_taken_all_foces = end_all_foces - start_all_foces

    return pot_e_total                                





# initial forces and model before simulation starts - need valid forces before starting
pot_e = calc_physics()
render_initial_model()

step = 0
relaxing = True 
relax_steps = 4000
E0 = None

# Main light (camera light)
cam_light = local_light(
    pos=scene.camera.pos - scene.camera.axis.norm()*4,
    color=color.gray(0.65)
)

# running simulation
while True:
    start_total_time = time.perf_counter()
    update_camera(TIME_STEP)
    
    # lighting position for camera
    cam_light.pos = scene.camera.pos - scene.camera.axis.norm()*3
    
    # Verlet integration method:
    # 1. half-step velocity update
    for i in range(no_balls):
        atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * TIME_STEP     # first half-step velocity update - uses current force to push velocity halfway forward

    # 2. position update
    for i in range(no_balls):
        atoms[i].pos += atoms[i].vel * TIME_STEP                      # update pos

    # 3. wrap positions back into box
    wrap_molecules()

    # 4. reset forces
    for i in range(no_balls):
        atoms[i].force = vector(0, 0, 0)                             # clear old forces before computing new ones

    # 5. compute new forces
    pot_e = calc_physics()                                           # now get new forces at the new positions

    # 6. second half-step velocity update
    for i in range(no_balls):
        atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * TIME_STEP     # compleate  the full velocity update using the new forces

    # 7. dampen starting strains - ensures the system is calm so it doesent blow up to begin with
    if step < relax_steps/5:
        damping = 0.99
    elif step < relax_steps/2:
        damping = 0.995
    elif step < relax_steps:
        damping = 0.999
    else:
        damping = 1                                             # no more dampening
        
    for i in range(no_balls):
        atoms[i].vel *= damping                                 # dampens some of the velocity at each step when begining

    k_e = calc_KE()                                             # energy tracking
    total_e = k_e + pot_e
    step += 1

    start_total_graphics = time.perf_counter() 
    # update graphics once per frame
    for i in range(no_balls):
        balls[i].pos = atoms[i].pos                                 # moves balls to current pos

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = atoms[bond.a1].pos
        bond_visuals[idx].axis = atoms[bond.a2].pos - atoms[bond.a1].pos

    end_total_graphics = time.perf_counter()
    time_taken_graphics = end_total_graphics - start_total_graphics


    end_total_time = time.perf_counter()
    time_taken_total = end_total_time - start_total_time

    time_other = time_taken_total - (time_taken_all_foces + time_taken_graphics)

    print(f"Time By % of Graphics: {(time_taken_graphics/time_taken_total)*100:.6f} %")
    print(f"Time By % of Exclusions: {(time_taken_exclusions/time_taken_total)*100:.6f} %")
    print(f"Time By % of Bonds: {(time_taken_bonds/time_taken_total)*100:.6f} %")
    print(f"Time By % of Bond Angles: {(time_taken_bond_angles/time_taken_total)*100:.6f} %")
    print(f"Time By % of Torsions: {(time_taken_torsions/time_taken_total)*100:.6f} %")
    print(f"Time By % of non_bonded: {(time_taken_non_bonded/time_taken_total)*100:.6f} %")
    print(f"Time By % of total forces: {(time_taken_all_foces/time_taken_total)*100:.6f} %")
    print(f"Time By % of other: {(time_other/time_taken_total)*100:.6f} %")








time_taken_all_foces = 0
time_taken_exclusions = 0
time_taken_bonds = 0
time_taken_bond_angles = 0
time_taken_torsions = 0
time_taken_non_bonded = 0
time_taken_graphics = 0
time_taken_total = 0
time_other = 0

