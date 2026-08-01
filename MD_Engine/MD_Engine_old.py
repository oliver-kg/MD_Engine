from vpython import *                   # imports all I use
import tkinter as tk
import math
import read_molecules
import time
from random import uniform
import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Length = nm
# Mass = amu
# Time = ps
# Charge = e
# Energy = KL/mol
#
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
    def __init__(self, index, element):

        props = ELEMENTS[element]

        self.index = index

        self.radius = props["radius"]
        self.element = element
        self.base_colour = props["color"]

        self.neighbours = []
        self.last_neighbour_reference = np.zeros(3)
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
    def __init__(self, ideal_ang_deg, a1, a2, a3, K):
        self.ideal_ang = math.radians(ideal_ang_deg)
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.K = K
    
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

class System:
    def __init__(self):

        # Dynamic state
        self.positions = None
        self.velocities = None
        self.forces = None

        # Constant per-atom data
        self.masses = None
        self.charges = None

        # Topology
        self.atoms = []
        self.molecules = []
        self.bonds = []
        self.bond_angles = []
        self.torsion_angles = []

        # Bonds
        self.bond_a = None
        self.bond_b = None
        self.bond_r0 = None
        self.bond_k = None

        # Angles

        self.angle_i = None
        self.angle_j = None
        self.angle_k = None
        self.angle_theta0 = None
        self.angle_kconst = None

        # Torsions

        self.torsion_i = None
        self.torsion_j = None
        self.torsion_k = None
        self.torsion_l = None
        self.torsion_terms = None


# molecule properties arrays
atoms = []
balls = []
bonds = []
bond_angles = []
torsion_angles = []
bond_visuals = []

real_space_PE = 0.0
reciprocal_PE = 0.0

molecules = []

# simulation pararmeters
pot_e = 0
k_e = 0
total_e = 0
TIME_STEP = 0.0004
LJ_SIGMA = 0.3
LJ_EPSILON = 0.02
LJ_CUTOFF = 2.5 * LJ_SIGMA                        # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
SKIN_CUTOFF = 0.5
REAL_CUTOFF = 9.0
NEIGHBOUR_CUTOFF = SKIN_CUTOFF+REAL_CUTOFF
COULOMB_CONSTANT = 138.935456
PBC_BOX_LENGTH = 8
MOLECULE_NUMBERS = 10
CELL_SIZE = NEIGHBOUR_CUTOFF
steps = []
avg_sys_energy_values = []

time_taken_all_foces = 0
time_taken_exclusions = 0
time_taken_bonds = 0
time_taken_bond_angles = 0
time_taken_torsions = 0
time_taken_non_bonded = 0
time_taken_LJ = 0
time_taken_coulombs = 0
time_taken_graphics = 0
time_taken_total = 0
time_other = 0

# ------------------------
# PME PARAMETERS
# ------------------------

PME_ALPHA = 0.35
PME_GRID = 32
BSPLINE_ORDER = 4

rho = np.zeros((PME_GRID, PME_GRID, PME_GRID), dtype=float)

# -- SETUP --

plt.ion()  # Interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlabel("Timestep")
ax.set_ylabel("Average Relative Error")
ax.set_title("Average Relative Error (100-step average)")
ax.grid(True)

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
    for i, a in enumerate(atoms):
        b = sphere(pos=vector(*positions[i]), radius=a.radius, color=a.base_colour, make_trail=False)
        balls.append(b)

    # Draw lines between bonds
    for bond in system.bonds:

        c = cylinder(
            pos=vector(*positions[bond.a1]),
            axis=vector(*positions[bond.a2]) - vector(*positions[bond.a1]),
            radius=0.01,
            color=vector(1,1,1)
        )
        bond_visuals.append(c)

def load_molecule_template():

    atom_pos = []
    atom_charge = []
    atom_type = []

    for ind in read_molecules.atom_index:
        atom_pos.append(read_molecules.atom_pos[ind])
        atom_charge.append(read_molecules.atom_charge[ind])
        atom_type.append(read_molecules.atom_type[ind])

    bonds = []

    for ind in range(len(read_molecules.a_a)):
        bonds.append((
            read_molecules.a_a[ind],
            read_molecules.a_b[ind],
            read_molecules.r_0[ind],
            read_molecules.k_engine[ind]
        ))

    angles = []

    for ind in range(len(read_molecules.a_i)):
        angles.append((
            read_molecules.a_i[ind],
            read_molecules.a_j[ind],
            read_molecules.a_k[ind],
            read_molecules.b_angle[ind],
            read_molecules.k_ang[ind]
        ))

    torsions = []

    for ind in range(len(read_molecules.d_i)):
        torsions.append((
            read_molecules.d_i[ind],
            read_molecules.d_j[ind],
            read_molecules.d_k[ind],
            read_molecules.d_l[ind],
            read_molecules.k_dih[ind],
            read_molecules.n[ind],
            read_molecules.ph[ind]
        ))

    return atom_pos, atom_charge, atom_type, bonds, angles, torsions

def insert_molecule(template):

    atom_pos, atom_charge, atom_type, bond_data, angle_data, torsion_data = template

    # Where this copy will go
    offset_position = vector(
        uniform(0, PBC_BOX_LENGTH/3),
        uniform(0, PBC_BOX_LENGTH/3),
        uniform(0, PBC_BOX_LENGTH/3)
    )

# First new atom index
    atom_offset = len(atoms)

    atom_indices = []

    # Create atoms
    for pos, charge, atom_type_name in zip(atom_pos,atom_charge,atom_type):

        index = len(atoms)

        atoms.append(
            Atom(
                index,
                atom_type_name
            )
        )

        new_pos = pos + offset_position

        position_list.append([
            new_pos.x,
            new_pos.y,
            new_pos.z
        ])

        velocity_list.append([0.0, 0.0, 0.0])
        force_list.append([0.0, 0.0, 0.0])

        charge_list.append(float(charge))
        mass_list.append(ELEMENTS[atom_type_name]["mass"])

        atom_indices.append(len(atoms)-1)

    molecules.append(Molecule(atom_indices))

    # Create bonds
    for a, b, r0, k in bond_data:

        bonds.append(
            Bond(
                a + atom_offset,
                b + atom_offset,
                r0,
                k
            )
        )

    # Create angles
    for i, j, k, theta0, k_theta in angle_data:

        bond_angles.append(
            BondAngle(
                theta0,
                i + atom_offset,
                j + atom_offset,
                k + atom_offset,
                k_theta
            )
        )

    # Create torsions
    for i, j, k, l, k_dih, n, phase in torsion_data:

        torsion_angles.append(
            TorsionAngle(
                i + atom_offset,
                j + atom_offset,
                k + atom_offset,
                l + atom_offset,
                [(k_dih, n, phase)]
            )
        )

position_list = []
velocity_list = []
force_list = []
charge_list = []
mass_list = []

#create test molecules

water = load_molecule_template()

for i in range(MOLECULE_NUMBERS):
    insert_molecule(water)

# create number of atoms after all atoms have been made
no_balls = len(atoms)
# create NumPy Vectors

positions = np.zeros((no_balls,3), dtype=np.float64)
velocities = np.zeros((no_balls,3), dtype=np.float64)
forces = np.zeros((no_balls,3), dtype=np.float64)

charges = np.zeros(no_balls)
masses = np.zeros(no_balls)

positions = np.asarray(position_list, dtype=np.float64)
velocities = np.asarray(velocity_list, dtype=np.float64)
forces = np.asarray(force_list, dtype=np.float64)

charges = np.asarray(charge_list, dtype=np.float64)
masses = np.asarray(mass_list, dtype=np.float64)



# Bonds
bond_a = np.asarray(
    [bond.a1 for bond in bonds],
    dtype=np.int32
)

bond_b = np.asarray(
    [bond.a2 for bond in bonds],
    dtype=np.int32
)

bond_r0 = np.asarray(
    [bond.ideal_dist for bond in bonds],
    dtype=np.float64
)

bond_k = np.asarray(
    [bond.K for bond in bonds],
    dtype=np.float64
)

# Bond Angles
angle_i = np.asarray(
    [bond_angle.a1 for bond_angle in bond_angles],
    dtype=np.int32
)

angle_j = np.asarray(
    [bond_angle.a2 for bond_angle in bond_angles],
    dtype=np.int32
)

angle_k = np.asarray(
    [bond_angle.a3 for bond_angle in bond_angles],
    dtype=np.int32
)

angle_theta0 = np.asarray(
    [bond_angle.ideal_ang for bond_angle in bond_angles],
    dtype=np.float64
)

angle_kconst = np.asarray(
    [bond_angle.K for bond_angle in bond_angles],
    dtype=np.float64
)


# Torsion Angles
torsion_i = np.asarray(
    [bond_angle.a1 for bond_angle in torsion_angles],
    dtype=np.int32
)

torsion_j = np.asarray(
    [bond_angle.a2 for bond_angle in torsion_angles],
    dtype=np.int32
)

torsion_k = np.asarray(
    [bond_angle.a3 for bond_angle in torsion_angles],
    dtype=np.int32
)

torsion_l = np.asarray(
    [bond_angle.a4 for bond_angle in torsion_angles],
    dtype=np.int32
)

torsion_psi = np.asarray(
    [bond_angle.psi for bond_angle in torsion_angles],
    dtype=np.float64
)

torsion_terms = np.asarray(
    [bond_angle.terms for bond_angle in torsion_angles],
    dtype=np.float64
)



system = System()

system.positions = positions
system.velocities = velocities
system.forces = forces

# Constant per-atom data
system.masses = masses
system.charges = charges

# Topology
system.atoms = atoms
system.molecules = molecules
system.bonds = bonds
system.bond_angles = bond_angles
system.torsion_angles = torsion_angles

# Bonds
system.bond_a = bond_a
system.bond_b = bond_b
system.bond_r0 = bond_r0
system.bond_k = bond_k

# Angles
system.angle_i = angle_i
system.angle_j = angle_j
system.angle_k = angle_k
system.angle_theta0 = angle_theta0
system.angle_kconst = angle_kconst

# Torsions

system.torsion_i = torsion_i
system.torsion_j = torsion_j
system.torsion_k = torsion_k
system.torsion_l = torsion_l
system.terms = torsion_terms

# Lennard Jones Force Equation (VDW's)
def calc_LJ_force(dist, direction):
    return  -((24*LJ_EPSILON*((2*((LJ_SIGMA/dist)**12)) - ((LJ_SIGMA/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def calc_KE(system):
    k_e = 0
    for i in range(no_balls): 
        k_e = k_e + (0.5*system.masses[i]*(np.dot(system.velocities[i], system.velocities[i]))) # use dot product as v^2 to calc ke
    return k_e

# calc Coulombs
def calc_coulombs(q1, q2, r, direction, k):

    if r == 0:
        return np.zeros(3)

    erfc_term = math.erfc(PME_ALPHA * r)

    exp_term = math.exp(-(PME_ALPHA * r) ** 2)

    force_mag = (
        k
        * q1
        * q2
        * (
            erfc_term / (r * r)
            +
            (2.0 * PME_ALPHA / math.sqrt(math.pi))
            * exp_term / r
        )
    )

    return -force_mag * direction

# builds the neibour list structure
def build_neighbours():
    neighbours = {i: set() for i in range(no_balls)}        # set makes sure no duplicates and fast to check
    for bond in system.bonds:
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
    time_taken_exclusions += exclusions_end - exclusions_start

    return bonded_12, bonded_13, bonded_14

bonded_12, bonded_13, bonded_14 = build_exclusion_sets()    # builds the sets to be used

# checks if a molecule has reached the box boundery and needs warping
def wrap_molecules(system):
    half = PBC_BOX_LENGTH / 2

    for mol in system.molecules:
        ref_pos = system.positions[mol.atom_indices[0]]
        shift = np.zeros(3)

        for d in range(3):                              # loop over x, y, z instead of one if/elif each
            while ref_pos[d] + shift[d] > half:
                shift[d] -= PBC_BOX_LENGTH
            while ref_pos[d] + shift[d] < -half:
                shift[d] += PBC_BOX_LENGTH

        if shift[0] != 0 or shift[1] != 0 or shift[2] != 0:
            for i in mol.atom_indices:
                system.positions[i] += shift

# minimum image periodic correction - turns the vector into the nearest image vector
def minimum_image(r_vector):                                  

    half = PBC_BOX_LENGTH * 0.5
    box = PBC_BOX_LENGTH

    if r_vector[0] > half:
        r_vector[0] -= box
    elif r_vector[0] < -half:
        r_vector[0] += box

    if r_vector[1] > half:
        r_vector[1] -= box
    elif r_vector[1] < -half:
        r_vector[1] += box

    if r_vector[2] > half:
        r_vector[2] -= box
    elif r_vector[2] < -half:
        r_vector[2] += box

    return r_vector

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

def coulombs_calculations(system, scale, q1, q2, r, r_hat, i, j):
    F_C = calc_coulombs(q1, q2, r, r_hat, COULOMB_CONSTANT)
    coulombs_force_calculated = scale*F_C
    system.forces[i] += coulombs_force_calculated                      # Newtons third law - apply to both atoms
    system.forces[j] -= coulombs_force_calculated
    global pot_e_total

    global real_space_PE
    energy = (
        scale
        * COULOMB_CONSTANT
        * q1
        * q2
        * math.erfc(PME_ALPHA * r)
        / r
    )

    real_space_PE += energy
    pot_e_total += energy

def LJ_calculations(system, r, r_hat, scale, LJ_energy_shift, i, j):
    LJ_force_calculated = scale * calc_LJ_force(r, r_hat)
    system.forces[i] += LJ_force_calculated                      # Newtons third law - apply to both atoms
    system.forces[j] -= LJ_force_calculated
    global pot_e_total
    pot_e_total += scale * (4 * LJ_EPSILON * ((LJ_SIGMA / r)**12 - (LJ_SIGMA / r)**6) - LJ_energy_shift) # add the LJ potential energy on

def build_neighbour_lists(system, neighbour_cutoff):
    cutoff2 = neighbour_cutoff**2
    
    for atom in system.atoms:
        atom.neighbours.clear()         # clear out all of the previous neighbour lists

    cells = build_cell_list(system)

    for i in range(no_balls):  
        my_cell = find_cell_of_atom(system.positions[i])

        for cell in neighbouring_cells(my_cell):
            for j in cells.get(cell, ()):
                if j <= i:
                    continue

                r_vec = minimum_image(system.positions[j] - system.positions[i])
                r2 = np.dot(r_vec, r_vec)                          # distance between pairs without square rooting - more efficient

                if r2 > cutoff2:
                    continue
                
                atoms[i].neighbours.append(j)
                atoms[j].neighbours.append(i)

    # update the new last neighbour refference for each atom
    for i, atom in enumerate(system.atoms):
        atom.last_neighbour_reference = system.positions[i].copy()

def find_cell_of_atom(atom_pos):
    return (
        math.floor(atom_pos[0] / CELL_SIZE),
        math.floor(atom_pos[1] / CELL_SIZE),
        math.floor(atom_pos[2] / CELL_SIZE)

    )

def build_cell_list(system):
    cells = {}

    for i in range(no_balls):
        cell = find_cell_of_atom(system.positions[i])

        if cell not in cells:
            cells[cell] = []

        cells[cell].append(i)

    return cells

def neighbouring_cells(cell):
    x,y,z = cell

    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield (x+dx, y+dy, z+dz)

# -------- PME --------

# Convert a VPython position into fractional mesh coordinates
def position_to_mesh(pos):
    half = PBC_BOX_LENGTH / 2
    x = (pos[0] + half) / PBC_BOX_LENGTH
    y = (pos[1] + half) / PBC_BOX_LENGTH
    z = (pos[2] + half) / PBC_BOX_LENGTH

    return (
        x * PME_GRID,
        y * PME_GRID,
        z * PME_GRID
    )

# Cardinal B-spline of arbitrary order
def bspline(n, u):

    # First-order spline (piecewise constant)
    if n == 1:
        if 0.0 <= u < 1.0:
            return 1.0
        return 0.0

    return (
        (u / (n - 1)) * bspline(n - 1, u)
        +
        ((n - u) / (n - 1)) * bspline(n - 1, u - 1)
    )

def compute_theta(f):

    theta = np.zeros(BSPLINE_ORDER)
    dtheta = np.zeros(BSPLINE_ORDER)

    for i in range(BSPLINE_ORDER):

        u = f + (BSPLINE_ORDER - 1 - i)

        theta[i] = bspline(BSPLINE_ORDER, u)
        dtheta[i] = bspline_derivative(BSPLINE_ORDER, u)

    return theta, dtheta

def build_charge_grid(system):

    rho.fill(0.0)
    spline_order = 4

    for i in range(no_balls):

        gx, gy, gz = position_to_mesh(system.positions[i])

        # Lower-left grid point
        ix = math.floor(gx)
        iy = math.floor(gy)
        iz = math.floor(gz)

        # Fractional offsets inside the cell
        fx = gx - ix
        fy = gy - iy
        fz = gz - iz

        theta_x, _ = compute_theta(fx)
        theta_y, _ = compute_theta(fy)
        theta_z, _ = compute_theta(fz)

        # Spread charge over a 4×4×4 cube
        for dx in range(spline_order):

            wx = theta_x[dx]

            for dy in range(spline_order):

                wy = theta_y[dy]

                for dz in range(spline_order):

                    wz = theta_z[dz]

                    rho[
                        (ix + dx) % PME_GRID,
                        (iy + dy) % PME_GRID,
                        (iz + dz) % PME_GRID
                    ] += system.charges[i] * wx * wy * wz

    return rho

def build_reciprocal_kernel():

    global C

    C = np.zeros(
        (PME_GRID, PME_GRID, PME_GRID),
        dtype=float
    )

    box = PBC_BOX_LENGTH

    for i in range(PME_GRID):

        if i <= PME_GRID//2:
            kx = i
        else:
            kx = i - PME_GRID

        for j in range(PME_GRID):

            if j <= PME_GRID//2:
                ky = j
            else:
                ky = j - PME_GRID

            for k in range(PME_GRID):

                if k <= PME_GRID//2:
                    kz = k
                else:
                    kz = k - PME_GRID

                if kx == ky == kz == 0:
                    continue

                k2 = (
                    (kx/box)**2 +
                    (ky/box)**2 +
                    (kz/box)**2
                )

                C[i,j,k] = (
                    np.exp(
                        -(math.pi**2)*k2/(PME_ALPHA**2)
                    )
                    /
                    (
                        math.pi
                        *
                        box**3
                        *
                        k2
                    )
                )

    return C

C = build_reciprocal_kernel()

def b_factor(m):

    total = 0j

    for k in range(BSPLINE_ORDER - 1):

        total += (
            bspline(BSPLINE_ORDER, k + 1)
            *
            np.exp(
                -2j * math.pi * m * k / PME_GRID
            )
        )

    if abs(total) < 1e-12:
        return 0.0

    phase = np.exp(
        -2j * math.pi * (BSPLINE_ORDER-1) * m / PME_GRID
    )

    b = phase / total

    return abs(b)**2

def build_B():

    global B

    B = np.zeros(
        (PME_GRID, PME_GRID, PME_GRID),
        dtype=float
    )

    for i in range(PME_GRID):

        if i <= PME_GRID//2:
            mx = i
        else:
            mx = i - PME_GRID

        bx = b_factor(mx)

        for j in range(PME_GRID):

            if j <= PME_GRID//2:
                my = j
            else:
                my = j - PME_GRID

            by = b_factor(my)

            for k in range(PME_GRID):

                if k <= PME_GRID//2:
                    mz = k
                else:
                    mz = k - PME_GRID

                bz = b_factor(mz)

                B[i,j,k] = bx * by * bz

    return B

B = build_B()

BC = B * C

def bspline_derivative(n, u):

    if n <= 1:
        return 0.0

    return (
        bspline(n - 1, u)
        -
        bspline(n - 1, u - 1)
    )

def reciprocal_interpolation(system, potential):

    scale = PME_GRID / PBC_BOX_LENGTH

    for i in range(no_balls):

        gx, gy, gz = position_to_mesh(system.positions[i])

        ix = math.floor(gx)
        iy = math.floor(gy)
        iz = math.floor(gz)

        fx = gx - ix
        fy = gy - iy
        fz = gz - iz

        theta_x, dtheta_x = compute_theta(fx)
        theta_y, dtheta_y = compute_theta(fy)
        theta_z, dtheta_z = compute_theta(fz)


        dPhidx = 0.0
        dPhidy = 0.0
        dPhidz = 0.0

        Phi = 0.0

        for dx in range(BSPLINE_ORDER):

            wx = theta_x[dx]
            dwx = dtheta_x[dx]

            for dy in range(BSPLINE_ORDER):

                wy = theta_y[dy]
                dwy = dtheta_y[dy]

                for dz in range(BSPLINE_ORDER):

                    wz = theta_z[dz]
                    dwz = dtheta_z[dz]

                    phi = potential[
                        (ix+dx) % PME_GRID,
                        (iy+dy) % PME_GRID,
                        (iz+dz) % PME_GRID
                    ].real

                    Phi += phi * wx * wy * wz

                    dPhidx += phi * dwx * wy * wz
                    dPhidy += phi * wx * dwy * wz
                    dPhidz += phi * wx * wy * dwz

        global pot_e_total
        global reciprocal_PE

        energy = 0.5 * COULOMB_CONSTANT * system.charges[i] * Phi

        reciprocal_PE += energy
        pot_e_total += energy

        Fx = COULOMB_CONSTANT * system.charges[i] * scale * dPhidx
        Fy = COULOMB_CONSTANT * system.charges[i] * scale * dPhidy
        Fz = COULOMB_CONSTANT * system.charges[i] * scale * dPhidz
        system.forces[i] -= np.array([Fx, Fy, Fz])

PME_SELF_ENERGY = 0.0
PME_SELF_ENERGY += np.sum(system.charges * system.charges)

PME_SELF_ENERGY *= (
    -COULOMB_CONSTANT
    * PME_ALPHA
    / math.sqrt(math.pi)
)

def coulomb_exclusion_correction(system, q1, q2, r, r_hat, i, j):
    global pot_e_total, real_space_PE
    if r < 1e-12:
        return
    erf_term = math.erf(PME_ALPHA * r)
    exp_term = math.exp(-(PME_ALPHA * r) ** 2)

    # cancels the erf(r)/r contribution the reciprocal grid implicitly added for this excluded pair
    force_mag = COULOMB_CONSTANT * q1 * q2 * (
        erf_term / (r * r) - (2.0 * PME_ALPHA / math.sqrt(math.pi)) * exp_term / r
    )
    F = force_mag * r_hat
    system.forces[i] += F
    system.forces[j] -= F

    energy = -COULOMB_CONSTANT * q1 * q2 * erf_term / r
    real_space_PE += energy
    pot_e_total += energy

def apply_exclusion_corrections(system):
    for (i, j) in bonded_12 | bonded_13:
        r_vec = minimum_image(system.positions[j] - system.positions[i])
        r = np.linalg.norm(r_vec)
        if r < 1e-12:
            continue
        r_hat = r_vec / r
        q1 = float(system.charges[i])
        q2 = float(system.charges[j])
        coulomb_exclusion_correction(system, q1, q2, r, r_hat, i, j)


# -------- Other Forces --------

def bond_forces(system):
    start_bonds = time.perf_counter()

    for i in range(len(system.bond_a)):
        a = system.bond_a[i]
        b = system.bond_b[i]
        r0 = system.bond_r0[i]
        k = system.bond_k[i]

        # compute the bond vector and bond lengths
        r_vec = minimum_image(system.positions[b] - system.positions[a])
        r = np.linalg.norm(r_vec)

        if r < 1e-12:                              # avoid devide by zero just in case
            continue

        r_hat = r_vec / r                           # bond direction

        F = 2 * k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds

        system.forces[a] += F                          # apply Newtons third law - equal and opposite forces
        system.forces[b] -= F

        global pot_e_total
        pot_e_total += k * (r - r0)**2                # calc the harmonic bond potential in function 2 type bond



    end_bonds = time.perf_counter()
    global time_taken_bonds
    time_taken_bonds += end_bonds - start_bonds

def angle_forces(system):
    start_bond_angles = time.perf_counter()
    
    for i in range(len(system.angle_i)):

        ideal_ang = system.angle_theta0[i]
        atomA = system.angle_i[i]
        atomB = system.angle_j[i]
        atomC = system.angle_k[i]
        k = system.angle_kconst[i]

        BA = system.positions[atomA] - system.positions[atomB]                                # find vector of B to A
        BC = system.positions[atomC] - system.positions[atomB]                                # find vector of B to C

        BA = minimum_image(BA)                                                 # update the coppy vectors
        BC = minimum_image(BC)
        
        r_BA = np.linalg.norm(BA)                                                          # get the lengths, as force depends on direction, and how long the arms are
        r_BC = np.linalg.norm(BC)
        
        if r_BA < 1e-12 or r_BC < 1e-12:
            continue
        
        theta_cos = (np.dot(BA, BC))/(r_BA*r_BC)                                           # calculate the dot product hrer
        theta_cos = max(-1.0, min(1.0, theta_cos))                                      # clamp for safety
        theta = math.acos(theta_cos)                                                    
        theta = theta                                                             # update theta

        sin_theta = math.sin(theta)                                                     # to not divide by zero
        if abs(sin_theta) < 1e-8:
            continue

        dU_dtheta = k * (theta - ideal_ang)                                             # the error of the current angle vs ideal bond angle

        f_a = -(dU_dtheta / sin_theta) * ((theta_cos / (r_BA * r_BA)) * BA - (1.0 / (r_BA * r_BC)) * BC) # build the part of the force on A that changes the angle without streaching BA

        f_c = -(dU_dtheta / sin_theta) * ((theta_cos / (r_BC * r_BC)) * BC - (1.0 / (r_BA * r_BC)) * BA) # build the part of the force on C that changes the angle without streaching BC

        f_b = -(f_a + f_c)                                                      # force to be applied on b to conserve energy
    
        system.forces[atomA] += f_a                                               # update the forces
        system.forces[atomB] += f_b
        system.forces[atomC] += f_c

        global pot_e_total
        pot_e_total += 0.5 * k * (theta - ideal_ang)**2                # update the potential energy change from the angle

    end_bond_angles = time.perf_counter()
    global time_taken_bond_angles
    time_taken_bond_angles += end_bond_angles - start_bond_angles

def torsion_forces(system):
    start_torsions = time.perf_counter()
    
    for i in range(len(system.torsion_i)):
        a1 = system.torsion_i[i]                           # 4 atoms involved in the torsion bond
        a2 = system.torsion_j[i] 
        a3 = system.torsion_k[i]
        a4 = system.torsion_l[i]
        psi = system.psi[i]


        b1 = system.positions[a2] - system.positions[a1]      # directions of the three bonds
        b2 = system.positions[a3] - system.positions[a2]
        b3 = system.positions[a4] - system.positions[a3]

        b1 = minimum_image(b1)                 # update the vectors of coppies
        b2 = minimum_image(b2)
        b3 = minimum_image(b3)

        n1 = np.cross(b1, b2)                      # normal to plane ABC
        n2 = np.cross(b2, b3)                      # normal to plane BCD

        eps = 1e-12
        n1_sq = np.dot(n1, n1)
        n2_sq = np.dot(n2, n2)
        b2_sq = np.dot(b2, b2)
        b2_mag = np.linalg.norm(b2)

        if n1_sq < eps or n2_sq < eps or b2_sq < eps:    # makes sure small values dont blow the system up
            continue
        
        x = np.dot(n1, n2)
        y = np.dot(np.cross(n1, n2),b2 / b2_mag)

        psi = math.atan2(y, x)                  # caluclate current psi angle 

        torsion_e = 0
        dV_dpsi = 0
        
        for k, n, delta in system.torsion_terms:                             # repeats for each term, updating the cos graph
            torsion_e += k * (1 + math.cos(n * psi - delta))   # calculate the potential energy
            dV_dpsi -= k * n * math.sin(n * psi - delta)            # how strongly the torsion wants to rotate
        
            
        fa_pref = dV_dpsi * (b2_mag / n1_sq)           # calculates the geometric scallings of the force
        fd_pref = -dV_dpsi * (b2_mag / n2_sq)

        f_a = fa_pref * n1                              # aligning the forces with the direction to the plane
        f_d = fd_pref * n2

        c1 = np.dot(b1, b2) / b2_sq                        # calculates how much b and c lean along the middle bond
        c2 = np.dot(b3, b2) / b2_sq

        f_b = -(1.0 + c1) * f_a + c2 * f_d                # calculate the final forces of b and c, by taking into account the forces of a and d
        f_c = f_c = -(f_a + f_b + f_d)

        system.forces[a1] += f_a                         # apply force to atom a
        system.forces[a2] += f_b                         # apply force to atom b
        system.forces[a3] += f_c                         # apply force to atom c
        system.forces[a4] += f_d                         # apply force to atom d

        global pot_e_total
        pot_e_total += torsion_e  # calculate the potential energy

    end_torsions = time.perf_counter()
    global time_taken_torsions
    time_taken_torsions += end_torsions - start_torsions

def non_bonded_forces(system):
    LJ_scale = 1.0
    coulomb_scale = 1.0
    LJ_energy_shift = 4 * LJ_EPSILON * ((LJ_SIGMA / LJ_CUTOFF)**12 - (LJ_SIGMA / LJ_CUTOFF)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

    for i, atom in enumerate(system.atoms):                   # loop over each unique atom pair once
        for j in atom.neighbours:

            if j <= i:                                  # makes sure not calculating the same force twice
                continue

            pair = tuple(sorted((i, j)))

            if pair in bonded_12 or pair in bonded_13:  # compleately skips 1-2 and 1-3 bonds
                continue
            
            if pair in bonded_14:                   # dapens 1-4 LJ forces
                LJ_scale = 0.5
            else:
                LJ_scale = 1                        # normal non-bonded forces on all other atoms
                coulomb_scale = 1

            r_vec = minimum_image(system.positions[j] - system.positions[i])     # find pair distances and update the "ghost molecules" 
            r = np.linalg.norm(r_vec)

            if r < 1e-12:
                continue

            r_hat = r_vec / r

            q1 = float(system.charges[i])             # get atom charges               
            q2 = float(system.charges[j])

            # -- add the coulombs back in when adding PME later

            if r > REAL_CUTOFF:             # skip invalid LJ, and only calculate coulomb forces
                coulombs_calculations(system, coulomb_scale, q1, q2, r, r_hat, i, j)
                continue

            else:                                   # dont skip any non bonded forces as pair is in the cutoff region 
                LJ_calculations(system, r, r_hat, LJ_scale, LJ_energy_shift, i, j)
                coulombs_calculations(system, coulomb_scale, q1, q2, r, r_hat, i, j)

# caculate physics stuff - the heart <3
def calc_physics():

    global pot_e_total
    pot_e_total = 0.0

    global real_space_PE
    global reciprocal_PE

    real_space_PE = 0.0
    reciprocal_PE = 0.0

    start_all_foces = time.perf_counter()

    # bond forces
    bond_forces(system)

    # bond angles
    angle_forces(system)

    # torsion angles
    torsion_forces(system)

    start_LJ = time.perf_counter()

    # #LJ forces and Coulombs (non bonded loop)
    non_bonded_forces(system)

    end_LJ = time.perf_counter()
    global time_taken_LJ
    time_taken_LJ += end_LJ - start_LJ

    # -----------------------
    # Reciprocal-space PME
    # -----------------------
    start_coulombs = time.perf_counter()

    apply_exclusion_corrections(system)
    rho = build_charge_grid(system)

    Qk = np.fft.fftn(rho)
    Qk *= BC

    potential = np.fft.ifftn(Qk).real
    potential *= PME_GRID**3

    reciprocal_interpolation(system, potential)

    # account for the self energy shift
    pot_e_total += PME_SELF_ENERGY

    end_coulombs = time.perf_counter()
    global time_taken_coulombs
    time_taken_coulombs += end_coulombs - start_coulombs

    global time_taken_non_bonded
    time_taken_non_bonded = time_taken_coulombs + time_taken_LJ
    
    end_all_foces = time.perf_counter()
    global time_taken_all_foces
    time_taken_all_foces += end_all_foces - start_all_foces

    return pot_e_total                                

# Main light (camera light)
cam_light = local_light(
    pos=scene.camera.pos - scene.camera.axis.norm()*4,
    color=color.gray(0.65)
)

# initial forces and model and neighbour list before simulation starts - need valid forces before starting
for i in range(no_balls):
    forces.fill(0.0)
pot_e = calc_physics()
render_initial_model()
build_neighbour_lists(system, NEIGHBOUR_CUTOFF)

step = 0
prev_step_count = 0
relax_steps = 2000
timestep_x = 100
average_total_energy = 0
max_displacement2 = 0
initial_total_energy = 0

# running simulation
while True:
    start_total_time = time.perf_counter()
    update_camera(TIME_STEP)
    
    # lighting position for camera
    cam_light.pos = scene.camera.pos - scene.camera.axis.norm()*3
    
    # Verlet integration method:
    # 1. half-step velocity update
    for i in range(no_balls):
        system.velocities[i] += 0.5 * (system.forces[i] /system.masses[i]) * TIME_STEP     # first half-step velocity update - uses current force to push velocity halfway forward

    # 2. position update
    for i in range(no_balls):
        system.positions[i] += system.velocities[i] * TIME_STEP                      # update pos

    # 3. wrap positions back into box
    wrap_molecules(system)

    # 4. reset forces
    for i in range(no_balls):
        forces.fill(0.0)                             # clear old forces before computing new ones

    # 5. compute new forces
    pot_e = calc_physics()                                           # now get new forces at the new positions

    # 6. second half-step velocity update
    for i in range(no_balls):
        system.velocities[i] += 0.5 * (system.forces[i] / system.masses[i]) * TIME_STEP     # compleate  the full velocity update using the new forces

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
        system.velocities[i] *= damping                                 # dampens some of the velocity at each step when begining

    k_e = calc_KE(system)                                             # energy tracking
    total_e = k_e + pot_e
    average_total_energy += total_e
    step += 1

    if step == relax_steps:
        average_total_energy = 0

    # find the displacement of the atom this timestep
    for i, atom in enumerate(system.atoms):
        disp2 = minimum_image(system.positions[i] - atom.last_neighbour_reference)

        displacement2 = r2 = np.dot(disp2,disp2)

        if displacement2 > max_displacement2:
            max_displacement2 = displacement2

    if max_displacement2 > (SKIN_CUTOFF / 2)**2:
        prev_step_count = step - prev_step_count
        print(f"New List After {prev_step_count} Steps")
        build_neighbour_lists(system, NEIGHBOUR_CUTOFF)
        max_displacement2 = 0
        
    start_total_graphics = time.perf_counter() 

    # ----- graphics render  ------

    # update graphics once per frame
    for i in range(no_balls):
        balls[i].pos = vector(*system.positions[i])                                 # moves balls to current pos

    for idx, bond in enumerate(system.bonds):                              # update the bond visuals
        bond_visuals[idx].pos = vector(*system.positions[bond.a1])

        bond_visuals[idx].axis = vector(
            *(system.positions[bond.a2] - system.positions[bond.a1])
        )

    # ----------------------------

    end_total_graphics = time.perf_counter()
    time_taken_graphics += end_total_graphics - start_total_graphics

    end_total_time = time.perf_counter()
    time_taken_total += end_total_time - start_total_time

    time_other = time_taken_total - (time_taken_graphics + time_taken_all_foces)



    # printing and debuging stuff that prints every x timesteps
    if step % timestep_x == 0:

        if step == (relax_steps+timestep_x):
            initial_total_energy = total_e

        if step >= (relax_steps+timestep_x):
            average_total_energy = average_total_energy/timestep_x
            steps.append(step)
            avg_sys_energy_values.append(average_total_energy-initial_total_energy)

            line.set_data(steps, avg_sys_energy_values)

            
            ax.relim()              # Recalculate limits
            ax.autoscale_view()     # Expand axes if needed
            plt.draw()
            plt.pause(0.001)
            average_total_energy = 0
        

        print()
        print(f"Step: {step}")
        print(f"KE: {k_e:.6f}  PE: {pot_e:.6f}  Total: {total_e:.6f}")
        print(f"Time Taken per cycle: {time_taken_total/timestep_x:.6f}s")
        print()
        print(f"Time By % of Graphics: {((time_taken_graphics/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bonds: {((time_taken_bonds/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bond Angles: {((time_taken_bond_angles/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Torsions: {((time_taken_torsions/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of LJ: {((time_taken_LJ/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Coulombs: {((time_taken_coulombs/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of non_bonded: {((time_taken_non_bonded/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of total forces: {((time_taken_all_foces/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of other: {((time_other/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        
        time_taken_all_foces = 0
        time_taken_exclusions = 0
        time_taken_bonds = 0
        time_taken_bond_angles = 0
        time_taken_torsions = 0
        time_taken_non_bonded = 0
        time_taken_coulombs = 0
        time_taken_LJ = 0
        time_taken_graphics = 0
        time_taken_total = 0
        time_other = 0