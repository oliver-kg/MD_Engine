from vpython import *                   # imports all I use
import tkinter as tk
import math
import read_molecules
import time
from random import uniform
import numpy as np
import matplotlib.pyplot as plt
import cupy as cp


# ----- Other Files -----
from bonded import bond_forces, angle_forces, torsion_forces
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

class Atom:
    def __init__(self, element):
        props = ELEMENTS[element]

        self.radius = props["radius"]
        self.element = element
        self.base_colour = props["color"]
        self.index = -1
        self.neighbours = []

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


class System:
    def __init__(self):

        self.positions = positions
        self.velocities = velocities
        self.forces = forces

        self.charges = charges
        self.masses = masses

        self.bonds = bonds
        self.bond_angles = bond_angles
        self.torsion_angles = torsion_angles
        self.molecules = molecules

        self.bond_a = bond_a
        self.bond_b = bond_b
        self.bond_r0 = bond_r0
        self.bond_k = bond_k

        self.PBC_BOX_LENGTH = PBC_BOX_LENGTH
   
    def minimum_image_many(self, r_vectors):
        half = self.PBC_BOX_LENGTH * 0.5
        box = self.PBC_BOX_LENGTH

        r_vectors[r_vectors > half] -= box
        r_vectors[r_vectors < -half] += box

        return r_vectors

    # minimum image periodic correction - turns the vector into the nearest image vector
    def minimum_image(self, r_vector):
        

        half = self.PBC_BOX_LENGTH * 0.5
        box = self.PBC_BOX_LENGTH

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

    # checks if a molecule has reached the box boundery and needs warping
    def wrap_molecules(self):
        positions = self.positions
        molecules = self.molecules
        PBC_BOX_LENGTH = self.PBC_BOX_LENGTH

        half = PBC_BOX_LENGTH / 2

        for mol in molecules:
            ref_pos = positions[mol.atom_indices[0]]
            shift = cp.zeros(3, dtype=cp.float64)

            if ref_pos[0] > half:
                shift[0] -= PBC_BOX_LENGTH
            elif ref_pos[0] < -half:
                shift[0] += PBC_BOX_LENGTH

            if ref_pos[1] > half:
                shift[1] -= PBC_BOX_LENGTH
            elif ref_pos[1] < -half:
                shift[1] += PBC_BOX_LENGTH

            if ref_pos[2] > half:
                shift[2] -= PBC_BOX_LENGTH
            elif ref_pos[2] < -half:
                shift[2] += PBC_BOX_LENGTH

            if shift[0] != 0 or shift[1] != 0 or shift[2] != 0:
                for i in mol.atom_indices:
                    positions[i] += shift



# molecule properties arrays
atoms = []
balls = []
bonds = []
bond_angles = []
torsion_angles = []
bond_visuals = []

initial_positions = []
initial_velocities = []
initial_charges = []
initial_masses = []

real_space_PE = 0.0
reciprocal_PE = 0.0

molecules = [
     Molecule([0])
    ]

# simulation pararmeters
pot_e = 0
k_e = 0
total_e = 0
TIME_STEP = 0.0001
LJ_SIGMA = 0.3
LJ_EPSILON = 0.02
LJ_CUTOFF = 2.5 * LJ_SIGMA                        # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
SKIN_CUTOFF = 0.5
REAL_CUTOFF = 1.1
NEIGHBOUR_CUTOFF = SKIN_CUTOFF+REAL_CUTOFF
COULOMB_CONSTANT = 138.935456
PBC_BOX_LENGTH = 5
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

PME_ALPHA = 3
PME_GRID = 32
BSPLINE_ORDER = 4

rho = cp.zeros((PME_GRID, PME_GRID, PME_GRID), dtype=float)

# ----- SETUP ------
'''
plt.ion()  # Interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlabel("Timestep")
ax.set_ylabel("Average System Energy")
ax.set_title("Average System Energy (100-step average)")
ax.grid(True)
'''
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
        b = sphere(pos=np_to_vec(positions[a.index]),radius=a.radius,color=a.base_colour,make_trail=False)
        balls.append(b)

    # Draw lines between bonds
    for bond in bonds:
        c = cylinder(pos=np_to_vec(positions[bond.a1]),axis=np_to_vec(positions[bond.a2] - positions[bond.a1]),radius=0.01,color=vector(1,1,1))
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
    for pos, charge, atom_type_name in zip(atom_pos,
                                           atom_charge,
                                           atom_type):

        atoms.append(Atom(atom_type_name))

        new_pos = pos + offset_position
        initial_positions.append(cp.asarray([new_pos.x, new_pos.y, new_pos.z], dtype=cp.float64))

        initial_velocities.append(cp.zeros(3))

        initial_charges.append(float(charge))

        initial_masses.append(ELEMENTS[atom_type_name]["mass"])

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
    for i, j, k_atom, theta0, k_theta in angle_data:

        bond_angles.append(
            BondAngle(
                theta0,
                0,
                i + atom_offset,
                j + atom_offset,
                k_atom + atom_offset,
                k_theta
            )
        )

    # Create torsions
    for i, j, k_atom, l, k_dih, n, phase in torsion_data:

        torsion_angles.append(
            TorsionAngle(
                i + atom_offset,
                j + atom_offset,
                k_atom + atom_offset,
                l + atom_offset,
                [(k_dih, n, phase)]
            )
        )

#create test molecules

molecule_template = load_molecule_template()

for i in range(2):
    insert_molecule(molecule_template)

def vec_to_np(v):
    return cp.array([v.x, v.y, v.z], dtype=cp.float64)

def np_to_vec(a):
    return vector(a[0], a[1], a[2])



# create number of atoms after all atoms have been made
no_balls = len(atoms)

positions = cp.asarray(initial_positions, dtype=cp.float64)
velocities = cp.asarray(initial_velocities, dtype=cp.float64)
forces = cp.zeros_like(positions)

charges = cp.asarray(initial_charges, dtype=cp.float64)
masses = cp.asarray(initial_masses, dtype=cp.float64)

last_neighbour_reference = cp.zeros((no_balls,3), dtype=cp.float64)

bond_a = cp.asarray(
    [bond.a1 for bond in bonds],
    dtype=cp.int32
)

bond_b = cp.asarray(
    [bond.a2 for bond in bonds],
    dtype=cp.int32
)

bond_r0 = cp.asarray(
    [bond.ideal_dist for bond in bonds],
    dtype=cp.float64
)

bond_k = cp.asarray(
    [bond.K for bond in bonds],
    dtype=cp.float64
)


for i, atom in enumerate(atoms):
    atom.index = i


system = System()
# ----- Components -----

# Lennard Jones Force Equation (VDW's)
def calc_LJ_force(dist, direction):
    return  -((24*LJ_EPSILON*((2*((LJ_SIGMA/dist)**12)) - ((LJ_SIGMA/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def calc_KE():
    k_e = 0
    for i in range(no_balls): 
        k_e = k_e + (0.5*masses[i]*(cp.dot(velocities[i], velocities[i]))) # use dot product as v^2 to calc ke
    return k_e

# calc Coulombs
def calc_coulombs(q1, q2, r, direction, k):

    if r < 1e-12:
        return cp.zeros(3, dtype=cp.float64)

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
    time_taken_exclusions += exclusions_end - exclusions_start

    return bonded_12, bonded_13, bonded_14

bonded_12, bonded_13, bonded_14 = build_exclusion_sets()    # builds the sets to be used

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

def coulombs_calculations(scale, q1, q2, r, r_hat, i, j):
    F_C = calc_coulombs(q1, q2, r, r_hat, COULOMB_CONSTANT)
    coulombs_force_calculated = scale*F_C
    forces[i] += coulombs_force_calculated                      # Newtons third law - apply to both atoms
    forces[j] -= coulombs_force_calculated
    global pot_e_total

    global real_space_PE
    #pot_e_total += (scale* COULOMB_CONSTANT*q1* q2*math.erfc(PME_ALPHA * r)/r)
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

def LJ_calculations(r, r_hat, scale, LJ_energy_shift, i, j):
    LJ_force_calculated = scale * calc_LJ_force(r, r_hat)
    forces[i] += LJ_force_calculated                      # Newtons third law - apply to both atoms
    forces[j] -= LJ_force_calculated
    global pot_e_total
    pot_e_total += scale * (4 * LJ_EPSILON * ((LJ_SIGMA / r)**12 - (LJ_SIGMA / r)**6) - LJ_energy_shift) # add the LJ potential energy on

def build_neighbour_lists(neighbour_cutoff):
    cutoff2 = neighbour_cutoff**2
    
    for atom in atoms:
        atom.neighbours.clear()         # clear out all of the previous neighbour lists

    cells = build_cell_list()

    for i, atom in enumerate(atoms):  
        my_cell = find_cell_of_atom(positions[i])

        for cell in neighbouring_cells(my_cell):
            for j in cells.get(cell, ()):
                if j <= i:
                    continue

                r_vec = system.minimum_image(positions[j] - positions[i])
                r2 = cp.dot(r_vec, r_vec)                          # distance between pairs without square rooting - more efficient

                if r2 > cutoff2:
                    continue
                
                atoms[i].neighbours.append(j)
                atoms[j].neighbours.append(i)

    # update the new last neighbour refference for each atom
    last_neighbour_reference[:] = positions

NUM_CELLS = max(3, int(PBC_BOX_LENGTH // CELL_SIZE))

def find_cell_of_atom(atom_pos):
    return (
        math.floor(atom_pos[0] / CELL_SIZE) % NUM_CELLS,
        math.floor(atom_pos[1] / CELL_SIZE) % NUM_CELLS,
        math.floor(atom_pos[2] / CELL_SIZE) % NUM_CELLS

    )

def build_cell_list():
    cells = {}

    for i, atom in enumerate(atoms):
        cell = find_cell_of_atom(positions[i])

        if cell not in cells:
            cells[cell] = []

        cells[cell].append(i)

    return cells

def neighbouring_cells(cell):
    x,y,z = cell

    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield ((x+dx) % NUM_CELLS, (y+dy) % NUM_CELLS, (z+dz) % NUM_CELLS)

# -------- PME --------

# Convert a VPython position into fractional mesh coordinates
def position_to_mesh(pos):
    half = PBC_BOX_LENGTH * 0.5
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

    theta = cp.zeros(BSPLINE_ORDER)
    dtheta = cp.zeros(BSPLINE_ORDER)

    for i in range(BSPLINE_ORDER):

        u = f + (BSPLINE_ORDER - 1 - i)

        theta[i] = bspline(BSPLINE_ORDER, u)
        dtheta[i] = bspline_derivative(BSPLINE_ORDER, u)

    return theta, dtheta

def build_charge_grid():

    rho.fill(0.0)
    spline_order = 4

    for i in range(no_balls):
        gx, gy, gz = position_to_mesh(positions[i])

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
                    ] += charges[i] * wx * wy * wz

    return rho

def build_reciprocal_kernel():

    global C

    C = cp.zeros(
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
                    cp.exp(
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
            cp.exp(
                -2j * math.pi * m * k / PME_GRID
            )
        )

    if abs(total) < 1e-12:
        return 0.0

    phase = cp.exp(
        -2j * math.pi * (BSPLINE_ORDER-1) * m / PME_GRID
    )

    b = phase / total

    return abs(b)**2

def build_B():

    global B

    B = cp.zeros(
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

def reciprocal_interpolation(potential):

    scale = PME_GRID / PBC_BOX_LENGTH

    for i in range(no_balls):
        gx, gy, gz = position_to_mesh(positions[i])

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

        energy = 0.5 * COULOMB_CONSTANT * charges[i] * Phi

        reciprocal_PE += energy
        pot_e_total += energy

        Fx = COULOMB_CONSTANT * charges[i] * scale * dPhidx
        Fy = COULOMB_CONSTANT * charges[i] * scale * dPhidy
        Fz = COULOMB_CONSTANT * charges[i] * scale * dPhidz
        forces[i] -= cp.array([Fx, Fy, Fz], dtype=cp.float64)

def coulomb_exclusion_correction(q1, q2, r, r_hat, i, j):
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
    forces[i] += F
    forces[j] -= F

    energy = -COULOMB_CONSTANT * q1 * q2 * erf_term / r
    real_space_PE += energy
    pot_e_total += energy

def apply_exclusion_corrections():
    for (i, j) in bonded_12 | bonded_13:
        r_vec = system.minimum_image(positions[j] - positions[i])
        r = cp.linalg.norm(r_vec)
        if r < 1e-12:
            continue
        r_hat = r_vec / r
        q1 = float(charges[i])
        q2 = float(charges[j])
        coulomb_exclusion_correction(q1, q2, r, r_hat, i, j)

PME_SELF_ENERGY = 0.0

for i in range(no_balls):
    PME_SELF_ENERGY += charges[i] * charges[i]

PME_SELF_ENERGY *= (
    -COULOMB_CONSTANT
    * PME_ALPHA
    / math.sqrt(math.pi)
)

# -------- Other Forces --------

def non_bonded_forces():
    LJ_scale = 1.0
    coulomb_scale = 1.0
    LJ_energy_shift = 4 * LJ_EPSILON * ((LJ_SIGMA / LJ_CUTOFF)**12 - (LJ_SIGMA / LJ_CUTOFF)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

    for i, atom in enumerate(atoms):                   # loop over each unique atom pair once
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

            r_vec = system.minimum_image(positions[j] - positions[i])     # find pair distances and update the "ghost molecules" 
            r = cp.linalg.norm(r_vec)
            r_hat = r_vec / cp.linalg.norm(r_vec)

            q1 = float(charges[i])             # get atom charges               
            q2 = float(charges[j])

            if r < 1e-12:
                continue

            #if r <= REAL_CUTOFF:             # skip invalid LJ, and only calculate coulomb forces
                #coulombs_calculations(coulomb_scale, q1, q2, r, r_hat, i, j)


            if r <= LJ_CUTOFF:                             
                LJ_calculations(r, r_hat, LJ_scale, LJ_energy_shift, i, j)

# caculate physics stuff - the heart <3
def calc_physics():

    global pot_e_total
    pot_e_total = 0.0
    bond_e = 0
    angle_e = 0
    torsion_e = 0

    global real_space_PE
    global reciprocal_PE

    real_space_PE = 0.0
    reciprocal_PE = 0.0

    start_all_foces = time.perf_counter()

    start_bonds = time.perf_counter()
    # bond forces - from py file
    bond_e = bond_forces(system)
    end_bonds = time.perf_counter()
    global time_taken_bonds
    time_taken_bonds += end_bonds - start_bonds

    start_bond_angles = time.perf_counter()
    # bond angles
    angle_e = angle_forces(system)
    end_bond_angles = time.perf_counter()
    global time_taken_bond_angles
    time_taken_bond_angles += end_bond_angles - start_bond_angles

    start_torsions = time.perf_counter()
    # torsion angles
    torsion_e = torsion_forces(system)
    end_torsions = time.perf_counter()
    global time_taken_torsions
    time_taken_torsions += end_torsions - start_torsions

    pot_e_total += bond_e + angle_e + torsion_e

    start_LJ = time.perf_counter()
    # #LJ forces and Coulombs (non bonded loop)
    non_bonded_forces()

    # ensure that coulombs also skip the excluded atoms
    #apply_exclusion_corrections()

    end_LJ = time.perf_counter()
    global time_taken_LJ
    time_taken_LJ += end_LJ - start_LJ

    # -----------------------
    # Reciprocal-space PME
    # -----------------------
    start_coulombs = time.perf_counter()

    rho = build_charge_grid()

    Qk = cp.fft.fftn(rho)
    Qk *= BC

    potential = cp.fft.ifftn(Qk).real
    potential *= PME_GRID**3

    reciprocal_interpolation(potential)

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
forces.fill(0)
pot_e = calc_physics()
render_initial_model()
build_neighbour_lists(NEIGHBOUR_CUTOFF)

step = 0
prev_step_count = 0
relax_steps = 2000
timestep_x = 1
average_total_energy = 0
max_displacement2 = 0


# running simulation
while True:
    start_total_time = time.perf_counter()
    update_camera(TIME_STEP)
    
    # lighting position for camera
    cam_light.index = scene.camera.pos - scene.camera.axis.norm()*3
    
    # Verlet integration method:
    # 1. half-step velocity update
    for i in range(no_balls):
        velocities[i] += 0.5 * (forces[i] / masses[i]) * TIME_STEP     # first half-step velocity update - uses current force to push velocity halfway forward


    # 2. position update
    for i in range(no_balls):
        positions[i] += velocities[i] * TIME_STEP                      # update pos

    # 3. wrap positions back into box
    system.wrap_molecules()

    # 4. reset forces
    for i in range(no_balls):
        forces[i] = cp.zeros(3, dtype=cp.float64)                             # clear old forces before computing new ones

    # 5. compute new forces
    pot_e = calc_physics()                                           # now get new forces at the new positions

    # 6. second half-step velocity update
    for i in range(no_balls):
        velocities[i] += 0.5 * (forces[i] / masses[i]) * TIME_STEP     # compleate  the full velocity update using the new forces

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
        velocities[i] *= damping                                 # dampens some of the velocity at each step when begining

    k_e = calc_KE()                                             # energy tracking
    total_e = k_e + pot_e
    average_total_energy += total_e
    step += 1




    # find the displacement of the atom this timestep
    for atom in atoms:
        disp2 = system.minimum_image(positions[atom.index] - last_neighbour_reference[atom.index])

        displacement2 = cp.dot(disp2,disp2)

        if displacement2 > max_displacement2:
            max_displacement2 = displacement2

    if max_displacement2 > (SKIN_CUTOFF / 2)**2:
        prev_step_count = step - prev_step_count
        print(f"New List After {prev_step_count} Steps")
        build_neighbour_lists(NEIGHBOUR_CUTOFF)
        max_displacement2 = 0
        
    start_total_graphics = time.perf_counter() 

    # ----- graphics render  ------

    # update graphics once per frame
    for i in range(no_balls):
        balls[i].pos = np_to_vec(positions[i])                                 # moves balls to current pos

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = np_to_vec(positions[bond.a1])
        bond_visuals[idx].axis = np_to_vec(positions[bond.a2] - positions[bond.a1])

    # ----------------------------

    end_total_graphics = time.perf_counter()
    time_taken_graphics += end_total_graphics - start_total_graphics

    end_total_time = time.perf_counter()
    time_taken_total += end_total_time - start_total_time

    time_other = time_taken_total - (time_taken_graphics + time_taken_all_foces)



    t0 = time.perf_counter()
    rho = build_charge_grid()
    t1 = time.perf_counter()

    Qk = cp.fft.fftn(rho)
    t2 = time.perf_counter()

    Qk *= BC
    t3 = time.perf_counter()

    potential = cp.fft.ifftn(Qk).real
    t4 = time.perf_counter()

    reciprocal_interpolation(potential)
    t5 = time.perf_counter()

    print("Grid:", t1-t0)
    print("FFT:", t2-t1)
    print("Multiply:", t3-t2)
    print("IFFT:", t4-t3)
    print("Interp:", t5-t4)



    # printing and debuging stuff that prints every x timesteps
    if step % timestep_x == 0:
        average_total_energy = average_total_energy/timestep_x
        steps.append(step)
        avg_sys_energy_values.append(average_total_energy)


        '''
        line.set_data(steps, avg_sys_energy_values)

        ax.relim()              # Recalculate limits
        ax.autoscale_view()     # Expand axes if needed
        plt.draw()
        plt.pause(0.001)
        average_total_energy = 0
        '''
        
        print()
        print(f"Step: {step}")
        print(f"KE: {k_e:.6f}  PE: {pot_e:.6f}  Total: {total_e:.6f}")
        print(f"Time Taken per cycle: {time_taken_total/timestep_x:.6f}s")
        print()

        print(f"Time By % of Graphics: {time_taken_graphics:.6f} s")
        print(f"Time By % of Bonds: {time_taken_bonds:.6f} s")
        print(f"Time By % of Bond Angles: {time_taken_bond_angles:.6f} s")
        print(f"Time By % of Torsions: {time_taken_torsions:.6f} s")
        print(f"Time By % of LJ: {time_taken_LJ:.6f} s")
        print(f"Time By % of Coulombs: {time_taken_coulombs:.6f} s")
        print(f"Time By % of non_bonded: {time_taken_non_bonded:.6f} s")
        print(f"Time By % of total forces: {time_taken_all_foces:.6f} s")
        print(f"Time By % of other: {time_other:.6f} s")

        '''
        print(f"Time By % of Graphics: {((time_taken_graphics/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bonds: {((time_taken_bonds/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bond Angles: {((time_taken_bond_angles/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Torsions: {((time_taken_torsions/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of LJ: {((time_taken_LJ/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Coulombs: {((time_taken_coulombs/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of non_bonded: {((time_taken_non_bonded/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of total forces: {((time_taken_all_foces/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of other: {((time_other/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        '''
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