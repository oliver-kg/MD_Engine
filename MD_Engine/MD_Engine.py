from vpython import *                   # imports all I use
import tkinter as tk
import math
import read_molecules
import time
from random import uniform
import numpy as np

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
    def __init__(self, pos, vel, force, charge, element):
        props = ELEMENTS[element]
        self.mass = props["mass"]
        self.pos = pos
        self.last_neighbour_reference = vector(0,0,0)
        self.vel = vel
        self.force = force
        self.charge = float(charge)
        self.radius = props["radius"]
        self.element = element
        self.base_colour = props["color"]
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


# molecule properties arrays
atoms = []
balls = []
bonds = []
bond_angles = []
torsion_angles = []
bond_visuals = []


real_space_PE = 0.0
reciprocal_PE = 0.0

molecules = [
     Molecule([0])
    ]

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
PBC_BOX_LENGTH = 32
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

        atoms.append(
            Atom(
                pos + offset_position,
                vector(0,0,0),
                vector(0,0,0),
                charge,
                atom_type_name
            )
        )
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

water = load_molecule_template()

for i in range(32):
    insert_molecule(water)

# create number of atoms after all atoms have been made
no_balls = len(atoms)

# create NumPy Vectors

positions = np.zeros((no_balls,3), dtype=np.float64)
velocities = np.zeros((no_balls,3), dtype=np.float64)
forces = np.zeros((no_balls,3), dtype=np.float64)

charges = np.zeros(no_balls)
masses = np.zeros(no_balls)

for i, atom in enumerate(atoms):
    positions[i] = [atom.pos.x, atom.pos.y, atom.pos.z]
    velocities[i] = [atom.vel.x, atom.vel.y, atom.vel.z]
    forces[i] = [0.0,0.0,0.0]
    charges[i] = atom.charge
    masses[i] = atom.mass

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

    if r == 0:
        return vector(0, 0, 0)

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
    half = PBC_BOX_LENGTH / 2
    if r_vector.x > half:
        r_vector.x = r_vector.x - PBC_BOX_LENGTH
    if r_vector.x < -half:
        r_vector.x = r_vector.x + PBC_BOX_LENGTH

    if r_vector.y > half:
        r_vector.y = r_vector.y - PBC_BOX_LENGTH
    if r_vector.y < -half:
        r_vector.y = r_vector.y + PBC_BOX_LENGTH

    if r_vector.z > half:
        r_vector.z = r_vector.z - PBC_BOX_LENGTH
    if r_vector.z < -half:
        r_vector.z = r_vector.z + PBC_BOX_LENGTH

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
    atoms[i].force += coulombs_force_calculated                      # Newtons third law - apply to both atoms
    atoms[j].force -= coulombs_force_calculated
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
    atoms[i].force += LJ_force_calculated                      # Newtons third law - apply to both atoms
    atoms[j].force -= LJ_force_calculated
    global pot_e_total
    pot_e_total += scale * (4 * LJ_EPSILON * ((LJ_SIGMA / r)**12 - (LJ_SIGMA / r)**6) - LJ_energy_shift) # add the LJ potential energy on

def build_neighbour_lists(neighbour_cutoff):
    cutoff2 = neighbour_cutoff**2
    
    for atom in atoms:
        atom.neighbours.clear()         # clear out all of the previous neighbour lists

    cells = build_cell_list()

    for i, atom in enumerate(atoms):  
        my_cell = find_cell_of_atom(atom.pos)

        for cell in neighbouring_cells(my_cell):
            for j in cells.get(cell, ()):
                if j <= i:
                    continue

                r_vec = PBC_box_for_vectors(atoms[j].pos - atoms[i].pos)
                r2 = mag2(r_vec)                          # distance between pairs without square rooting - more efficient

                if r2 > cutoff2:
                    continue
                
                atoms[i].neighbours.append(j)
                atoms[j].neighbours.append(i)

    # update the new last neighbour refference for each atom
    for atom in atoms:
        atom.last_neighbour_reference = atom.pos

def find_cell_of_atom(atom_pos):
    return (
        math.floor(atom_pos.x / CELL_SIZE),
        math.floor(atom_pos.y / CELL_SIZE),
        math.floor(atom_pos.z / CELL_SIZE)

    )

def build_cell_list():
    cells = {}

    for i, atom in enumerate(atoms):
        cell = find_cell_of_atom(atom.pos)

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
    x = (pos.x + half) / PBC_BOX_LENGTH
    y = (pos.y + half) / PBC_BOX_LENGTH
    z = (pos.z + half) / PBC_BOX_LENGTH

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

def build_charge_grid():

    rho.fill(0.0)
    spline_order = 4

    for atom in atoms:

        gx, gy, gz = position_to_mesh(atom.pos)

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
                    ] += atom.charge * wx * wy * wz

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

def reciprocal_interpolation(potential):

    scale = PME_GRID / PBC_BOX_LENGTH

    for atom in atoms:

        gx, gy, gz = position_to_mesh(atom.pos)

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

        energy = 0.5 * COULOMB_CONSTANT * atom.charge * Phi

        reciprocal_PE += energy
        pot_e_total += energy

        Fx = COULOMB_CONSTANT * atom.charge * scale * dPhidx
        Fy = COULOMB_CONSTANT * atom.charge * scale * dPhidy
        Fz = COULOMB_CONSTANT * atom.charge * scale * dPhidz
        atom.force -= vector(Fx, Fy, Fz)

PME_SELF_ENERGY = 0.0

for atom in atoms:
    PME_SELF_ENERGY += atom.charge * atom.charge

PME_SELF_ENERGY *= (
    -COULOMB_CONSTANT
    * PME_ALPHA
    / math.sqrt(math.pi)
)

# -------- Other Forces --------

def bond_forces():
    start_bonds = time.perf_counter()

    for bond in bonds:
        a = bond.a1
        b = bond.a2
        r0 = bond.ideal_dist
        k = bond.K

        # compute the bond vector and bond lengths
        r_vec = PBC_box_for_vectors(atoms[b].pos - atoms[a].pos)
        r = mag(r_vec)

        if r < 1e-12:                              # avoid devide by zero just in case
            continue

        r_hat = norm(r_vec)                     # bond direction

        F = 2 * k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds

        atoms[a].force += F                          # apply Newtons third law - equal and opposite forces
        atoms[b].force -= F

        global pot_e_total
        pot_e_total += k * (r - r0)**2                # calc the harmonic bond potential in function 2 type bond



    end_bonds = time.perf_counter()
    global time_taken_bonds
    time_taken_bonds += end_bonds - start_bonds

def angle_forces():
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

        BA = PBC_box_for_vectors(BA)                                                 # update the coppy vectors
        BC = PBC_box_for_vectors(BC)
        
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

        global pot_e_total
        pot_e_total += 0.5 * k * (theta - ideal_ang)**2                # update the potential energy change from the angle

    end_bond_angles = time.perf_counter()
    global time_taken_bond_angles
    time_taken_bond_angles += end_bond_angles - start_bond_angles

def torsion_forces():
    start_torsions = time.perf_counter()
    
    for angle in torsion_angles:
        a1 = angle.a1                           # 4 atoms involved in the torsion bond
        a2 = angle.a2
        a3 = angle.a3
        a4 = angle.a4


        b1 = atoms[a2].pos - atoms[a1].pos      # directions of the three bonds
        b2 = atoms[a3].pos - atoms[a2].pos
        b3 = atoms[a4].pos - atoms[a3].pos

        b1 = PBC_box_for_vectors(b1)                 # update the vectors of coppies
        b2 = PBC_box_for_vectors(b2)
        b3 = PBC_box_for_vectors(b3)

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

        global pot_e_total
        pot_e_total += torsion_e  # calculate the potential energy

    end_torsions = time.perf_counter()
    global time_taken_torsions
    time_taken_torsions += end_torsions - start_torsions

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

            r_vec = PBC_box_for_vectors(atoms[j].pos - atoms[i].pos)     # find pair distances and update the "ghost molecules" 
            r = mag(r_vec)
            r_hat = norm(r_vec)

            q1 = float(atoms[i].charge)             # get atom charges               
            q2 = float(atoms[j].charge)

            # -- add the coulombs back in when adding PME later

            if r < 1e-12 or r > REAL_CUTOFF:             # skip invalid LJ, and only calculate coulomb forces
                coulombs_calculations(coulomb_scale, q1, q2, r, r_hat, i, j)
                continue

            else:                                   # dont skip any non bonded forces as pair is in the cutoff region 
                LJ_calculations(r, r_hat, LJ_scale, LJ_energy_shift, i, j)
                coulombs_calculations(coulomb_scale, q1, q2, r, r_hat, i, j)

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
    bond_forces()

    # bond angles
    angle_forces()

    # torsion angles
    torsion_forces()

    start_LJ = time.perf_counter()

    # #LJ forces and Coulombs (non bonded loop)
    non_bonded_forces()

    end_LJ = time.perf_counter()
    global time_taken_LJ
    time_taken_LJ += end_LJ - start_LJ

    # -----------------------
    # Reciprocal-space PME
    # -----------------------
    start_coulombs = time.perf_counter()

    rho = build_charge_grid()

    Qk = np.fft.fftn(rho)
    Qk *= BC

    potential = np.fft.ifftn(Qk).real
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
for i in range(no_balls):
    atoms[i].force = vector(0, 0, 0)
pot_e = calc_physics()
render_initial_model()
build_neighbour_lists(NEIGHBOUR_CUTOFF)

step = 0
prev_step_count = 0
relax_steps = 2000
timestep_x = 50
average_total_energy = 0
max_displacement2 = 0

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
    average_total_energy += total_e
    step += 1

    # find the displacement of the atom this timestep
    for atom in atoms:
        disp2 = PBC_box_for_vectors(atom.pos - atom.last_neighbour_reference)

        displacement2 = mag2(disp2)

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
        balls[i].pos = atoms[i].pos                                 # moves balls to current pos

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = atoms[bond.a1].pos
        bond_visuals[idx].axis = atoms[bond.a2].pos - atoms[bond.a1].pos

    # ----------------------------

    end_total_graphics = time.perf_counter()
    time_taken_graphics += end_total_graphics - start_total_graphics

    end_total_time = time.perf_counter()
    time_taken_total += end_total_time - start_total_time

    time_other = time_taken_total - (time_taken_graphics + time_taken_all_foces)



    # printing and debuging stuff that prints every x timesteps
    if step % timestep_x == 0:
        average_total_energy = average_total_energy/timestep_x
        steps.append(step)
        avg_sys_energy_values.append(average_total_energy)

        # line.set_data(steps, avg_sys_energy_values)

        '''
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