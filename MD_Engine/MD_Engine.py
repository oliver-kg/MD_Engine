from vpython import *                   # imports all I use
import tkinter as tk
import math
import time
from random import uniform
import numpy as np
import cupy as cp
from cupyx.scipy.special import erfc, erf
import matplotlib.pyplot as plt

from system import System, Atom, Molecule, ELEMENTS
from constants import *
import read_molecules


# ------------------------
# Length = nm
# Mass = amu
# Time = ps
# Charge = e
# Energy = KL/mol
# ------------------------



system = System()

time_taken_all_foces = 0
time_taken_exclusions = 0
time_taken_bonds = 0
time_taken_bond_angles = 0
time_taken_torsions = 0
time_taken_non_bonded = 0
time_taken_LJ_and_coulombs = 0
time_taken_PME = 0
time_taken_graphics = 0
time_taken_total = 0
time_other = 0

# ---------------------------------------- SETUP ----------------------------------------

plt.ion()  # Interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlabel("Timestep")
ax.set_ylabel("Average Total System Energy")
ax.set_title("Average Total System Energy (x-step average)")
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

# Main light (camera light)
cam_light = local_light(
    pos=scene.camera.pos - scene.camera.axis.norm()*4,
    color=color.gray(0.65)
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
    for i, a in enumerate(system.atoms):
        b = sphere(pos=vector(*system.positions[i]), radius=a.radius, color=a.base_colour, make_trail=False)
        system.balls.append(b)

    # Draw lines between bonds
    for bond_index in range(len(system.bond_a)):
        a1 = system.bond_a[bond_index]
        a2 = system.bond_b[bond_index]
        c = cylinder(
            pos=vector(*system.positions[a1]),
            axis=vector(*system.positions[a2]) - vector(*system.positions[a1]),
            radius=0.01,
            color=vector(1,1,1)
        )
        system.bond_visuals.append(c)

def load_molecule_template():

    atom_pos = []
    atom_charge = []
    atom_type = []
    ff_atom_type = []

    for ind in read_molecules.atom_index:
        atom_pos.append(read_molecules.atom_pos[ind])
        atom_charge.append(read_molecules.atom_charge[ind])
        atom_type.append(read_molecules.atom_type[ind])
        ff_atom_type.append(read_molecules.ff_atom_type[ind])

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

    return atom_pos, atom_charge, atom_type, ff_atom_type, bonds, angles, torsions

def insert_molecule(template):

    atom_pos, atom_charge, atom_type, ff_atom_type, bond_data, angle_data, torsion_data = template

    # Where this copy will go
    offset_position = vector(
        uniform(0, PBC_BOX_LENGTH),
        uniform(0, PBC_BOX_LENGTH),
        uniform(0, PBC_BOX_LENGTH)
    )

# First new atom index
    atom_offset = len(system.atoms)

    atom_indices = []

    # Create atoms
    for pos, charge, atom_type_name, ff_type_name in zip(atom_pos, atom_charge, atom_type, ff_atom_type):

        index = len(system.atoms)

        system.atoms.append(
            Atom(
                index,
                atom_type_name
            )
        )

        new_pos = pos + offset_position

        system.positions.append([
        new_pos.x,
        new_pos.y,
        new_pos.z
    ])

        system.velocities.append([
            0.0,
            0.0,
            0.0
        ])

        system.forces.append([
            0.0,
            0.0,
            0.0
        ])

        system.charges.append(float(charge))

        system.ff_atom_types.append(ff_type_name)

        system.ff_type_ids.append(system.ff_type_to_id[ff_type_name])

        system.masses.append(
            ELEMENTS[atom_type_name]["mass"]
        )

        atom_indices.append(len(system.atoms)-1)

    system.molecules.append(Molecule(atom_indices))

    # Create bonds
    for a, b, r0, k in bond_data:

        global_a = a + atom_offset
        global_b = b + atom_offset

        system.bond_a.append(global_a)
        system.bond_b.append(global_b)
        system.bond_r0.append(r0)
        system.bond_k.append(k)

        # Determine whether this bond should be constrained
        atom_a_type = atom_type[a]
        atom_b_type = atom_type[b]

        is_h_bond = (
            atom_a_type == "H" or
            atom_b_type == "H"
        )

        system.bond_constrained.append(
            CONSTRAIN_H_BONDS and is_h_bond
        )

        if CONSTRAIN_H_BONDS and is_h_bond:

            system.constraint_a.append(global_a)
            system.constraint_b.append(global_b)
            system.constraint_r0.append(r0)

    # Create angles
    for i, j, k, theta0, k_theta in angle_data:

        system.angle_i.append(i + atom_offset)
        system.angle_j.append(j + atom_offset)
        system.angle_k.append(k + atom_offset)

        system.angle_theta0.append(math.radians(theta0))

        system.angle_kconst.append(k_theta)

    # Create torsions
    for i, j, k, l, k_dih, n, phase in torsion_data:

        system.torsion_i.append(i + atom_offset)
        system.torsion_j.append(j + atom_offset)
        system.torsion_k.append(k + atom_offset)
        system.torsion_l.append(l + atom_offset)

        system.torsion_kterm.append(k_dih)
        system.torsion_n.append(n)
        system.torsion_delta.append(math.radians(phase))


system.ff_type_to_id = read_molecules.ff_type_to_id
system.num_ff_types = read_molecules.num_ff_types

system.lj_c6_matrix = cp.asarray(
    read_molecules.lj_c6_matrix,
    dtype=cp.float64
)

system.lj_c12_matrix = cp.asarray(
    read_molecules.lj_c12_matrix,
    dtype=cp.float64
)


# Create test molecules
water = load_molecule_template()

for i in range(MOLECULE_NUMBERS):
    insert_molecule(water)


# create Vectors
system.positions = cp.asarray(
    system.positions,
    dtype=cp.float64
)

system.velocities = cp.asarray(
    system.velocities,
    dtype=cp.float64
)

system.forces = cp.asarray(
    system.forces,
    dtype=cp.float64
)

system.charges = cp.asarray(
    system.charges,
    dtype=cp.float64
)

system.masses = cp.asarray(
    system.masses,
    dtype=cp.float64
)


lj_c6_cpu = read_molecules.lj_c6_matrix
lj_c12_cpu = read_molecules.lj_c12_matrix

lj_shift_matrix = np.zeros(
    (read_molecules.num_ff_types, read_molecules.num_ff_types),
    dtype=np.float64
)

for i in range(system.num_ff_types):
    for j in range(system.num_ff_types):

        C6 = lj_c6_cpu[i, j]
        C12 = lj_c12_cpu[i, j]

        lj_shift_matrix[i, j] = (
            C12 / LJ_CUTOFF**12
            - C6 / LJ_CUTOFF**6
        )

system.lj_shift_matrix = cp.asarray(
    lj_shift_matrix,
    dtype=cp.float64
)

system.ff_type_ids = cp.asarray(
    system.ff_type_ids,
    dtype=cp.int32
)
# Topology
system.n_atoms = len(system.positions)

# Bonds
system.bond_a = cp.asarray(system.bond_a, dtype=cp.int32)
system.bond_b = cp.asarray(system.bond_b, dtype=cp.int32)

system.bond_r0 = cp.asarray(system.bond_r0, dtype=cp.float64)
system.bond_k = cp.asarray(system.bond_k, dtype=cp.float64)

system.bond_constrained = cp.asarray(system.bond_constrained, dtype=cp.bool_)
system.constraint_a = cp.asarray(system.constraint_a, dtype=cp.int32)
system.constraint_b = cp.asarray(system.constraint_b, dtype=cp.int32)
system.constraint_r0 = cp.asarray(system.constraint_r0, dtype=cp.float64)

# Angles
system.angle_i = cp.asarray(system.angle_i, dtype=cp.int32)
system.angle_j = cp.asarray(system.angle_j, dtype=cp.int32)
system.angle_k = cp.asarray(system.angle_k, dtype=cp.int32)
system.angle_theta0 = cp.asarray(system.angle_theta0, dtype=cp.float64)
system.angle_kconst = cp.asarray(system.angle_kconst, dtype=cp.float64)

# Torsions
system.torsion_i = cp.asarray(system.torsion_i, dtype=cp.int32)
system.torsion_j = cp.asarray(system.torsion_j, dtype=cp.int32)
system.torsion_k = cp.asarray(system.torsion_k, dtype=cp.int32)
system.torsion_l = cp.asarray(system.torsion_l, dtype=cp.int32)
system.torsion_kterm = cp.asarray(system.torsion_kterm,dtype=cp.float64)
system.torsion_n = cp.asarray(system.torsion_n,dtype=cp.int32)
system.torsion_delta = cp.asarray(system.torsion_delta,dtype=cp.float64)

system.n_bonds = len(system.bond_a)
system.n_angles = len(system.angle_i)
system.n_torsions = len(system.torsion_i)



# Checks if a molecule has reached the box boundery and needs warping
def wrap_molecules(system):
    half = PBC_BOX_LENGTH / 2

    for mol in system.molecules:
        ref_pos = system.positions[mol.atom_indices[0]]
        shift = cp.zeros(3)

        for d in range(3):                              # loop over x, y, z instead of one if/elif each
            while ref_pos[d] + shift[d] > half:
                shift[d] -= PBC_BOX_LENGTH
            while ref_pos[d] + shift[d] < -half:
                shift[d] += PBC_BOX_LENGTH

        if shift[0] != 0 or shift[1] != 0 or shift[2] != 0:
            for i in mol.atom_indices:
                system.positions[i] += shift

# Simple camera movement via key press
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

# SHAKE-like constraints on H-X bonds to help increase dt
def rattle_position_constraints(system, old_positions, free_half_velocities):

    if system.constraint_a.size == 0:
        return free_half_velocities

    dt = TIME_STEP

    # Start from the unconstrained Verlet trial position
    trial_positions = old_positions + free_half_velocities * dt
    system.positions[:] = trial_positions

    # Bond vectors at q_n - these are the gradients used by the RATTLE position correction.
    old_bond_vectors = []

    for n in range(system.constraint_a.size):
        a = int(system.constraint_a[n])
        b = int(system.constraint_b[n])

        r_vec = system.minimum_image(
            old_positions[b] - old_positions[a]
        )

        old_bond_vectors.append(r_vec)

    for _ in range(CONSTRAINT_MAX_ITERATIONS):

        max_error = 0.0

        for n in range(system.constraint_a.size):

            a = int(system.constraint_a[n])
            b = int(system.constraint_b[n])
            r0 = system.constraint_r0[n]

            r_old = old_bond_vectors[n]

            # Current trial bond vector
            r = system.minimum_image(
                system.positions[b] - system.positions[a]
            )

            r2 = cp.dot(r, r)
            target2 = r0 * r0

            # Constraint:
            # g(q) = 1/2 (r^2 - r0^2)
            g = 0.5 * (r2 - target2)

            error = abs(float(g))
            max_error = max(max_error, error)

            if error <= CONSTRAINT_TOLERANCE:
                continue

            ma = system.masses[a]
            mb = system.masses[b]

            inv_mass_sum = (
                1.0 / ma +
                1.0 / mb
            )

            # Derivative of constraint wrt lambda
            denominator = (
                dt * dt
                * inv_mass_sum
                * cp.dot(r, r_old)
            )

            if abs(float(denominator)) < 1e-20:
                raise RuntimeError(
                    f"RATTLE position solve became singular for "
                    f"bond {a}-{b}."
                )

            # Newton correction to the constraint multiplier
            delta_lambda = -g / denominator

            # Position corrections
            delta_a = (
                -dt * dt
                * delta_lambda
                * r_old
                / ma
            )

            delta_b = (
                dt * dt
                * delta_lambda
                * r_old
                / mb
            )

            system.positions[a] += delta_a
            system.positions[b] += delta_b

        if max_error <= CONSTRAINT_TOLERANCE:
            break

    else:
        raise RuntimeError(
            "RATTLE position constraints failed to converge. "
            f"Maximum squared-distance constraint error = "
            f"{max_error:.3e}"
        )

    constrained_half_velocities = (
        system.positions - old_positions
    ) / dt

    return constrained_half_velocities

def rattle_velocity_constraints(system):

    if system.constraint_a.size == 0:
        return

    dt = TIME_STEP

    for _ in range(CONSTRAINT_MAX_ITERATIONS):

        max_error = 0.0

        for n in range(system.constraint_a.size):

            a = int(system.constraint_a[n])
            b = int(system.constraint_b[n])

            r_vec = system.minimum_image(
                system.positions[b] - system.positions[a]
            )

            r2 = cp.dot(r_vec, r_vec)

            if float(r2) < 1e-20:
                raise RuntimeError(
                    f"RATTLE velocity bond {a}-{b} has zero length."
                )

            v_rel = (
                system.velocities[b]
                - system.velocities[a]
            )

            # Velocity constraint:
            # r · v_rel = 0
            constraint = cp.dot(r_vec, v_rel)

            error = abs(float(constraint))
            max_error = max(max_error, error)

            if error <= CONSTRAINT_VELOCITY_TOLERANCE:
                continue

            ma = system.masses[a]
            mb = system.masses[b]

            inv_mass_sum = (
                1.0 / ma +
                1.0 / mb
            )

            lambda_value = (
                -constraint
                / (dt * inv_mass_sum * r2)
            )

            delta_va = (
                -dt
                * lambda_value
                * r_vec
                / ma
            )

            delta_vb = (
                dt
                * lambda_value
                * r_vec
                / mb
            )

            system.velocities[a] += delta_va
            system.velocities[b] += delta_vb

        if max_error <= CONSTRAINT_VELOCITY_TOLERANCE:
            break

    else:
        raise RuntimeError(
            "RATTLE velocity constraints failed to converge. "
            f"Maximum r·v constraint error = "
            f"{max_error:.3e}"
        )
    
# ---------------------------------------- NEIGHBOURS ----------------------------------------


# Builds the neibour list structure
def build_neighbours():
    bond_a_cpu = cp.asnumpy(system.bond_a)
    bond_b_cpu = cp.asnumpy(system.bond_b)

    neighbours = {i: set() for i in range(system.n_atoms)}        # set makes sure no duplicates and fast to check

    for bond_index in range(len(system.bond_a)):
        a1 = bond_a_cpu[bond_index]
        a2 = bond_b_cpu[bond_index]
        neighbours[a1].add(a2)                    # bonding works both ways so adds the bond to the specific atom you are looking at at the moment
        neighbours[a2].add(a1)

    return neighbours

# Filters out all the bonds into 1-2, 1-3 or 1-4 bonds
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

def build_neighbour_lists(system, neighbour_cutoff):

    cutoff2 = neighbour_cutoff**2
    positions_cpu = cp.asnumpy(system.positions)
    # clear out all of the previous neighbour lists

    pair_i = []
    pair_j = []

    pair_lj_scale = []
    pair_coulomb_scale = []

    cells = build_cell_list(system, positions_cpu)

    old_pairs = set(
        zip(
            cp.asnumpy(system.pair_i).tolist(),
            cp.asnumpy(system.pair_j).tolist()
        )
    )



    for i in range(system.n_atoms):  
        my_cell = find_cell_of_atom(positions_cpu[i])

        for cell in neighbouring_cells(my_cell):
            for j in cells.get(cell, ()):
                if j <= i:
                    continue

                r_vec = system.minimum_image_cpu(positions_cpu[j] - positions_cpu[i])
                r2 = np.dot(r_vec, r_vec)                          # distance between pairs without square rooting - more efficient

                if r2 > cutoff2:
                    continue
                
                pair = (i, j)

                if pair in bonded_12 or pair in bonded_13:
                    continue

                pair_i.append(i)
                pair_j.append(j)

                if pair in bonded_14:
                    pair_lj_scale.append(0.5)
                else:
                    pair_lj_scale.append(1.0)

                pair_coulomb_scale.append(1.0)

    # update the new last neighbour refference for each atom
    for i, atom in enumerate(system.atoms):
        atom.last_neighbour_reference = cp.asnumpy(system.positions[i]).copy()

    system.pair_i = cp.asarray(pair_i, dtype=cp.int32)
    system.pair_j = cp.asarray(pair_j, dtype=cp.int32)

    system.pair_lj_scale = cp.asarray(pair_lj_scale, dtype=np.float64)
    system.pair_coulomb_scale = cp.asarray(pair_coulomb_scale, dtype=np.float64)


    new_pairs = set(
        zip(
            cp.asnumpy(pair_i).tolist(),
            cp.asnumpy(pair_j).tolist()
        )
    )

    added_pairs = sorted(new_pairs - old_pairs)
    removed_pairs = sorted(old_pairs - new_pairs)

    with open(DEBUG_LOG_FILE, "a", encoding="utf-8") as f:

        f.write("\n")
        f.write(f"NEIGHBOUR LIST REBUILD AT STEP {system.step}\n")
        f.write(f"OLD PAIRS: {len(old_pairs)}\n")
        f.write(f"NEW PAIRS: {len(new_pairs)}\n")
        f.write(f"ADDED: {len(added_pairs)}\n")
        f.write(f"REMOVED: {len(removed_pairs)}\n")

        for i, j in added_pairs:

            rij = system.minimum_image(
                system.positions[j] - system.positions[i]
            )

            r = float(cp.linalg.norm(rij).get())

            f.write(
                f"  ADDED {i}-{j} r={r:.9f} nm\n"
            )

        for i, j in removed_pairs:

            rij = system.minimum_image(
                system.positions[j] - system.positions[i]
            )

            r = float(cp.linalg.norm(rij).get())

            f.write(
                f"  REMOVED {i}-{j} r={r:.9f} nm\n"
            )

def find_cell_of_atom(atom_pos):

    # Convert coordinates from [-L/2, L/2)
    # into [0, L)
    x = atom_pos[0] + HALF_BOX
    y = atom_pos[1] + HALF_BOX
    z = atom_pos[2] + HALF_BOX

    return (
        int(math.floor(x / CELL_SIZE)) % NUM_CELLS,
        int(math.floor(y / CELL_SIZE)) % NUM_CELLS,
        int(math.floor(z / CELL_SIZE)) % NUM_CELLS
    )


def build_cell_list(system, positions_cpu):

    cells = {}

    for i in range(system.n_atoms):

        cell = find_cell_of_atom(positions_cpu[i])

        if cell not in cells:
            cells[cell] = []

        cells[cell].append(i)

    return cells

def neighbouring_cells(cell):

    x, y, z = cell

    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):

                yield (
                    (x + dx) % NUM_CELLS,
                    (y + dy) % NUM_CELLS,
                    (z + dz) % NUM_CELLS
                )


# ---------------------------------------- NON-BONDED CALCS ----------------------------------------


# Lennard Jones Force Equation (VDW's)
def calc_LJ_force(dist, direction):
    return  -((24*LJ_EPSILON*((2*((LJ_SIGMA/dist)**12)) - ((LJ_SIGMA/dist)**6)))*(1/dist))*direction

# Calc Kintetic energy of particles
def calc_KE_old(system):
    system.kinetic_energy = 0
    for i in range(system.n_atoms): 
        system.kinetic_energy = system.kinetic_energy + (0.5*system.masses[i]*(cp.dot(system.velocities[i], system.velocities[i]))) # use dot product as v^2 to calc ke
    return system.kinetic_energy

def calc_KE(system):
    speed2 = cp.sum(system.velocities**2, axis=1)

    ke_per_atom = 0.5 * system.masses * speed2

    system.kinetic_energy = cp.sum(ke_per_atom)

    return system.kinetic_energy 

# Calc Coulombs
def calc_coulombs(q1, q2, r, direction, k):

    if r == 0:
        return cp.zeros(3)

    erfc_term = erfc(PME_ALPHA * r)

    exp_term = cp.exp(-(PME_ALPHA * r) ** 2)

    force_mag = (
        k
        * q1
        * q2
        * (
            erfc_term / (r * r)
            +
            (2.0 * PME_ALPHA / cp.sqrt(cp.pi))
            * exp_term / r
        )
    )

    return -force_mag * direction

def coulombs_calculations(system, scale, q1, q2, r, r_hat, i, j):
    F_C = calc_coulombs(q1, q2, r, r_hat, COULOMB_CONSTANT)
    coulombs_force_calculated = scale*F_C
    system.forces[i] += coulombs_force_calculated                      # Newtons third law - apply to both atoms
    system.forces[j] -= coulombs_force_calculated

    energy = (
        scale
        * COULOMB_CONSTANT
        * q1
        * q2
        * math.erfc(PME_ALPHA * r)
        / r
    )

    system.real_space_PE += energy
    system.potential_energy += energy

def LJ_calculations(system, r, r_hat, scale, LJ_energy_shift, i, j):
    LJ_force_calculated = scale * calc_LJ_force(r, r_hat)
    print(i, j, r, LJ_force_calculated)
    system.forces[i] += LJ_force_calculated                      # Newtons third law - apply to both atoms
    system.forces[j] -= LJ_force_calculated
    system.potential_energy += scale * (4 * LJ_EPSILON * ((LJ_SIGMA / r)**12 - (LJ_SIGMA / r)**6) - LJ_energy_shift) # add the LJ potential energy on


# ---------------------------------------- PME ----------------------------------------


# Convert all particle positions into fractional mesh coordinates
def position_to_mesh(system):

    half = PBC_BOX_LENGTH * 0.5

    gx = ((system.positions[:,0] + half) / PBC_BOX_LENGTH) * PME_GRID
    gy = ((system.positions[:,1] + half) / PBC_BOX_LENGTH) * PME_GRID
    gz = ((system.positions[:,2] + half) / PBC_BOX_LENGTH) * PME_GRID

    return gx, gy, gz

def bspline(n, u):

    u = cp.atleast_1d(cp.asarray(u, dtype=cp.float64))

    spline = cp.zeros_like(u)

    factorial = math.factorial(n - 1)

    for k in range(n + 1):

        coeff = ((-1) ** k) * math.comb(n, k)

        spline += (
            coeff *
            cp.maximum(u - k, 0.0) ** (n - 1)
        )

    spline /= factorial

    return spline

def bspline_cpu(n, u):

    u = np.atleast_1d(np.asarray(u, dtype=np.float64))

    spline = np.zeros_like(u)

    factorial = math.factorial(n - 1)

    for k in range(n + 1):

        coeff = ((-1) ** k) * math.comb(n, k)

        spline += (
            coeff *
            np.maximum(u - k, 0.0) ** (n - 1)
        )

    spline /= factorial

    return spline

def bspline_derivative(n, u):

    u = cp.atleast_1d(cp.asarray(u, dtype=cp.float64))

    if n <= 1:
        return cp.zeros_like(u)

    return (
        bspline(n - 1, u)
        - bspline(n - 1, u - 1)
    )

def compute_theta(f):

    f = cp.atleast_1d(cp.asarray(f, dtype=cp.float64))

    theta = cp.empty((f.size, BSPLINE_ORDER), dtype=cp.float64)
    dtheta = cp.empty((f.size, BSPLINE_ORDER), dtype=cp.float64)

    for i in range(BSPLINE_ORDER):

        u = f + (BSPLINE_ORDER - 1 - i)

        theta[:, i] = bspline(BSPLINE_ORDER, u)
        dtheta[:, i] = bspline_derivative(BSPLINE_ORDER, u)

    return theta, dtheta

def build_charge_grid(system):

    system.rho.fill(0.0)

    gx, gy, gz = position_to_mesh(system)

    # Lower-left grid point
    ix = cp.floor(gx).astype(cp.int32)
    iy = cp.floor(gy).astype(cp.int32)
    iz = cp.floor(gz).astype(cp.int32)

    # Fractional offsets inside the cell
    fx = gx - ix
    fy = gy - iy
    fz = gz - iz

    theta_x, _ = compute_theta(fx)
    theta_y, _ = compute_theta(fy)
    theta_z, _ = compute_theta(fz)

    # Spread charge over a 4×4×4 cube
    for dx in range(BSPLINE_ORDER):

        wx = theta_x[:, dx]

        for dy in range(BSPLINE_ORDER):

            wy = theta_y[:, dy]

            for dz in range(BSPLINE_ORDER):

                wz = theta_z[:, dz]

                cp.add.at(
                    system.rho,
                    (
                        (ix + dx) % PME_GRID,
                        (iy + dy) % PME_GRID,
                        (iz + dz) % PME_GRID
                    ),
                    system.charges * wx * wy * wz
                )

    return system.rho

def build_reciprocal_kernel():

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

def b_factor(m):

    total = 0j

    for k in range(BSPLINE_ORDER - 1):

        total += (
            bspline_cpu(BSPLINE_ORDER, k + 1).item()
            *
            np.exp(
                -2j * math.pi * m * k / PME_GRID
            )
        )

    if np.abs(total) < 1e-12:
        return 0.0

    phase = np.exp(
        -2j * math.pi * (BSPLINE_ORDER-1) * m / PME_GRID
    )

    b = phase / total

    return abs(b)**2

def build_B():

    Bx = np.empty(PME_GRID)

    for i in range(PME_GRID):

        m = i if i <= PME_GRID//2 else i-PME_GRID
        Bx[i] = b_factor(m)

    B = (
        Bx[:,None,None]
        * Bx[None,:,None]
        * Bx[None,None,:]
    )

    return B

system.BC = cp.asarray(
    build_B()
    *
    build_reciprocal_kernel()
)

def reciprocal_interpolation(system, potential):

    scale = PME_GRID / PBC_BOX_LENGTH

    # Compute mesh coordinates for every atom
    gx, gy, gz = position_to_mesh(system)

    ix = cp.floor(gx).astype(cp.int32)
    iy = cp.floor(gy).astype(cp.int32)
    iz = cp.floor(gz).astype(cp.int32)

    fx = gx - ix
    fy = gy - iy
    fz = gz - iz

    theta_x, dtheta_x = compute_theta(fx)
    theta_y, dtheta_y = compute_theta(fy)
    theta_z, dtheta_z = compute_theta(fz)



    Phi = cp.zeros(system.n_atoms, dtype=cp.float64)

    dPhidx = cp.zeros(system.n_atoms, dtype=cp.float64)
    dPhidy = cp.zeros(system.n_atoms, dtype=cp.float64)
    dPhidz = cp.zeros(system.n_atoms, dtype=cp.float64)

    for dx in range(BSPLINE_ORDER):

        wx = theta_x[:, dx]
        dwx = dtheta_x[:, dx]

        for dy in range(BSPLINE_ORDER):

            wy = theta_y[:, dy]
            dwy = dtheta_y[:, dy]

            for dz in range(BSPLINE_ORDER):

                wz = theta_z[:, dz]
                dwz = dtheta_z[:, dz]

                phi = cp.real(
                    potential[
                        (ix + dx) % PME_GRID,
                        (iy + dy) % PME_GRID,
                        (iz + dz) % PME_GRID
                    ]
                )

                Phi += phi * wx * wy * wz

                dPhidx += phi * dwx * wy * wz
                dPhidy += phi * wx * dwy * wz
                dPhidz += phi * wx * wy * dwz


    energy = 0.5 * COULOMB_CONSTANT * system.charges * Phi

    system.reciprocal_PE += cp.sum(energy)
    system.potential_energy += cp.sum(energy)

    Fx = COULOMB_CONSTANT * system.charges * scale * dPhidx
    Fy = COULOMB_CONSTANT * system.charges * scale * dPhidy
    Fz = COULOMB_CONSTANT * system.charges * scale * dPhidz

    system.forces -= cp.stack((Fx, Fy, Fz), axis=1)

PME_SELF_ENERGY = 0.0
PME_SELF_ENERGY += cp.sum(system.charges * system.charges)

PME_SELF_ENERGY *= (
    -COULOMB_CONSTANT
    * PME_ALPHA
    / math.sqrt(math.pi)
)

def coulomb_exclusion_correction(system, q1, q2, r, r_hat, i, j):
    if r < 1e-12:
        return
    erf_term = erf(PME_ALPHA * r)
    exp_term = cp.exp(-(PME_ALPHA * r) ** 2)

    # cancels the erf(r)/r contribution the reciprocal grid implicitly added for this excluded pair
    force_mag = COULOMB_CONSTANT * q1 * q2 * (
        erf_term / (r * r) - (2.0 * PME_ALPHA / cp.sqrt(math.pi)) * exp_term / r
    )
    F = force_mag * r_hat
    system.forces[i] += F
    system.forces[j] -= F

    energy = -COULOMB_CONSTANT * q1 * q2 * erf_term / r
    system.real_space_PE += energy
    system.potential_energy += energy

def apply_exclusion_corrections(system):
    for (i, j) in bonded_12 | bonded_13:
        r_vec = system.minimum_image(system.positions[j] - system.positions[i])
        r = cp.linalg.norm(r_vec)
        if r < 1e-12:
            continue
        r_hat = r_vec / r
        q1 = system.charges[i]
        q2 = system.charges[j]
        coulomb_exclusion_correction(system, q1, q2, r, r_hat, i, j)


# ---------------------------------------- OTHER FORCES ----------------------------------------


def bond_forces(system):

    start_bonds = time.perf_counter()

    # Only flexible bonds receive harmonic forces
    flexible = ~system.bond_constrained

    if not cp.any(flexible):
        return

    a = system.bond_a[flexible]
    b = system.bond_b[flexible]
    r0 = system.bond_r0[flexible]
    k = system.bond_k[flexible]


    # Compute the bond vector and bond lengths
    r_vec = system.minimum_image((system.positions[b] - system.positions[a]))
    r = cp.linalg.norm(r_vec, axis=1)

    valid = r > 1e-12                              # avoid devide by zero just in case

    a = a[valid]
    b = b[valid]

    r = r[valid]
    r0 = r0[valid]
    k = k[valid]

    r_vec = r_vec[valid]

    r_hat = r_vec / r[:, None]                           # bond direction

    F = 2 * k[:, None] * (r - r0)[:, None] * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds

    cp.add.at(system.forces, a,  F)                          # apply Newtons third law - equal and opposite forces
    cp.add.at(system.forces, b, -F)
    system.potential_energy += cp.sum(k * (r - r0)**2)                # calc the harmonic bond potential in function 2 type bond



    end_bonds = time.perf_counter()
    global time_taken_bonds
    time_taken_bonds += end_bonds - start_bonds

def angle_forces(system):
    start_bond_angles = time.perf_counter()

    with open("cuda/angle_kernel.cu") as f:
        code = f.read()

    angle_kernel = cp.RawKernel(
        code,
        "angle_kernel"
    )

    threads = 256
    blocks = (system.n_angles + threads - 1) // threads

    system.potential_energy_gpu.fill(0)

    angle_kernel(
    (blocks,),
    (threads,),
    (
        system.positions,
        system.forces,
        system.angle_i,
        system.angle_j,
        system.angle_k,
        system.angle_theta0,
        system.angle_kconst,
        cp.float64(PBC_BOX_LENGTH),
        system.potential_energy_gpu,
        system.n_angles
    )
)
    cp.cuda.runtime.deviceSynchronize()

    system.potential_energy += float(system.potential_energy_gpu[0])

    end_bond_angles = time.perf_counter()
    global time_taken_bond_angles
    time_taken_bond_angles += end_bond_angles - start_bond_angles

def torsion_forces(system):
    start_torsions = time.perf_counter()
    
    with open("cuda/torsion_kernel.cu") as f:
        code = f.read()

    torsion_kernel = cp.RawKernel(
        code,
        "torsion_kernel"
    )

    threads = 256
    blocks = (system.n_torsions + threads - 1) // threads

    system.potential_energy_gpu.fill(0)

    torsion_kernel(
    (blocks,),
    (threads,),
    (
        system.positions,
        system.forces,
        system.torsion_i,
        system.torsion_j,
        system.torsion_k,
        system.torsion_l,
        system.torsion_kterm,
        system.torsion_n,
        system.torsion_delta,
        cp.float64(PBC_BOX_LENGTH),
        system.potential_energy_gpu,
        system.n_torsions
    )
)
    cp.cuda.runtime.deviceSynchronize()

    system.potential_energy += float(system.potential_energy_gpu[0])

    end_torsions = time.perf_counter()
    global time_taken_torsions
    time_taken_torsions += end_torsions - start_torsions

def non_bonded_forces(system):

    closest_r = float("inf")

    system.debug_max_lj_force = 0.0
    system.debug_max_lj_pair = None
    system.debug_max_lj_r = 0.0

    system.debug_max_coul_force = 0.0
    system.debug_max_coul_pair = None
    system.debug_max_coul_r = 0.0


    lj_PE = 0.0

    pair_i = system.pair_i
    pair_j = system.pair_j

    # Get force-field type IDs for every neighbour pair
    type_i = system.ff_type_ids[pair_i]
    type_j = system.ff_type_ids[pair_j]

    # Get LJ parameters for every neighbour pair
    pair_C6 = system.lj_c6_matrix[type_i, type_j]
    pair_C12 = system.lj_c12_matrix[type_i, type_j]
    pair_shift = system.lj_shift_matrix[type_i, type_j]

    ri = system.positions[pair_i]
    rj = system.positions[pair_j]

    qi = system.charges[pair_i]
    qj = system.charges[pair_j]

    lj_scale = system.pair_lj_scale
    coulomb_scale = system.pair_coulomb_scale

    r_vec = system.minimum_image(rj - ri)

    r = cp.linalg.norm(r_vec, axis=1)

    lj_candidate = (pair_C6 > 0) | (pair_C12 > 0)

    if cp.any(lj_candidate):

        lj_indices = cp.where(lj_candidate)[0]

        closest_lj_local = cp.argmin(r[lj_candidate])
        closest = lj_indices[closest_lj_local]

        closest_i = int(pair_i[closest])
        closest_j = int(pair_j[closest])
        closest_r = float(r[closest])

        closest_C6 = float(pair_C6[closest])
        closest_C12 = float(pair_C12[closest])
        closest_shift = float(pair_shift[closest])

    if system.step % 10 == 0 and closest_r < LJ_CUTOFF:

        lj_force_debug = (
            12.0 * closest_C12 / closest_r**13
            - 6.0 * closest_C6 / closest_r**7
        )

        lj_energy_debug = (
            closest_C12 / closest_r**12
            - closest_C6 / closest_r**6
            - closest_shift
        )

        print(
            f"STEP {system.step} | "
            f"PAIR {closest_i}-{closest_j} | "
            f"TYPES {system.ff_atom_types[closest_i]}-"
            f"{system.ff_atom_types[closest_j]} | "
            f"r={closest_r:.6f} | "
            f"C6={closest_C6:.6e} | "
            f"C12={closest_C12:.6e} | "
            f"SHIFT={closest_shift:.6e} | "
            f"FLJ={lj_force_debug:.3e} | "
            f"LJ_E={lj_energy_debug:.3e}"
        )

    valid = r > 1e-12

    pair_i = pair_i[valid]
    pair_j = pair_j[valid]

    pair_C6 = pair_C6[valid]
    pair_C12 = pair_C12[valid]
    pair_shift = pair_shift[valid]

    qi = qi[valid]
    qj = qj[valid]

    lj_scale = lj_scale[valid]
    coulomb_scale = coulomb_scale[valid]

    # LJ Force
    r_vec = r_vec[valid]
    r = r[valid]

    r_hat = r_vec / r[:, None]

    lj_mask = r <= LJ_CUTOFF

    if cp.any(lj_mask):

        lj_r = r[lj_mask]
        lj_hat = r_hat[lj_mask]

        lj_i = pair_i[lj_mask]
        lj_j = pair_j[lj_mask]

        lj_scale_local = lj_scale[lj_mask]
        shift_local = pair_shift[lj_mask]

        C6_local = pair_C6[lj_mask]
        C12_local = pair_C12[lj_mask]

        inv_r = 1.0 / lj_r

        inv_r2 = inv_r * inv_r
        inv_r6 = inv_r2 * inv_r2 * inv_r2
        inv_r7 = inv_r6 * inv_r

        inv_r12 = inv_r6 * inv_r6
        inv_r13 = inv_r12 * inv_r

        force_mag = (
            12.0 * C12_local * inv_r13
            - 6.0 * C6_local * inv_r7
        )

        force_mag *= lj_scale_local

        F_lj = -force_mag[:, None] * lj_hat

        # ---------------------------------------------------------
        # DEBUG: largest LJ pair force
        # ---------------------------------------------------------
        lj_force_abs = cp.abs(force_mag)

        if cp.any(lj_force_abs > 0):
            idx = int(cp.argmax(lj_force_abs))

            max_lj_force = float(lj_force_abs[idx].get())

            if max_lj_force > system.debug_max_lj_force:
                system.debug_max_lj_force = max_lj_force
                system.debug_max_lj_pair = (
                    int(lj_i[idx]),
                    int(lj_j[idx])
                )
                system.debug_max_lj_r = float(lj_r[idx].get())

        cp.add.at(system.forces, lj_i,  F_lj)
        cp.add.at(system.forces, lj_j, -F_lj)

        U_lj = (
            C12_local * inv_r12
            - C6_local * inv_r6
            - shift_local
        )

        current_lj_PE = cp.sum(
            lj_scale_local * U_lj
        )

        system.potential_energy += current_lj_PE
        lj_PE += float(current_lj_PE.get())
    
    # Coulomb Force
    
    real_mask = r <= REAL_CUTOFF

    if cp.any(real_mask):

        real_r = r[real_mask]
        real_hat = r_hat[real_mask]

        real_i = pair_i[real_mask]
        real_j = pair_j[real_mask]

        q1 = qi[real_mask]
        q2 = qj[real_mask]

        scale = coulomb_scale[real_mask]

        alpha_r = PME_ALPHA * real_r

        erfc_term = erfc(alpha_r)
        exp_term = cp.exp(-(alpha_r**2))

        force_mag = (
            COULOMB_CONSTANT
            * q1
            * q2
            * (
                erfc_term / (real_r**2)
                +
                (2.0 * PME_ALPHA / cp.sqrt(cp.pi))
                * exp_term
                / real_r
            )
        )

        coulomb_strength = cp.abs(q1 * q2)

        closest_coulomb = cp.argmin(
            cp.where(coulomb_strength > 0, real_r, cp.inf)
        )

        force_mag *= scale

        F_coul = -force_mag[:, None] * real_hat

        # ---------------------------------------------------------
        # DEBUG: largest real-space Coulomb pair force
        # ---------------------------------------------------------
        coul_force_abs = cp.abs(force_mag)

        if cp.any(coul_force_abs > 0):
            idx = int(cp.argmax(coul_force_abs))

            max_coul_force = float(coul_force_abs[idx].get())

            if max_coul_force > system.debug_max_coul_force:
                system.debug_max_coul_force = max_coul_force
                system.debug_max_coul_pair = (
                    int(real_i[idx]),
                    int(real_j[idx])
                )
                system.debug_max_coul_r = float(real_r[idx].get())

        cp.add.at(system.forces, real_i, F_coul)
        cp.add.at(system.forces, real_j, -F_coul)

        U_coul = (
            scale
            * COULOMB_CONSTANT
            * q1
            * q2
            * erfc_term
            / real_r
        )

        energy = cp.sum(U_coul)

        system.lj_PE = lj_PE
        system.real_space_PE += energy
        system.potential_energy += energy
    


# -------------------------------------
# ENERGY / COLLISION DEBUG LOGGER
# -------------------------------------


DEBUG_LOG_FILE = "MD_debug_log.txt"

def initialise_debug_log():

    with open(DEBUG_LOG_FILE, "w", encoding="utf-8") as f:

        f.write("=" * 120 + "\n")
        f.write("MD ENGINE DEBUG LOG\n")
        f.write("=" * 120 + "\n")

        f.write(f"TIMESTEP (ps): {TIME_STEP}\n")
        f.write(f"LJ CUTOFF (nm): {LJ_CUTOFF}\n")
        f.write(f"REAL CUTOFF (nm): {REAL_CUTOFF}\n")
        f.write(f"PME ALPHA: {PME_ALPHA}\n")
        f.write(f"PBC BOX LENGTH (nm): {PBC_BOX_LENGTH}\n")

        f.write("\n")
        f.write("=" * 120 + "\n")
        f.write("STEP-BY-STEP DATA\n")
        f.write("=" * 120 + "\n\n")

        f.write(
            "STEP\t"
            "TIME_PS\t"
            "KE\t"
            "LJ_PE\t"
            "REAL_PE\t"
            "RECIP_PE\t"
            "SELF_PE\t"
            "EXCL_PE\t"
            "TOTAL_PE\t"
            "TOTAL_E\t"
            "MIN_H_OW\t"
            "MIN_O_O\t"
            "MAX_FORCE\t"
            "MAX_FORCE_ATOM\t"
            "MAX_VELOCITY\t"
            "MAX_LJ_FORCE\t"
            "LJ_PAIR\t"
            "LJ_R\t"
            "MAX_COUL_FORCE\t"
            "COUL_PAIR\t"
            "COUL_R\t"
            "MAX_PME_FORCE\t"
            "PME_ATOM\t"
            "NEIGHBOUR_PAIRS\n"
        )

def debug_log_step(system):

    # ---------------------------------
    # Basic energy values
    # ---------------------------------

    step = int(system.step)
    time_ps = step * TIME_STEP

    KE = float(system.kinetic_energy)

    LJ = float(getattr(system, "lj_PE", 0.0))
    REAL = float(getattr(system, "real_space_PE", 0.0))
    RECIP = float(getattr(system, "reciprocal_PE", 0.0))
    SELF = float(getattr(system, "self_PE", 0.0))
    EXCL = float(getattr(system, "exclusion_PE", 0.0))

    TOTAL_PE = float(system.potential_energy)
    TOTAL_E = KE + TOTAL_PE

    # ---------------------------------
    # Atom information
    # ---------------------------------

    positions = cp.asarray(system.positions)
    velocities = cp.asarray(system.velocities)
    forces = cp.asarray(system.forces)

    # maximum force
    force_magnitudes = cp.linalg.norm(forces, axis=1)
    max_force = float(cp.max(force_magnitudes).get())
    max_force_atom = int(cp.argmax(force_magnitudes))

    # maximum velocity
    velocity_magnitudes = cp.linalg.norm(velocities, axis=1)
    max_velocity = float(cp.max(velocity_magnitudes).get())

    # ---------------------------------
    # Force-source diagnostics
    # ---------------------------------

    lj_pair = getattr(system, "debug_max_lj_pair", None)
    coul_pair = getattr(system, "debug_max_coul_pair", None)

    if lj_pair is not None:
        lj_i, lj_j = lj_pair
        lj_i_type = system.ff_atom_types[lj_i]
        lj_j_type = system.ff_atom_types[lj_j]
    else:
        lj_i = lj_j = -1
        lj_i_type = lj_j_type = "NONE"

    if coul_pair is not None:
        coul_i, coul_j = coul_pair
        coul_i_type = system.ff_atom_types[coul_i]
        coul_j_type = system.ff_atom_types[coul_j]
    else:
        coul_i = coul_j = -1
        coul_i_type = coul_j_type = "NONE"

    reciprocal_force = getattr(
        system,
        "debug_max_reciprocal_force",
        0.0
    )

    reciprocal_atom = getattr(
        system,
        "debug_max_reciprocal_atom",
        -1
    )


    # ---------------------------------
    # Minimum H-OW and O-O distances
    #
    # This uses atom types from system.atom_types
    # ---------------------------------

    min_h_ow = float("inf")
    min_o_o = float("inf")

    h_indices = []
    ow_indices = []
    o_indices = []

    for idx, atom_type in enumerate(system.ff_atom_types):

        if atom_type == "H":
            h_indices.append(idx)

        if atom_type == "OW":
            ow_indices.append(idx)

        if atom_type == "O":
            o_indices.append(idx)

    # H-OW
    for h in h_indices:

        for ow in ow_indices:

            if h == ow:
                continue

            r_vec = system.minimum_image(
                positions[ow] - positions[h]
            )

            r = float(cp.linalg.norm(r_vec).get())

            if r < min_h_ow:
                min_h_ow = r

    # O-O
    oxygen_indices = ow_indices + o_indices

    for a in range(len(oxygen_indices)):

        for b in range(a + 1, len(oxygen_indices)):

            i = oxygen_indices[a]
            j = oxygen_indices[b]

            r_vec = system.minimum_image(
                positions[j] - positions[i]
            )

            r = float(cp.linalg.norm(r_vec).get())

            if r < min_o_o:
                min_o_o = r

    # -------------------------------------------------------------------------
    # Neighbour list size
    # -------------------------------------------------------------------------

    try:
        neighbour_pairs = len(system.pair_i)
    except Exception:
        neighbour_pairs = -1

    # -------------------------------------------------------------------------
    # Write to log
    # -------------------------------------------------------------------------

    with open(DEBUG_LOG_FILE, "a", encoding="utf-8") as f:

        f.write(
            f"{step}\t"
            f"{time_ps:.8f}\t"
            f"{KE:.9f}\t"
            f"{LJ:.9f}\t"
            f"{REAL:.9f}\t"
            f"{RECIP:.9f}\t"
            f"{SELF:.9f}\t"
            f"{EXCL:.9f}\t"
            f"{TOTAL_PE:.9f}\t"
            f"{TOTAL_E:.9f}\t"
            f"{min_h_ow:.9f}\t"
            f"{min_o_o:.9f}\t"
            f"{max_force:.9f}\t"
            f"{max_force_atom}\t"
            f"{max_velocity:.9f}\t"
            f"{system.debug_max_lj_force:.9f}\t"
            f"{lj_i}-{lj_j} ({lj_i_type}-{lj_j_type})\t"
            f"{system.debug_max_lj_r:.9f}\t"
            f"{system.debug_max_coul_force:.9f}\t"
            f"{coul_i}-{coul_j} ({coul_i_type}-{coul_j_type})\t"
            f"{system.debug_max_coul_r:.9f}\t"
            f"{reciprocal_force:.9f}\t"
            f"{reciprocal_atom}\t"
            f"{neighbour_pairs}\n"
        )

def debug_log_closest_pairs(system):

    positions = cp.asarray(system.positions)
    charges = cp.asarray(system.charges)

    h_indices = []
    ow_indices = []

    for idx, atom_type in enumerate(system.ff_atom_types):

        if atom_type == "H":
            h_indices.append(idx)

        elif atom_type == "OW":
            ow_indices.append(idx)

    closest_h_ow = None

    # ---------------------------------
    # Closest H-OW
    # ---------------------------------

    for h in h_indices:

        for ow in ow_indices:

            r_vec = system.minimum_image(
                positions[ow] - positions[h]
            )

            r = float(cp.linalg.norm(r_vec).get())

            if closest_h_ow is None or r < closest_h_ow["r"]:

                closest_h_ow = {
                    "i": h,
                    "j": ow,
                    "r": r
                }

    # ---------------------------------
    # Closest O-O
    # ---------------------------------

    oxygen_indices = ow_indices

    closest_o_o = None

    for a in range(len(oxygen_indices)):

        for b in range(a + 1, len(oxygen_indices)):

            i = oxygen_indices[a]
            j = oxygen_indices[b]

            r_vec = system.minimum_image(
                positions[j] - positions[i]
            )

            r = float(cp.linalg.norm(r_vec).get())

            if closest_o_o is None or r < closest_o_o["r"]:

                closest_o_o = {
                    "i": i,
                    "j": j,
                    "r": r
                }

    with open(DEBUG_LOG_FILE, "a", encoding="utf-8") as f:

        f.write("\n")
        f.write(f"STEP {system.step} PAIR DETAILS\n")

        # H-OW
        if closest_h_ow is not None:

            i = closest_h_ow["i"]
            j = closest_h_ow["j"]
            r = closest_h_ow["r"]

            qi = float(charges[i].get())
            qj = float(charges[j].get())

            f.write(
                f"  H-OW CLOSEST: "
                f"{i}-{j} "
                f"r={r:.9f} nm "
                f"qH={qi:.6f} "
                f"qOW={qj:.6f}\n"
            )

        # O-O
        if closest_o_o is not None:

            i = closest_o_o["i"]
            j = closest_o_o["j"]
            r = closest_o_o["r"]

            f.write(
                f"  O-O CLOSEST: "
                f"{i}-{j} "
                f"r={r:.9f} nm\n"
            )

        f.write("\n")


initialise_debug_log()


# ---------------------------------------- PHYSICS LOOP ----------------------------------------


# caculate physics stuff - the heart <3
def calc_physics():

    system.potential_energy = 0.0
    system.potential_energy_gpu = cp.array([0.0], dtype=cp.float64)
    system.real_space_PE = 0.0
    system.reciprocal_PE = 0.0

    start_all_foces = time.perf_counter()

    # bond forces
    bond_forces(system)

    # bond angles
    angle_forces(system)

    # torsion angles
    torsion_forces(system)

    start_LJ_and_coulombs = time.perf_counter()

    # #LJ forces and Coulombs (non bonded loop)
    non_bonded_forces(system)

    end_LJ_and_coulombs = time.perf_counter()
    global time_taken_LJ_and_coulombs
    time_taken_LJ_and_coulombs += end_LJ_and_coulombs - start_LJ_and_coulombs


    # ----- Reciprocal-space PME -----
    
    
    start_PME = time.perf_counter()

    apply_exclusion_corrections(system)
    build_charge_grid(system)

    Qk = cp.fft.fftn(system.rho)
    Qk *= system.BC

    potential = cp.fft.ifftn(Qk).real
    potential *= PME_GRID**3

    forces_before = system.forces.copy()
    reciprocal_interpolation(system, potential)
    reciprocal_forces = system.forces - forces_before

    # ---------------------------------------------------------
    # DEBUG: largest reciprocal PME force
    # ---------------------------------------------------------
    reciprocal_force_mag = cp.linalg.norm(
        reciprocal_forces,
        axis=1
    )

    max_recip_idx = int(cp.argmax(reciprocal_force_mag))

    system.debug_max_reciprocal_force = float(
        reciprocal_force_mag[max_recip_idx].get()
    )

    system.debug_max_reciprocal_atom = max_recip_idx


    # account for the self energy shift
    system.self_PE = PME_SELF_ENERGY
    system.potential_energy += system.self_PE

    end_PME = time.perf_counter()
    global time_taken_PME
    time_taken_PME += end_PME - start_PME

    global time_taken_non_bonded
    time_taken_non_bonded = time_taken_PME + time_taken_LJ_and_coulombs
    
    end_all_foces = time.perf_counter()
    global time_taken_all_foces
    time_taken_all_foces += end_all_foces - start_all_foces
                                
    return system.potential_energy


TRAJECTORY_FILE = "MD_trajectory.xyz"

# Save one frame every N MD steps
TRAJECTORY_INTERVAL = 100


def initialise_trajectory():
    """
    Create an empty multi-frame XYZ trajectory file.

    XYZ format per frame:

        N
        comment
        Element x y z
        Element x y z
        ...

    Coordinates are in nm, as used internally by the engine.
    """

    with open(TRAJECTORY_FILE, "w", encoding="utf-8"):
        pass

def write_trajectory_frame(system):
    """
    Append one complete multi-frame XYZ frame.

    The complete frame is assembled in memory first so that the
    trajectory file is never intentionally left halfway through
    a frame by this function.
    """

    element_map = {
        "H": "H",
        "OW": "O",
        "O": "O",
        "HW": "H",
        "C": "C",
        "N": "N",
        "S": "S",
    }

    positions_cpu = cp.asnumpy(system.positions)
    n_atoms = system.n_atoms
    time_ps = system.step * TIME_STEP

    lines = [
        f"{n_atoms}\n",
        f"Step {system.step} Time {time_ps:.8f} ps "
        f"L={PBC_BOX_LENGTH:.8f} nm\n"
    ]

    for i in range(n_atoms):
        ff_type = system.ff_atom_types[i]

        if ff_type not in element_map:
            raise ValueError(
                f"No XYZ element mapping for atom type "
                f"'{ff_type}' at atom {i}"
            )

        element = element_map[ff_type]
        x, y, z = positions_cpu[i]

        lines.append(
            f"{element} "
            f"{x:.8f} "
            f"{y:.8f} "
            f"{z:.8f}\n"
        )

    # Sanity check before touching the file
    if len(lines) != n_atoms + 2:
        raise RuntimeError(
            "XYZ frame construction failed: wrong atom count."
        )

    frame = "".join(lines)

    with open(TRAJECTORY_FILE, "a", encoding="utf-8") as f:
        f.write(frame)
        f.flush()


initialise_trajectory()

# ---------------------------------------- SIMULATION LOOP ----------------------------------------


# initial forces and model and neighbour list before simulation starts - need valid forces before starting
system.forces.fill(0.0)
build_neighbour_lists(system, NEIGHBOUR_CUTOFF)
system.potential_energy = calc_physics()
render_initial_model()


last_neighbour_build_step = 0
relax_steps = 0
timestep_x = 1
max_displacement2 = 0

# running simulation
while True:
    start_total_time = time.perf_counter()
    update_camera(TIME_STEP)
    old_positions = system.positions.copy()
    # lighting position for camera
    cam_light.pos = scene.camera.pos - scene.camera.axis.norm()*3
    
    # Verlet integration method:

    # 1. First physical half-step
    free_half_velocities = (system.velocities + 0.5 * (system.forces / system.masses[:, None]) * TIME_STEP)

    # 2. Position RATTLE
    constrained_half_velocities = rattle_position_constraints(system, old_positions, free_half_velocities)

    # Replace velocity with the constrained half-step velocity
    system.velocities[:] = constrained_half_velocities

    # 3. Now wrap molecules
    wrap_molecules(system)

    # 4. Forces at q_(n+1)
    system.forces.fill(0.0)

    system.potential_energy = calc_physics()

    # 5. Second physical half-step
    system.velocities += (0.5 * (system.forces / system.masses[:, None]) * TIME_STEP)

    # 6. Relaxation damping
    if system.step < relax_steps / 5:
        damping = 0.99
    elif system.step < relax_steps / 2:
        damping = 0.995
    elif system.step < relax_steps:
        damping = 0.999
    else:
        damping = 1.0

    system.velocities *= damping

    # 7. Final RATTLE velocity constraint
    rattle_velocity_constraints(system)

    system.kinetic_energy = calc_KE(system)                                             # energy tracking
    system.total_energy = system.kinetic_energy + system.potential_energy
    system.average_total_energy += system.total_energy

    debug_log_step(system)
    debug_log_closest_pairs(system)


    system.step += 1

    if system.step == relax_steps:
        system.average_total_energy = 0


    # Get positions onto np CPU arrays
    positions_cpu = cp.asnumpy(system.positions)

    # find the displacement of the atom this timestep
    for i, atom in enumerate(system.atoms):
        disp2 = system.minimum_image_cpu(positions_cpu[i] - atom.last_neighbour_reference)

        displacement2 = r2 = cp.dot(disp2,disp2)

        if displacement2 > max_displacement2:
            max_displacement2 = displacement2

    if max_displacement2 > (SKIN_CUTOFF / 2)**2:
        print(f"New List After "f"{system.step - last_neighbour_build_step} Steps")
        last_neighbour_build_step = system.step
        build_neighbour_lists(system, NEIGHBOUR_CUTOFF)
        max_displacement2 = 0
        

    start_total_graphics = time.perf_counter() 

    # ----- graphics render  ------

    # update graphics once per frame
    for i in range(system.n_atoms):
        system.balls[i].pos = vector(*positions_cpu[i])                                 # moves balls to current pos

    bond_a_cpu = cp.asnumpy(system.bond_a)
    bond_b_cpu = cp.asnumpy(system.bond_b)
    for bond_index in range(len(system.bond_a)):                              # update the bond visuals

        a1 = bond_a_cpu[bond_index]
        a2 = bond_b_cpu[bond_index]
        system.bond_visuals[bond_index].pos = vector(*positions_cpu[a1])

        system.bond_visuals[bond_index].axis = vector(
            *(positions_cpu[a2] - positions_cpu[a1])
        )

    end_total_graphics = time.perf_counter()
    time_taken_graphics += end_total_graphics - start_total_graphics

    end_total_time = time.perf_counter()
    time_taken_total += end_total_time - start_total_time

    time_other = time_taken_total - (time_taken_graphics + time_taken_all_foces)


    if system.step % TRAJECTORY_INTERVAL == 0:
        write_trajectory_frame(system)


    # printing and debuging stuff that prints every x timesteps
    if system.step % timestep_x == 0:

        if system.step == (relax_steps+timestep_x):
            system.initial_total_energy = float(cp.asnumpy(system.total_energy))

        if system.step >= (relax_steps+timestep_x):
            system.average_total_energy = system.average_total_energy/timestep_x
            system.steps.append(system.step)
            system.avg_sys_energy_values.append(float(cp.asnumpy(system.total_energy)))

            line.set_data(system.steps, system.avg_sys_energy_values)

            ax.relim()
            ax.autoscale_view()
            plt.draw()
            plt.pause(0.001)
            system.average_total_energy = 0

            #print(f"KE: {system.kinetic_energy:.12f}  PE: {system.potential_energy:.12f}  Total: {system.total_energy:.12f}")
        '''
        print()
        print(f"Step: {system.step}")
        print(f"KE: {system.kinetic_energy:.6f}  PE: {system.potential_energy:.6f}  Total: {system.total_energy:.6f}")
        print(f"Time Taken per cycle: {time_taken_total/timestep_x:.6f}s")
        print()
        print(f"Time By % of Graphics: {((time_taken_graphics/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bonds: {((time_taken_bonds/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Bond Angles: {((time_taken_bond_angles/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Torsions: {((time_taken_torsions/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of LJ: {((time_taken_LJ_and_coulombs/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
        print(f"Time By % of Coulombs: {((time_taken_PME/timestep_x)/time_taken_total)*100*timestep_x:.6f} %")
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
        time_taken_PME = 0
        time_taken_LJ_and_coulombs = 0
        time_taken_graphics = 0
        time_taken_total = 0
        time_other = 0