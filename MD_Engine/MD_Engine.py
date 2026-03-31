from vpython import canvas, rate, sphere, vector,mag, norm,dot, cylinder, box
import tkinter as tk
import math
import random

# Get screen dimensions
root = tk.Tk()
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
root.destroy()

# Create a canvas that fills the screen 
scene = canvas(title="Full Screen VPython",
               width=screen_width,
               height=screen_height)

# atom class
class Atom:
    def __init__(self, mass, pos, vel, force, charge):
        self.mass = mass
        self.pos = pos
        self.vel = vel
        self.force = force
        self.charge = charge

# molecule class
class Molecule:
    def __init__(self, atom_indices):
        self.atom_indices = atom_indices

# bonds class
class Bond:
    def __init__(self, a1, a2, Ideal_Dist, K):
        self.a1 = a1
        self.a2 = a2
        self.Ideal_Dist = Ideal_Dist
        self.K = K

atoms = []
balls = []
bonds = []
bond_visuals = []

# create test bonds
bonds = [
    Bond(0, 1, 1, 200),
    Bond(1, 2, 1, 200),
    Bond(1, 3, 1, 200),
    Bond(4, 5, 1, 200),
    Bond(5, 6, 1, 200),
    Bond(6, 7, 1, 200)
]

#create test molecules
molecules = [
    Molecule([0,1,2,3]),
    Molecule([4,5,6,7]),
    Molecule([8]),
    Molecule([9]),
    Molecule([10]),
    Molecule([11]),
    Molecule([12]),
    Molecule([13]),
    Molecule([14])
    ]

                # simulation pararmeters
E_pot = 0
K_E = 0
E_total = 0
No_Balls = 15
dt = 0.0008
sigma = 0.8
epsilon = 0.2
PBC_Box_Length = 10


# draws a faint box of the simulation area
L = PBC_Box_Length

box_visual = box(
    pos=vector(0,0,0),
    size=vector(L, L, L),
    opacity=0.08,        # faint
    color=vector(1,1,1)
)


# create pos

count = 0                                   # chooses a rough 3d grid size and spacing
n_per_axis = math.ceil(No_Balls ** (1/3))
spacing = 1 * sigma
jitter_amount = 0.2 * sigma                 # adds randomness


for i in range(n_per_axis):                 # create rough grid of balls
    for j in range(n_per_axis):
        for k in range(n_per_axis):
            if count >= No_Balls:
                break
            x = i * spacing + random.uniform(-jitter_amount, jitter_amount)
            y = j * spacing + random.uniform(-jitter_amount, jitter_amount)
            z = k * spacing + random.uniform(-jitter_amount, jitter_amount)
            
            atoms.append(Atom(
                mass=1,                     # create empty masses
                pos=vector(x, y, z),        # create positions
                vel=vector(0, 0, 0),        # create empty velocities
                force=vector(0,0,0),        # create empty forces
                charge=0
            ))
            count += 1

# re centres balls from 0,0,0 centre
for a in atoms:
    a.pos.x -= (n_per_axis * spacing) / 2
    a.pos.y -= (n_per_axis * spacing) / 2
    a.pos.z -= (n_per_axis * spacing) / 2


# create balls
for a in atoms:
    b = sphere(pos=a.pos, radius=0.1, color=vector(0,0,0), make_trail=False)
    balls.append(b)

# Lennard Jones Force Equation (VDW's)
def Calc_LJ(dist, direction):
    return  -((24*epsilon*((2*((sigma/dist)**12)) - ((sigma/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def Calc_KE():
    K_E = 0
    for i in range(No_Balls): 
        K_E = K_E + (0.5*atoms[i].mass*(dot(atoms[i].vel, atoms[i].vel))) # use dot product as v^2 to calc ke
    return K_E

# Draw lines between bonds
for bond in bonds:
    c = cylinder(
        pos=atoms[bond.a1].pos,
        axis=atoms[bond.a2].pos - atoms[bond.a1].pos,
        radius=0.02,
        color=vector(1,1,1)
    )
    bond_visuals.append(c)

# checks if the two atoms are bonded, called a 1-2 exclusion in MD
# this is used to skip directly bonded atoms so their LJ forces dont get calculated
def Are_Bonded(i, j):                                                                   
    for bond in bonds:
        if (bond.a1 == i and bond.a2 == j) or (bond.a1 == j and bond.a2 == i):
            return True
    return False

# checks if a molecule has reached the box boundery and needs warping
def Wrap_Molecules():
    half = PBC_Box_Length / 2

    for mol in molecules:
        ref_pos = atoms[mol.atom_indices[0]].pos
        shift = vector(0, 0, 0)

        if ref_pos.x > half:
            shift.x -= PBC_Box_Length
        elif ref_pos.x < -half:
            shift.x += PBC_Box_Length

        if ref_pos.y > half:
            shift.y -= PBC_Box_Length
        elif ref_pos.y < -half:
            shift.y += PBC_Box_Length

        if ref_pos.z > half:
            shift.z -= PBC_Box_Length
        elif ref_pos.z < -half:
            shift.z += PBC_Box_Length

        if shift.x != 0 or shift.y != 0 or shift.z != 0:
            for i in mol.atom_indices:
                atoms[i].pos += shift

# minimum image periodic correction - turns the vector into the nearest image vector
def PBC_Box_For_Vectors(r_vector):                                  
    if r_vector.x > PBC_Box_Length/2:
        r_vector.x = r_vector.x - PBC_Box_Length
    if r_vector.x < -PBC_Box_Length/2:
        r_vector.x = r_vector.x + PBC_Box_Length

    if r_vector.y > PBC_Box_Length/2:
        r_vector.y = r_vector.y - PBC_Box_Length
    if r_vector.y < -PBC_Box_Length/2:
        r_vector.y = r_vector.y + PBC_Box_Length

    if r_vector.z > PBC_Box_Length/2:
        r_vector.z = r_vector.z - PBC_Box_Length
    if r_vector.z < -PBC_Box_Length/2:
        r_vector.z = r_vector.z + PBC_Box_Length

# caculate physics stuff - the heart <3
def Calc_Physics():
    E_pot = 0.0

    # bond forces
    for bond in bonds:
        a = bond.a1
        b = bond.a2
        r0 = bond.Ideal_Dist
        k = bond.K

        # compute the bond vector and bond lengths
        r_vec = atoms[b].pos - atoms[a].pos
        PBC_Box_For_Vectors(r_vec)
        r = mag(r_vec)


        if r == 0:                              # avoid devide by zero just in case
            continue

        r_hat = norm(r_vec)                     # bond direction

        F = k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector

        atoms[a].force += F                          # apply Newtons third law - equal and opposite forces
        atoms[b].force -= F

        E_pot += 0.5 * k * (r - r0)**2          # calc the harmonic bond potential

    # LJ for nonbonded pairs only
    cutoff = 2.5 * sigma                        # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
    U_shift = 4 * epsilon * ((sigma / cutoff)**12 - (sigma / cutoff)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

    for i in range(No_Balls):                   # loop over each unique atom pair once
        for j in range(i + 1, No_Balls):

            if Are_Bonded(i, j):                # skips bonded pairs
                continue

            r_vec = atoms[j].pos - atoms[i].pos # find pair distances ect
            PBC_Box_For_Vectors(r_vec)
            r = mag(r_vec)
            
            if r == 0 or r > cutoff:            # skip invalid or too distant pairs
                continue

            r_hat = norm(r_vec)

            F = Calc_LJ(r, r_hat)               # calculate LJ force

            atoms[i].force += F                      # Newtons third law - apply to both atoms
            atoms[j].force -= F

            E_pot += 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) - U_shift # add the LJ potential energy on

    return E_pot                                


# initial forces before simulation starts - need valid forces before starting
for i in range(No_Balls):
    atoms[i].force = vector(0, 0, 0)
E_pot = Calc_Physics()

step = 0
relaxing = True 
relax_steps = 1000
while True:
    rate(120)                                                       # capped at 120 render

    for _ in range(10):                                             # run physics faster

        # 1. half-step velocity update
        for i in range(No_Balls):
            atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * dt     # first half-step velocity update - uses current force to push velocity halfway forward

        # 2. position update
        for i in range(No_Balls):
            atoms[i].pos += atoms[i].vel * dt                      # update pos

        # 3. wrap positions back into box
        Wrap_Molecules()

        # 4. reset forces
        for i in range(No_Balls):
            atoms[i].force = vector(0, 0, 0)                             # clear old forces before computing new ones

        # 5. compute new forces
        E_pot = Calc_Physics()                                      # now get new forces at the new positions

        # 6. second half-step velocity update
        for i in range(No_Balls):
            atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * dt     # compleate  the full velocity update using the new forces

        # 7. dampen starting strains - ensures the system is calm so it doesent blow up
        if step < relax_steps/5:
            
            damping = 0.99
        elif step < relax_steps/2:
            damping = 0.995
        elif step < relax_steps:
            damping = 0.999
        else:
            damping = 1.0                                           # no more dampening
            

        for i in range(No_Balls):
            atoms[i].vel *= damping                                # dampens some of the velocity at each step when begining


            
        K_E = Calc_KE()                                             # energy tracking
        E_total = K_E + E_pot
        step += 1

        if step % 100 == 0:                                         # debug printout stuff
            print(f"KE: {K_E:.6f}  PE: {E_pot:.6f}  Total: {E_total:.6f}")
            for bond in bonds:
                a = bond.a1
                b = bond.a2
                r_vec = atoms[b].pos - atoms[a].pos
                PBC_Box_For_Vectors(r_vec)
                print(f"Bond {a}-{b}: {mag(r_vec):.6f}")

    # update graphics once per frame
    for i in range(No_Balls):
        balls[i].pos = atoms[i].pos                                 # moves balls to current pos
        speed = mag(atoms[i].vel)                                  # ball colour stuff
        ColourPS = min(speed / 2.5, 1)
        balls[i].color = vector(ColourPS, 0, 1 - ColourPS)
        balls[i].trail_color = vector(ColourPS, 0, 1 - ColourPS)

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = atoms[bond.a1].pos
        bond_visuals[idx].axis = atoms[bond.a2].pos - atoms[bond.a1].pos