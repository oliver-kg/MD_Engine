from vpython import canvas, rate, sphere, vector,mag, norm,dot, cross, cylinder, box
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

class Bond_Angle:
    def __init__(self, ideal_ang_deg, theta, first, second, third, k):
        self.ideal_ang = math.radians(ideal_ang_deg)
        self.theta = theta
        self.first = first
        self.second = second
        self.third = third
        self.k = k
    

atoms = []
balls = []
bonds = []
bond_angles = []
bond_visuals = []


# create test bonds
bonds = [
    #methane
    Bond(0, 1, 1, 500),
    Bond(0, 2, 1, 500),
    Bond(0, 3, 1, 500),
    Bond(0, 4, 1, 500),

    #methane
    Bond(5, 6, 1, 500),
    Bond(5, 7, 1, 500),
    Bond(5, 8, 1, 500),
    Bond(5, 9, 1, 500),

    #ethane
    Bond(30, 31, 1, 500),
    Bond(30, 32, 1, 500),
    Bond(30, 33, 1, 500),
    Bond(30, 34, 1, 500),
    Bond(31, 35, 1, 500),
    Bond(31, 36, 1, 500),
    Bond(31, 37, 1, 500),

    #water
    Bond(38, 39, 1, 500),
    Bond(39, 40, 1, 500)

]

#create test angles
bond_angles = [
    
    #methane
    Bond_Angle(109.5, 0, 1, 0, 2, 200),
    Bond_Angle(109.5, 0, 1, 0, 3, 200),
    Bond_Angle(109.5, 0, 1, 0, 4, 200),
    Bond_Angle(109.5, 0, 2, 0, 3, 200),
    Bond_Angle(109.5, 0, 2, 0, 4, 200),
    Bond_Angle(109.5, 0, 3, 0, 4, 200),

    #methane
    Bond_Angle(109.5, 0, 6, 5, 7, 200),
    Bond_Angle(109.5, 0, 6, 5, 8, 200),
    Bond_Angle(109.5, 0, 6, 5, 9, 200),
    Bond_Angle(109.5, 0, 7, 5, 8, 200),
    Bond_Angle(109.5, 0, 7, 5, 9, 200),
    Bond_Angle(109.5, 0, 8, 5, 9, 200),

    #ethane
    Bond_Angle(107.4, 0, 32, 30, 33, 200),
    Bond_Angle(107.4, 0, 32, 30, 34, 200),
    Bond_Angle(107.4, 0, 33, 30, 34, 200),
    Bond_Angle(111, 0, 32, 30, 31, 200),
    Bond_Angle(111, 0, 33, 30, 31, 200),
    Bond_Angle(111, 0, 34, 30, 31, 200),
    Bond_Angle(107.4, 0, 36, 31, 35, 200),
    Bond_Angle(107.4, 0, 36, 31, 37, 200),
    Bond_Angle(107.4, 0, 37, 31, 35, 200),
    Bond_Angle(111, 0, 35, 31, 30, 200),
    Bond_Angle(111, 0, 36, 31, 30, 200),
    Bond_Angle(111, 0, 37, 31, 30, 200),

    #water
    Bond_Angle(104.5, 0, 38, 39, 40, 200),
    ]

#create test molecules
molecules = [
    #methane
    Molecule([0,1,2,3,4]),
    
    #methane
    Molecule([5,6,7,8,9]),
    
    #ethane
    Molecule([30,31,32,33,34,35,36,37]),

    #water
    Molecule([38,39,40]),

    #singles
    Molecule([10]),
    Molecule([11]),
    Molecule([12]),
    Molecule([13]),
    Molecule([14]),
    Molecule([15]),
    Molecule([16]),
    Molecule([17]),
    Molecule([18]),
    Molecule([19]),
    Molecule([20]),
    Molecule([21]),
    Molecule([22]),
    Molecule([23]),
    Molecule([24]),
    Molecule([25]),
    Molecule([26]),
    Molecule([27]),
    Molecule([28]),
    Molecule([29])
    ]

                # simulation pararmeters
e_pot = 0
k_e = 0
e_total = 0
no_balls = 41
dt = 0.001
sigma = 0.8
epsilon = 0.2
PBC_box_length = 10

# draws a faint box of the simulation area
L = PBC_box_length

box_visual = box(
    pos=vector(0,0,0),
    size=vector(L, L, L),
    opacity=0.08,        # faint
    color=vector(1,1,1)
)

# create pos
count = 0                                   # chooses a rough 3d grid size and spacing
n_per_axis = math.ceil(no_balls ** (1/3))
spacing = 1 * sigma
jitter_amount = 0.2 * sigma                 # adds randomness

# create rough grid of balls
for i in range(n_per_axis):
    for j in range(n_per_axis):
        for k in range(n_per_axis):
            if count >= no_balls:
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


# change masses test
atoms[0].mass = 12
atoms[5].mass = 12
atoms[30].mass = 12
atoms[31].mass = 12
atoms[39].mass = 16

# create balls
for a in atoms:
    b = sphere(pos=a.pos, radius=0.1, color=vector(0,0,0), make_trail=False)
    balls.append(b)

# change radius test
balls[0].radius = 0.2
balls[5].radius = 0.2
balls[30].radius = 0.2
balls[31].radius = 0.2
balls[39].radius = 0.22

# Lennard Jones Force Equation (VDW's)
def Calc_LJ(dist, direction):
    return  -((24*epsilon*((2*((sigma/dist)**12)) - ((sigma/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def Calc_KE():
    k_e = 0
    for i in range(no_balls): 
        k_e = k_e + (0.5*atoms[i].mass*(dot(atoms[i].vel, atoms[i].vel))) # use dot product as v^2 to calc ke
    return k_e

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
    half = PBC_box_length / 2

    for mol in molecules:
        ref_pos = atoms[mol.atom_indices[0]].pos
        shift = vector(0, 0, 0)

        if ref_pos.x > half:
            shift.x -= PBC_box_length
        elif ref_pos.x < -half:
            shift.x += PBC_box_length

        if ref_pos.y > half:
            shift.y -= PBC_box_length
        elif ref_pos.y < -half:
            shift.y += PBC_box_length

        if ref_pos.z > half:
            shift.z -= PBC_box_length
        elif ref_pos.z < -half:
            shift.z += PBC_box_length

        if shift.x != 0 or shift.y != 0 or shift.z != 0:
            for i in mol.atom_indices:
                atoms[i].pos += shift

# minimum image periodic correction - turns the vector into the nearest image vector
def PBC_Box_For_Vectors(r_vector):                                  
    if r_vector.x > PBC_box_length/2:
        r_vector.x = r_vector.x - PBC_box_length
    if r_vector.x < -PBC_box_length/2:
        r_vector.x = r_vector.x + PBC_box_length

    if r_vector.y > PBC_box_length/2:
        r_vector.y = r_vector.y - PBC_box_length
    if r_vector.y < -PBC_box_length/2:
        r_vector.y = r_vector.y + PBC_box_length

    if r_vector.z > PBC_box_length/2:
        r_vector.z = r_vector.z - PBC_box_length
    if r_vector.z < -PBC_box_length/2:
        r_vector.z = r_vector.z + PBC_box_length

# caculate physics stuff - the heart <3
def Calc_Physics():
    e_pot = 0.0

    # bond angles
    for angle in bond_angles:
        theta = angle.theta
        ideal_ang = angle.ideal_ang
        atomA = angle.first
        atomB = angle.second
        atomC = angle.third
        k = angle.k

        BA = atoms[atomA].pos - atoms[atomB].pos                                # find vector of B to A
        BC = atoms[atomC].pos - atoms[atomB].pos                                # find vector of B to C
        
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

        e_pot += 0.5 * k * (theta - ideal_ang)**2                               # update the potential energy change from the angle

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

        e_pot += 0.5 * k * (r - r0)**2          # calc the harmonic bond potential

    # LJ for nonbonded pairs only
    cutoff = 2.5 * sigma                        # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
    U_shift = 4 * epsilon * ((sigma / cutoff)**12 - (sigma / cutoff)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

    for i in range(no_balls):                   # loop over each unique atom pair once
        for j in range(i + 1, no_balls):

            if Are_Bonded(i, j):                # skips bonded pairs
                continue

            r_vec = atoms[j].pos - atoms[i].pos # find pair distances ect
            PBC_Box_For_Vectors(r_vec)          # update "copy vectors"
            r = mag(r_vec)
            
            if r == 0 or r > cutoff:            # skip invalid or too distant pairs
                continue

            r_hat = norm(r_vec)

            F = Calc_LJ(r, r_hat)               # calculate LJ force
            atoms[i].force += F                      # Newtons third law - apply to both atoms
            atoms[j].force -= F

            e_pot += 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) - U_shift # add the LJ potential energy on

    return e_pot                                


# initial forces before simulation starts - need valid forces before starting
for i in range(no_balls):
    atoms[i].force = vector(0, 0, 0)
e_pot = Calc_Physics()

step = 0
relaxing = True 
relax_steps = 3500
E0 = None
while True:
    rate(120)                                                       # capped at 120 render

    for _ in range(10):                                             # run physics faster

        # 1. half-step velocity update
        for i in range(no_balls):
            atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * dt     # first half-step velocity update - uses current force to push velocity halfway forward

        # 2. position update
        for i in range(no_balls):
            atoms[i].pos += atoms[i].vel * dt                      # update pos

        # 3. wrap positions back into box
        Wrap_Molecules()

        # 4. reset forces
        for i in range(no_balls):
            atoms[i].force = vector(0, 0, 0)                             # clear old forces before computing new ones

        # 5. compute new forces
        e_pot = Calc_Physics()                                           # now get new forces at the new positions

        # 6. second half-step velocity update
        for i in range(no_balls):
            atoms[i].vel += 0.5 * (atoms[i].force / atoms[i].mass) * dt     # compleate  the full velocity update using the new forces

        # 7. dampen starting strains - ensures the system is calm so it doesent blow up
        if step < relax_steps/5:
            damping = 0.99
        elif step < relax_steps/2:
            damping = 0.995
        elif step < relax_steps:
            damping = 0.999
        else:
            damping = 1                                           # no more dampening
            
        for i in range(no_balls):
            atoms[i].vel *= damping                                 # dampens some of the velocity at each step when begining


        k_e = Calc_KE()                                             # energy tracking
        e_total = k_e + e_pot
        step += 1

        if step % 100 == 0:                                         # debug printout stuff
            print(f"KE: {k_e:.6f}  PE: {e_pot:.6f}  Total: {e_total:.6f}")
            for bond in bonds:
                a = bond.a1
                b = bond.a2
                r_vec = atoms[b].pos - atoms[a].pos
                PBC_Box_For_Vectors(r_vec)
                print(f"Bond {a}-{b}: {mag(r_vec):.6f}")
                
            for angles in bond_angles:
                print(f"Angle: {math.degrees(angles.theta)}")

    # update graphics once per frame
    for i in range(no_balls):
        balls[i].pos = atoms[i].pos                                 # moves balls to current pos
        speed = mag(atoms[i].vel)                                   # ball colour stuff
        ColourPS = min(speed / 2, 1)
        balls[i].color = vector(ColourPS, 0, 1 - ColourPS)
        balls[i].trail_color = vector(ColourPS, 0, 1 - ColourPS)

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = atoms[bond.a1].pos
        bond_visuals[idx].axis = atoms[bond.a2].pos - atoms[bond.a1].pos