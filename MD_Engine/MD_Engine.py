from vpython import canvas, rate, sphere, vector,mag, norm,dot, cylinder
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

                  # atom variables
positions = []
velocities = []
forces = []
balls = []
masses = []
bonds = [[0, 1, 1, 200], [1, 2, 1, 200], [1, 3, 1, 200], [4, 5, 1, 200], [5, 6, 1, 200], [6, 7, 1, 200]]
bond_visuals = []

                # simulation pararmeters
E_pot = 0
K_E = 0
E_total = 0
No_Balls = 9
dt = 0.0005
sigma = 0.8
epsilon = 0.2
PBC_Box_Length = 10

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
            positions.append(vector(x, y, z))
            count += 1

# re centres balls from 0,0,0 centre
for p in positions:
    p.x -= (n_per_axis * spacing) / 2
    p.y -= (n_per_axis * spacing) / 2
    p.z -= (n_per_axis * spacing) / 2

# create empty velocities
for i in range(No_Balls):
    vel = vector(0, 0, 0)
    velocities.append(vel)

# create empty forces
for i in range(No_Balls):
    force = vector(0, 0, 0)
    forces.append(force)

# create empty masses
for i in range(No_Balls):
    mass = 1
    masses.append(mass)

# create balls
for i in range(No_Balls):
    b = sphere(pos=positions[i], radius=0.1, color=vector(0,0,0), make_trail=False)
    balls.append(b)

# Lennard Jones Force Equation (VDW's)
def Calc_LJ(dist, direction):
    return  -((24*epsilon*((2*((sigma/dist)**12)) - ((sigma/dist)**6)))*(1/dist))*direction

# calc Kintetic energy of particles
def Calc_KE():
    K_E = 0
    for i in range(No_Balls): 
        K_E = K_E + (0.5*masses[i]*(dot(velocities[i], velocities[i]))) # use dot product as v^2 to calc ke
    return K_E

# Draw lines between bonds
for bond in bonds:
    c = cylinder(
        pos=positions[bond[0]],
        axis=positions[bond[1]] - positions[bond[0]],
        radius=0.02,
        color=vector(1,1,1)
    )
    bond_visuals.append(c)

# checks if the two atoms are bonded, called a 1-2 exclusion in MD
# this is used to skip directly bonded atoms so their LJ forces dont get calculated
def Are_Bonded(i, j):                                                                   
    for bond in bonds:
        if (bond[0] == i and bond[1] == j) or (bond[0] == j and bond[1] == i):
            return True
    return False

# checks if a particle has reached the box boundery and needs warping
def PBC_Box_For_Pos():
    for i in range(No_Balls):                       
            if positions[i].x > PBC_Box_Length/2:
                positions[i].x -= PBC_Box_Length
            if positions[i].x < -PBC_Box_Length/2:
                positions[i].x += PBC_Box_Length

            if positions[i].y > PBC_Box_Length/2:
                positions[i].y -= PBC_Box_Length
            if positions[i].y < -PBC_Box_Length/2:
                positions[i].y += PBC_Box_Length

            if positions[i].z > PBC_Box_Length/2:
                positions[i].z -= PBC_Box_Length
            if positions[i].z < -PBC_Box_Length/2:
                positions[i].z += PBC_Box_Length

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
        a = bond[0]
        b = bond[1]
        r0 = bond[2]
        k = bond[3]

        # compute the bond vector and bond lengths
        r_vec = positions[b] - positions[a]
        PBC_Box_For_Vectors(r_vec)
        r = mag(r_vec)


        if r == 0:                              # avoid devide by zero just in case
            continue

        r_hat = norm(r_vec)                     # bond direction

        F = k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector

        forces[a] += F                          # apply Newtons third law - equal and opposite forces
        forces[b] -= F

        E_pot += 0.5 * k * (r - r0)**2          # calc the harmonic bond potential

    # LJ for nonbonded pairs only
    cutoff = 2.5 * sigma                        # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
    U_shift = 4 * epsilon * ((sigma / cutoff)**12 - (sigma / cutoff)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

    for i in range(No_Balls):                   # loop over each unique atom pair once
        for j in range(i + 1, No_Balls):

            if Are_Bonded(i, j):                # skips bonded pairs
                continue

            r_vec = positions[j] - positions[i] # find pair distances ect
            PBC_Box_For_Vectors(r_vec)
            r = mag(r_vec)
            
            if r == 0 or r > cutoff:            # skip invalid or too distant pairs
                continue

            r_hat = norm(r_vec)

            F = Calc_LJ(r, r_hat)               # calculate LJ force

            forces[i] += F                      # Newtons third law - apply to both atoms
            forces[j] -= F

            E_pot += 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) - U_shift # add the LJ potential energy on

    return E_pot                                


# initial forces before simulation starts - need valid forces before starting
for i in range(No_Balls):
    forces[i] = vector(0, 0, 0)
E_pot = Calc_Physics()

step = 0
relaxing = True 
relax_steps = 1000
while True:
    rate(120)                                                       # capped at 120 render

    for _ in range(10):                                             # run physics faster

        # 1. half-step velocity update
        for i in range(No_Balls):
            velocities[i] += 0.5 * (forces[i] / masses[i]) * dt     # first half-step velocity update - uses current force to push velocity halfway forward

        # 2. position update
        for i in range(No_Balls):
            positions[i] += velocities[i] * dt                      # update pos

        # 3. wrap positions back into box
        PBC_Box_For_Pos()

        # 4. reset forces
        for i in range(No_Balls):
            forces[i] = vector(0, 0, 0)                             # clear old forces before computing new ones

        # 5. compute new forces
        E_pot = Calc_Physics()                                      # now get new forces at the new positions

        # 6. second half-step velocity update
        for i in range(No_Balls):
            velocities[i] += 0.5 * (forces[i] / masses[i]) * dt     # compleate  the full velocity update using the new forces

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
            velocities[i] *= damping                                # dampens some of the velocity at each step when begining


            
        K_E = Calc_KE()                                             # energy tracking
        E_total = K_E + E_pot
        step += 1

        if step % 100 == 0:                                         # debug printout stuff
            print(f"KE: {K_E:.6f}  PE: {E_pot:.6f}  Total: {E_total:.6f}")
            for bond in bonds:
                a = bond[0]
                b = bond[1]
                r_vec = positions[b] - positions[a]
                PBC_Box_For_Vectors(r_vec)
                print(f"Bond {a}-{b}: {mag(r_vec):.6f}")

    # update graphics once per frame
    for i in range(No_Balls):
        balls[i].pos = positions[i]                                 # moves balls to current pos
        speed = mag(velocities[i])                                  # ball colour stuff
        ColourPS = min(speed / 2.5, 1)
        balls[i].color = vector(ColourPS, 0, 1 - ColourPS)
        balls[i].trail_color = vector(ColourPS, 0, 1 - ColourPS)

    for idx, bond in enumerate(bonds):                              # update the bond visuals
        bond_visuals[idx].pos = positions[bond[0]]
        bond_visuals[idx].axis = positions[bond[1]] - positions[bond[0]]