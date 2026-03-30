from vpython import canvas, dist, rate, sphere, vector, color,mag, norm,dot, cylinder
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
bonds = [[0, 1, 1, 50], [1, 2, 1, 50]]
bond_visuals = []
E_pot = 0
K_E = 0
E_total = 0
No_Balls = 3
dt = 0.0001
sigma = 1
epsilon = 0.5
PBC_Box_Length = 100

# create pos
count = 0
n_per_axis = math.ceil(No_Balls ** (1/3))
spacing = 1 * sigma
jitter_amount = 0.2 * sigma  # how much randomness to add

# create rough grid of balls
for i in range(n_per_axis):
    for j in range(n_per_axis):
        for k in range(n_per_axis):
            if count >= No_Balls:
                break
            x = i * spacing + random.uniform(-jitter_amount, jitter_amount)
            y = j * spacing + random.uniform(-jitter_amount, jitter_amount)
            z = k * spacing + random.uniform(-jitter_amount, jitter_amount)
            positions.append(vector(x, y, z))
            count += 1

# space balls from 0,0,0 centre
for p in positions:
    p.x -= (n_per_axis * spacing) / 2
    p.y -= (n_per_axis * spacing) / 2
    p.z -= (n_per_axis * spacing) / 2

positions = [
    vector(0,0,0),
    vector(1,0,0),
    vector(2,0,0)
]

# create empty velocities
for i in range(No_Balls):
    vel = vector(0, 0, 0)
    velocities.append(vel)
velocities = [
    vector(0.01,0,0),
    vector(0,0,0),
    vector(-0.01,0,0)
]

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

def PBC_Box(i):
    if positions[i].x > PBC_Box_Length/2:                               # checks if the x,y or z position are out of bounds
        positions[i].x = positions[i].x - PBC_Box_Length                # adjust the new current position

    if positions[i].x < -PBC_Box_Length/2:                              # repeat for the two ways the x axis can move, eg, y and z direction
        positions[i].x = positions[i].x + PBC_Box_Length

    if positions[i].y > PBC_Box_Length/2:                               # repeat for the y axis
        positions[i].y = positions[i].y - PBC_Box_Length

    if positions[i].y < -PBC_Box_Length/2:
        positions[i].y = positions[i].y + PBC_Box_Length

    if positions[i].z > PBC_Box_Length/2:                               # repeat for the z axis
        positions[i].z = positions[i].z - PBC_Box_Length

    if positions[i].z < -PBC_Box_Length/2:
        positions[i].z = positions[i].z + PBC_Box_Length


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



def Are_Bonded(i, j):                                                                   # checks if the two atoms are a bond, called a 1-2 exclusion in MD
    for bond in bonds:
        if (bond[0] == i and bond[1] == j) or (bond[0] == j and bond[1] == i):
            return True
    return False

# caculate physics stuff
def Calc_Physics():
    E_pot = 0
    for i in range(No_Balls):                                               # calculates distance between all pairs of points, avoiding common pairs
        for j in range(i+1, No_Balls):           
            
            if Are_Bonded(i, j):                                            # skip bonded pairs
                continue
            
            cutoff = 2.5 * sigma                           # cut off value to not calculate particles far away.
            r_vector = positions[j] - positions[i]                          # calc resultant vector between two points
            PBC_Box_For_Vectors(r_vector)                                   # Checks the copies around the real sim
            dist = mag(r_vector)                                            # calc distance
            
            if dist > cutoff:                                               # if dist > than cutoff, skip that pair calc
                continue
            
            direction = norm(r_vector)                                      # unit vector
            #print("dist between: "+str(i)+"  and "+str(j)+": "+str(dist))
            
            # Calculate the lennard jones force
            #F = Calc_LJ(dist, direction)                                    # Direction needed to apply direction to the scalar Force
            #forces[i] = forces[i] + F                                       # store forces on particles
            #forces[j] = forces[j] - F

            # calculate total energy
            U_shift = 4*epsilon*(((sigma/cutoff)**12)-((sigma/cutoff)**6))  # remove the small changes in pot e when you cut off particles from a small pot e to 0
            E_pot = E_pot + 4*epsilon*(((sigma/dist)**12)-((sigma/dist)**6)) - U_shift # calc the potential energy of LJ forces

    for bond in bonds:
        i = bond[0]
        j = bond[1]
        r0 = bond[2]
        k = bond[3]

        r_vector = positions[j] - positions[i]
        PBC_Box_For_Vectors(r_vector)
        dist = mag(r_vector)

        direction = norm(r_vector)

        F = -k * (dist - r0) * direction

        forces[i] += F
        forces[j] -= F

        E_pot += 0.5 * k * (dist - r0)**2
    
    return E_pot

# Initialise forces before simulation starts
for i in range(No_Balls):
    forces[i] = vector(0,0,0)

E_pot = Calc_Physics()
 
step=0   
# sim loop
while True:
    rate(120)
    for _ in range(10):                                                     # Run physics faster
        # calc the next particle position with verlet equations
        # first half-step velocity + position + colours 
        
        for i in range(No_Balls):
            A = forces[i] / masses[i]                                       # calculate acceleration
            velocities[i] += 0.5 * A * dt                                   # half-step velocities
            
        for i in range(No_Balls):
            positions[i] += velocities[i] * dt                              # update positions: move

        
        for i in range(No_Balls):
            PBC_Box(i)                                                      # update copy boxes
        
        
        for i in range(No_Balls):                                           # reset forces
            forces[i] = vector(0,0,0)

        E_pot = Calc_Physics()                                              # run the physics (Calc distances, Forces, and Potential energy) 

        for i in range(No_Balls):                                           # update second half of velocities, required for verlet - this is after the "initial push" of the first velocity after calculating the forces
            A_new = forces[i] / masses[i]
            velocities[i] += 0.5 * A_new * dt
 
            
        K_E = Calc_KE()
        E_total = K_E + E_pot                                            
        
        step += 1                                                           # print every so often

        if step % 100 == 0:
            # DEBUG: check bond length
            for bond in bonds:
                i = bond[0]
                j = bond[1]
                r = mag(positions[j] - positions[i])
                print(f"Bond {i}-{j}: {r:.4f}")
            

        
        E_total = 0                                                         # reset energies after each timestep compleated
        K_E = 0
        E_pot = 0

    for i in range(No_Balls):                                               # update balls
            balls[i].pos = positions[i]
            speed = mag(velocities[i])                                      # calc speed
            ColourPS = speed/2.5                                            # colour %
            ColourPS = min(ColourPS, 1)
            balls[i].color = vector((ColourPS),0,(1-(ColourPS)))            # change ball colour
            balls[i].trail_color = vector((ColourPS),0,(1-(ColourPS)))      # change trail colour

            for idx, bond in enumerate(bonds):                              # update the bond sticks positions
                bond_visuals[idx].pos = positions[bond[0]]
                bond_visuals[idx].axis = positions[bond[1]] - positions[bond[0]]
    