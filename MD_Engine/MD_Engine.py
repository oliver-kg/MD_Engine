from vpython import canvas, dist, rate, sphere, vector, color,mag, norm,dot
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
E_pot = 0
K_E = 0
E_total = 0
No_Balls = 27
dt = 0.002
sigma = 1.2
epsilon = 0.5
PBC_Box_Length = math.sqrt(No_Balls)

# create pos
count = 0
n_per_axis = math.ceil(No_Balls ** (1/3))
spacing = 1.2 * sigma
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

def PBC_Box(i):
    wrapped = False                                                     # wrapped to stop trails from forming when warping
    if positions[i].x > PBC_Box_Length/2:                               # checks if the x,y or z position are out of bounds
        positions[i].x = positions[i].x - PBC_Box_Length                # adjust the new current position
        wrapped = True
    if positions[i].x < -PBC_Box_Length/2:                              # repeat for the two ways the x axis can move, eg, y and z direction
        positions[i].x = positions[i].x + PBC_Box_Length
        wrapped = True

    if positions[i].y > PBC_Box_Length/2:                               # repeat for the y axis
        positions[i].y = positions[i].y - PBC_Box_Length
        wrapped = True
    if positions[i].y < -PBC_Box_Length/2:
        positions[i].y = positions[i].y + PBC_Box_Length
        wrapped = True

    if positions[i].z > PBC_Box_Length/2:                               # repeat for the z axis
        positions[i].z = positions[i].z - PBC_Box_Length
        wrapped = True
    if positions[i].z < -PBC_Box_Length/2:
        positions[i].z = positions[i].z + PBC_Box_Length
        wrapped = True

    if wrapped:
        balls[i].clear_trail()                                          # clears the line

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

def Calc_Physics():
    E_pot = 0
    for i in range(No_Balls):                                               # calculates distance between all pairs of points, avoiding common pairs
        for j in range(i+1, No_Balls):           
            cutoff = (PBC_Box_Length-1)/2 * sigma                                            # cut off value to not calculate particles far away.
            r_vector = positions[j] - positions[i]                          # calc resultant vector between two points
            PBC_Box_For_Vectors(r_vector)                              # Checks the copies around the real sim
            dist = mag(r_vector)                                            # calc distance
            
            if dist > cutoff:                                               # if dist > than cutoff, skip that pair calc
                continue
            
            direction = norm(r_vector)                                      # unit vector
            #print("dist between: "+str(i)+"  and "+str(j)+": "+str(dist))

            # Calculate the lennard jones force
            F = Calc_LJ(dist, direction)                                    # Direction needed to apply direction to the scalar Force

            forces[i] = forces[i] + F                                       # store forces on particles
            forces[j] = forces[j] - F

            # calculate total energy
            U_shift = 4*epsilon*(((sigma/cutoff)**12)-((sigma/cutoff)**6))  # remove the small changes in pot e when you cut off particles from a small pot e to 0
            E_pot = E_pot + 4*epsilon*(((sigma/dist)**12)-((sigma/dist)**6)) - U_shift # calc the potential energy of LJ forces
    
    return E_pot
 
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
            positions[i] += velocities[i] * dt                              # update positions
            PBC_Box(i)

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
            print(E_total)
        
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
    