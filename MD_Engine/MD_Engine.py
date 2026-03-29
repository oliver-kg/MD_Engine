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

                  #atom variables
positions = []
velocities = []
forces = []
balls = []
masses = []
E_pot = 0
K_E = 0
E_total = 0
No_Balls = 10
dt = 0.00002
Ideal_Dist = 5       #ideal distance from F=0
K = 0.2             #spring constant (strngth of bond)

# create pos
for i in range(No_Balls):
    rand_x = random.uniform(-2, 2)
    rand_y = random.uniform(-2, 2)
    rand_z = random.uniform(-2, 2)
    pos = vector(rand_x, rand_y, rand_z)
    positions.append(pos)

# create empty velocities
for i in range(No_Balls):
    vel = vector(0, 0, 0)
    velocities.append(vel)

#create empty forces
for i in range(No_Balls):
    force = vector(0, 0, 0)
    forces.append(force)

#create empty masses
for i in range(No_Balls):
    mass = 1
    masses.append(mass)

# create balls
for i in range(No_Balls):
    b = sphere(pos=positions[i], radius=0.1, color=vector(0,0,0), make_trail=True, interval = 1, retain = 500)
    balls.append(b)

def Calc_Diff():
    E_pot = 0
    for i in range(No_Balls-1):                                             #calculates distance between all pairs of points, avoiding common pairs
        for j in range(i, No_Balls-1):
            j=j+1
            r_vector = positions[j] - positions[i]                          # calc resultant vector
            dist = mag(r_vector)                                            # calc distance
            direction = norm(r_vector)                                      # unit vector
            #print("dist between: "+str(i)+"  and "+str(j)+": "+str(dist))

            #Calculate the Spring Force
            F = -K*(dist-Ideal_Dist)*direction                              #Direction needed to apply direction to the scalar Force

            forces[i] = forces[i] + F                                       #store forces on particles
            forces[j] = forces[j] - F

            #calculate total energy
            E_pot = E_pot+(0.5*K*((dist-Ideal_Dist)**2))
    
    return E_pot

for i in range(No_Balls):                                               # initial force calculation
    forces[i] = vector(0,0,0)

E_pot = Calc_Diff()

step=0   
# sim loop
while True:
    rate(60)
    
    for _ in range(200):                                                    #Run physics faster
        # first half-step velocity + position + colours 
        for i in range(No_Balls):
            A = forces[i] / masses[i]                                       #calculate acceleration
            velocities[i] += 0.5 * A * dt                                   #half-step velocities
            positions[i] += velocities[i] * dt                              # update positions
            balls[i].pos = positions[i]

            speed = mag(velocities[i])                                      #calc speed
            ColourPS = speed/40                                             #colour  %
            ColourPS = min(ColourPS, 1)
            balls[i].color = vector((ColourPS),0,(1-(ColourPS)))            #change ball colour
            balls[i].trail_color = vector((ColourPS),0,(1-(ColourPS)))      #change trail colour

    
        for i in range(No_Balls):                                           #reset forces
            forces[i] = vector(0,0,0)

        E_pot = Calc_Diff()                                                 #calc forces

        for i in range(No_Balls):                                           #update second half of velocities - this is after the "initial push" of the first velocity after calculating the forces
            A_new = forces[i] / masses[i]
            velocities[i] += 0.5 * A_new * dt

        for i in range(No_Balls): 
            K_E = K_E + (0.5*masses[i]*(dot(velocities[i], velocities[i])))#use dot product as v^2 to calc ke


        E_total = K_E + E_pot
        step += 1                                                         #print every so often
        if step % 100 == 0:
            print(E_total)
        E_total = 0
        K_E = 0
        E_pot = 0
    