import math
import time
import numpy as np
from periodic_boundary import minimum_image


def bond_forces(system):

    atoms = system.atoms
    bonds = system.bonds
    global pot_e_tota
    pot_e_total = system.pot_e_total
    PBC_BOX_LENGTH = system.PBC_BOX_LENGTH


    for bond in bonds:
        a = bond.a1
        b = bond.a2
        r0 = bond.ideal_dist
        k = bond.K

        # compute the bond vector and bond lengths
        r_vec = minimum_image(atoms[b].pos - atoms[a].pos, PBC_BOX_LENGTH)
        r = np.linalg.norm(r_vec)

        if r < 1e-12:                              # avoid devide by zero just in case
            continue

        r_hat = r_vec / r                     # bond direction

        F = 2 * k * (r - r0) * r_hat                # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds



        atoms[a].force += F                          # apply Newtons third law - equal and opposite forces
        atoms[b].force -= F

        pot_e_total += k * (r - r0)**2                # calc the harmonic bond potential in function 2 type bond

def angle_forces(system):

    atoms = system.atoms
    bond_angles = system.bond_angles
    global pot_e_total
    pot_e_total = system.pot_e_total
    PBC_BOX_LENGTH = system.PBC_BOX_LENGTH
    
    for angle in bond_angles:

        theta = angle.theta
        ideal_ang = angle.ideal_ang
        atomA = angle.first
        atomB = angle.second
        atomC = angle.third
        k = angle.k

        BA = atoms[atomA].pos - atoms[atomB].pos                                # find vector of B to A
        BC = atoms[atomC].pos - atoms[atomB].pos                                # find vector of B to C

        BA = minimum_image(BA,PBC_BOX_LENGTH)                                                 # update the coppy vectors
        BC = minimum_image(BC,PBC_BOX_LENGTH)
        
        r_BA = np.linalg.norm(BA)                                                          # get the lengths, as force depends on direction, and how long the arms are
        r_BC = np.linalg.norm(BC)
        
        if r_BA < 1e-12 or r_BC < 1e-12:
            continue
        
        theta_cos = np.dot(BA, BC) / (r_BA * r_BC)                                          # calculate the dot product hrer
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

        pot_e_total += 0.5 * k * (theta - ideal_ang)**2                # update the potential energy change from the angle

def torsion_forces(system):

    atoms = system.atoms
    torsion_angles = system.torsion_angles
    global pot_e_total
    pot_e_total = system.pot_e_total
    PBC_BOX_LENGTH = system.PBC_BOX_LENGTH
    
    for angle in torsion_angles:
        a1 = angle.a1                           # 4 atoms involved in the torsion bond
        a2 = angle.a2
        a3 = angle.a3
        a4 = angle.a4


        b1 = atoms[a2].pos - atoms[a1].pos      # directions of the three bonds
        b2 = atoms[a3].pos - atoms[a2].pos
        b3 = atoms[a4].pos - atoms[a3].pos

        b1 = minimum_image(b1,PBC_BOX_LENGTH)                 # update the vectors of coppies
        b2 = minimum_image(b2,PBC_BOX_LENGTH)
        b3 = minimum_image(b3,PBC_BOX_LENGTH)

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
        y = np.dot(np.cross(n1, n2), b2 / b2_mag)

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

        c1 = np.dot(b1, b2) / b2_sq                        # calculates how much b and c lean along the middle bond
        c2 = np.dot(b3, b2) / b2_sq

        f_b = -(1.0 + c1) * f_a + c2 * f_d                # calculate the final forces of b and c, by taking into account the forces of a and d
        f_c = f_c = -(f_a + f_b + f_d)

        atoms[a1].force += f_a                         # apply force to atom a
        atoms[a2].force += f_b                         # apply force to atom b
        atoms[a3].force += f_c                         # apply force to atom c
        atoms[a4].force += f_d                         # apply force to atom d

        pot_e_total += torsion_e  # calculate the potential energy
