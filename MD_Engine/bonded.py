import math
import cupy as cp


def bond_forces(system):
    positions = system.positions
    forces = system.forces
    bond_a = system.bond_a
    bond_b = system.bond_b
    bond_r0 = system.bond_r0
    bond_k = system.bond_k
    bond_e = 0

    # compute the bond vector and bond lengths
    r_vec = system.minimum_image_many(positions[bond_b] - positions[bond_a])
    r = cp.linalg.norm(r_vec, axis=1)

    # avoid dividing by 0
    valid = r > 1e-12

    bond_a = bond_a[valid]
    bond_b = bond_b[valid]
    bond_k = bond_k[valid]
    bond_r0 = bond_r0[valid]

    r = r[valid]
    r_vec = r_vec[valid]

    r_hat = r_vec / r[:,None]               # calc direction

    F_mag = 2 * bond_k * (r - bond_r0)      # harmonic restoring force (bonds) - direction is used to turn scalar into vector. Switched to function 2 type bonds
    F_vec = F_mag[:, None] * r_hat               

    cp.add.at(forces, bond_a,  F_vec)     # apply Newtons third law - equal and opposite forces
    cp.add.at(forces, bond_b, -F_vec)

    bond_e = cp.sum(bond_k * (r - bond_r0)**2)                # calc the harmonic bond potential in function 2 type bond

    return bond_e

def angle_forces(system):

    positions = system.positions
    forces = system.forces
    bond_angles = system.bond_angles
    angle_e = 0
    
    for angle in bond_angles:

        theta = angle.theta
        ideal_ang = angle.ideal_ang
        atomA = angle.first
        atomB = angle.second
        atomC = angle.third
        k = angle.k

        BA = positions[atomA] - positions[atomB]                                # find vector of B to A
        BC = positions[atomC] - positions[atomB]                                # find vector of B to C

        BA = system.minimum_image(BA)                                                 # update the coppy vectors
        BC = system.minimum_image(BC)
        
        r_BA = cp.linalg.norm(BA)                                                          # get the lengths, as force depends on direction, and how long the arms are
        r_BC = cp.linalg.norm(BC)
        
        if r_BA < 1e-12 or r_BC < 1e-12:
            continue
        
        theta_cos = cp.dot(BA, BC) / (r_BA * r_BC)                                          # calculate the dot product hrer
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
    
        forces[atomA] += f_a                                               # update the forces
        forces[atomB] += f_b
        forces[atomC] += f_c

        angle_e += 0.5 * k * (theta - ideal_ang)**2                # update the potential energy change from the angle
    return angle_e

def torsion_forces(system):

    positions = system.positions
    forces = system.forces
    torsion_angles = system.torsion_angles
    torsion_e = 0
    
    for angle in torsion_angles:
        a1 = angle.a1                           # 4 atoms involved in the torsion bond
        a2 = angle.a2
        a3 = angle.a3
        a4 = angle.a4


        b1 = positions[a2] - positions[a1]      # directions of the three bonds
        b2 = positions[a3] - positions[a2]
        b3 = positions[a4] - positions[a3]

        b1 = system.minimum_image(b1)                 # update the vectors of coppies
        b2 = system.minimum_image(b2)
        b3 = system.minimum_image(b3)

        n1 = cp.cross(b1, b2)                      # normal to plane ABC
        n2 = cp.cross(b2, b3)                      # normal to plane BCD

        eps = 1e-12
        n1_sq = cp.dot(n1, n1)
        n2_sq = cp.dot(n2, n2)
        b2_sq = cp.dot(b2, b2)
        b2_mag = cp.linalg.norm(b2)

        if n1_sq < eps or n2_sq < eps or b2_sq < eps:    # makes sure small values dont blow the system up
            continue
        
        x = cp.dot(n1, n2)
        y = cp.dot(cp.cross(n1, n2), b2 / b2_mag)

        psi = math.atan2(y, x)                  # caluclate current psi angle
        angle.psi = psi 

        dV_dpsi = 0
        
        for k, n, delta in angle.terms:                             # repeats for each term, updating the cos graph
            torsion_e += k * (1 + math.cos(n * psi - delta))   # calculate the potential energy
            dV_dpsi -= k * n * math.sin(n * psi - delta)            # how strongly the torsion wants to rotate
        
        fa_pref = dV_dpsi * (b2_mag / n1_sq)           # calculates the geometric scallings of the force
        fd_pref = -dV_dpsi * (b2_mag / n2_sq)

        f_a = fa_pref * n1                              # aligning the forces with the direction to the plane
        f_d = fd_pref * n2

        c1 = cp.dot(b1, b2) / b2_sq                        # calculates how much b and c lean along the middle bond
        c2 = cp.dot(b3, b2) / b2_sq

        f_b = -(1.0 + c1) * f_a + c2 * f_d                # calculate the final forces of b and c, by taking into account the forces of a and d
        f_c = f_c = -(f_a + f_b + f_d)

        forces[a1] += f_a                         # apply force to atom a
        forces[a2] += f_b                         # apply force to atom b
        forces[a3] += f_c                         # apply force to atom c
        forces[a4] += f_d                         # apply force to atom d

    return torsion_e