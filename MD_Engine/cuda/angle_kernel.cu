extern "C" __global__

void angle_kernel(
    const double *positions,
    double *forces,
    const int *angle_i,
    const int *angle_j,
    const int *angle_k,
    const double *angle_theta0,
    const double *angle_kconst,
    double box_length,
    double * potential_energy,
    int n_angles
)

{
    // -------------------------------------------------
    // 1. Which angle does this thread own?
    // -------------------------------------------------

    int angle = blockIdx.x * blockDim.x + threadIdx.x;

    if (angle >= n_angles)
        return;

    // -------------------------------------------------
    // 2. Read angle information
    // -------------------------------------------------

    int a = angle_i[angle];
    int b = angle_j[angle];
    int c = angle_k[angle];
    double k = angle_kconst[angle];
    double ideal_ang = angle_theta0[angle];

    // -------------------------------------------------
    // 3. Read atom coordinates
    // -------------------------------------------------
 
    double x1 = positions[3 * a + 0];
    double y1 = positions[3 * a + 1];
    double z1 = positions[3 * a + 2];

    double x2 = positions[3 * b + 0];
    double y2 = positions[3 * b + 1];
    double z2 = positions[3 * b + 2];

    double x3 = positions[3 * c + 0];
    double y3 = positions[3 * c + 1];
    double z3 = positions[3 * c + 2];

    // -------------------------------------------------
    // 4. Bond vector - same as (system.positions[b] - system.positions[a]) ...
    // -------------------------------------------------

    double BAx = x1 - x2;
    double BAy = y1 - y2;
    double BAz = z1 - z2;

    double BCx = x3 - x2;
    double BCy = y3 - y2;
    double BCz = z3 - z2;

    // -------------------------------------------------
    // 5. Minimum image - same as system.minimum_image()
    // -------------------------------------------------

    BAx -= box_length * round(BAx / box_length);
    BAy -= box_length * round(BAy / box_length);
    BAz -= box_length * round(BAz / box_length);

    BCx -= box_length * round(BCx / box_length);
    BCy -= box_length * round(BCy / box_length);
    BCz -= box_length * round(BCz / box_length);

    // -------------------------------------------------
    // 6. Distance - same as cp.linalg.norm(BA, axis=1)...
    // -------------------------------------------------

    double r_BA2 = BAx * BAx + BAy * BAy + BAz * BAz;
    double r_BC2 = BCx * BCx + BCy * BCy + BCz * BCz;

    if (r_BA2 < 1e-24f || r_BC2 <  1e-24f)
        return;

    // -------------------------------------------------
    // 7. Calculate theta_cos stuff - same as theta_cos /= r_BA * r_BC ...
    // -------------------------------------------------

    double dot = BAx * BCx + BAy * BCy + BAz * BCz;

    double r_BA = sqrt(r_BA2);                          // bond lenths not squared
    double r_BC = sqrt(r_BC2);

    double theta_cos = dot / (r_BA * r_BC);

    // -------------------------------------------------
    // 8. Clip the values to make sure they arent extreme - same as theta_cos = cp.clip(theta_cos,-1,1)
    // -------------------------------------------------

    theta_cos = min(1.0f, theta_cos);
    theta_cos = max(-1.0f, theta_cos);

    // -------------------------------------------------
    // 9. Calc theta and sin_theta - same as theta = cp.arccos(theta_cos) ...
    // -------------------------------------------------
    double theta = acos(theta_cos);
    double sin_theta = sin(theta);

    if (abs(sin_theta) < 1e-8f)                         // makes sure that sin_theta isnt near 0 when dividing
        return;

    // -------------------------------------------------
    // 10. Calculating the forces - sames as f_a ...
    // -------------------------------------------------

    double dU_dtheta = k * (theta - ideal_ang);
    double common = -(dU_dtheta / sin_theta);
    double inv_rBA_rBC = (1.0f / (r_BA * r_BC));
    double inv_rBA2 = (theta_cos / r_BA2);
    double inv_rBC2 = (theta_cos / r_BC2);

    double FAx = common * (inv_rBA2 * BAx - inv_rBA_rBC * BCx);
    double FAy = common * (inv_rBA2 * BAy - inv_rBA_rBC * BCy);
    double FAz = common * (inv_rBA2 * BAz - inv_rBA_rBC * BCz);

    double FCx = common * (inv_rBC2 * BCx - inv_rBA_rBC * BAx);
    double FCy = common * (inv_rBC2 * BCy - inv_rBA_rBC * BAy);
    double FCz = common * (inv_rBC2 * BCz - inv_rBA_rBC * BAz);

    double FBx = -(FAx + FCx);
    double FBy = -(FAy + FCy);
    double FBz = -(FAz + FCz);

    // -------------------------------------------------
    // 11. Apply Newton's Third Law - same as cp.add.at(system.forces, atomA, f_a) ...
    // -------------------------------------------------

    atomicAdd(&forces[3 * a + 0], FAx);
    atomicAdd(&forces[3 * a + 1], FAy);
    atomicAdd(&forces[3 * a + 2], FAz);

    atomicAdd(&forces[3 * b + 0], FBx);
    atomicAdd(&forces[3 * b + 1], FBy);
    atomicAdd(&forces[3 * b + 2], FBz);

    atomicAdd(&forces[3 * c + 0], FCx);
    atomicAdd(&forces[3 * c + 1], FCy);
    atomicAdd(&forces[3 * c + 2], FCz);

    // -------------------------------------------------
    // 12. Bond potential energy - same as pot_e_total += cp.sum(0.5 * k * (theta - ideal_ang)**2) 
    // -------------------------------------------------

    double energy = 0.5 * k * ((theta - ideal_ang) * (theta - ideal_ang));

    atomicAdd(potential_energy, energy);
}