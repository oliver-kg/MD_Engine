extern "C" __global__

void angle_kernel(
    const float *positions,
    float *forces,
    const int *angle_i,
    const int *angle_j,
    const int *angle_k,
    const float *angle_theta0,
    const float *angle_kconst,
    float box_length,
    float * potential_energy,
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
    float k = angle_kconst[angle];
    float ideal_ang = angle_theta0[angle];

    // -------------------------------------------------
    // 3. Read atom coordinates
    // -------------------------------------------------
 
    float x1 = positions[3 * a + 0];
    float y1 = positions[3 * a + 1];
    float z1 = positions[3 * a + 2];

    float x2 = positions[3 * b + 0];
    float y2 = positions[3 * b + 1];
    float z2 = positions[3 * b + 2];

    float x3 = positions[3 * c + 0];
    float y3 = positions[3 * c + 1];
    float z3 = positions[3 * c + 2];

    // -------------------------------------------------
    // 4. Bond vector - same as (system.positions[b] - system.positions[a]) ...
    // -------------------------------------------------

    float BAx = x1 - x2;
    float BAy = y1 - y2;
    float BAz = z1 - z2;

    float BCx = x3 - x2;
    float BCy = y3 - y2;
    float BCz = z3 - z2;

    // -------------------------------------------------
    // 5. Minimum image - same as system.minimum_image()
    // -------------------------------------------------

    BAx -= box_length * roundf(BAx / box_length);
    BAy -= box_length * roundf(BAy / box_length);
    BAz -= box_length * roundf(BAz / box_length);

    BCx -= box_length * roundf(BCx / box_length);
    BCy -= box_length * roundf(BCy / box_length);
    BCz -= box_length * roundf(BCz / box_length);

    // -------------------------------------------------
    // 6. Distance - same as cp.linalg.norm(BA, axis=1)...
    // -------------------------------------------------

    float r_BA2 = BAx * BAx + BAy * BAy + BAz * BAz;
    float r_BC2 = BCx * BCx + BCy * BCy + BCz * BCz;

    if (r_BA2 < 1e-24f || r_BC2 <  1e-24f)
        return;

    // -------------------------------------------------
    // 7. Calculate theta_cos stuff - same as theta_cos /= r_BA * r_BC ...
    // -------------------------------------------------

    float dot = BAx * BCx + BAy * BCy + BAz * BCz;

    float r_BA = sqrtf(r_BA2);                          // bond lenths not squared
    float r_BC = sqrtf(r_BC2);

    float theta_cos = dot / (r_BA * r_BC);

    // -------------------------------------------------
    // 8. Clip the values to make sure they arent extreme - same as theta_cos = cp.clip(theta_cos,-1,1)
    // -------------------------------------------------

    theta_cos = fminf(1.0f, theta_cos);
    theta_cos = fmaxf(-1.0f, theta_cos);

    // -------------------------------------------------
    // 9. Calc theta and sin_theta - same as theta = cp.arccos(theta_cos) ...
    // -------------------------------------------------
    float theta = acosf(theta_cos);
    float sin_theta = sinf(theta);

    if (abs(sin_theta) < 1e-8f)                         // makes sure that sin_theta isnt near 0 when dividing
        return;

    // -------------------------------------------------
    // 10. Calculating the forces - sames as f_a ...
    // -------------------------------------------------

    float dU_dtheta = k * (theta - ideal_ang);
    float common = -(dU_dtheta / sin_theta);
    float inv_rBA_rBC = (1.0f / (r_BA * r_BC));
    float inv_rBA2 = (theta_cos / r_BA2);
    float inv_rBC2 = (theta_cos / r_BC2);

    float FAx = common * (inv_rBA2 * BAx - inv_rBA_rBC * BCx);
    float FAy = common * (inv_rBA2 * BAy - inv_rBA_rBC * BCy);
    float FAz = common * (inv_rBA2 * BAz - inv_rBA_rBC * BCz);

    float FCx = common * (inv_rBC2 * BCx - inv_rBA_rBC * BAx);
    float FCy = common * (inv_rBC2 * BCy - inv_rBA_rBC * BAy);
    float FCz = common * (inv_rBC2 * BCz - inv_rBA_rBC * BAz);

    float FBx = -(FAx + FCx);
    float FBy = -(FAy + FCy);
    float FBz = -(FAz + FCz);

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

    float energy = 0.5 * k * ((theta - ideal_ang) * (theta - ideal_ang));

    atomicAdd(potential_energy, energy);
}