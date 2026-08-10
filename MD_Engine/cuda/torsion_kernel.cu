extern "C" __global__

void torsion_kernel(
    const double *positions,
    double *forces,
    const int *torsion_i,
    const int *torsion_j,
    const int *torsion_k,
    const int *torsion_l,
    const double *torsion_kterm,
    const int *torsion_n,
    const double *torsion_delta,
    double box_length,
    double * potential_energy,
    int n_torsions
)

{
    // -------------------------------------------------
    // 1. Which bond does this thread own?
    // -------------------------------------------------

    int torsion = blockIdx.x * blockDim.x + threadIdx.x; // this gives each individual thread a bond to calculate

    if (torsion >= n_torsions)
        return;

    // -------------------------------------------------
    // 2. Read bond information
    // -------------------------------------------------

    int a1 = torsion_i[torsion];
    int a2 = torsion_j[torsion];
    int a3 = torsion_k[torsion];
    int a4 = torsion_l[torsion];

    double k = torsion_kterm[torsion];
    int n = torsion_n[torsion];
    double delta = torsion_delta[torsion];

    // -------------------------------------------------
    // 3. Read atom coordinates
    // -------------------------------------------------
    
    double x1 = positions[3 * a1 + 0];
    double y1 = positions[3 * a1 + 1];
    double z1 = positions[3 * a1 + 2];

    double x2 = positions[3 * a2 + 0];
    double y2 = positions[3 * a2 + 1];
    double z2 = positions[3 * a2 + 2];

    double x3 = positions[3 * a3 + 0];
    double y3 = positions[3 * a3 + 1];
    double z3 = positions[3 * a3 + 2];

    double x4 = positions[3 * a4 + 0];
    double y4 = positions[3 * a4 + 1];
    double z4 = positions[3 * a4 + 2];

    // -------------------------------------------------
    // 4. Bond vector - same as (system.positions[b] - system.positions[a]) ...
    // -------------------------------------------------

    double b1x = x2 - x1;
    double b1y = y2 - y1;
    double b1z = z2 - z1;

    double b2x = x3 - x2;
    double b2y = y3 - y2;
    double b2z = z3 - z2;

    double b3x = x4 - x3;
    double b3y = y4 - y3;
    double b3z = z4 - z3;

    // -------------------------------------------------
    // 5. Minimum image - same as system.minimum_image()
    // -------------------------------------------------

    b1x -=  box_length * round(b1x / box_length);
    b1y -=  box_length * round(b1y / box_length);
    b1z -=  box_length * round(b1z / box_length);

    b2x -=  box_length * round(b2x / box_length);
    b2y -=  box_length * round(b2y / box_length);
    b2z -=  box_length * round(b2z / box_length);

    b3x -=  box_length * round(b3x / box_length);
    b3y -=  box_length * round(b3y / box_length);
    b3z -=  box_length * round(b3z / box_length);

    // -------------------------------------------------
    // 6. Calc normals to planes - same as n1 = cp.cross(b1, b2) ...
    // -------------------------------------------------

    double n1x = b1y*b2z - b1z*b2y;
    double n1y = b1z*b2x - b1x*b2z;
    double n1z = b1x*b2y - b1y*b2x;

    double n2x = b2y*b3z - b2z*b3y;
    double n2y = b2z*b3x - b2x*b3z;
    double n2z = b2x*b3y - b2y*b3x;

    // -------------------------------------------------
    // 7. Calc the magnitudes - same as n1_sq = cp.sum(n1*n1, axis=1) ...
    // -------------------------------------------------

    double n1_sq = n1x*n1x + n1y*n1y + n1z*n1z;
    double n2_sq = n2x*n2x + n2y*n2y + n2z*n2z;

    double b2_sq = b2x*b2x + b2y*b2y + b2z*b2z;
    double b2_mag = sqrt(b2_sq);

    if (n1_sq < 1e-12f || n2_sq < 1e-12f || b2_sq < 1e-12f)
        return;

    // -------------------------------------------------
    // 8. Calc x - same as x = cp.sum(n1*n2, axis=1)...
    // -------------------------------------------------
    
        double x = n1x*n2x + n1y*n2y + n1z*n2z;

    // -------------------------------------------------
    // 9. Calc b2_hat - same as b2_hat = b2 / b2_mag ...
    // -------------------------------------------------

    double inv_b2 = 1.0f / b2_mag;

    double b2hatx = b2x * inv_b2;
    double b2haty = b2y * inv_b2;
    double b2hatz = b2z * inv_b2;

    // -------------------------------------------------
    // 10. Cross product - same as cross12 = cp.cross(n1,n2) ...
    // -------------------------------------------------

    double cross12x = n1y*n2z - n1z*n2y;
    double cross12y = n1z*n2x - n1x*n2z;
    double cross12z = n1x*n2y - n1y*n2x;

    // -------------------------------------------------
    // 11. Calc y - same as y = cp.sum(cross12 * b2_hat) ...
    // -------------------------------------------------

    double y = cross12x*b2hatx + cross12y*b2haty + cross12z*b2hatz;

    // -------------------------------------------------
    // 12. Calc psi - same as psi = cp.arctan2(y, x) ...
    // -------------------------------------------------

    double psi = atan2(y, x);

    // -------------------------------------------------
    // 13. Calc torsion energy and derivative - same as torsion_e = k * (1 + cp.cos(n * psi - delta)) ...
    // -------------------------------------------------

    double energy = k * (1.0f + cos(n * psi - delta));
    double dV_dpsi = -k * n * sin(n * psi - delta);

    // -------------------------------------------------
    // 13. Calc prefactors - same as fa_pref = dV_dpsi * b2_mag / n1_sq ...
    // ------------------------------------------------- 

    double inv_n1_sq = 1.0f / n1_sq;
    double inv_n2_sq = 1.0f / n2_sq;

    double fa_pref = dV_dpsi * b2_mag * inv_n1_sq;
    double fd_pref = -dV_dpsi * b2_mag * inv_n2_sq;

    // -------------------------------------------------
    // 14. Calc forces - same as f_a = fa_pref[:,None] * n1 ...
    // ------------------------------------------------- 

    double FAx = fa_pref * n1x;
    double FAy = fa_pref * n1y;
    double FAz = fa_pref * n1z;    

    double FDx = fd_pref * n2x;
    double FDy = fd_pref * n2y;
    double FDz = fd_pref * n2z;

    double inv_b2_sq = 1.0f / b2_sq;
    double c1 = (b1x*b2x + b1y*b2y + b1z*b2z) * inv_b2_sq;
    double c2 = (b3x*b2x + b3y*b2y + b3z*b2z) * inv_b2_sq;

    double FBx = -(1.0f + c1) * FAx + c2 * FDx;
    double FBy = -(1.0f + c1) * FAy + c2 * FDy;
    double FBz = -(1.0f + c1) * FAz + c2 * FDz;

    double FCx = -(FAx + FBx + FDx);
    double FCy = -(FAy + FBy + FDy);
    double FCz = -(FAz + FBz + FDz);

    // -------------------------------------------------
    // 15. Apply the forces - same as cp.add.at(system.forces, a1,  f_a) ...
    // ------------------------------------------------- 

    atomicAdd(&forces[3 * a1 + 0], FAx);
    atomicAdd(&forces[3 * a1 + 1], FAy);
    atomicAdd(&forces[3 * a1 + 2], FAz);

    atomicAdd(&forces[3 * a2 + 0], FBx);
    atomicAdd(&forces[3 * a2 + 1], FBy);
    atomicAdd(&forces[3 * a2 + 2], FBz);

    atomicAdd(&forces[3 * a3 + 0], FCx);
    atomicAdd(&forces[3 * a3 + 1], FCy);
    atomicAdd(&forces[3 * a3 + 2], FCz);

    atomicAdd(&forces[3 * a4 + 0], FDx);
    atomicAdd(&forces[3 * a4 + 1], FDy);
    atomicAdd(&forces[3 * a4 + 2], FDz);

    // -------------------------------------------------
    // 16. Apply the potential energy - same as pot_e_total += cp.sum(torsion_e) ...
    // ------------------------------------------------- 

    atomicAdd(potential_energy, energy);
}