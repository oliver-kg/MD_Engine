extern "C" __global__
void bond_kernel(
    const float *positions, // the * is reciving the adress of the array itself on the GPU instead of an actual array
    float *forces,          // this makes the variable a pointer, pointing to the position in memory
    const int *bond_a,
    const int *bond_b,
    const float *bond_r0,
    const float *bond_k,
    float box_length,
    float *potential_energy,
    int n_bonds)

{
    // -------------------------------------------------
    // 1. Which bond does this thread own?
    // -------------------------------------------------

    int bond = blockIdx.x * blockDim.x + threadIdx.x; // this gives each individual thread a bond to calculate

    if (bond >= n_bonds)
        return;

    // -------------------------------------------------
    // 2. Read bond information
    // -------------------------------------------------

    int a = bond_a[bond]; // gets the values for one specific bond
    int b = bond_b[bond];

    float r0 = bond_r0[bond];
    float k = bond_k[bond];

    // -------------------------------------------------
    // 3. Read atom coordinates
    // The data is stored as a long list of numbers, and every three numbers is equiviant to the x,y and z values
    // -------------------------------------------------

    float x1 = positions[3 * a + 0]; // gets the x pos value of atom a, reads it, and store in the register x1
    float y1 = positions[3 * a + 1]; // gets the y pos value of atom a, reads it, and store in the register y1
    float z1 = positions[3 * a + 2]; // gets the z pos value of atom a, reads it, and store in the register z1

    float x2 = positions[3 * b + 0]; // gets the x pos value of atom b, reads it, and store in the register x2
    float y2 = positions[3 * b + 1]; // gets the y pos value of atom b, reads it, and store in the register y2
    float z2 = positions[3 * b + 2]; // gets the z pos value of atom b, reads it, and store in the register z2

    // -------------------------------------------------
    // 4. Bond vector - same as (system.positions[b] - system.positions[a])
    // -------------------------------------------------

    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    // -------------------------------------------------
    // 5. Minimum image - same as system.minimum_image()
    // -------------------------------------------------

    dx -= box_length * roundf(dx / box_length);
    dy -= box_length * roundf(dy / box_length);
    dz -= box_length * roundf(dz / box_length);

    // -------------------------------------------------
    // 6. Distance - same as cp.linalg.norm(r_vec, axis=1)
    // -------------------------------------------------

    float r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < 1e-24f) // same as valid = r > 1e-12, without taking the sqrt
        return;

    float r = sqrtf(r2);

    float inv_r = 1.0f / r;

    // -------------------------------------------------
    // 7. Harmonic bond force - same as F = 2 * k[:, None] * (r - r0)[:, None] * r_hat ...
    // -------------------------------------------------

    float force_mag = 2.0f * k * (r - r0);

    float Fx = force_mag * dx * inv_r;
    float Fy = force_mag * dy * inv_r;
    float Fz = force_mag * dz * inv_r;

    // -------------------------------------------------
    // 8. Apply Newton's Third Law - same as cp.add.at(system.forces, a,  F) ...
    // -------------------------------------------------

    atomicAdd(&forces[3 * a + 0], Fx); // the & means "give me the adress" of this, and then now it knows where in the memory to write the new force
    atomicAdd(&forces[3 * a + 1], Fy);
    atomicAdd(&forces[3 * a + 2], Fz);

    atomicAdd(&forces[3 * b + 0], -Fx);
    atomicAdd(&forces[3 * b + 1], -Fy);
    atomicAdd(&forces[3 * b + 2], -Fz);

    // -------------------------------------------------
    // 9. Bond energy - same as system.potential_energy += cp.sum(k * (r - r0)**2)
    // -------------------------------------------------

    float energy = k * (r - r0) * (r - r0);

    atomicAdd(potential_energy, energy);
}