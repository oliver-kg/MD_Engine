extern "C" __global__

void torsion_kernel(
    const float *positions,
    float *forces,
    const int *torsion_i,
    const int *torsion_j,
    const int *torsion_k,
    const int *torsion_l,
    float box_length,
    float * potential_energy,
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

    int a1 = torsion_i[torsion]
    int a2 = torsion_j[torsion]
    int a3 = torsion_k[torsion]
    int a4 = torsion_l[torsion]

    // -------------------------------------------------
    // 3. Read atom coordinates
    // -------------------------------------------------
    
    float x1 = positions[3 * a1 + 0];
    float y1 = positions[3 * a1 + 1];
    float z1 = positions[3 * a1 + 2];

    float x2 = positions[3 * a2 + 0];
    float y2 = positions[3 * a2 + 1];
    float z2 = positions[3 * a2 + 2];

    float x3 = positions[3 * a3 + 0];
    float y3 = positions[3 * a3 + 1];
    float z3 = positions[3 * a3 + 2];

    float x4 = positions[3 * a4 + 0];
    float y4 = positions[3 * a4 + 1];
    float z4 = positions[3 * a4 + 2];

    // -------------------------------------------------
    // 4. Bond vector - same as (system.positions[b] - system.positions[a]) ...
    // -------------------------------------------------

    float ABx = x2 - x1;
    float ABy = y2 - y1;
    float ABz = z2 - z1;

    float BCx = x3 - x2;
    float BCy = y3 - y2;
    float BCz = z3 - z2;

    float CDx = x4 - x3;
    float CDy = y4 - y3;
    float CDz = z4 - z3;

    // -------------------------------------------------
    // 5. Minimum image - same as system.minimum_image()
    // -------------------------------------------------

    ABx -=  box_length * roundf(BAx / box_length);
    ABy -=  box_length * roundf(BAy / box_length);
    ABz -=  box_length * roundf(BAz / box_length);

    BCx -=  box_length * roundf(BCx / box_length);
    BCy -=  box_length * roundf(BCy / box_length);
    BCz -=  box_length * roundf(BCz / box_length);

    CDx -=  box_length * roundf(CDx / box_length);
    CDy -=  box_length * roundf(CDy / box_length);
    CDz -=  box_length * roundf(CDz / box_length);


    
}