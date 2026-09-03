from vpython import vector
from scipy.special import erfcinv
import math

# ------------------------
# NON BONDED PARAMETERS
# ------------------------
# LJ
LJ_SIGMA = 0.3
LJ_EPSILON = 0.08
LJ_CUTOFF = 2.5 * LJ_SIGMA                                                                  # compute LJ forces only within a cutoff reigon, as forces are too weak anyway

# ------------------------
# COULOMBS
# ------------------------

COULOMB_CONSTANT = 138.935456
REAL_CUTOFF = 1.2

# ------------------------
# SIMULATION PARAMETERS
# ------------------------

TIME_STEP = 0.00005
PBC_BOX_LENGTH = 4
MOLECULE_NUMBERS = 3

# ------------------------
# NEIGHBOURS
# ------------------------

SKIN_CUTOFF = 0.1

NONBONDED_CUTOFF = max(LJ_CUTOFF, REAL_CUTOFF)

NEIGHBOUR_CUTOFF = NONBONDED_CUTOFF + SKIN_CUTOFF

# ------------------------
# PME PARAMETERS
# ------------------------

def next_good_fft_size(n):

    good = [2, 3, 4, 5, 6, 8, 9, 10, 12,
            15, 16, 18, 20, 24, 25, 27,
            30, 32, 36, 40, 45, 48, 50,
            54, 60, 64, 72, 75, 80, 81, 90,
            96, 100, 108, 120, 125, 128]

    for x in good:
        if x >= n:
            return x

    return 2 ** math.ceil(math.log2(n))


PME_REAL_TOL = 1e-8
PME_ALPHA = erfcinv(PME_REAL_TOL) / REAL_CUTOFF

PME_GRID_SPACING = 0.10

PME_GRID = next_good_fft_size(
    math.ceil(PBC_BOX_LENGTH / PME_GRID_SPACING)
)

BSPLINE_ORDER = 4

# ------------------------
# CELL LIST
# ------------------------

NUM_CELLS = max(1, int(PBC_BOX_LENGTH // NEIGHBOUR_CUTOFF))

CELL_SIZE = PBC_BOX_LENGTH / NUM_CELLS
HALF_BOX = PBC_BOX_LENGTH * 0.5



ELEMENTS = {
    "H": {"mass": 1.0080,  "radius": 0.03, "color": vector(1, 1, 1)},
    "C": {"mass": 12.0110, "radius": 0.05, "color": vector(0.2, 0.2, 0.2)},
    "O": {"mass": 15.9994, "radius": 0.04, "color": vector(1, 0, 0)},
    "N": {"mass": 14.0067, "radius": 0.047, "color": vector(0, 0, 1)},
    "S": {"mass": 32.0, "radius": 0.15, "color": vector(1, 1, 0)},
    "Cl": {"mass": 35.453, "radius": 0.05, "color": vector(0, 1, 0)}
}
