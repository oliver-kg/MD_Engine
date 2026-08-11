from vpython import vector
# ------------------------
# PME PARAMETERS
# ------------------------

PME_ALPHA = 3.5
PME_GRID = 32
BSPLINE_ORDER = 4

# ------------------------
# NON BONDED PARAMETERS
# ------------------------
# LJ
LJ_SIGMA = 0.3
LJ_EPSILON = 0.08
LJ_CUTOFF = 2.5 * LJ_SIGMA                                                                  # compute LJ forces only within a cutoff reigon, as forces are too weak anyway
LJ_ENERGY_SHIFT = 4 * LJ_EPSILON * ((LJ_SIGMA / LJ_CUTOFF)**12 - (LJ_SIGMA / LJ_CUTOFF)**6) # used so the LJ potential smoothly becomes 0 instead of cutting off

# COULOMBS

COULOMB_CONSTANT = 138.935456
REAL_CUTOFF = 1.2
SKIN_CUTOFF = 0.1

# ------------------------
# SIMULATION PARAMETERS
# ------------------------

TIME_STEP = 0.00005
PBC_BOX_LENGTH = 3
MOLECULE_NUMBERS = 3
NEIGHBOUR_CUTOFF = SKIN_CUTOFF + REAL_CUTOFF
CELL_SIZE = NEIGHBOUR_CUTOFF

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
