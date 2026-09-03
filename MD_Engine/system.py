import cupy as cp
import numpy as np
from constants import *

# atom class
class Atom:
    def __init__(self, index, element):

        props = ELEMENTS[element]

        self.radius = props["radius"]
        self.element = element
        self.base_colour = props["color"]

        self.neighbours = []
        self.last_neighbour_reference = np.zeros(3)
# molecule class
class Molecule:
    def __init__(self, atom_indices):
        self.atom_indices = atom_indices


class System:
    def __init__(self):

        # Dynamic state
        self.positions = []
        self.velocities = []
        self.forces = []

        # Constant per-atom data
        self.masses = []
        self.charges = []
        self.ff_atom_types = []
        self.ff_type_ids = []
        self.num_ff_types = []
        self.lj_c6_matrix = []
        self.lj_c12_matrix = []
        self.ff_type_to_id = []
        self.lj_shift_matrix = []

        # Topology
        self.atoms = []
        self.molecules = []
        self.n_atoms = 0

        # Bonds
        self.bond_a = []
        self.bond_b = []
        self.bond_r0 = []
        self.bond_k = []
        self.step = 0

        # Angles

        self.angle_i = []
        self.angle_j = []
        self.angle_k = []
        self.angle_theta0 = []
        self.angle_kconst = []
        self.n_angles = 0

        # Torsions

        self.torsion_i = []
        self.torsion_j = []
        self.torsion_k = []
        self.torsion_l = []
        self.torsion_kterm  = []
        self.torsion_n = []
        self.torsion_delta =[]
        self.n_torsions = 0

        # Non Bonded loop
        self.pair_i = cp.empty(0, dtype=cp.int32)
        self.pair_j = cp.empty(0, dtype=cp.int32)

        self.pair_lj_scale = cp.empty(0, dtype=cp.float64)
        self.pair_coulomb_scale = cp.empty(0, dtype=cp.float64)

        self.rho = cp.zeros(
            (PME_GRID, PME_GRID, PME_GRID),
            dtype=cp.float64
        )

        self.BC = None

        # Energies
        self.potential_energy = cp.array([0.0], dtype=cp.float64)
        self.potential_energy_gpu = cp.array([0.0], dtype=cp.float64)
        self.kinetic_energy = 0.0
        self.real_space_PE = 0.0
        self.reciprocal_PE = 0.0
        self.total_energy = 0.0
        self.self_PE = 0.0
        self.exclusion_PE = 0.0


        self.force_bond = None
        self.force_angle = None
        self.force_torsion = None
        self.force_lj = None
        self.force_real_coulomb = None
        self.force_exclusion = None
        self.force_reciprocal = None




        # Graph Stuff
        self.step = 0
        self.initial_total_energy = 0.0
        self.average_total_energy = 0.0
        self.steps = []
        self.avg_sys_energy_values = []

        # Visuals
        self.balls = []
        self.bond_visuals = []
        self.n_bonds = 0

    def minimum_image(self, r_vector):
        return r_vector - PBC_BOX_LENGTH * cp.round(r_vector / PBC_BOX_LENGTH)

    def minimum_image_cpu(self, r_vector):
        return r_vector - PBC_BOX_LENGTH * np.round(r_vector / PBC_BOX_LENGTH)