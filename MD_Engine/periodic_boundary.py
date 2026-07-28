import numpy as np


# minimum image periodic correction - turns the vector into the nearest image vector
def minimum_image(r_vector, PBC_BOX_LENGTH):

    half = PBC_BOX_LENGTH * 0.5
    box = PBC_BOX_LENGTH

    if r_vector[0] > half:
        r_vector[0] -= box
    elif r_vector[0] < -half:
        r_vector[0] += box

    if r_vector[1] > half:
        r_vector[1] -= box
    elif r_vector[1] < -half:
        r_vector[1] += box

    if r_vector[2] > half:
        r_vector[2] -= box
    elif r_vector[2] < -half:
        r_vector[2] += box

    return r_vector


# checks if a molecule has reached the box boundery and needs warping
def wrap_molecules(system):
    atoms = system.atoms
    molecules = system.molecules
    PBC_BOX_LENGTH = system.PBC_BOX_LENGTH

    half = PBC_BOX_LENGTH / 2

    for mol in molecules:
        ref_pos = atoms[mol.atom_indices[0]].pos
        shift = np.zeros(3, dtype=np.float64)

        if ref_pos[0] > half:
            shift[0] -= PBC_BOX_LENGTH
        elif ref_pos[0] < -half:
            shift[0] += PBC_BOX_LENGTH

        if ref_pos[1] > half:
            shift[1] -= PBC_BOX_LENGTH
        elif ref_pos[1] < -half:
            shift[1] += PBC_BOX_LENGTH

        if ref_pos[2] > half:
            shift[2] -= PBC_BOX_LENGTH
        elif ref_pos[2] < -half:
            shift[2] += PBC_BOX_LENGTH

        if shift[0] != 0 or shift[1] != 0 or shift[2] != 0:
            for i in mol.atom_indices:
                atoms[i].pos += shift

