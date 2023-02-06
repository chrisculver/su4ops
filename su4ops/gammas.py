import numpy as np

NS=4

def get_basis(name):
    if name=='dirac-pauli':
        return get_dirac_pauli_gammas()
    elif name=='weyl-chiral':
        return get_weyl_chiral_gammas()
    elif name=='degrand-rossi':
        return get_degrand_rossi_gammas()
    elif name=='chroma':
        return get_chroma_gammas()
    else:
        raise ValueError('Unknown basis, you requested {} but it is not implemented'.format(name))


def get_dirac_pauli_gammas():
    if NS!=4:
        raise ValueError('hardcoded for NS=4...')

    gammas=[]
    for i in range(6):
        gammas.append(np.zeros((NS,NS),dtype=complex))

    gammas[1][0,3] = -1j
    gammas[1][1,2] = -1j
    gammas[1][2,1] = 1j
    gammas[1][3,0] = 1j

    gammas[2][0,3] = -1
    gammas[2][1,2] = 1
    gammas[2][2,1] = 1
    gammas[2][3,0] = -1

    gammas[3][0,2] = -1j
    gammas[3][1,3] = 1j
    gammas[3][2,0] = 1j
    gammas[3][3,1] = -1j

    gammas[4][0,0] = 1
    gammas[4][1,1] = 1
    gammas[4][2,2] = -1
    gammas[4][3,3] = -1

    gammas[5][0,2] = 1
    gammas[5][1,3] = 1
    gammas[5][2,0] = 1
    gammas[5][3,1] = 1

    return gammas



def get_weyl_chiral_gammas():
    if NS!=4:
        raise ValueError('hardcoded for NS=4...')

    gammas=[]
    for i in range(6):
        gammas.append(np.zeros((NS,NS),dtype=complex))

    gammas[1][0,3] = -1j
    gammas[1][1,2] = -1j
    gammas[1][2,1] = 1j
    gammas[1][3,0] = 1j

    gammas[2][0,3] = -1
    gammas[2][1,2] = 1
    gammas[2][2,1] = 1
    gammas[2][3,0] = -1

    gammas[3][0,2] = -1j
    gammas[3][1,3] = 1j
    gammas[3][2,0] = 1j
    gammas[3][3,1] = -1j

    gammas[4][0,2] = 1
    gammas[4][1,3] = 1
    gammas[4][2,0] = 1
    gammas[4][3,1] = 1

    gammas[5][0,0] = -1
    gammas[5][1,1] = -1
    gammas[5][2,2] = 1
    gammas[5][3,3] = 1

    return gammas

def get_degrand_rossi_gammas():
    if NS!=4:
        raise ValueError('hardcoded for NS=4...')

    gammas=[]
    for i in range(6):
        gammas.append(np.zeros((NS,NS),dtype=complex))

    gammas[1][0,3] = 1j
    gammas[1][1,2] = 1j
    gammas[1][2,1] = -1j
    gammas[1][3,0] = -1j

    gammas[2][0,3] = -1
    gammas[2][1,2] = 1
    gammas[2][2,1] = 1
    gammas[2][3,0] = -1

    gammas[3][0,2] = 1j
    gammas[3][1,3] = -1j
    gammas[3][2,0] = -1j
    gammas[3][3,1] = 1j

    gammas[4][0,2] = 1
    gammas[4][1,3] = 1
    gammas[4][2,0] = 1
    gammas[4][3,1] = 1

    gammas[5][0,0] = -1
    gammas[5][1,1] = -1
    gammas[5][2,2] = 1
    gammas[5][3,3] = 1

    return gammas


def get_chroma_gammas():
    if NS!=4:
        raise ValueError('hardcoded for NS=4...')

    gammas=[]
    for i in range(6):
        gammas.append(np.zeros((NS,NS),dtype=complex))

    gammas[1][0,3] = -1j
    gammas[1][1,2] = -1j
    gammas[1][2,1] = 1j
    gammas[1][3,0] = 1j

    gammas[2][0,3] = 1
    gammas[2][1,2] = -1
    gammas[2][2,1] = -1
    gammas[2][3,0] = 1

    gammas[3][0,2] = -1j
    gammas[3][1,3] = 1j
    gammas[3][2,0] = 1j
    gammas[3][3,1] = -1j

    gammas[4][0,2] = -1
    gammas[4][1,3] = -1
    gammas[4][2,0] = -1
    gammas[4][3,1] = -1

    gammas[5][0,0] = 1
    gammas[5][1,1] = 1
    gammas[5][2,2] = -1
    gammas[5][3,3] = -1

    return gammas