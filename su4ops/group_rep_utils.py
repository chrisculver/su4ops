from su4ops.utils import fullVec_to_reduced
import numpy as np
import sympy as sp

def makeRepMat(gElem, basis, mapping):
    rotMat = []
    for b in basis:
        rotMat.append(fullVec_to_reduced(
            b.spatial_rotate(gElem), mapping, len(basis)))

    return np.transpose(np.array(rotMat))

def projectorMat(irrep, rep, group, row=0):
    dim = len(rep[0])
    res = np.zeros((dim, dim), dtype=complex)
    for i, elem in enumerate(group.elements):
        #print(type(float(elem.irreps[irrep][row,row])))
        res += complex(elem.irreps[irrep][row, row])*np.transpose(rep[i])
    return (res*len(irrep)/len(group.elements)).round(8)


def operators(irrep, rep, group, row=0):
    proj = projectorMat(irrep, rep, group, row)
    rrefMat = sp.Matrix(proj).rref()[0].tolist()
    for row in rrefMat[:]:
        if (row == [0 for i in range(len(row))]):
            rrefMat.remove(row)
    rrefMat = np.array(rrefMat).astype(complex)
    return rrefMat