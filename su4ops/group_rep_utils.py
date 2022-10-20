from utils import fullVec_to_reduced
import numpy as np
import sympy as sp

def makeRepMat(basis, extraBasis, gElem, id, f1=0, f2=0):
    rotMat = []
    refMat = []
    for b in basis:
        rotMat.append(fullVec_to_reduced(
        b.spatial_rotate(gElem), basis, extraBasis, f1, f2))
        refMat.append(fullVec_to_reduced(
            b.spatial_rotate(id), basis, extraBasis, f1, f2))
    #res = np.zeros((len(basis), len(basis)), dtype=complex)

    #for r in range(len(basis)):
    #  for c in range(len(basis)):
    #    res[r, c] = np.dot(rotMat[r], refMat[c])

    #return np.transpose(res)
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