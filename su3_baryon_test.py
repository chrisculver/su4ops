import quark as Q
import numpy as np
import sympy as sp
import FiniteVolumeGroups as fvg
from constants import *

def makeRepMat(basis, gElem):
    res = []
    for b in basis:
        res.append(b.spatial_rotate(gElem))
    return res

oh = fvg.cubic.Oh()

quark = lambda spin : Q.Quark({
    'bar': False,
    'flavor': 1,
    'color': 0,
    'spin': spin,
})

basis = []
for s0 in range(0,NS):
    for s1 in range(0,NS):
        for s2 in range(0,NS):
            basis.append(Q.Elemental(1,[quark(s0),quark(s1),quark(s2)]))

len(basis)

coefMat=[]
for b in basis:
    for g in oh.elements:
        coefMat.append(b.spatial_rotate(g))
len(coefMat)
len(coefMat[0])



def get_vectors(proj, basis):
    rrefMat = np.matrix(sp.Matrix(proj).rref()[0],dtype=float)
    for i in range(0,30):
        print(rrefMat[i])

    #print(np.matmul(rrefMat,basis))

get_vectors(coefMat, basis)
