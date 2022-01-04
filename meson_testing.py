# Note that we are not using the group OhD....

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
quarkBar = lambda spin : Q.Quark({
    'bar': True,
    'flavor': 1,
    'color': 0,
    'spin': spin,
})

basis = []
for s0 in range(0,NS):
    for s1 in range(0,NS):
        basis.append(Q.Elemental(1,[quarkBar(s0),quark(s1)]))

print(len(basis))

idRepMat = []
for b in basis:
    idRepMat.append(b.spatial_rotate(oh.elements[0]))

(idRepMat==np.identity(16)).all()


# What about C4z

expectedC4z = np.diag([1,1j,1,1j,
                     -1j,1,-1j,1,
                       1,1j,1,1j,
                     -1j,1,-1j,1])

import math
for i,g in enumerate(oh.elements):
    if g.identifier['direction']==[0,0,1] and g.identifier['angle']==math.pi/2. and g.identifier['parity']==1:
        print("Rz rotation index is {}".format(i))
    if g.identifier['direction']==[0,1,0] and g.identifier['angle']==math.pi/2. and g.identifier['parity']==1:
        print("Ry rotation index is {}".format(i))

c4zRepMat = []
for b in basis:
    c4zRepMat.append(b.spatial_rotate(oh.elements[22]))
c4yRepMat = []
for b in basis:
    c4yRepMat.append(b.spatial_rotate(oh.elements[20]))

np.allclose(c4zRepMat,expectedC4z)

np.diag(c4zRepMat)

c4yRepMat


mesonRep = []
for g in oh.elements:
    mesonRep.append(makeRepMat(basis, g))

def projectorMat(irrep,rep,group,row=0):
    dim = len(rep[0])
    res = np.zeros((dim,dim),dtype=complex)
    for i,elem in enumerate(group.elements):
        #print(type(float(elem.irreps[irrep][row,row])))
        res += complex(elem.irreps[irrep][row,row])*np.transpose(rep[i])
    return res*len(irrep)/len(group.elements)


tmp=projectorMat('A1g',mesonRep,oh)
tmp=tmp.round(8)

print(str(basis[0]+basis[1]))
print(str(1.0*basis[0]+1.0*basis[1]))

def get_vectors(proj, basis):
    rrefMat = sp.Matrix(proj).rref()[0]
    print(rrefMat)
    print(np.matmul(np.matrix(rrefMat,dtype=float),basis))

get_vectors(tmp,basis)

print(basis[0])
print(basis[0] + basis[5])
