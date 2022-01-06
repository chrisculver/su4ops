# Note that we are not using the group OhG....
# so we don't have G-parity...

import math
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


def quark(spin): return Q.Quark({
    'bar': False,
    'flavor': 0,
    'color': 0,
    'spin': spin,
})


def quarkBar(spin): return Q.Quark({
    'bar': True,
    'flavor': 1,
    'color': 0,
    'spin': spin,
})


quark(0).spatial_rotate(oh.elements[0])

quark(0).g_parity_rotate()
print(quark(0))

basis = []
for s0 in range(0, NS):
    for s1 in range(0, NS):
        basis.append(Q.Elemental(1, [quarkBar(s1), quark(s0)]))

print(len(basis))
for b in basis:
  print(b)

idRepMat = []
for b in basis:
    idRepMat.append(b.spatial_rotate(oh.elements[0]))

(idRepMat == np.identity(16)).all()

print(idRepMat)
# What about C4z
expectedIs = np.diag([1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])

expectedC4z = np.diag([1, 1j, 1, 1j,
                       -1j, 1, -1j, 1,
                       1, 1j, 1, 1j,
                       -1j, 1, -1j, 1])

expectedC4y = 0.5*np.array([[1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, -1, -1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1],
                            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, -1, -1],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1],
                            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1]])

for i, g in enumerate(oh.elements):
    if g.identifier['direction'] == [0, 0, 1] and g.identifier['angle'] == math.pi/2. and g.identifier['parity'] == 1:
        print("Rz rotation index is {}".format(i))
    if g.identifier['direction'] == [0, 1, 0] and g.identifier['angle'] == math.pi/2. and g.identifier['parity'] == 1:
        print("Ry rotation index is {}".format(i))
    if g.identifier['name'] == "E" and g.identifier['parity'] == -1:
      print("I_s rotation index is {}".format(i))

c4zRepMat = []
for b in basis:
    c4zRepMat.append(b.spatial_rotate(oh.elements[22]))
c4yRepMat = []
for b in basis:
    c4yRepMat.append(b.spatial_rotate(oh.elements[20]))

isRepMat = []
for b in basis:
  isRepMat.append(b.spatial_rotate(oh.elements[24]))

np.allclose(c4zRepMat, expectedC4z)
np.allclose(c4yRepMat, expectedC4y)

np.allclose(isRepMat, expectedIs)

mesonRep = []
for g in oh.elements:
    mesonRep.append(makeRepMat(basis, g))

wG = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      ]

gparityRep = []
for rep in mesonRep:
  gparityRep.append(np.matmul(rep, wG))


def projectorMat(irrep, rep, gRep, group, row=0):
    dim = len(rep[0])
    res = np.zeros((dim, dim), dtype=complex)
    for i, elem in enumerate(group.elements):
        #print(type(float(elem.irreps[irrep][row,row])))
        res += complex(elem.irreps[irrep][row, row])*np.transpose(rep[i])
        res += complex(elem.irreps[irrep][row, row])*np.transpose(gRep[i])
    return res*len(irrep)/len(group.elements)


tmp = projectorMat('A1g', mesonRep, gparityRep, oh)
tmp = tmp.round(8)


def get_vectors(proj, basis):
  rrefMat = sp.Matrix(proj).rref()[0].tolist()
  for row in rrefMat[:]:
    if (row == [0 for i in range(len(row))]):
      rrefMat.remove(row)
  rrefMat = np.array(rrefMat).astype(float)

  return rrefMat


vs = get_vectors(tmp, basis)
for v in vs:
  print(v)

print("{}+{}+{}+{}".format(basis[0], basis[5], basis[10], basis[15]))


tmp = projectorMat('T1g', mesonRep, gparityRep, oh)
tmp = tmp.round(8)
vs = get_vectors(tmp, basis)
for v in vs:
  print(v)

print("{}+{}-{}-{}".format(basis[1], basis[4], basis[11], basis[14]))


tmp = projectorMat('T1u', mesonRep, gparityRep, oh)
tmp = tmp.round(8)
vs = get_vectors(tmp, basis)
for v in vs:
  print(v)
print("{}+{}".format(basis[3], basis[6]))
print("{}+{}".format(basis[9], basis[12]))
