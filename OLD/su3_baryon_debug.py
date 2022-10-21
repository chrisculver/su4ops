import su4ops.utils as utils
import quark as Q
import numpy as np
import sympy as sp
import FiniteVolumeGroups as fvg
from constants import NS

oh = fvg.cubic.Oh()


def fullVec_to_reduced(vec, basis, extraBasis):

  newVec = [0 for b in basis]
  for i, val in enumerate(vec):
    s0 = i % NS
    tmp = i//NS
    s1 = tmp % NS
    tmp = tmp//NS
    s2 = tmp % NS
    tmp = tmp//NS
    s3 = tmp
    elemental = Q.Elemental(1, [quark(s0), quark(s1), quark(s2), quark(s3)])
    newIdx = 0
    if elemental in extraBasis:
      newIdx = basis.index(extraBasis[elemental])
    else:
      newIdx = basis.index(elemental)

    newVec[newIdx] += val
  return np.array(newVec).round(8)


def makeRepMat(basis, extraBasis, gElem, id):
  rotMat = []
  refMat = []
  for b in basis:
    rotMat.append(fullVec_to_reduced(
      b.spatial_rotate(gElem), basis, extraBasis))
    refMat.append(fullVec_to_reduced(b.spatial_rotate(id), basis, extraBasis))

  res = np.zeros((len(basis), len(basis)), dtype=complex)
  for r in range(len(basis)):
    for c in range(len(basis)):
      res[r, c] = np.dot(rotMat[r], refMat[c])

  return np.transpose(res)


def quark(spin): return Q.Quark({
    'bar': False,
    'flavor': 0,
    'color': 0,
    'spin': spin,
})


basis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
        newOp = Q.Elemental(1, [quark(s0), quark(s1), quark(s2)])
        foundRelated = False
        for b in basis:
          if newOp.quarks in utils.permutations(b.quarks):
            foundRelated = True
            extraBasis[newOp] = Q.Elemental(1, b.quarks)
        if not foundRelated:
          basis.append(newOp)

len(basis)
len(extraBasis)

coefMat = []
for b in basis:
    for g in oh.elements:
        coefMat.append(b.spatial_rotate(g))
len(coefMat)
len(coefMat[0])


def get_vectors(proj, basis):
    rrefMat = np.matrix(sp.Matrix(proj).rref()[0], dtype=float)
    for i in range(0, 30):
        print(rrefMat[i])

    #print(np.matmul(rrefMat,basis))


get_vectors(coefMat, basis)
