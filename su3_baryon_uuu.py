import utils
import su4ops.quark as Q
import su4ops.elemental as E
import numpy as np
import FiniteVolumeGroups as fvg
from su4ops.constants import NS,gammas

def matPrint(mat):
  dim = mat.shape[0]
  s = ""
  for i in range(dim):
    for j in range(dim):
      s += "{}+{}i ".format(mat[i][j].real, mat[i][j].imag)
    s += "\n"
  print(s)


o2h = fvg.cubic.O2h()
len(o2h.elements)


basis = []
fullbasis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
      newOp = E.Elemental(
        1, [utils.quark(s0), utils.quark(s1), utils.quark(s2)])
      fullbasis.append(newOp)
      foundRelated = False
      for b in basis:
        if newOp.quarks in utils.permutations(b.quarks):
          foundRelated = True
          sign = 1
          if not utils.arePermsEqualParity(newOp.quarks, b.quarks):
            sign = 1
          extraBasis[newOp] = {
            'elemental': E.Elemental(1, b.quarks), 'sign': sign}
      if not foundRelated:
        basis.append(newOp)

len(basis)
len(extraBasis)
len(fullbasis)
type(extraBasis)
for b in basis:
  print(b)
for e, val in extraBasis.items():
  print("{} = {} * {}".format(e, val['sign'], val['elemental']))

rep = []
for g in o2h.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, o2h.elements[0]))


metric = np.matrix([[0. for j in range(len(rep[0]))] for i in range(len(rep[0]))], dtype='complex')
for g in rep:
  metric += np.matmul(np.matrix(g).getH(), np.matrix(g))
print(3.0*np.diag(metric)/len(rep))



fvg.representation_checks.is_valid_rep(rep)


import math
matPrint((rep[20]*2*math.sqrt(2)).round(4))

utils.quark(0).spatial_rotate(o2h.elements[20])
utils.quark(1).spatial_rotate(o2h.elements[20])
np.matmul((np.identity(4)+np.matmul(gammas[3],gammas[1])),np.array([1,0,0,0]))







for bi,b in enumerate(basis):
  for gi,elem in enumerate(o2h.elements):
    dim=len(basis)
    lhs=utils.su3_fullVec_to_reduced(b.spatial_rotate(elem),basis,extraBasis)
    rhs=np.zeros(dim,dtype=complex)
    bVec=np.ones(dim)
    for j in range(dim):
      rhs[j]+=bVec[j]*rep[gi][j,bi]
    if not np.allclose(lhs,rhs):
      print(lhs)
      print(rhs)
      raise AssertionError("not a valid rep matrix for basis element {} and group element {}".format(bi,gi))
matPrint(rep[1])





o2h.elements.index(o2h.get_element('C4y'))
o2h.elements.index(o2h.get_element('C4z'))
o2h.elements.index(o2h.get_element('E',parity=-1))

import math
Wc4z=np.diag(np.matrix((rep[22]*math.sqrt(2)).round(4)))

expC4z = np.array([-1+1j,1+1j,-1+1j,1+1j,1-1j,1+1j,1-1j,-1+1j,1+1j,1-1j,-1-1j,1-1j,-1-1j,1+1j,1-1j,-1-1j,-1+1j,1+1j,1-1j,-1-1j])
np.allclose(Wc4z,expC4z)

np.allclose(np.diag(np.matrix(rep[48])),np.array([1,1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,-1,-1,-1]))


matPrint((rep[20]*2*math.sqrt(2)).round(4))




tot = 0
for irrep in o2h.elements[0].irreps:
  ops = utils.operators(irrep, rep, o2h)
  tot += len(ops)*len(o2h.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))

print("{} operators across all irreps".format(tot))

for op in utils.operators('G1g', rep, o2h):
  utils.print_vec(op, basis)


utils.operators('Hg', rep, o2h)

utils.operators('Hu', rep, o2h)
