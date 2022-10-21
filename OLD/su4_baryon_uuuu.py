import math
import su4ops.utils as utils
import su4ops.quark as Q
import su4ops.elemental as E
import numpy as np
import FiniteVolumeGroups as fvg
from su4ops.constants import NS, gammas


def matPrint(mat):
  dim = mat.shape[0]
  s = ""
  for i in range(dim):
    for j in range(dim):
      s += "{}+{}i ".format(mat[i][j].real, mat[i][j].imag)
    s += "\n"
  print(s)


oh = fvg.cubic.Oh()
len(oh.elements)


basis = []
fullbasis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
      for s3 in range(0, NS):
        newOp = E.Elemental(
          1, [utils.quark(s0), utils.quark(s1), utils.quark(s2), utils.quark(s3)])
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
#for b in basis:
#  print(b)
#for e, val in extraBasis.items():
#  print("{} = {} * {}".format(e, val['sign'], val['elemental']))

rep = []
for g in oh.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, oh.elements[0]))


metric = np.matrix([[0. for j in range(len(rep[0]))] for i in range(len(rep[0]))], dtype='complex')
for g in rep:
  metric += np.matmul(np.matrix(g).getH(), np.matrix(g))

for i in range(len(rep[0])):
  for j in range(len(rep[0])):
    if i != j and not np.isclose(metric[i, j], 0):
      print("Off diagonal element is non-zero")
      break
matPrint(np.array(metric))


fvg.representation_checks.is_valid_rep(rep)


for bi, b in enumerate(basis):
  for gi, elem in enumerate(oh.elements):
    dim = len(basis)
    lhs = utils.su4_fullVec_to_reduced(
      b.spatial_rotate(elem), basis, extraBasis)
    rhs = np.zeros(dim, dtype=complex)
    bVec = np.ones(dim)
    for j in range(dim):
      rhs[j] += bVec[j]*rep[gi][j, bi]
    if not np.allclose(lhs, rhs):
      print(lhs)
      print(rhs)
      raise AssertionError(
        "not a valid rep matrix for basis element {} and group element {}".format(bi, gi))


tot = 0
for irrep in oh.elements[0].irreps:
  ops = utils.operators(irrep, rep, oh)
  tot += len(ops)*len(oh.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))

print("{} operators across all irreps".format(tot))


for op in utils.operators('A1g', rep, oh):
  utils.print_vec(op, basis)
  print()

# same for ubar operator
# +1*uuuu(0,0,3,3) -2*uuuu(0,1,2,3) + 1*uuuu(1,1,2,2)

#O_1(t)= + 1*uuuu(0,0,3,3) -2*uuuu(0,1,2,3) + 1*uuuu(1,1,2,2)
#
#O_2(t)=O_1(t) + 1*uuuu(0303) + uuuu(3300) + ... uuuu(2211) + ....

for op in utils.operators('Eg', rep, oh):
  utils.print_vec(op, basis)
  print()

for op in utils.operators('Eu', rep, oh):
  utils.print_vec(op, basis)
  print()
