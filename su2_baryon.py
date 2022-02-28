import utils
import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
from constants import NS

oh = fvg.cubic.Oh()

basis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
      newOp = Q.Elemental(1, [utils.quark(s0), utils.quark(s1)])
      foundRelated = False
      for b in basis:
        if newOp.quarks in utils.permutations(b.quarks):
          foundRelated = True
          extraBasis[newOp] = Q.Elemental(1, b.quarks)
      if not foundRelated:
        basis.append(newOp)

len(basis)

for b in basis:
  print(b)
for b in extraBasis:
  print(b)

tVec = basis[0].spatial_rotate(oh.elements[1])

tVec
print(basis[1])
print(basis[4])

utils.fullVec_to_reduced(tVec, basis, extraBasis)


rep = []
for g in oh.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, oh.elements[0]))

tot = 0
for irrep in oh.elements[0].irreps:
  ops = utils.operators(irrep, rep, oh)
  tot += len(ops)*len(oh.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))

print("{} operators across all irreps".format(tot))

utils.operators('A1u', rep, oh)
print("{}-{}".format(basis[3], basis[5]))
