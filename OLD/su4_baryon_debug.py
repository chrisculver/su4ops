

import su4ops.utils as utils
import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
from constants import NS

oh = fvg.cubic.Oh()


basis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
      for s3 in range(0, NS):
        newOp = Q.Elemental(1, [utils.quark(s0), utils.quark(
          s1), utils.quark(s2), utils.quark(s3)])
        foundRelated = False
        for b in basis:
          if newOp.quarks in utils.permutations(b.quarks):
            foundRelated = True
            extraBasis[newOp] = Q.Elemental(1, b.quarks)
        if not foundRelated:
          basis.append(newOp)

len(basis)
len(extraBasis)

rep = []
for g in oh.elements:
    rep.append(utils.makeRepMat(basis, extraBasis, g, oh.elements[0]))


tot = 0
for irrep in oh.elements[0].irreps:
  ops = utils.operators(irrep, rep, oh)
  tot += len(ops)*len(oh.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))
print("{} operators across all irreps".format(tot))

utils.operators('A1g', rep, oh)

print("{}-2*{}+{}".format(basis[9], basis[14], basis[23]))

#TODO: is the full Gamma_{abcd} symmetric or anti-symmetric.


utils.operators('Eg', rep, oh)
print("{}+{}+{}+{}".format(basis[0], basis[20], basis[30], basis[34]))


utils.operators('T2g', rep, oh)
