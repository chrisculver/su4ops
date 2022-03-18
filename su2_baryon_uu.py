import utils
import su4ops.quark as Q
import su4ops.elemental as E
import numpy as np
import FiniteVolumeGroups as fvg
from su4ops.constants import NS

oh = fvg.cubic.Oh()

basis = []
fullbasis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
      newOp = E.Elemental(1, [utils.quark(s0), utils.quark(s1)])
      fullbasis.append(newOp)
      foundRelated = False
      for b in basis:
        if newOp.quarks in utils.permutations(b.quarks):
          foundRelated = True
          sign = 1
          if not utils.arePermsEqualParity(newOp.quarks, b.quarks):
            sign = -1
          extraBasis[newOp] = {
            'elemental': E.Elemental(1, b.quarks), 'sign': sign}
      if not foundRelated:
        basis.append(newOp)

len(basis)

for b in fullbasis:
  print(b)


for b in basis:
  print(b)
for e, val in extraBasis.items():
  print("{} = {} * {}".format(e, val['sign'], val['elemental']))


print(basis[1])
tst = E.Elemental(1, [utils.quark(1), utils.quark(0)])


rot01 = basis[1].spatial_rotate(oh.elements[1])
rot10 = tst.spatial_rotate(oh.elements[1])


utils.print_vec(rot01.round(4), fullbasis)
utils.print_vec(rot10.round(4), fullbasis)


rep = []
for g in oh.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, oh.elements[0]))


fvg.representation_checks.is_valid_rep(rep)


tot = 0
for irrep in oh.elements[0].irreps:
  ops = utils.operators(irrep, rep, oh)
  tot += len(ops)*len(oh.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))

print("{} operators across all irreps".format(tot))

utils.operators('A1u', rep, oh)
print("{}-{}".format(basis[3], basis[5]))

#[[0,0,0,1]]
#[0,0,-1,0]
#[0,0,0,0]
#[0,0,0,0]


#[[0,0,0,1],
#[0,0,-1,0],
#[0,1,0,0],
#[-1,0,0,0]]
