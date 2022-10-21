import su4ops.utils as utils
import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
from constants import NS

o2h = fvg.cubic.O2h()

basis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
      newOp = Q.Elemental(
        1, [utils.quark(s0), utils.quark(s1), utils.quark(s2)])
      foundRelated = False
      for b in basis:
        if newOp.quarks in utils.permutations(b.quarks):
          foundRelated = True
          sign = 1
          if not utils.arePermsEqualParity(newOp.quarks, b.quarks):
            sign = -1
          extraBasis[newOp] = {
            'elemental': Q.Elemental(1, b.quarks), 'sign': sign}
      if not foundRelated:
        basis.append(newOp)

len(basis)
len(extraBasis)
type(extraBasis)
for b in basis:
  print(b)
for e, val in extraBasis.items():
  print("{} = {} * {}".format(e, val['sign'], val['elemental']))


rep = []
for g in o2h.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, o2h.elements[0]))


#check that the rep is closed
for i, g1 in enumerate(rep):
    for j, g2 in enumerate(rep):
        prod = np.matmul(g1, g2)
        if not any(np.allclose(prod, g) for g in rep):
            print("not closed for {}, {}".format(i, j))
            break

for op in utils.operators('G1g', rep, o2h):
  utils.print_vec(op, basis)
