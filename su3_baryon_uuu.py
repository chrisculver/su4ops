import utils
import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
from constants import NS

o2h = fvg.cubic.O2h()
len(o2h.elements)

basis = []
fullbasis = []
extraBasis = {}
for s0 in range(0, NS):
  for s1 in range(0, NS):
    for s2 in range(0, NS):
      newOp = Q.Elemental(
        1, [utils.quark(s0), utils.quark(s1), utils.quark(s2)])
      fullbasis.append(newOp)
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
len(fullbasis)
type(extraBasis)
for b in basis:
  print(b)
for e, val in extraBasis.items():
  print("{} = {} * {}".format(e, val['sign'], val['elemental']))

rep = []
for g in o2h.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, o2h.elements[0]))


metric = np.matrix([[0. for j in range(len(rep[0]))]
                   for i in range(len(rep[0]))], dtype='complex')
for g in rep:
  metric += np.matmul(np.matrix(g).getH(), np.matrix(g))
print(3.0*np.diag(metric)/len(rep))

o2h.elements[48].rotation
utils.print_vec(np.diag(rep[48]), basis)
np.diag(rep[49])

print(basis[0])
utils.print_vec(basis[0].spatial_rotate(o2h.elements[48]), fullbasis)
utils.print_vec(utils.fullVec_to_reduced(
  basis[0].spatial_rotate(o2h.elements[48]), basis, extraBasis), basis)


#check that the rep is closed
for i, g1 in enumerate(rep):
    for j, g2 in enumerate(rep):
        prod = np.matmul(g1, g2)
        if not any(np.allclose(prod, g) for g in rep):
            print("not closed for {}, {}".format(i, j))
            break

# check associativity
for a in rep:
    for b in rep:
        for c in rep:
            lhs = np.matmul(np.matmul(a, b), c)
            rhs = np.matmul(a, np.matmul(b, c))
            if not np.allclose(lhs, rhs):
                print("associativity fails")
                break

# identity
len(rep[0])
np.allclose(rep[0], np.identity(len(rep[0])))

# inverse
for g1 in rep:
    has_inverse = False
    for g2 in rep:
        if np.allclose(np.matmul(g1, g2), np.identity(len(rep[0]))):
            has_inverse = True
    if not has_inverse:
        print("g1 doesn't have an inverse")

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
