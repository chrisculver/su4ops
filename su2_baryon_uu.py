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


#check that the rep is closed
for g1 in rep:
    for g2 in rep:
        prod = np.matmul(g1, g2)
        if not any(np.allclose(prod, g) for g in rep):
            print("not closed")

# check associativity
for a in rep:
    for b in rep:
        for c in rep:
            lhs = np.matmul(np.matmul(a, b), c)
            rhs = np.matmul(a, np.matmul(b, c))
            if not np.allclose(lhs, rhs):
                print("associativity fails")

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
