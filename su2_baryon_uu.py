from su4ops.constants import gammas
import utils
import su4ops.quark as Q
import su4ops.elemental as E
import numpy as np
import FiniteVolumeGroups as fvg
from su4ops.constants import NS


def matPrint(mat):
  dim = mat.shape[0]
  s = ""
  for i in range(dim):
    for j in range(dim):
      s += "{}+{}i ".format(mat[i][j].real, mat[i][j].imag)
    s += "\n"
  print(s)


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

print(fullbasis[0])
print(fullbasis[1])
print(fullbasis[4])
print(fullbasis[5])

rot00 = fullbasis[0].spatial_rotate(oh.elements[1])
rot01 = fullbasis[1].spatial_rotate(oh.elements[1])
rot10 = fullbasis[4].spatial_rotate(oh.elements[1])
rot11 = fullbasis[5].spatial_rotate(oh.elements[1])


print(rot00.round(4))
print(utils.su2_fullVec_to_reduced(rot00, basis, extraBasis))


oh.elements[1].identifier['name']
id = np.identity(4)
g2g3 = np.matmul(gammas[2], gammas[3])
g3g1 = np.matmul(gammas[3], gammas[1])
g1g2 = np.matmul(gammas[1], gammas[2])


utils.print_vec(rot00.round(4), fullbasis)
utils.print_vec(rot01.round(4), fullbasis)
utils.print_vec(rot10.round(4), fullbasis)
utils.print_vec(rot11.round(4), fullbasis)

utils.su2_fullVec_to_reduced(rot00,basis,extraBasis)
utils.su2_fullVec_to_reduced(rot01,basis,extraBasis)
utils.su2_fullVec_to_reduced(rot10,basis,extraBasis)
utils.su2_fullVec_to_reduced(rot11,basis,extraBasis)


(id+g2g3+g3g1+g1g2)/2.


utils.quark(0).spatial_rotate(oh.elements[1])
utils.quark(1).spatial_rotate(oh.elements[1])


rep = []
for g in oh.elements:
  rep.append(utils.makeRepMat(basis, extraBasis, g, oh.elements[0]))


fvg.representation_checks.is_valid_rep(rep)


r1=rep[1]
r1adj=np.array(np.matrix(rep[1]).getH())


matPrint(np.matmul(r1,r1adj))

matPrint(np.matmul(rep[1], rep[1]))

matPrint(rep[24])
matPrint(rep[34])


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
