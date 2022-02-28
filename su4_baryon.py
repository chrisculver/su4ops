import utils
import quark as Q
import numpy as np
import sympy as sp
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
len(basis)+len(extraBasis)

rep = []
for g in oh.elements:
  rep.append(makeRepMat(basis, extraBasis, g, oh.elements[0]))


np.allclose(rep[0], np.eye(len(basis), len(basis)))

quark(0).spatial_rotate(oh.elements[1])
print(basis[0])
basis[0].spatial_rotate(oh.elements[1]).round(8)
fullVec_to_reduced(basis[0].spatial_rotate(oh.elements[1]), basis, extraBasis)


for r in range(len(basis)):
  line = ""
  for c in range(len(basis)):
    line += "{} ".format(rep[1][r, c].round(4))
  line += "\n"
  print(line)


def projectorMat(irrep, rep, group, row=0):
    dim = len(rep[0])
    res = np.zeros((dim, dim), dtype=complex)
    for i, elem in enumerate(group.elements):
        #print(type(float(elem.irreps[irrep][row,row])))
        res += complex(elem.irreps[irrep][row, row])*np.transpose(rep[i])
    return (res*len(irrep)/len(group.elements)).round(8)


def operators(irrep, rep, group, row=0):
  proj = projectorMat(irrep, rep, group, row)
  rrefMat = sp.Matrix(proj).rref()[0].tolist()
  for row in rrefMat[:]:
    if (row == [0 for i in range(len(row))]):
      rrefMat.remove(row)
  rrefMat = np.array(rrefMat).astype(complex)
  return rrefMat


def op_basis_map(op):
  basis_map = {}
  for i, val in enumerate(op):
    if not np.isclose(val, 0):
      basis_map[i] = val
  return basis_map


tot = 0
for irrep in oh.elements[0].irreps:
  ops = operators(irrep, rep, oh)
  tot += len(ops)*len(oh.elements[0].irreps[irrep])
  print("{} ops in {}".format(len(ops), irrep))
print("{} ops across all irreps".format(tot))

op_basis_map(operators('A1g', rep, oh)[0])


### OLD STUFF FOR FULL 256 BASIS
#print("{}-{}-{}+{}\n+{}-{}-{}+{}".format(basis[5],basis[20],basis[65],basis[80],basis[175],basis[190],basis[235],basis[250]))


print("{} operators across all irreps".format(tot))

for irrep in ['Eg', 'Eu', 'T2g', 'T2u']:
  ops = operators(irrep, rep, oh, row=0)
  print("Ops in {}:".format(irrep))
  for op in ops:
    print("  {}".format(op_basis_map(op)))
  print("")

for i in [0, 10, 34, 40, 85, 95, 119, 125, 130, 136, 160, 170, 215, 221, 245, 255]:
  print(basis[i])
op_map = op_basis_map(operators('T2g', rep, oh, row=0)[0])

file = open("chris_baryon.txt", "w")

for s0 in range(NS):
  for s1 in range(NS):
    for s2 in range(NS):
      for s3 in range(NS):
        idx = s0*NS*NS*NS + s1*NS*NS + s2*NS + s3
        if idx in op_map.keys():
          file.write("{} {} {} {} {}\n".format(s0, s1, s2, s3, op_map[idx]))
        else:
          file.write("{} {} {} {} {}\n".format(s0, s1, s2, s3, 0))

file.close()
