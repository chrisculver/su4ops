import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
from constants import *

oh = fvg.cubic.Oh()

quark = lambda spin : Q.Quark({
    'bar': False,
    'flavor': 0,
    'color': 0,
    'spin': spin,
})

basis = []

for s0 in range(0,NS):
    for s1 in range(0,NS):
        for s2 in range(0,NS):
            for s3 in range(0,NS):
                basis.append(Q.Elemental([quark(s0),quark(s1),quark(s2),quark(s3)]))



rep = []
for g in oh.elements:
    rep.append(makeRepMat(basis,g))

np.allclose(rep[0],np.identity(4**4))
