from constants import *
import quark as Q
import numpy as np
import FiniteVolumeGroups as fvg
import math

o = fvg.cubic.O()
oh = fvg.cubic.Oh()

id = np.identity(4)
g2g3=np.matmul(gammas[2],gammas[3])
g3g1=np.matmul(gammas[3],gammas[1])
g1g2=np.matmul(gammas[1],gammas[2])

GammaExpected = {
    "E": id,
    "C4x": (id+g2g3)/math.sqrt(2),
    "C4y": (id+g3g1)/math.sqrt(2),
    "C4z": (id+g1g2)/math.sqrt(2),
    "C4x-1": (id-g2g3)/math.sqrt(2),
    "C4y-1": (id-g3g1)/math.sqrt(2),
    "C4z-1": (id-g1g2)/math.sqrt(2),
    "C2x": g2g3,
    "C2y": g3g1,
    "C2z": g1g2,
    "C3alpha": (id-g2g3-g3g1+g1g2)/2.,
    "C3beta": (id-g2g3+g3g1-g1g2)/2.,
    "C3gamma": (id+g2g3-g3g1-g1g2)/2.,
    "C3delta": (id+g2g3+g3g1+g1g2)/2.,
    "C3alpha-1": (id+g2g3+g3g1-g1g2)/2.,
    "C3beta-1": (id+g2g3-g3g1+g1g2)/2.,
    "C3gamma-1": (id-g2g3+g3g1+g1g2)/2.,
    "C3delta-1": (id-g2g3-g3g1-g1g2)/2.,
    "C2a": (g2g3+g3g1)/math.sqrt(2),
    "C2b": (g2g3-g3g1)/math.sqrt(2),
    "C2c": (g2g3+g1g2)/math.sqrt(2),
    "C2d": -(g2g3-g1g2)/math.sqrt(2),
    "C2e": (g3g1+g1g2)/math.sqrt(2),
    "C2f": (g3g1-g1g2)/math.sqrt(2)
}

GammaBarExpected = {
    "E": id,
    "C4x": (id-g2g3)/math.sqrt(2),
    "C4y": (id-g3g1)/math.sqrt(2),
    "C4z": (id-g1g2)/math.sqrt(2),
    "C4x-1": (id+g2g3)/math.sqrt(2),
    "C4y-1": (id+g3g1)/math.sqrt(2),
    "C4z-1": (id+g1g2)/math.sqrt(2),
    "C2x": -g2g3,
    "C2y": -g3g1,
    "C2z": -g1g2,
    "C3alpha": (id+g2g3+g3g1-g1g2)/2.,
    "C3beta": (id+g2g3-g3g1+g1g2)/2.,
    "C3gamma": (id-g2g3+g3g1+g1g2)/2.,
    "C3delta": (id-g2g3-g3g1-g1g2)/2.,
    "C3alpha-1": (id-g2g3-g3g1+g1g2)/2.,
    "C3beta-1": (id-g2g3+g3g1-g1g2)/2.,
    "C3gamma-1": (id+g2g3-g3g1-g1g2)/2.,
    "C3delta-1": (id+g2g3+g3g1+g1g2)/2.,
    "C2a": -(g2g3+g3g1)/math.sqrt(2),
    "C2b": -(g2g3-g3g1)/math.sqrt(2),
    "C2c": -(g2g3+g1g2)/math.sqrt(2),
    "C2d": (g2g3-g1g2)/math.sqrt(2),
    "C2e": -(g3g1+g1g2)/math.sqrt(2),
    "C2f": -(g3g1-g1g2)/math.sqrt(2)
}

for i,g in enumerate(o.elements):
    if not np.allclose(Q.Gamma(g),GammaExpected[g.identifier["name"]]):
        print("Constructing gamma for {} failed...".format(g.identifier['name']))
    if not np.allclose(Q.Gamma(g,True),GammaBarExpected[g.identifier["name"]]):
        print("Constructing gammaBar for {} failed...".format(g.identifier['name']))

#this is a trivial check
inv = gammas[4]

for i,g in enumerate(oh.elements):
    if g.identifier["parity"]==-1:
        if not np.allclose(Q.Gamma(g),np.matmul(inv,GammaExpected[g.identifier["name"]])):
            print("Constructing gamma for {} failed...".format(g.identifier['name']))
        if not np.allclose(Q.Gamma(g,True),np.matmul(inv,GammaBarExpected[g.identifier["name"]])):
            print("Constructing gammaBar for {} failed...".format(g.identifier['name']))
