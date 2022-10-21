import su4ops.gammas as G

NC = 0
NS = G.NS

gammas=G.get_basis('dirac-pauli')

flavorLabels = {0: 'u', 1: 'd'}
NF = len(flavorLabels)
colorLabels = {0: 'a', 1: 'b', 2: 'c', 3: 'd'}
