import su4ops.gammas as G

NC = 2
NS = G.NS
gamma_basis = 'degrand-rossi'  # 'dirac-pauli' #weyl-chiral #degrand-rossi

gammas = G.get_basis(gamma_basis)

flavorLabels = {0: 'u', 1: 'd'}
NF = len(flavorLabels)
colorLabels = {0: 'a', 1: 'b', 2: 'c', 3: 'd'}
