import copy

# from https://www.bernardosulzbach.com/heaps-algorithm/
def _swap(elements, i, j):
    elements[i], elements[j] = elements[j], elements[i]

# from https://www.bernardosulzbach.com/heaps-algorithm/
def _generate_permutations(elements, n):
    # As by Robert Sedgewick in Permutation Generation Methods
    c = [0] * n
    yield elements
    i = 0
    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                _swap(elements, 0, i)
            else:
                _swap(elements, c[i], i)
            yield elements
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1


def permutations(elems):
    elements = copy.deepcopy(elems) #leave original list in order
    return _generate_permutations(elements, len(elements))
