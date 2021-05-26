from collections import defaultdict
from PhenotypeGraphFun import *
from copy import deepcopy

def stronglycc_iterative(vertices, edges):
    """
    This is a non-recursive version of strongly_connected_components_path.
    See the docstring of that function for more details.

    Examples
    --------
    Example from Gabow's paper [1]_.

    >>> vertices = [1, 2, 3, 4, 5, 6]
    >>> edges = {1: [2, 3], 2: [3, 4], 3: [], 4: [3, 5], 5: [2, 6], 6: [3, 4]}
    >>> for scc in strongly_connected_components_iterative(vertices, edges):
    ...     print(scc)
    ...
    set([3])
    set([2, 4, 5, 6])
    set([1])

    Example from Tarjan's paper [2]_.

    >>> vertices = [1, 2, 3, 4, 5, 6, 7, 8]
    >>> edges = {1: [2], 2: [3, 8], 3: [4, 7], 4: [5],
    ...          5: [3, 6], 6: [], 7: [4, 6], 8: [1, 7]}
    >>> for scc in  strongly_connected_components_iterative(vertices, edges):
    ...     print(scc)
    ...
    set([6])
    set([3, 4, 5, 7])
    set([8, 1, 2])

    """
    identified = set()
    stack = []
    index = {}
    boundaries = []

    for v in vertices:
        if v not in index:
            to_do = [('VISIT', v)]
            while to_do:
                operation_type, v = to_do.pop()
                if operation_type == 'VISIT':
                    index[v] = len(stack)
                    stack.append(v)
                    boundaries.append(index[v])
                    to_do.append(('POSTVISIT', v))
                    # We reverse to keep the search order identical to that of
                    # the recursive code;  the reversal is not necessary for
                    # correctness, and can be omitted.
                    to_do.extend(
                        reversed([('VISITEDGE', w) for w in edges[v]]))
                elif operation_type == 'VISITEDGE':
                    if v not in index:
                        to_do.append(('VISIT', v))
                    elif v not in identified:
                        while index[v] < boundaries[-1]:
                            boundaries.pop()
                else:
                    # operation_type == 'POSTVISIT'
                    if boundaries[-1] == index[v]:
                        boundaries.pop()
                        scc = list(stack[index[v]:])
                        del stack[index[v]:]
                        identified.update(scc)
                        yield scc


def condensation_graph(edges, paramslist):
    vertices = [p for p in edges]
    #compute strongly connected components
    stronglycc = [scc for scc in stronglycc_iterative(vertices, edges) if scc !=[]]

    #collapse phenotype graph key values to first item in strongly connected groups
    condensation = {}
    for s in stronglycc:
        condensation[s[0]] = list(set([i for j in s for i in edges[j] if i != s[0]]))

    for i in condensation:
        condensation[i] = list(set([s[0] for j in range(len(condensation[i])) for s in stronglycc if condensation[i][j] in s and i != s[0]]))

    reduced_paramslist = []
    for layer in paramslist:
        p = [s[0] for s in stronglycc if s[0] in layer]
        reduced_paramslist.append(p)

    return condensation, reduced_paramslist, stronglycc
