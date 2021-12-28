import itertools
import networkx as nx

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
    '''
    Computes condensation of edges into its strongly connected components.
    :param edges: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    '''
    vertices = [p for p in edges]
    # compute strongly connected components
    stronglycc = [scc for scc in stronglycc_iterative(vertices, edges) if scc != []]

    # collapse phenotype graph dict key's to first item in strongly connected groups
    condensation = {}
    for s in stronglycc:
        condensation[s[0]] = list(set([i for j in s for i in edges[j] if i != s[0]]))

    # collapse phenotype graph dict list_element to first item in strongly connected groups
    for i in condensation:
        condensation[i] = list(set(
            [s[0] for j in range(len(condensation[i])) for s in stronglycc if condensation[i][j] in s and i != s[0]]))

    # also give a reduction of the phenotype pattern
    reduced_paramslist = []
    for layer in paramslist:
        p = [s[0] for s in stronglycc if s[0] in layer]
        reduced_paramslist.append(p)

    return condensation, reduced_paramslist, stronglycc


def condensation_graph_optimized(edges):
    '''
    Computes condensation of edges into its strongly connected components.
    :param edges: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    '''
    key = lambda x : x[0]
    vertices_by_layer = sorted([sorted(list(g[1])) for g in itertools.groupby(edges,key)],key=key)

    # compute strongly connected components
    # stronglycc has layer integers for keys
    stronglycc = {}
    for vertices in vertices_by_layer:
        layer_edges = { node : [e for e in edges[node] if e[0] == node[0]] for node in vertices}
        stronglycc[vertices[0][0]] = [sorted(scc) for scc in stronglycc_iterative(vertices, layer_edges) if scc != []]

    # make condensation graph
    condensation = {}
    for layer,scc_k in stronglycc.items():
        # scc_k = list of scc's for layer k
        if layer < len(vertices_by_layer)-1:
            # scc_kp1 = list of scc's for layer k+1
            scc_kp1 = stronglycc[layer+1]
        else:
            scc_kp1 = []
        for scc in scc_k:
            # scc = list of (layer, param)
            if scc[0] not in condensation:
                condensation[scc[0]] = []
            for tcc in scc_k + scc_kp1:
                # tcc = list of (layer, param)
                if tcc != scc:
                    for s in scc:
                        if set(edges[s]).intersection(tcc):
                            condensation[scc[0]].append(tcc[0])
                            break

    # also give a reduction of the phenotype pattern
    reduced_paramslist = sorted([sorted(list(g[1])) for g in itertools.groupby(condensation,key)],key=key)

    return condensation, reduced_paramslist, stronglycc

def condensation_graph_gradient(database, edges, FP_list):
    '''
    Computes condensation of edges into its strongly connected components.
    :param edges: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    '''
    c = database.conn.cursor()

    key = lambda x : x[0]
    vertices_by_layer = sorted([sorted(list(g[1])) for g in itertools.groupby(edges,key)],key=key)

    # compute strongly connected components
    # stronglycc has layer integers for keys
    strongcc = {}
    for vertices in vertices_by_layer:
        layer_edges = { node : [e for e in edges[node] if e[0] == node[0]] for node in vertices}
        strongcc[vertices[0][0]] = [sorted(scc) for scc in stronglycc_iterative(vertices, layer_edges) if scc != []]
    
    stronglycc = {}
    for key in strongcc:
        stronglycc[key] = []

    for i in strongcc:
        for scc in strongcc[i]:
            MGI = []
            for node in scc:
                    p = node[2]
                    MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
                    MGI.append((node, MGI_result.fetchone()[0]))

                    data = {}
                    for mgi in MGI:
                        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(mgi[-1]))]
                        inter = tuple(set(x for x in FP_result if x in FP_list))
                        if inter not in data:
                            data[inter] = [mgi[0]]
                        else:
                            data[inter] += [mgi[0]]
            for fp in data:
                stronglycc[i].append([node for node in data[fp]])   

    # make condensation graph
    condensation = {}
    for layer,scc_k in stronglycc.items():
        # scc_k = list of scc's for layer k
        if layer < len(vertices_by_layer)-1:
            # scc_kp1 = list of scc's for layer k+1
            scc_kp1 = stronglycc[layer+1]
        else:
            scc_kp1 = []
        for scc in scc_k:
            # scc = list of (layer, param)
            if scc[0] not in condensation:
                condensation[scc[0]] = []
            for tcc in scc_k + scc_kp1:
                # tcc = list of (layer, param)
                if tcc != scc:
                    for s in scc:
                        if set(edges[s]).intersection(tcc):
                            condensation[scc[0]].append(tcc[0])
                            break

    return condensation, stronglycc

def all_cond_paths(cond, stop_set_num, path_length):
    G = nx.DiGraph() #building networkx graph grom grad_graph
    for node in cond:
        for edge in cond[node]:
            G.add_edge(node, edge)

    start_set = []
    for i in cond:
        a, b, c = i
        if a==0 and b==0:
            start_set.append(i)

    stop_set = []
    for i in cond:
        a, b, c = i
        if a==stop_set_num and b==stop_set_num:
            stop_set.append(i)

    Paths = []
    for s in start_set:
        for t in stop_set:
            for path in nx.all_simple_paths(G, s, t, cutoff = path_length):
                Paths.append(path)
    return Paths