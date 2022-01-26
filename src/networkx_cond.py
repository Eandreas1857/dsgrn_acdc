import networkx as nx

def condensation(G, scc=None):
    """Returns the condensation of G.
 
    The condensation of G is the graph with each of the strongly connected
    components contracted into a single node.
 
    Parameters
    ----------
    G : NetworkX DiGraph
       A directed graph.
 
    scc:  list or generator (optional, default=None)
       Strongly connected components. If provided, the elements in
       `scc` must partition the nodes in `G`. If not provided, it will be
       calculated as scc=nx.strongly_connected_components(G).
 
    Returns
    -------
    C : NetworkX DiGraph
       The condensation graph C of G. The node labels are integers
       corresponding to the index of the component in the list of
       strongly connected components of G. C has a graph attribute named
       'mapping' with a dictionary mapping the original nodes to the
       nodes in C to which they belong. Each node in C also has a node
       attribute 'members' with the set of original nodes in G that
       form the SCC that the node in C represents.
    
    members : dictionary of scc of G associated to each node in C
 
    Raises
    ------
    NetworkXNotImplemented:
        If G is not directed
 
    Notes
    -----
    After contracting all strongly connected components to a single node,
    the resulting graph is a directed acyclic graph.
 
    """
    if scc is None:
        scc = nx.strongly_connected_components(G)
    mapping = {}
    members = {}
    C = nx.DiGraph()
    # Add mapping dict as graph attribute
    C.graph["mapping"] = mapping
    if len(G) == 0:
        return C
    for i, component in enumerate(scc):
        members[i] = component
        mapping.update((n, i) for n in component)
    number_of_components = i + 1
    C.add_nodes_from(range(number_of_components))
    C.add_edges_from(
        (mapping[u], mapping[v]) for u, v in G.edges() if mapping[u] != mapping[v]
    )
    # Add a list of members (ie original nodes) to each node (ie scc) in C.
    nx.set_node_attributes(C, members, "members")
    return C, members

def strongly_connected_components_by_MGI(G, database):
    """Returns strongly connected components of G refined by morse graph index 
    """
    c = database.conn.cursor()

    strongcc = []
    for scc in nx.strongly_connected_components(G):
        sub_scc = {}
        for node in scc:
            p = node[2]
            MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
            MGI = MGI_result.fetchone()[0]
            if MGI not in sub_scc:
                sub_scc[MGI] = [node]
            else:
                sub_scc[MGI] += [node]
        for i in sub_scc:
            strongcc.append(sub_scc[i])
    return strongcc