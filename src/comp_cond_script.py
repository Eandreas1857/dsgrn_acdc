import networkx as nx
import DSGRN
from DSGRN import *
import itertools
import sys
sys.setrecursionlimit(10**8)
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/notebooks')
from PhenotypeGraphFun import *
from PhenotypeGraphviz import *
from GradientFun import *

from save_files import *

def load_json(file_to_load):
    with open(file_to_load) as f:
        data = json.load(f)
        
    new_data = {}
    for key in data:
        new_data[ast.literal_eval(key)] = [tuple(i) for i in data[key] ]
    return new_data

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

def condensation_graph_gradient(database, edges, FP_list):
    '''
    Computes condensation of edges into its strongly connected components.
    :param edges: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    '''
    c = database.conn.cursor()

    vertices_by_layer = []
    for kni in range(13):
        vert = []
        for node in grad_graph:
            if node[1] == kni:
                vert.append(node)
        if vert != []:
            vertices_by_layer.append(vert)

    # compute strongly connected components
    # stronglycc has layer integers for keys
    strongcc = {}
    for vertices in vertices_by_layer:
        layer_edges = { node : [e for e in edges[node] if e[1] == node[1]] for node in vertices}
        strongcc[vertices[0][1]] = [sorted(scc) for scc in stronglycc_iterative(vertices, layer_edges) if scc != []]
    
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

def Fullconn():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_Fullconn.db")

    Kni_list = {0: ['00'], 1:['40'], 2:['44','50','C0' ], 3: ['54','C4','D0'], 4: ['55','D4','CC','F0'], 5: ['D5','DC','F4'], 6: ['DD','F5','FC'], 7: ['FD'], 8:['FF']}
    Hb_list = {0:['FF'], 1: ['FD'], 2: ['DD','F5','FC'], 3: ['D5','DC','F4'], 4: ['55','D4','CC','F0'], 5: ['54','C4','D0'], 6:['44','50','C0' ], 7:['40'], 8: ['00']}
    
    FP_Poset = {
        'FP { 1, 2, 0, 0 }': [ 'FP { 1, 1, 0, 0 }',  'FP { 2, 2, 0, 0 }'],
        'FP { 1, 1, 0, 0 }': ['FP { 2, 1, 0, 0 }'],
        'FP { 2, 2, 0, 0 }': ['FP { 2, 1, 0, 0 }'],

        'FP { 2, 1, 0, 0 }': [ 'FP { 2, 0, 0, 0 }',  'FP { 2, 1, 1, 0 }'],
        'FP { 2, 0, 0, 0 }': ['FP { 2, 0, 1, 0 }'],
        'FP { 2, 1, 1, 0 }': ['FP { 2, 0, 1, 0 }'],

        'FP { 2, 0, 1, 0 }': [ 'FP { 1, 0, 1, 0 }',  'FP { 2, 0, 2, 0 }'],
        'FP { 1, 0, 1, 0 }': ['FP { 1, 0, 2, 0 }'],
        'FP { 2, 0, 2, 0 }': ['FP { 1, 0, 2, 0 }'],

        'FP { 1, 0, 2, 0 }': [ 'FP { 0, 0, 2, 0 }',  'FP { 1, 0, 2, 1 }'],
        'FP { 0, 0, 2, 0 }': ['FP { 0, 0, 2, 1 }'],
        'FP { 1, 0, 2, 1 }': ['FP { 0, 0, 2, 1 }'],

        'FP { 0, 0, 2, 1 }': [ 'FP { 0, 0, 1, 1 }',  'FP { 0, 0, 2, 2 }'],
        'FP { 0, 0, 1, 1 }': ['FP { 0, 0, 1, 2 }'],
        'FP { 0, 0, 2, 2 }': ['FP { 0, 0, 1, 2 }'],

        'FP { 0, 0, 1, 2 }': [ 'FP { 0, 1, 1, 2 }',  'FP { 0, 0, 0, 2 }'],
        'FP { 0, 1, 1, 2 }': ['FP { 0, 1, 0, 2 }'],
        'FP { 0, 0, 0, 2 }': ['FP { 0, 1, 0, 2 }'],

        'FP { 0, 1, 0, 2 }': [ 'FP { 0, 2, 0, 2 }',  'FP { 0, 1, 0, 1 }'],
        'FP { 0, 2, 0, 2 }': ['FP { 0, 2, 0, 1 }'],
        'FP { 0, 1, 0, 1 }': ['FP { 0, 2, 0, 1 }'],
        'FP { 0, 2, 0, 1 }': []
    }
    FP_keep = [node for node in FP_Poset.keys()]

    return database, Hb_list, Kni_list, FP_Poset, FP_keep
   

def StrongEdges():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_StrongEdges.db") 

    Hb_list = {0: ['F'], 1: ['E'], 2: ['A', 'C'], 3: ['8'], 4: ['0']}
    Kni_list = {0:['000'], 1: ['200'], 2: ['208', '240', '600'], 3: ['248','608','640','E00'], 4: ['249','648','618','6C0','E08','E40'], 5:['649','658','6C8','E18','E48','EC0'], 6: ['659','6C9','6D8', 'E49','E58','EC8','E38','FC0'], 7: ['6D9','E59','EC9','ED8','E78','FC8'], 8: ['6DB','ED9','E79','EF8','FC9','FD8'], 9: ['EDB','EF9','FD9','FF8'], 10: ['EFB','FDB','FF9'], 11:['FFB'], 12: ['FFF']}
    
    FP_Poset = {
     'FP { 0, 3, 0, 0 }': [ 'FP { 1, 3, 0, 0 }',  'FP { 0, 2, 0, 0 }'],
     'FP { 1, 3, 0, 0 }': ['FP { 1, 2, 0, 0 }'],
     'FP { 0, 2, 0, 0 }': ['FP { 1, 2, 0, 0 }', 'FP { 0, 1, 0, 0 }'],
     'FP { 0, 1, 0, 0 }': ['FP { 1, 1, 0, 0 }'],


     'FP { 1, 2, 0, 0 }': ['FP { 1, 1, 0, 0 }'],
     'FP { 1, 1, 0, 0 }': ['FP { 1, 0, 0, 0 }'],

     'FP { 1, 0, 0, 0 }': ['FP { 0, 0, 0, 0 }', 'FP { 1, 0, 1, 0 }'],
     'FP { 0, 0, 0, 0 }': ['FP { 0, 0, 1, 0 }'],
     'FP { 1, 0, 1, 0 }': ['FP { 0, 0, 1, 0 }'],

     'FP { 0, 0, 1, 0 }': ['FP { 0, 0, 1, 1 }'],

     'FP { 0, 0, 1, 1 }': ['FP { 0, 0, 1, 2 }', 'FP { 0, 0, 0, 1 }'],
     'FP { 0, 0, 0, 1 }': ['FP { 0, 0, 0, 2 }'],
     'FP { 0, 0, 1, 2 }': ['FP { 0, 0, 1, 3 }', 'FP { 0, 0, 0, 2 }'],
     'FP { 0, 0, 0, 2 }': ['FP { 0, 0, 0, 3 }'],
     'FP { 0, 0, 1, 3 }': ['FP { 0, 0, 0, 3 }'],

     'FP { 0, 0, 0, 3 }': ['FP { 0, 1, 0, 3 }'],

     'FP { 0, 1, 0, 3 }': ['FP { 0, 2, 0, 3 }', 'FP { 0, 1, 0, 2 }'],
     'FP { 0, 1, 0, 2 }': ['FP { 0, 2, 0, 2 }', 'FP { 0, 1, 0, 1 }'],
     'FP { 0, 2, 0, 3 }': ['FP { 0, 2, 0, 2 }', 'FP { 0, 3, 0, 3 }'],
     'FP { 0, 2, 0, 2 }': ['FP { 0, 2, 0, 1 }', 'FP { 0, 3, 0, 2 }'],
     'FP { 0, 3, 0, 3 }': ['FP { 0, 3, 0, 2 }'],
     'FP { 0, 1, 0, 1 }': ['FP { 0, 2, 0, 1 }'],
     'FP { 0, 2, 0, 1 }': ['FP { 0, 3, 0, 1 }'],

     'FP { 0, 3, 0, 2 }': ['FP { 0, 3, 0, 1 }'], 
     'FP { 0, 3, 0, 1 }': []
    }
    FP_keep = [node for node in FP_Poset.keys()]

    def load_json(file_to_load):
        with open(file_to_load) as f:
            data = json.load(f)
            
        new_data = {}
        for key in data:
            new_data[ast.literal_eval(key)] = [tuple(i) for i in data[key] ]
        return new_data

    grad_graph = load_json('StrongEdges_grad_graph_new') #hb,kni <9 factor graph steps

    
    return database, Hb_list, Kni_list, FP_Poset, FP_keep, grad_graph

def run(net):
    if net == 'StrongEdges':
        return StrongEdges()
    if net == 'Fullconn':
        return Fullconn()
        
if __name__=="__main__":
    database, Hb_list, Kni_list, FP_Poset, FP_keep, grad_graph = run('StrongEdges')
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)


    G = nx.DiGraph() #building networkx graph grom grad_graph
    for node in grad_graph:
        G.add_node(node)
        for edge in grad_graph[node]:
            G.add_edge(node, edge)

    del_list = [] #list we want to delete from grad_graph

    for node in grad_graph: #delete node if fp not in fp_keep
        p = node[2]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if not set(FP_result).intersection(set(FP_keep)):
            del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]

    del_list = [] 
    for node in G: #want close to hb = 4 kni in (kni, hb) "graph"
        if node[0] == 0:
            if node[1]>3:
                del_list.append(node)
        elif node[0] == 1:
            if node[1]>6:
                del_list.append(node)
        elif node[0] == 2:
            if node[1] not in [3, 4, 5, 6, 7, 8, 9]:
                del_list.append(node)
        elif node[0] == 3:
            if node[1]<6:
                del_list.append(node)
        elif node[0] == 4:
            if node[1]<3:
                del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]

    nodelist = grad_graph.keys()
    new_grad_graph = {}
    for n in nodelist:
        new_grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    cond, scc = condensation_graph_gradient(database, new_grad_graph, FP_keep)
    
    save_json(cond, 'StrongEdges_grad_graph_cond_new')


#softwrap is alt Z
## how to run in terminal: 
# python comp_cond_script.py </dev/null >comp_cond_script.log 2>&1 &

# ps -a (shows al process that are running on machine)

# cat comp_cond_script.log (prints log)   head does first 10 lines or tail does last 10 lines


## how to use megplux

# ssh elizabethandreas@msu.montana.edu

# pip install DSGRN