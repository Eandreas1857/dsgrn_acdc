import itertools
import networkx as nx
import DSGRN
from DSGRN import *

# Function to generate all binary strings
def generateAllBinaryStrings(binary_list, n, arr, i):
    if i == n:
        a = []
        for i in range(0, n):
            a.append(arr[i])
        binary_list.append(a)
        return
     
    # First assign "0" at ith position
    # and try for all other permutations
    # for remaining positions
    arr[i] = 0
    generateAllBinaryStrings(binary_list, n, arr, i + 1)
 
    # And then assign "1" at ith position
    # and try for all other permutations
    # for remaining positions
    arr[i] = 1
    generateAllBinaryStrings(binary_list, n, arr, i + 1)

def get_network_string(net):
    """
    net: tuple, with first item a set of edges, and second item with a list of 0,1's depicting if edge in first item
    is repressing (0) or activating (1). 
    
    Example:  net = ((('Hb', 'Gt'),
                ('Hb', 'Kr'),
                ('Hb', 'Kni'),
                ('Gt', 'Hb'),
                ('Gt', 'Kr'),
                ('Gt', 'Kni'),
                ('Kr', 'Hb'),
                ('Kr', 'Gt')),
                [0, 0, 1, 1, 0, 0, 1, 0])
    then the edge from 'Hb' to 'Gt' is repressing while the 'Kr' to 'Hb' edge is activating.

    returns: string for use with DSGRN network input.
    """

    net_dict = {'Hb': [], 'Gt': [], 'Kr': [], 'Kni':[]}

    for i in net[0]:
        index = net[0].index(i)
        net_dict[i[1]] += [(net[1][index], i[0])]
        
    new = {}
    for node in net_dict:
        d = defaultdict(list)
        act_str = ''
        rep_str = ''
        for k, *v in net_dict[node]:
            d[k].append(v)
        
        for edge_type in list(d.items()):
            if edge_type[0] == 0:
                rep_str = ''
                for i in edge_type[1]:
                    rep_str += '(~' + i[0] + ')'

            if edge_type[0] == 1:
                act_str = '('
                for i in edge_type[1]:
                    if edge_type[1].index(i) == 0:
                        act_str +=  i[0]
                    else:
                        act_str +=  '+' + i[0]
                act_str += ')'
        new[node] = act_str + rep_str

    return '"""Hb : ' + new['Hb'] + '\n' + 'Gt : ' + new['Gt'] + '\n' + 'Kr : ' + new['Kr'] + '\n' + 'Kni : ' + new['Kni'] + '"""'


def get_all_networks(node_list, n):
    """
    node_list: set of nodes wanting in network.
    n: number of edges wanting in network. 
    returns: list of tuples. First element in tuple is set of edges, second is a binary list 
             where 0 depictis edge is rep and 1 depicts if edge is act.
    """
    binary_list = []
    arr = [None]*n
    generateAllBinaryStrings(binary_list, n, arr, 0)

    edge_list = [(a, b) for a in node_list for b in node_list if a != b]
    all_edge_comb = list(itertools.combinations(edge_list, 8))
    all_network_comb =  [(a, b) for a in all_edge_comb for b in binary_list]

    return all_network_comb

def computable_networks(all_network_comb):  

    computable = []
    for net in all_network_comb:
        try:
            string = get_network_string(net)
            network = DSGRN.Network(string)
            pg = ParameterGraph(network)
            p = pg.size()
            computable.append((all_network_comb.index(net), p))
    
        except:
            continue
    return computable

def return_computable_net_w_limited_PG_size(computable, size_limit = 3240000):
    d = defaultdict(list)

    for k, v in computable:
        d[v].append(k)

    allowed = []
    for i in sorted(list(d.items())):
        if i[0] <= size_limit:
            for j in i[1]:
                allowed.append(j) 
    return allowed

def convert_edges_to_networkx(edges):
    H = nx.DiGraph() 
    for edge in edges:
        H.add_edge(edge[0], edge[1])
    return H

def convert_dict_to_networkx(dict):
    H = nx.DiGraph() 
    for s in dict:
        for t in dict[s]:
            H.add_edge(s,t)
    return H

def save_networkx_as_png(G, filename):
    g = nx.drawing.nx_pydot.to_pydot(G)
    g.write_png(filename)

def save_networkx_as_png_w_edge_weight(G, data, filename):
    pdot = nx.drawing.nx_pydot.to_pydot(G)
    for i, edge in enumerate(pdot.get_edges()):
        edge.set_label(edge.get_attributes()[data])
    pdot.write_png(filename)