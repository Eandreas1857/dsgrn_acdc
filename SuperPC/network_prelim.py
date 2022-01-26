import DSGRN
from DSGRN import *
import collections
import json, os
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import math

from get_FG import *
from get_FP_Poset import *
from networkx_cond import *

def save_json(edges, filename):
    mydict = {}
    for key in edges.keys():
        if type(key) is not str:
            mydict[str(key)] = edges[key]
            
    with open(filename, 'w') as fp:
        json.dump(mydict, fp)

def plot_FG_layer_comb_in_G(G, scc, Hb_max, network_filename):
    plt.figure(figsize=(5,5))
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    x = []
    y = []

    for s in G:
        x.append(Hb_max-scc[s][0][0])
        y.append(scc[s][0][1])
    
    title = 'Kept FG layer combos ' + network_filename
    plt.scatter(y, x)
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    plt.title(title)
    plt.savefig(network_filename + '_diagP.png')   # save the figure to file
    plt.close() 

def get_network_string(edges, bool):
    """
    edges: a tuple of edges
    bool: a list of 0,1's depicting if edge in first item is repressing (0) or activating (1). 
    
    Example:  edges = (('Hb', 'Gt'),
                ('Hb', 'Kr'),
                ('Hb', 'Kni'),
                ('Gt', 'Hb'),
                ('Gt', 'Kr'),
                ('Gt', 'Kni'),
                ('Kr', 'Hb'),
                ('Kr', 'Gt'))

              bool = [0, 0, 1, 1, 0, 0, 1, 0]
    then the edge from 'Hb' to 'Gt' is repressing while the 'Kr' to 'Hb' edge is activating.

    returns: string for use with DSGRN network input.
    """

    net_dict = {'Hb': [], 'Gt': [], 'Kr': [], 'Kni':[]}

    for i in edges:
        index = edges.index(i)
        net_dict[i[1]] += [(bool[index], i[0])]
        
    new = {}
    for node in net_dict:
        d = collections.defaultdict(list)
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

def get_grad_graph_strict_bagged(database, network_string):
    pg = DSGRN.ParameterGraph(database.network)
    c = database.conn.cursor()

    out_edges = get_number_out_edges_from_string(network_string)
    FP_Poset = get_FP_Poset(out_edges)[0]
    FP_keep = [node for node in FP_Poset.keys()]

    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb = {}
    for i in Hb_list:
        for j in Hb_list[i]:
            Hb[j] = i
    Kni = {}
    for i in Kni_list:
        for j in Kni_list[i]:
            Kni[j] = i

    G = nx.DiGraph()
    for s in range(pg.size()):
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(s))
        MGI = MGI_result.fetchone()[0]
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if len(FP_result) == 1 and FP_result[0][0:2] == 'FP':
            if set(FP_result).intersection(set(FP_keep)):
                sHb = Hb[((((pg.parameter(s)).logic())[0]).stringify())[6:-2]]
                sKni = Kni[((((pg.parameter(s)).logic())[3]).stringify())[6:-2]]
                for t in list(pg.adjacencies(s, 'codim1')):
                    MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(t))
                    MGI = MGI_result.fetchone()[0]
                    FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
                    if len(FP_result) == 1 and FP_result[0][0:2] == 'FP':
                        if set(FP_result).intersection(set(FP_keep)):
                            tHb = Hb[((((pg.parameter(t)).logic())[0]).stringify())[6:-2]]
                            tKni = Kni[((((pg.parameter(t)).logic())[3]).stringify())[6:-2]] 
                            if sHb+1 == tHb and sKni == tKni:
                                G.add_edge((sHb, sKni, s), (tHb, tKni, t))  
                            elif sHb == tHb and sKni+1 == tKni:
                                G.add_edge((sHb, sKni, s), (tHb, tKni, t))
                            elif sHb == tHb and sKni == tKni:
                                G.add_edge((sHb, sKni, s), (tHb, tKni, t))  
    return G

def get_product_graph(database, cG, scc, FP_Poset):
    '''
    cG: condensation graph of gradient graph with only monostable fixed points and nodes in FP_Poset, expects networkx object.
    scc: dictonary of stongly connected components, where keys are node labels in given graph 
         and values are nodes in original graph the condensation is derived from. 

    returns: Product Graph, i.e., reduces cG to having only edges that appear in FP_Poset. 
             Then removes all parts of graph not connected to a node in the start set. 
    '''
    c = database.conn.cursor()
    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    H = nx.DiGraph() #building networkx graph from FP_poset
    for node in FP_Poset:
        for edge in FP_Poset[node]:
            H.add_edge(node, edge)

    del_list = [] #currently written with FP repeats

    P = deepcopy(cG)

    for edge in P.edges():
        s = scc[edge[0]][0][-1]
        t = scc[edge[1]][0][-1]
        sMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(s))
        MGI = sMGI.fetchone()[0]
        sFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

        tMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(t))
        MGI = tMGI.fetchone()[0]
        tFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

        keep = False
        if (sFP[0],tFP[0]) in H.edges():
            keep = True
        if sFP[0] == tFP[0]:
            keep = True
        if keep == False:
            del_list.append(edge)

    for edge in del_list:   
        P.remove_edge(edge[0],edge[1])

    P.remove_nodes_from(list(nx.isolates(cG)))

    start_set = []
    for node in P:
        p = scc[node][0]
        if p[0] == 0 and p[1] == 0:
            start_set.append(node)

    del_list = []
    for node in P.nodes():
        for i in start_set:
            if node not in start_set:
                try:
                    nx.shortest_path(P, i, node)
                    break
                except:
                    if i == start_set[-1]:
                        del_list.append(node)
                        break
                    else:
                        continue

    for node in del_list:   
        P.remove_node(node)

    stop_set = []
    for node in P:
        p = scc[node][0]
        if p[0] == Hb_max and p[1] == Kni_max:
            stop_set.append(node)

    del_list = []
    for node in P.nodes():
        for i in stop_set:
            if node not in stop_set:
                try:
                    nx.shortest_path(P, node, i)
                    break
                except:
                    if i == stop_set[-1]:
                        del_list.append(node)
                        break
                    else:
                        continue

    for node in del_list:   
        P.remove_node(node)

    return P

def return_start_stop_set(database, graph, scc, Hb_max, Kni_max, FP_Regions):
    '''
    graph: can be in dictinary or networkx form, function expects a condensation graph.
    scc: dictonary of stongly connected components, where keys are node labels in given graph 
         and values are nodes in original graph the condensation is derived from. 
    Hb_max, Kni_max: Highest factor graph layers.
    start_FP_list, stop_FP_list: list of fixed points wanting to constrain the starting and stoping sets to.
    returns: set of nodes considered starting nodes for a path and stoping nodes.  
    '''
    c = database.conn.cursor()
    
    start_set = []
    stop_set = []

    for node in graph:
        n = scc[node][0]
        #print(node, p)
        if n[0] == 0 and n[1] == 0:

            p = scc[node][0][-1]
            MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
            MGI = MGI_result.fetchone()[0]
            FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
            if FP_result[0] in FP_Regions[1] or FP_result[0] in FP_Regions[2]:
                start_set.append(node)


        if n[0] == Hb_max and n[1] == Kni_max:

            p = scc[node][0][-1]
            MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
            MGI = MGI_result.fetchone()[0]
            FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
            if FP_result[0] in FP_Regions[7] or FP_result[0] in FP_Regions[8]:
                stop_set.append(node)

    return start_set, stop_set

def find_breaks_in_FG_comb(database, network_filename, P, scc, Hb_max, Kni_max, FP_Regions):

    breaks = [(0,0), (Hb_max, Kni_max)]
    for h in range(Hb_max+1):
        for k in range(Kni_max+1):
            if (h,k) != (0,0) and (h,k) != (Hb_max, Kni_max):
                
                remove = (h,k)

                T = deepcopy(P)
                for node in P.nodes():
                    if scc[node][0][0:2] == remove:
                        T.remove_node(node)

                start_set, stop_set = return_start_stop_set(database, T, scc, Hb_max, Kni_max, FP_Regions)

                if start_set == []:
                    result = False
                if stop_set == []:
                    result = False
                            
                else:
                    for s in start_set:
                        for t in stop_set:
                            try:
                                nx.shortest_path(T, s, t)
                                result = True
                                break
                            except:
                                if s == start_set[-1]:
                                    if t == stop_set[-1]:
                                        result = False
                                        break
                                    else:
                                        continue
                        else:
                            continue        
                        break

                if result == False:
                    breaks.append((h,k))
    plt.figure(figsize=(5,5)) 
    ax = plt.gca()
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    

    x = []
    y = []
    for s in scc:
        x.append(Hb_max-scc[s][0][0])
        y.append(scc[s][0][1])

    
    plt.scatter(y, x)
    for i in breaks:
        plt.scatter([i[1]],[Hb_max-i[0]], color = 'r')
    
    title = 'Required (red) FG layer combos in ' + network_filename
    plt.title(title)
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    
    plt.savefig(network_filename + '_breaks.png')   # save the figure to file
    plt.close() 

    return breaks

def build_diag(Hb_max, Kni_max, breaks):
    keep = set()

    rk = math.ceil((Kni_max)/(Hb_max))
    rh = math.ceil((Hb_max)/(Kni_max))
    r = max(rk,rh)
    keep.add((0,0))

    if rk >= rh:
        h = 0
        k = 0
        
        while k<Kni_max:
            for i in range(r):
                if k + 1 <= Kni_max:
                    k += 1 
                    keep.add((h,k))
                    if i == r-1:
                        if k + 1 <= Kni_max:
                            keep.add((h,k+1))
                            keep.add((h+2, k-1))
                        h+=1
                        keep.add((h,k))
                        
        h = 0
        k = 0
        while h<Hb_max:
            h+=1
            keep.add((h,k))
            for i in range(r):
                if k + 1 <= Kni_max:
                    k += 1 
                    keep.add((h,k))

    if rk < rh:
        h = 0
        k = 0
        
        while h<Hb_max:
            for i in range(r):
                if h + 1 <= Hb_max:
                    h += 1 
                    keep.add((h,k))
                    if i == r-1:
                        if h + 1 <= Hb_max:
                            keep.add((h+1,k))
                            keep.add((h-1, k+2))
                        k+=1
                        keep.add((h,k))
                        
        h = 0
        k = 0
        while k<Kni_max:
            k+=1
            keep.add((h,k))
            for i in range(r):
                if h + 1 <= Hb_max:
                    h += 1 
                    keep.add((h,k))

    for i in keep.copy():
        keep.add((Hb_max - i[0], Kni_max - i[1]))       

    for i in breaks:
        if i not in keep:
            h = i[0]
            k = i[1]
            keep.add((h,k))
            if k != 0 and k > Kni_max/2:
                while (h,k-1) not in keep:
                    k -= 1
                    keep.add((h,k))
            if k != Kni_max and k < Kni_max/2:
                while (h,k+1) not in keep:
                    k += 1
                    keep.add((h,k))
            if k == Kni_max/2:
                if h < Hb_max/2:
                    while (h,k-1) not in keep:
                        k -= 1
                        keep.add((h,k))

                if h > Hb_max/2:
                    while (h,k+1) not in keep:
                        k += 1
                        keep.add((h,k))
            k = i[1]
            if h != 0 and h > Hb_max/2:
                while (h-1,k) not in keep:
                    h -= 1
                    keep.add((h,k))
            if h != Hb_max and h < Hb_max/2:
                while (h+1,k) not in keep:
                    h += 1
                    keep.add((h,k))
            if h == Hb_max/2:
                if k < Kni_max/2:
                    while (h-1,k) not in keep:
                        h -= 1
                        keep.add((h,k))
                if k > Kni_max/2:
                    while (h+1,k) not in keep:
                        h += 1
                        keep.add((h,k))      
    return keep

def remove_unnecessary_nodes_in_P(P, breaks, keep, scc, Kni_max):
    diagP = deepcopy(P)
    for b in breaks:
        if (b[0]-1,b[1]) in breaks:
            for k in range(b[1]+1, Kni_max+1):
                try:
                    keep.remove((b[0]-1, k))
                except:
                    KeyError
        if b[1]<Kni_max:
            for h in range(b[0]):
                try:
                    keep.remove((h, k+1))
                except:
                    KeyError

    for node in P:
        if scc[node][0][0:2] not in keep:
            diagP.remove_node(node)

    return diagP

def test_any_path_exists_in_product(string, network_filename):
    '''
    string: network string.
    network_filename: Name wanting to save network text file as, expects that .txt not at end.
    returns: True if paths exists, False if no paths exists in product graph.
    '''
    
    # Make DSGRN database and network txt file
    
    txt_filename = network_filename + ".txt"
    f = open(txt_filename,"w") # Make txt file for network, needed to build DSGRN database
    f.write(string)
    f.close()

    net = DSGRN.Network(string)
    pg = ParameterGraph(net)

    db_filename = network_filename + ".db"
    os.system("mpiexec -n 2 Signatures "+ txt_filename + ' ' + db_filename)
    database = Database(db_filename)

    out_edges = get_number_out_edges_from_string(string)
    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    FP_Poset, FP_Regions = get_FP_Poset(out_edges)

    # Compute G (grad_graph) and save.

    G = get_grad_graph_strict_bagged(database, string)
    grad_graph_filename = "grad_graph_strict_" + network_filename

    grad_graph = {}
    for n in G:
        grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    save_json(grad_graph, grad_graph_filename)

    # Compute condensation cG of G

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    # Compute Product graph P and the restricted diagonal Product diagP

    P = get_product_graph(database, cG, scc, FP_Poset)
    lenP_nodes = len(P.nodes())
    lenP_edges = len(P.edges())

    breaks = find_breaks_in_FG_comb(database, network_filename, P, scc, Hb_max, Kni_max, FP_Regions) 
    keep = build_diag(Hb_max, Kni_max, breaks)
    diagP = remove_unnecessary_nodes_in_P(P, breaks, keep, scc, Kni_max)

    plot_FG_layer_comb_in_G(diagP, scc, Hb_max, 'scc FG layer combos ' + network_filename)

    # Test is any paths exist
    start_set, stop_set = return_start_stop_set(database, diagP, scc, Hb_max, Kni_max, FP_Regions)

    if start_set == []:
        #print("Empty start set")
        result = False
    if stop_set == []:
        #print("Empty stop set")
        result = False
                
    else:
        for s in start_set:
            for t in stop_set:
                try:
                    nx.shortest_path(cG, s, t)
                    #print('Path exists from ' + str(s) + ' to '+ str(t))
                    result = True
                    break
                except:
                    if s == start_set[-1]:
                        if t == stop_set[-1]:
                            #print('No Path Exists')
                            result = False
                            break
                    else:
                        continue
            else:
                continue        
            break
        
    return network_filename, (pg.size(), len(G.nodes()), len(G.edges()), len(cG.nodes()), len(cG.edges()), lenP_nodes, lenP_edges, len(diagP.nodes()), len(diagP.edges), result)


def main(network_tup):

    network_filename = 'network' + str(network_tup[0])
    network = get_network_string(network_tup[1], network_tup[-1])

    results = test_any_path_exists_in_product(network, network_filename)

    return results