
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import math
from copy import deepcopy
import os
from all_networks_with_n_nodes_e_edges import *
from save_files import *
from GradientFun import *
from get_FG import *
from get_FP_Poset import *
from networkx_cond import *
from metrics import *

import sys
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/Optimized_Ncut_Directed_and_Undirected/src')
from Clustering_by_weighted_cuts_in_directed_graphs import *

import collections

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

    monostable_query_object = MonostableQuery(database)
    matches = monostable_query_object.matches()

    mgs = {}
    for row in c.execute('select MorseGraphIndex, Label from MorseGraphAnnotations'):
        if row[1] in FP_keep:
            if row[0] in matches:
                if row[0] in mgs:
                    mgs[row[0]].append(row[1][0:2]) 
                else:
                    mgs[row[0]] = [row[1][0:2]] 

    set_of_MGIm = list(i for i in mgs if mgs[i].count('FP') == 1)

    string = 'select * from Signatures where MorseGraphIndex in ({seq})'.format(seq=','.join(['?'] * len(set_of_MGIm)))
    PGIset = [row[0] for row in c.execute(string, set_of_MGIm)]

    adj = {}
    S = set(PGIset) 
    for s in PGIset:
        adj[s] = set(pg.adjacencies(s, 'codim1')).intersection(S)

    layers_PGI = {}
    for s in PGIset:
        logic = (pg.parameter(s)).logic()
        sHb = Hb[((logic[0]).stringify())[6:-2]]
        sKni = Kni[((logic[3]).stringify())[6:-2]]
        layers_PGI[s] = (sHb, sKni)

    G = nx.DiGraph()
    for s in PGIset:
        for t in adj[s]:
            sHb = layers_PGI[s][0]
            sKni = layers_PGI[s][1]
            tHb = layers_PGI[t][0]
            tKni = layers_PGI[t][1] 
            if sHb+1 == tHb and sKni == tKni:
                G.add_edge((sHb, sKni, s), (tHb, tKni, t))  
            elif sHb == tHb and sKni+1 == tKni:
                G.add_edge((sHb, sKni, s), (tHb, tKni, t))
            elif sHb == tHb and sKni == tKni:
                G.add_edge((sHb, sKni, s), (tHb, tKni, t))
    return G

def reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset):
    c = database.conn.cursor()

    FP_keep = [node for node in FP_Poset.keys()]

    G = nx.DiGraph() #building networkx graph
    for node in grad_graph:
        G.add_node(node)
        for edge in grad_graph[node]:
            G.add_edge(node, edge)

    del_list = []

    for node in grad_graph:
        p = node[-1]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if not set(FP_result).intersection(set(FP_keep)):
            del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]

    return G, grad_graph

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

def test_any_path_exists_in_product(string, network_filename, database = None, grad_graph = None, reduce = True):
    '''
    string: network string.
    network_filename: Name wanting to save network text file as, expects that .txt not at end.
    returns: True if paths exists, False if no paths exists in product graph.
    '''
    
    # Make DSGRN database
    if database == None:
        txt_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".txt"
        f = open(txt_filename,"w") # Make txt file for network, needed to build DSGRN database
        f.write(string)
        f.close()
        db_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".db"
        os.system("mpiexec -n 2 Signatures "+ txt_filename + ' ' + db_filename)
        database = Database(db_filename)

    out_edges = get_number_out_edges_from_string(string)
    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    FP_Poset, FP_Regions = get_FP_Poset(out_edges)

    # If grad_graph has not already been computed for this network, compute it and save.
    if grad_graph == None:
        G = get_grad_graph_strict_bagged(database, string)
        grad_graph_filename = "grad_graph_strict_"+network_filename

        grad_graph = {}
        for n in G:
            grad_graph[n] = [nbr for nbr in G.neighbors(n)]

        save_json(grad_graph, grad_graph_filename)

    else:
        G = nx.DiGraph() #building networkx graph
        for node in grad_graph:
            G.add_node(node)
            for edge in grad_graph[node]:
                G.add_edge(node, edge)

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    P = get_product_graph(database, cG, scc, FP_Poset)

    start_set, stop_set = return_start_stop_set(database, P, scc, Hb_max, Kni_max, FP_Regions)

    if start_set == []:
        print("Empty start set")
        result = False
    if stop_set == []:
        print("Empty stop set")
        result = False
                
    else:
        for s in start_set:
            for t in stop_set:
                try:
                    nx.shortest_path(cG, s, t)
                    print('Path exists from ' + str(s) + ' to '+ str(t))
                    result = True
                    break
                except:
                    if s == start_set[-1]:
                        if t == stop_set[-1]:
                            print('No Path Exists')
                            result = False
                            break
                    else:
                        continue
            else:
                continue        
            break

    if result == True:
        find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max, FP_Regions)
        
    return result

def find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max, FP_Regions):

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
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    
    plt.show()

    return breaks

def add_source_weight_to_cond(G, cond, scc):

    for node in cond:
        count = 0
        for edge in cond[node]:
            yes_count = 0
            for s in scc[node]:
                for t in scc[edge]:
                    if G.has_edge(s,t) == True:
                        yes_count += 1
                        count +=1
            cond[node][edge]['weight'] = yes_count
            
        for edge in cond[node]:
            cond[node][edge]['weight'] = cond[node][edge]['weight']/count
    return cond

def avg_two_way_edge_weight_for_OpenOrd(cG, data = 'weight'):
    weight = nx.get_edge_attributes(cG, data)
    edges = cG.edges()
    for s,t in edges:
        if (t,s) in edges:
            w1 = weight[s,t]
            w2 = weight[t,s]
            avg_w = (w1+w2)/2
            cG[s][t]['weight'] = avg_w
            cG[t][s]['weight'] = avg_w
    return cG

def load_G_cG_P_weighted(network_filename, grad_graph_filename, return_product = True):
    database_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/Saved_Files/MegaPlex_results/Analysing_random_networks_20220212/" + network_filename + ".db"
    database = Database(database_filename) 
    network_txt_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/Saved_Files/MegaPlex_results/Analysing_random_networks_20220212/" + network_filename + ".txt"

    with open(network_txt_filename,"r") as f:
        network = f.read()

    grad_graph = load_json(grad_graph_filename)

    out_edges = get_number_out_edges_from_string(network)
    FP_Poset, FP_Regions = get_FP_Poset(out_edges)
    G = reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset)[0]

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)

    if return_product == True:
        Hb_list, Kni_list = get_Hb_Kni_list(database)
        Hb_max = len(Hb_list)-1
        Kni_max = len(Kni_list)-1
        P = get_product_graph(database, N, scc, FP_Poset)

        breaks = find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max, FP_Regions) 
        keep = build_diag(Hb_max, Kni_max, breaks)
        diagP = remove_unnecessary_nodes_in_P(P, breaks, keep, scc, Kni_max)

    return (G, cG, diagP, scc) if return_product == True else (G, cG)

def create_cond_subgraphs_graphml(database, cond, k, P, path_nodes, scc, FP_Region, start_set, stop_set, Filename, in_out_degree):
    ''' graphml filetype '''
    c = database.conn.cursor()
    
    Kni_att = {}
    Hb_att = {}
    MGI_att = {}
    Region_att = {}
    scc_size_att = {}
    graph = {}
    s_t = {}
    C = {}
    PC = {}
    for node in cond:
            
        p = scc[node][0][-1]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]

        MGI_att[node] = MGI
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if len(FP_result) == 1:
            for r in FP_Region:
                if FP_result[0] in FP_Region[r]:
                    Region_att[node] = r
        else:
            Region_att[node] = 'not mono-stable'
        
        Hb_att[node] = scc[node][0][0]
        Kni_att[node] = scc[node][0][1]

        if node in path_nodes:
            graph[node] = 'path'
        elif node in P:
            graph[node] = 'product'
        else:
            graph[node] = 'cond' 

        scc_size_att[node] = len(scc[node])
    for node in start_set:
        s_t[node] = 'starting'
    for node in stop_set:
        s_t[node] = 'stoping'
    
    nx.set_node_attributes(cond, Hb_att, 'Hb_FG_layer')
    nx.set_node_attributes(cond, Kni_att, 'Kni_FG_layer')
    nx.set_node_attributes(cond, MGI_att, 'MGI')
    nx.set_node_attributes(cond, Region_att, 'Region')
    nx.set_node_attributes(cond, graph, 'Group')
    nx.set_node_attributes(cond, scc_size_att, 'scc size')
    nx.set_node_attributes(cond, s_t, 'start_stop')
    
    N = deepcopy(cond)

    undir = deepcopy(cond).to_undirected()          
    for component in list(nx.connected_components(undir)):
        for node in start_set:
            if node in component:
                break
            elif node == start_set[-1]:
                for node in component:
                    N.remove_node(node)

    Y, cut, clusters = asym_optimized_normalized_cut(N, k, nodelist = None, data = 'weight', in_out_degree = in_out_degree)
    print('WCut_cG', cut)
    for c in clusters:
        for n in clusters[c]:
            C[n] = c
    nx.set_node_attributes(cond, C, 'cG_Clusters')

    Y, cut, clusters = asym_optimized_normalized_cut(P, k, nodelist = None, data = 'weight', in_out_degree = in_out_degree)
    print('WCut_P', cut)

    for c in clusters:
        for n in clusters[c]:
            PC[n] = c
    nx.set_node_attributes(cond, PC, 'P_Clusters')

    nx.write_graphml(cond, Filename)

def get_gephi_graph_for_cond(database, network, G, graphml_filename, path_nodes = []):
    '''
    grad_graph: expects graph as dictinary
    network_txt_filename: filename and place where txt file format of the network string is saved.
    graphml_filename: name wanting for graphml file, will add location automatically. Expects .graphml at end.
    '''

    out_edges = get_number_out_edges_from_string(network)
    FP_Poset, FP_Regions = get_FP_Poset(out_edges)

    #G = reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset)[0]

    strongcc = strongly_connected_components_by_MGI(G, database)

    cG, scc = condensation(G, strongcc)
    P = get_product_graph(database, cG, scc, FP_Poset)
    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    breaks = find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max, FP_Regions) 
    keep = build_diag(Hb_max, Kni_max, breaks)
    diagP = remove_unnecessary_nodes_in_P(P, breaks, keep, Kni_max)
    
    start_set, stop_set = return_start_stop_set(database, diagP, scc, Hb_max, Kni_max, FP_Regions)

    filename = '/home/elizabeth/Desktop/GIT/dsgrn_acdc/Saved_Files/Graphml/' + graphml_filename
    
    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)
    avg_two_way_edge_weight_for_OpenOrd(N, data = 'weight')

    ### Notice graph only has bagged FP in it, the condensation of the gradient graph only, without removing all nodes not in bag is much larger.
    create_cond_subgraphs_graphml(database, N, k, diagP, path_nodes, scc, FP_Regions, start_set, stop_set, filename)

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

def get_info(network_filename, graphml_filename, grad_graph_filename, k):
    database_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".db"
    database = Database(database_filename) 
    network_txt_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".txt"

    with open(network_txt_filename,"r") as f:
        network = f.read()

    net = DSGRN.Network(network)
    pg = ParameterGraph(net)
    print('pg size', pg.size())

    grad_graph = load_json(grad_graph_filename)
    print('grad size', len(grad_graph))

    out_edges = get_number_out_edges_from_string(network)
    FP_Poset, FP_Regions = get_FP_Poset(out_edges)
    
    print('r grad size', len(G), len(G.edges()))
    avg, prob, percentile = score_region_transitions(database, G, FP_Regions)
    print("avg, prob, percentile: ", avg, prob, percentile)
    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)
    print('cond size', len(cG), len(cG.edges()))

    undir = deepcopy(cG).to_undirected()          
    for component in list(nx.connected_components(undir)):
        print(len(component))

    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    P = get_product_graph(database, cG, scc, FP_Poset)
    print('prod size', len(P), len(P.edges()))

    breaks = find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max, FP_Regions) 
    keep = build_diag(Hb_max, Kni_max, breaks)
    diagP = remove_unnecessary_nodes_in_P(P, breaks, keep, scc, Kni_max)
    print('diag prod size', len(diagP), len(diagP.edges()))
    plot_FG_layer_comb_in_G(diagP, scc, Hb_max, 'scc FG layer combos ' + network_filename)

    start_set, stop_set = return_start_stop_set(database, diagP, scc, Hb_max, Kni_max, FP_Regions)

    filename = '/home/elizabeth/Desktop/GIT/dsgrn_acdc/Saved_Files/Graphml/' + graphml_filename

    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)
    avg_two_way_edge_weight_for_OpenOrd(N, data = 'weight')

    ### Notice graph only has bagged FP in it, the condensation of the gradient graph only, without removing all nodes not in bag is much larger.
    create_cond_subgraphs_graphml(database, N, k, diagP, [], scc, FP_Regions, start_set, stop_set, filename)

def plot_FG_layer_comb(keep, Hb_max, breaks, title):
    plt.figure(figsize=(5,5))
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    x = []
    y = []

    for s in keep:
        x.append(Hb_max - s[0])
        y.append(s[1])

    plt.scatter(y, x)
    for i in breaks:
        if i not in keep:
            plt.scatter([i[1]],[Hb_max - i[0]], color = 'r')
        else:
            plt.scatter([i[1]],[Hb_max - i[0]], color = 'g')
    plt.title(title)
    
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    plt.show()

def plot_FG_layer_comb_in_G(G, scc, Hb_max, title):
    plt.figure(figsize=(5,5))
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    x = []
    y = []

    for s in G:
        x.append(Hb_max-scc[s][0][0])
        y.append(scc[s][0][1])

    plt.scatter(y, x)
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    plt.title(title)
    plt.show()





