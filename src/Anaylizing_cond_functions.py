import DSGRN
from DSGRN import *
import networkx as nx
import matplotlib.pyplot as plt
from copy import deepcopy
import os
from all_networks_with_n_nodes_e_edges import *
from save_files import *
from GradientFun import *
from get_FG import *
from get_FP_Poset import *
from networkx_cond import *

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
            if i != node:
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

    return P

def return_start_stop_set(database, graph, scc, Hb_max, Kni_max, start_FP_list = None, stop_FP_list = None):
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
            if start_FP_list != None:
                p = scc[node][0][-1]
                MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
                MGI = MGI_result.fetchone()[0]
                FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
                if FP_result in start_FP_list:
                    start_set.append(node)
            else:
                start_set.append(node)

        if n[0] == Hb_max and n[1] == Kni_max:
            if stop_FP_list != None:
                p = scc[node][0][-1]
                MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
                MGI = MGI_result.fetchone()[0]
                FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
                if FP_result in stop_FP_list:
                    stop_set.append(node)
            else:
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

    FP_Poset = get_FP_Poset(out_edges)[0]

    # If grad_graph has not already been computed for this network, compute it and save.
    if grad_graph == None:
        gradlist = get_gradlist_strict(database, Hb_list, Kni_list)
        grad_graph = get_gradient_graph_parallel(database, gradlist, 7, Hb_list, Kni_list)
        grad_graph_filename = "grad_graph_strict_"+network_filename
        save_json(grad_graph, grad_graph_filename)

        if reduce == True:
            G, ngg = reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset)
            ngg_filename = "reduced_grad_graph_strict_"+network_filename
            save_json(ngg, ngg_filename)

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    P = get_product_graph(database, cG, scc, FP_Poset)

    start_set, stop_set = return_start_stop_set(database, P, scc, Hb_max, Kni_max)

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

    return result

def find_breaks_in_FG_comb(database, P, scc, Hb_max, Kni_max):

    breaks = []
    for h in range(Hb_max+1):
        for k in range(Kni_max+1):
            if (h,k) != (0,0) and (h,k) != (Hb_max, Kni_max):
                
                remove = (h,k)

                T = deepcopy(P)
                for node in P.nodes():
                    if scc[node][0][0:2] == remove:
                        T.remove_node(node)

                start_set, stop_set = return_start_stop_set(database, T, scc, Hb_max, Kni_max)

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
    x = []
    y = []

    for s in scc:
        x.append(Hb_max-scc[s][0][0])
        y.append(scc[s][0][1])

    plt.scatter(y, x)
    for i in breaks:
        plt.scatter([Hb_max-i[0]],[i[1]], color = 'r')
    plt.xlabel('Kni Facter Graph Layer')
    plt.ylabel('Hb Facter Graph Layer')
    plt.show()

    return breaks

def create_cond_subgraphs_graphml(database, grad_graph, cond, prod_graph_nodes, path_nodes, scc, FP_Region, start_set, stop_set, Filename):
    ''' graphml filetype '''
    c = database.conn.cursor()

    N = nx.DiGraph()
    for node in grad_graph:
        N.add_node(node)
        for edge in grad_graph[node]:
            N.add_edge(node, edge)

    G = nx.DiGraph()
    
    Kni_att = {}
    Hb_att = {}
    MGI_att = {}
    Region_att = {}
    scc_size_att = {}
    graph = {}
    s_t = {}

    for node in cond:
        G.add_node(node)
        count = 0
        for edge in cond[node]:
            G.add_edge(node, edge)
            yes_count = 0
            for s in scc[node]:
                for t in scc[edge]:
                    if N.has_edge(s,t) == True:
                        yes_count += 1
                        count +=1
            G[node][edge]['weight'] = yes_count
            
        for edge in cond[node]:
            G[node][edge]['weight'] = G[node][edge]['weight']/count

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
        elif node in prod_graph_nodes:
            graph[node] = 'product'
        else:
            graph[node] = 'cond' 

        scc_size_att[node] = len(scc[node])
    for node in start_set:
        s_t[node] = 'starting'
    for node in stop_set:
        s_t[node] = 'stoping'
    
    nx.set_node_attributes(G, 'Hb_FG_layer', Hb_att)
    nx.set_node_attributes(G, 'Kni_FG_layer', Kni_att)
    nx.set_node_attributes(G, 'MGI', MGI_att)
    nx.set_node_attributes(G, 'Region', Region_att)
    nx.set_node_attributes(G, 'group', graph)
    nx.set_node_attributes(G, 'scc size', scc_size_att)
    nx.set_node_attributes(G, 'start_stop', s_t)
    
    group=nx.get_node_attributes(G,'group')
    att = {}
    for edge in G.edges():
        s = edge[0]
        t = edge[1]
        if group[s] == 'path':
            if group[t] != 'path':
                att[s] = 'leave path'
    nx.set_node_attributes(G, 'leaving', att)

    nx.write_graphml(G, Filename)

def get_gephi_graph_for_cond(database, network, grad_graph, graphml_filename, path_nodes = []):
    '''
    grad_graph: expects graph as dictinary
    network_txt_filename: filename and place where txt file format of the network string is saved.
    graphml_filename: name wanting for graphml file, will add location automatically. Expects .graphml at end.
    '''

    out_edges = get_number_out_edges_from_string(network)
    FP_Poset, FP_Region = get_FP_Poset(out_edges)

    G = reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset)[0]

    strongcc = strongly_connected_components_by_MGI(G, database)

    cG, scc = condensation(G, strongcc)
    P = get_product_graph(database, cG, scc, FP_Poset)

    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1
    start_set, stop_set = return_start_stop_set(database, P, scc, Hb_max, Kni_max)

    filename = '/home/elizabeth/Desktop/GIT/dsgrn_acdc/Saved_Files/Graphml/' + graphml_filename
    
    ### Notice graph only has bagged FP in it, the condensation of the gradient graph only, without removing all nodes not in bag is much larger.
    create_cond_subgraphs_graphml(database, grad_graph, cG, P, path_nodes, scc, FP_Region, start_set, stop_set, filename)