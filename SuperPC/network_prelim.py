import DSGRN
from DSGRN import *
import collections
import json, os
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import psutil
import math
from copy import deepcopy

from get_FG import *
from get_FP_Poset import *
from networkx_cond import *
from Cut import *

def save_json(edges, filename):
    mydict = {}
    for key in edges.keys():
        if type(key) is not str:
            mydict[str(key)] = edges[key]
            
    with open(filename, 'w') as fp:
        json.dump(mydict, fp)

def check_memory(network_filename, Max_mem):
    virtualMemoryInfo = psutil.virtual_memory()
    MemoryUsage = virtualMemoryInfo.percent

    if MemoryUsage>Max_mem:
        print(MemoryUsage, "WARNING! Memory usage high. Stopping thread for " + network_filename, flush=True)

    return MemoryUsage

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
    pg = ParameterGraph(database.network)
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

def check(database, s, t, scc, FP_Regions):
    ps = scc[s][0][-1]
    pt = scc[t][0][-1]  
    c = database.conn.cursor() 
    sMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(ps))
    MGI = sMGI.fetchone()[0]
    sFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

    tMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(pt))
    MGI = tMGI.fetchone()[0]
    tFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

    keep = False
    for r in FP_Regions:
        if r < 8:
            if sFP[0] in FP_Regions[r] and tFP[0] in FP_Regions[r+1]:
                keep = True
            if sFP[0] == tFP[0]:
                keep = 'same'
    return keep

def P_with_absorbing_nodes(database, N, diagP, scc, FP_Regions, stop_set):
    w = nx.get_edge_attributes(N, 'weight')

    mP = diagP.copy()
    mP.add_edge('leak', 'leak', weight = 1)
    mP.add_edge('rpert', 'rpert', weight = 1)
    mP.add_edge('pert', 'pert', weight = 1)
    mP.add_edge('skip', 'skip', weight = 1)

    for s in stop_set:
        for t in [n for n in mP.neighbors(s)]:
            mP.remove_edge(s,t)
        mP.add_edge(s,s, weight = 1)

    pV = diagP.nodes()
    pE = diagP.edges()

    for n in pV:
        if n not in stop_set:
            leak_sum = 0
            rpert_sum = 0
            skip_sum = 0
            pert_sum = 0
            for (s,t) in N.edges(n):
                if (s,t) not in pE:
                    if t not in pV:
                        leak_sum += w[s,t]
                    else:
                        c = check(database, n, t, scc, FP_Regions)
                        if c == True:
                            rpert_sum += w[s,t]
                            #mP.add_edge(s,t, weight = w[s,t])
                        if c == 'same':
                            pert_sum += w[s,t]
                        if c == False:
                            skip_sum += w[s,t]

            mP.add_edge(n,'leak', weight =leak_sum) #leaving P
            mP.add_edge(n,'rpert', weight =rpert_sum) #pert to next region
            mP.add_edge(n,'skip', weight =skip_sum) #skipping full regions
            mP.add_edge(n,'pert', weight =pert_sum) #pert in region
    return mP

def absorbing_Markov_prob(mP, scc, start_set):
    weight = nx.get_edge_attributes(mP, 'weight')
    nlen = len(mP)
    nodelist = [n for n in mP]

    index = dict(zip(nodelist, range(nlen)))

    M = np.full((nlen, nlen), np.nan, order=None)
    for s in mP:
        for t in mP.neighbors(s):
            M[index[s], index[t]] = weight[s,t]

    M[np.isnan(M)] = 0
    M = np.asarray(M, dtype=None)

    absorb = []
    for s in mP:
        if  M[index[s], index[s]] == 1:
            absorb.append(s)
    column_names = [r for r in mP]
    row_names = [r for r in mP]

    df = pd.DataFrame(M, columns=column_names, index=row_names)

    #Need to put all absorbing node rows to the bottom of the matrix
    new_order = [r for r in mP]

    for n in absorb:
        new_order.remove(n)
        new_order.append(n)
            
    #Now we can get an array reordered to obtain Q and R, new_order is the new ordering of rows/columns
    arr = df[new_order].loc[new_order].to_numpy()

    #Obtain Q!
    temp = []
    for row in range(len(arr)-len(absorb)):
        temp.append(list(arr[row])[:-len(absorb)])
    Q = np.array([np.array(xi) for xi in temp])

    #Obtain R!
    temp = []
    for row in range(len(arr)-len(absorb)):
        temp.append(list(arr[row])[-len(absorb):])
    R = np.array([np.array(xi) for xi in temp])

    n = len(new_order)-len(absorb)
    a = np.zeros((n,))
    size = 0 
    for s in start_set:
        size += len(scc[s])
    for s in start_set:
        i = new_order.index(s)
        a[i] = len(scc[s])/size

    I = np.identity(len(Q))
    N = (np.linalg.inv((I-Q)))
    B = np.dot(N,R)
    V = np.dot(a,B)

    results = {'B sum':sum(V)}
    print(f"{sum(V):.9f}", flush = True)
    l_s = 0
    for i in range(len(V)):
        #if round(B[new_order.index(starting_state)][i]*100,2) != 0.0:
        results[absorb[i]] = V[i]
        print(f"{V[i]:.9f}", '\tending on node', absorb[i], flush = True)
        if absorb[i] == 'skip' or absorb[i] == 'leak':
            l_s += V[i]
    print(f"{l_s:.9f}", flush = True)
    results['leak + skip'] = l_s
    return  results

def test_any_path_exists_in_product(string, network_filename):
    '''
    string: network string.
    network_filename: Name wanting to save network text file as, expects that .txt not at end.
    returns: True if paths exists, False if no paths exists in product graph.
    '''
    Max_mem = 90
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

    if check_memory(network_filename, Max_mem)>Max_mem:
        return network_filename, {'PG size': pg.size(), 'G size': len(G.nodes()), 'G edges': len(G.edges()), 'mem': 'process killed due to memory error'} 

    grad_graph = {}
    for n in G:
        grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    save_json(grad_graph, grad_graph_filename)

    # Compute condensation cG of G

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)

    # Compute Product graph P and the restricted diagonal Product diagP

    P = get_product_graph(database, N, scc, FP_Poset)
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

    if result == True:
        c = find_best_clustering(diagP, network_filename, 15, nodelist = None, data = 'weight', in_out_degree = 'out', save_file = True)[0]
    
    mP =  P_with_absorbing_nodes(database, N, diagP, scc, FP_Regions, stop_set)
    markov_results = absorbing_Markov_prob(mP, scc, start_set)

    results = network_filename, {'PG size': pg.size(), 'G size': len(G.nodes()), 'G edges': len(G.edges()), 'cG size': len(cG.nodes()), 'cG edges': len(cG.edges()), 'P size' : lenP_nodes, 'P edges': lenP_edges, 'diagP size': len(diagP.nodes()), 'diagP edges': len(diagP.edges), 'path exists': result, 'WCut': c, 'markov_results': markov_results}

    print(results, flush=True)

    return results


def main(network_tup):

    network_filename = 'network' + str(network_tup[0])
    network = get_network_string(network_tup[1], network_tup[-1])

    results = test_any_path_exists_in_product(network, network_filename)
    print(results, flush=True)
    return results