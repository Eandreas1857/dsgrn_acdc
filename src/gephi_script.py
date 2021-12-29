import networkx as nx
import DSGRN
from DSGRN import *
import pandas as pd

import sys
sys.setrecursionlimit(10**8)
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/src')
from PhenotypeGraphFun import *
from PhenotypeGraphviz import *
from networkx_cond import *
from GradientFun import *
from CondensationGraph_iter import *
from save_files import *

def create_cond_subgraphs_graphml(database, grad_graph, cond, prod_graph_nodes, path_nodes, scc, start_set, stop_set, Filename):
    ''' graphml filetype '''
    c = database.conn.cursor()

    N = nx.DiGraph() #building networkx graph from cond
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
        for edge in cond[node]:
            G.add_edge(node, edge)
            yes_count = 0
            for s in scc[node]:
                for t in scc[edge]:
                    if N.has_edge(s,t) == True:
                        yes_count += 1
            G[node][edge]['weight'] = yes_count#/(len(scc[node])*len(scc[edge]))

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

    nx.write_graphml(G, Filename)

def create_grad_graph_w_subgraphs_graphml(database, G, cond, prod_graph_nodes, path_nodes, start_set, stop_set, Filename):
    ''' graphml filetype '''
    c = database.conn.cursor()
    
    Kni_att = {}
    Hb_att = {}
    MGI_att = {}
    Region_att = {}
    
    graph = {}
    s_t = {}


    for node in G:
       
            
        p = node[-1]
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
        
        Hb_att[node] = node[0]
        Kni_att[node] = node[1]

        if node in path_nodes:
            graph[node] = 'path'
        elif node in prod_graph_nodes:
            graph[node] = 'product'
        elif node in cond:
            graph[node] = 'cond'
        else:
            graph[node] = 'gradient' 

    for node in start_set:
        s_t[node] = 'starting'
    for node in stop_set:
        s_t[node] = 'stoping'
    
    nx.set_node_attributes(G, 'Hb_FG_layer', Hb_att)
    nx.set_node_attributes(G, 'Kni_FG_layer', Kni_att)
    nx.set_node_attributes(G, 'MGI', MGI_att)
    nx.set_node_attributes(G, 'Region', Region_att)
    nx.set_node_attributes(G, 'group', graph)
    
    nx.set_node_attributes(G, 'start_stop', s_t)

    nx.write_graphml(G, Filename)

def Fullconn():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_Fullconn.db")

    Kni_list = {0: ['00'], 1:['40'], 2:['44','50','C0' ], 3: ['54','C4','D0'], 4: ['55','D4','CC','F0'], 5: ['D5','DC','F4'], 6: ['DD','F5','FC'], 7: ['FD'], 8:['FF']}
    Hb_list = {0:['FF'], 1: ['FD'], 2: ['DD','F5','FC'], 3: ['D5','DC','F4'], 4: ['55','D4','CC','F0'], 5: ['54','C4','D0'], 6:['44','50','C0' ], 7:['40'], 8: ['00']}
    Hb_max = 8
    Kni_max = 8

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

    grad_graph = load_json('Fullconn_grad_graph_Hb_Kni_expanded')

    c = database.conn.cursor()
    G = nx.DiGraph() #building networkx graph from cond
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
        if len(FP_result) != 1:
            del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]
    
    nodelist = grad_graph.keys()
    new_grad_graph = {}
    for n in nodelist:
        new_grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    scc = load_json('Fullconn_scc_correct')
    cond = load_json_cond('Fullconn_cond_correct')
    prod_graph = load_json_cond('Fullconn_prod_graph_correct')

    Paths = []
    for df in pd.read_csv('Fullconn_Prod_G_Paths_wo_multistab_correct.csv',header = None, chunksize=1000000):
        df_list = df.values.tolist()
        Paths += [path for path in df_list]

    starting = ['FP { 2, 1, 0, 0 }', 'FP { 2, 2, 0, 0 }']
    stoping = ['FP { 0, 2, 0, 2 }']

    FP_Region = {1: ['FP { 1, 2, 0, 0 }', 'FP { 1, 1, 0, 0 }', 'FP { 2, 2, 0, 0 }'],
                2: ['FP { 2, 1, 0, 0 }', 'FP { 2, 0, 0, 0 }', 'FP { 2, 1, 1, 0 }'],
                3: ['FP { 2, 0, 1, 0 }', 'FP { 1, 0, 1, 0 }', 'FP { 2, 0, 2, 0 }'],
                4: ['FP { 1, 0, 2, 0 }', 'FP { 0, 0, 2, 0 }', 'FP { 1, 0, 2, 1 }'],
                5: ['FP { 0, 0, 2, 1 }', 'FP { 0, 0, 1, 1 }', 'FP { 0, 0, 2, 2 }'],
                6: ['FP { 0, 0, 1, 2 }', 'FP { 0, 1, 1, 2 }', 'FP { 0, 0, 0, 2 }'],
                7: ['FP { 0, 1, 0, 2 }', 'FP { 0, 2, 0, 2 }', 'FP { 0, 1, 0, 1 }'],
                8: ['FP { 0, 2, 0, 1 }']
                }

    return database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, new_grad_graph, scc, cond, prod_graph, Paths, starting, stoping, FP_Region
   

def StrongEdges():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_StrongEdges.db") 

    Hb_list = {0: ['F'], 1: ['E'], 2: ['A', 'C'], 3: ['8'], 4: ['0']}
    Kni_list = {0:['000'], 1: ['200'], 2: ['208', '240', '600'], 3: ['248','608','640','E00'], 4: ['249','648','618','6C0','E08','E40'], 5:['649','658','6C8','E18','E48','EC0'], 6: ['659','6C9','6D8', 'E49','E58','EC8','E38','FC0'], 7: ['6D9','E59','EC9','ED8','E78','FC8'], 8: ['6DB','ED9','E79','EF8','FC9','FD8'], 9: ['EDB','EF9','FD9','FF8'], 10: ['EFB','FDB','FF9'], 11:['FFB'], 12: ['FFF']}
    Hb_max = 4
    Kni_max = 12
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

    grad_graph = load_json('StrongEdges_grad_graph_new')
 
    c = database.conn.cursor()

    G = nx.DiGraph() #building networkx graph from cond
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
        if len(FP_result) != 1:
            del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]

    del_list = [] 

    for node in G: #want close to hb = 4 kni in (kni, hb) "graph"
        if node[0] == 0:
            if node[1]>4:
                del_list.append(node)
        elif node[0] == 1:
            if node[1]>6:
                del_list.append(node)
        elif node[0] == 3:
            if node[1]<6:
                del_list.append(node)
        elif node[0] == 4:
            if node[1]<8:
                del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]
    
    nodelist = grad_graph.keys()
    new_grad_graph = {}
    for n in nodelist:
        new_grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    scc = load_json('StrongEdges_scc_correct')
    cond = load_json_cond('StrongEdges_cond_correct')
    prod_graph = load_json_cond('StrongEdges_prod_graph_correct')

    Paths = []
    for df in pd.read_csv('StrongEdges_Prod_G_Paths_wo_multistab_correct (copy).csv',chunksize=1000000, engine='python'):
        df_list = df.values.tolist()
        Paths += [path for path in df_list]
    
    new_paths = []
    for path in Paths:
        new_paths += [[int(i) for i in path if str(i) != 'nan']]
    
    starting = ['FP { 1, 1, 0, 0 }', 'FP { 1, 2, 0, 0 }']
    stoping = ['FP { 0, 3, 0, 3 }', 'FP { 0, 2, 0, 3 }', 'FP { 0, 1, 0, 3 }']

    FP_Region = {1: ['FP { 0, 3, 0, 0 }', 'FP { 1, 3, 0, 0 }',
                'FP { 0, 2, 0, 0 }',
                'FP { 0, 1, 0, 0 }'],
                            2: [ 'FP { 1, 2, 0, 0 }',
                'FP { 1, 1, 0, 0 }'],
                            3: ['FP { 1, 0, 0, 0 }', 'FP { 0, 0, 0, 0 }',
                'FP { 1, 0, 1, 0 }'],
                            4: [ 'FP { 0, 0, 1, 0 }'],
                            5: ['FP { 0, 0, 1, 1 }', 'FP { 0, 0, 0, 1 }',
                'FP { 0, 0, 1, 2 }',
                'FP { 0, 0, 0, 2 }',
                'FP { 0, 0, 1, 3 }'],
                            6: ['FP { 0, 0, 0, 3 }'],
                            7: ['FP { 0, 1, 0, 3 }', 'FP { 0, 1, 0, 2 }',
                'FP { 0, 2, 0, 3 }',
                'FP { 0, 2, 0, 2 }',
                'FP { 0, 3, 0, 3 }',
                'FP { 0, 1, 0, 1 }',
                'FP { 0, 2, 0, 1 }'],
                            8: ['FP { 0, 3, 0, 2 }',
                'FP { 0, 3, 0, 1 }']}

    return database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, new_grad_graph, scc, cond, prod_graph, new_paths, starting, stoping, FP_Region


def run(net):
    if net == 'StrongEdges':
        return StrongEdges()
    if net == 'Fullconn':
        return Fullconn()

if __name__=="__main__":
    '''Change comp_grad_gephi to True if computing Gephi file for gradient graph, False if computing Gephi file for condensation'''
    
    comp_grad_gephi = False

    net = 'StrongEdges'
    Filename = 'Testing_script_again.graphml'# MUST HAVE .graphml AT END!!!

    '''Change how edges are weighted for condensation in the function create_cond_subgraphs_graphml'''

    database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, grad_graph, scc, cond, prod_graph, Paths, starting, stoping, FP_Region = run(net)
    
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)

    start_set = []
    stop_set = []
    c = database.conn.cursor()
    for node in prod_graph:
        n = scc[node][0]
        #print(node, p)
        if n[0] == 0 and n[1] == 0:
            p = scc[node][0][-1]
            MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
            MGI = MGI_result.fetchone()[0]
            FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
            if FP_result[0] in starting:
                start_set.append(node)

        if n[0] == Hb_max and n[1] == Kni_max:
            p = scc[node][0][-1]
            MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
            MGI = MGI_result.fetchone()[0]
            FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
            if FP_result[0] in stoping:
                stop_set.append(node)

    path_nodes = []
    for path in Paths:
        for node in path:
            if node not in path_nodes:
                path_nodes.append(int(node))

    if comp_grad_gephi == False:

        prod_graph_nodes = []
        for node in prod_graph:
            prod_graph_nodes.append(node)

        create_cond_subgraphs_graphml(database, grad_graph, cond, prod_graph_nodes, path_nodes, scc, start_set, stop_set, Filename)

    if comp_grad_gephi == True:

        G = nx.DiGraph() #building networkx graph from cond
        for node in grad_graph:
            G.add_node(node)
            for edge in grad_graph[node]:
                G.add_edge(node, edge)

        cond_node_list = []
        for node in cond:
            for s in scc[node]:
                cond_node_list.append(s)
        
        prod_grad_node_list = []
        for node in prod_graph:
            for s in scc[node]:
                prod_grad_node_list.append(s)

        path_grad_node_list = []
        for node in path_nodes:
            for s in scc[node]:
                path_grad_node_list.append(s)

        path_start = []
        for node in start_set:
            for s in scc[node]:
                path_start.append(s)

        path_stop = []
        for node in stop_set:
            for s in scc[node]:
                path_stop.append(s)

        create_grad_graph_w_subgraphs_graphml(database, G, cond_node_list, prod_grad_node_list, path_grad_node_list, path_start, path_stop, Filename)

## how to run in terminal: 
# python gephi_script.py </dev/null >gephi_script.log 2>&1 &

# ps -a (shows al process that are running on machine)

# cat cond_paths_script.log (prints log)   head does first 10 lines or tail does last 10 lines