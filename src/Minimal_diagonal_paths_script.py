import networkx as nx
import DSGRN
from DSGRN import *
import psutil
import random
import csv

import sys
sys.setrecursionlimit(10**8)
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/src')
from PhenotypeGraphFun import *
from PhenotypeGraphviz import *
from networkx_cond import *
from GradientFun import *
from CondensationGraph_iter import *
from save_files import *

def load_json(file_to_load):
    with open(file_to_load) as f:
        data = json.load(f)
        
    new_data = {}
    for key in data:
        new_data[ast.literal_eval(key)] = [tuple(i) for i in data[key] ]
    return new_data



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

    return database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, grad_graph
   

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

    
    return database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, grad_graph

def run(net):
    if net == 'StrongEdges':
        return StrongEdges()
    if net == 'Fullconn':
        return Fullconn()
        
if __name__=="__main__":
    net = 'StrongEdges'

    database, Hb_list, Kni_list, Hb_max, Kni_max, FP_Poset, FP_keep, grad_graph = run(net)

    paths_filename = "Fullconn_Prod_G_Paths_wo_multistab_correct.csv"

    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)

    G = nx.DiGraph() #building networkx graph from cond
    for node in grad_graph:
        G.add_node(node)
        for edge in grad_graph[node]:
            G.add_edge(node, edge)

    #### STRONGEDGES ####

    if net == 'StrongEdges':
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

        strongcc = strongly_connected_components_by_MGI(G, database)

        cG, scc = condensation(G, strongcc)

        remove = [(3,0), (4,0),
            (1,1), (3,1), (4,1),
            (3,2), (4,2),
            (3,3), (4,3),
            (0,4), (3,4), (4,4),
            (0,5), (3,5), (4,5),
            (0,6), (1,6), (4,6),
            (0,7), (1,7), (2,7), (4,7),
            (0,8), (1,8), (2,8), (4,8),
            (0,9), (1,9), (2,9),
            (0,10), (1,10), (2,10),
            (0,11), (1,11), (2,11), (3,11),
            (0,12), (1,12), (2,12), (3,12), (0,1), (0,2), (0,3), (1,2), (1,3), (1,4), (1,5)]

        T = deepcopy(G)     
        for node in G.nodes():
            if node[0:2] in remove:
                T.remove_node(node)
        
        strongcc = strongly_connected_components_by_MGI(T, database)

        cG, scc = condensation(T, strongcc)

        H = nx.DiGraph() #building networkx graph from FP_poset
        for node in FP_Poset:
            for edge in FP_Poset[node]:
                H.add_edge(node, edge)

        del_list = [] #currently written with FP repeats
        
        for edge in cG.edges():
            s = scc[edge[0]][0][-1]
            t = scc[edge[1]][0][-1]
            sMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(s))
            MGI = sMGI.fetchone()[0]
            sFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

            tMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(t))
            MGI = tMGI.fetchone()[0]
            tFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
            #print(node, sFP, tFP)

            keep = False
            if (sFP[0],tFP[0]) in H.edges():
                keep = True
            if sFP[0] == tFP[0]:
                keep = True
            if keep == False:
                del_list.append(edge)

        for edge in del_list:   
            cG.remove_edge(edge[0],edge[1])

        cG.remove_nodes_from(list(nx.isolates(cG)))

        start_set = []
        for node in cG:
            p = scc[node][0]
            #print(node, p)
            if p[0] == 0 and p[1] == 0:
                start_set.append(node)

        del_list = []
        for node in cG.nodes():
            for i in start_set:
                if i != node:
                    try:
                        nx.shortest_path(cG, i, node)
                        break
                    except:
                        if i == start_set[-1]:
                            del_list.append(node)
                            break
                        else:
                            continue

        for node in del_list:   
            cG.remove_node(node)
    
        starting = ['FP { 1, 1, 0, 0 }', 'FP { 1, 2, 0, 0 }']
        stoping = ['FP { 0, 3, 0, 3 }', 'FP { 0, 2, 0, 3 }', 'FP { 0, 1, 0, 3 }']

    # Starting/Stoping for Fullconn
    #starting = [ 'FP { 2, 2, 0, 0 }', 'FP { 2, 1, 0, 0 }']
    #stoping = ['FP { 0, 2, 0, 2 }']

    start_set = []
    stop_set = []

    for node in cG:
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


    Paths = []
    count = 0

    for s in start_set:
        for t in stop_set:
            for path in nx.all_simple_paths(cG, s, t, cutoff=25):
                    count +=1
        
    
    with open(paths_filename, "w") as f:
        wr = csv.writer(f)
        wr.writerows(Paths)


    
## how to run in terminal: 
# python Minimal_diagonal_paths_script.py </dev/null >cond_paths_script.log 2>&1 &

# ps -a (shows al process that are running on machine)

# cat Minimal_diagonal_paths_script.log (prints log)   head does first 10 lines or tail does last 10 lines
