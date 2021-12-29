import networkx as nx
import DSGRN
from DSGRN import *

import sys
sys.setrecursionlimit(10**8)
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/notebooks')
from PhenotypeGraphFun import *
from PhenotypeGraphviz import *
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

    
    return database, Hb_list, Kni_list, FP_Poset, FP_keep

def run(net):
    if net == 'StrongEdges':
        return StrongEdges()
    if net == 'Fullconn':
        return Fullconn()
        
if __name__=="__main__":
    database, Hb_list, Kni_list, FP_Poset, FP_keep = run('StrongEdges')
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)

    gradlist = get_gradlist(database, Hb_list, Kni_list)
    grad_graph = get_gradient_graph_parallel(database, gradlist, 4, Hb_list, Kni_list)
    
    save_json(grad_graph, 'StrongEdges_grad_graph')

    G = nx.DiGraph() #building networkx graph grom grad_graph
    for node in grad_graph:
        G.add_node(node)
        for edge in grad_graph[node]:
            G.add_edge(node, edge)

    del_list = [] #list we want to delete from grad_graph

    for node in grad_graph:
        p = node[2]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if not set(FP_result).intersection(set(FP_keep)):
            del_list.append(node)

    for n in del_list: #removes del_list nodes in networkx graph and grad_graph keys
        G.remove_node(n)
        del grad_graph[n]

    nodelist = grad_graph.keys()
    new_grad_graph = {}
    for n in nodelist:
        new_grad_graph[n] = [nbr for nbr in G.neighbors(n)]

    cond, scc = condensation_graph_gradient(database, new_grad_graph, FP_keep)
    
    save_json(cond, 'StrongEdges_grad_graph_cond_test')


#softwrap is alt Z
## how to run in terminal: 
# python script_run.py </dev/null >script_run.log 2>&1 &

# ps -a (shows al process that are running on machine)

# cat script_run.log (prints log)   head does first 10 lines or tail does last 10 lines


## how to use megplux

# ssh elizabethandreas@msu.montana.edu

# pip install DSGRN