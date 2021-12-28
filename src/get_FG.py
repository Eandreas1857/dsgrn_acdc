import DSGRN
from DSGRN import *
import sys, os

# Disable printing
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def get_FG_layer(graph, current_layer):
    new_layer = []
    while current_layer != []:
        node = current_layer.pop()
        for edge in graph.edges:
            if node == edge[0]:
                new_layer.append(edge[-1])
    return list(set(new_layer))

def get_FG_layer_hex(graph, current_layer):

    FG_layer = [graph.data[node][0] for node in current_layer]
    
    return FG_layer

def get_hex_FG(database, gene):
    
    single_gene_query = SingleGeneQuery(database, gene)
    graph = single_gene_query(2)

    FG = {}
    layer = 0
    current_layer = [0]
    FG[layer] = get_FG_layer_hex(graph, current_layer)
    new_layer = get_FG_layer(graph,current_layer)

    while new_layer !=[]:
        layer +=1
        current_layer = new_layer
        FG[layer] = get_FG_layer_hex(graph, current_layer)
        new_layer = get_FG_layer(graph,current_layer)

    return FG

def get_Hb_Kni_list(database):

    with HiddenPrints():
        Hb_list = {}
        Hb = get_hex_FG(database, 'Hb')
        for key in reversed(list(Hb.keys())):
            Hb_list[len(Hb)-1-key] = Hb[key]  #all code requires Hb_list to be in reverse order, as it simplfies code.

        Kni_list = get_hex_FG(database, 'Kni')

    return Hb_list, Kni_list


