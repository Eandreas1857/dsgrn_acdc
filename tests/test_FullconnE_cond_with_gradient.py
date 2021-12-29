# test_FullconnE_cond_with_gradient.py

import pytest
import sys
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/src')
from PhenotypeGraphFun import *
from PhenotypeGraphFun_Gradient import *

def FullconnE_new_edges_24735():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_FullconnE.db") 
    AP35 = {"Hb":[0,2], "Gt":2, "Kr":0, "Kni":0}
    AP37 = {"Hb":2, "Gt":[0,2], "Kr":0, "Kni":0}
    AP40 = {"Hb":2, "Gt":0, "Kr":[0,2], "Kni":0} 
    AP45 = {"Hb":[0,2], "Gt":0, "Kr":2, "Kni":0} 
    AP47 = {"Hb":[0,2], "Gt":0, "Kr":2, "Kni":0} 
    AP51 = {"Hb":0, "Gt":0, "Kr":2, "Kni":[0,2]} 
    AP57 = {"Hb":0, "Gt":0, "Kr":[0,2], "Kni":2} 
    AP61 = {"Hb":0, "Gt":0, "Kr":[0,2], "Kni":2}
    AP63 = {"Hb":0, "Gt":[0,2], "Kr":0, "Kni":2} 
    AP67 = {"Hb":0, "Gt":2, "Kr":0, "Kni":[0,2]}

    D = [[AP35], [AP37], [AP40], [AP45], [AP51], [AP57], [AP63], [AP67]]
    paramslist = get_paramslist_optimized(database, D, '=')

    FullE_hex = {0: ['C0'], 1:['C4','D0'], 2:['D4'], 3: ['DC','F4'], 4: ['FC']}
    edges = get_phenotype_graph_optimized(database, paramslist, 1)
    
    new_edges = add_phenotype_gradient(database, edges, FullE_hex, FullE_hex)

    return new_edges[(0, 24735)]

def FullconnE_new_edges_1372():
    database = Database("/home/elizabeth/Desktop/ACDC/ACDC_FullconnE.db") 
    AP35 = {"Hb":[0,2], "Gt":2, "Kr":0, "Kni":0}
    AP37 = {"Hb":2, "Gt":[0,2], "Kr":0, "Kni":0}
    AP40 = {"Hb":2, "Gt":0, "Kr":[0,2], "Kni":0} 
    AP45 = {"Hb":[0,2], "Gt":0, "Kr":2, "Kni":0} 
    AP47 = {"Hb":[0,2], "Gt":0, "Kr":2, "Kni":0} 
    AP51 = {"Hb":0, "Gt":0, "Kr":2, "Kni":[0,2]} 
    AP57 = {"Hb":0, "Gt":0, "Kr":[0,2], "Kni":2} 
    AP61 = {"Hb":0, "Gt":0, "Kr":[0,2], "Kni":2}
    AP63 = {"Hb":0, "Gt":[0,2], "Kr":0, "Kni":2} 
    AP67 = {"Hb":0, "Gt":2, "Kr":0, "Kni":[0,2]}

    D = [[AP35], [AP37], [AP40], [AP45], [AP51], [AP57], [AP63], [AP67]]
    paramslist = get_paramslist_optimized(database, D, '=')

    FullE_hex = {0: ['C0'], 1:['C4','D0'], 2:['D4'], 3: ['DC','F4'], 4: ['FC']}
    edges = get_phenotype_graph_optimized(database, paramslist, 1)
    
    new_edges = add_phenotype_gradient(database, edges, FullE_hex, FullE_hex)

    return new_edges[(6, 1372)]

def test_FullconnE_new_edges():
    assert FullconnE_new_edges_24735() == [ (0, 25078), (1, 34339), (1, 25078), (1, 24784), (0, 24742), (0, 24784), (1, 24742), (0, 24734), (0, 34339), (1, 24734), (1, 24735)]

    assert FullconnE_new_edges_1372() == [ (7, 10976), (6, 3773), (7, 1379), (6, 2058), (6, 1470), (7, 2058), (6, 1386), (6, 6174), (6, 1379), (7, 1470), (7, 1386), (7, 1421), (6, 10976), (7, 6174), (6, 1421), (7, 1372), (7, 3773)]
