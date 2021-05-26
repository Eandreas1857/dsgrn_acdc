import DSGRN, itertools
from DSGRN import *
from NFixedPointQuery import NFixedPointQuery

from copy import deepcopy

import multiprocessing
from functools import partial
import itertools


def stronglycc_iterative(vertices, edges):
    """
    This is a non-recursive version of strongly_connected_components_path.
    See the docstring of that function for more details.
    Examples
    --------
    Example from Gabow's paper [1]_.
    >>> vertices = [1, 2, 3, 4, 5, 6]
    >>> edges = {1: [2, 3], 2: [3, 4], 3: [], 4: [3, 5], 5: [2, 6], 6: [3, 4]}
    >>> for scc in strongly_connected_components_iterative(vertices, edges):
    ...     print(scc)
    ...
    set([3])
    set([2, 4, 5, 6])
    set([1])
    Example from Tarjan's paper [2]_.
    >>> vertices = [1, 2, 3, 4, 5, 6, 7, 8]
    >>> edges = {1: [2], 2: [3, 8], 3: [4, 7], 4: [5],
    ...          5: [3, 6], 6: [], 7: [4, 6], 8: [1, 7]}
    >>> for scc in  strongly_connected_components_iterative(vertices, edges):
    ...     print(scc)
    ...
    set([6])
    set([3, 4, 5, 7])
    set([8, 1, 2])
    """
    identified = set()
    stack = []
    index = {}
    boundaries = []

    for v in vertices:
        if v not in index:
            to_do = [('VISIT', v)]
            while to_do:
                operation_type, v = to_do.pop()
                if operation_type == 'VISIT':
                    index[v] = len(stack)
                    stack.append(v)
                    boundaries.append(index[v])
                    to_do.append(('POSTVISIT', v))
                    # We reverse to keep the search order identical to that of
                    # the recursive code;  the reversal is not necessary for
                    # correctness, and can be omitted.
                    to_do.extend(
                        reversed([('VISITEDGE', w) for w in edges[v]]))
                elif operation_type == 'VISITEDGE':
                    if v not in index:
                        to_do.append(('VISIT', v))
                    elif v not in identified:
                        while index[v] < boundaries[-1]:
                            boundaries.pop()
                else:
                    # operation_type == 'POSTVISIT'
                    if boundaries[-1] == index[v]:
                        boundaries.pop()
                        scc = list(stack[index[v]:])
                        del stack[index[v]:]
                        identified.update(scc)
                        yield scc

def get_phenotype_graph_gradient(database, paramslist, num_processes, Hb_list, Kni_list, repeat_layer=True):
    '''
    Perform successive intersections of allowable parameter transitions with codimension 1 parameter adjacencies. This results in a list of pairs of neighboring parameters with the desired bounds properties. Allows repeats of bound matches except for the last list. The result is saved as a graph, making the end result amenable to graph searches.
    :param database: DSGRN.Database object
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    :param Hb_list: A dict of the hexcode poset for Hb where the values are hexcodes and the keys are the hexcode orders, 0 is low.
    :param Kni_list: A dict of the hexcode poset for Kni where the values are hexcodes and the keys are the hexcode orders, 0 is low.
    :param repeat_layer: allow repeated bounds matches except at the last layer.
    :return: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    '''

    pg = DSGRN.ParameterGraph(database.network)
    todo = list(itertools.chain.from_iterable(paramslist))
    
    pool = multiprocessing.Pool(processes=num_processes)
    work = partial(todo_node_gradient, pg = pg, paramslist = paramslist, Hb_list = Hb_list, Kni_list = Kni_list)
    results = pool.map_async(work, todo)

    edges = {}
    for m,p in results.get():
        edges[m] = p

    return edges

def todo_node_gradient(node, pg, paramslist, Hb_list, Kni_list, repeat_layer = True):
    (k,p) = node
    
    params_p = pg.parameter(p)
    logic_p = params_p.logic()
    Hb_logic_p = logic_p[0]
    Kni_logic_p = logic_p[3]
    
    Hb_hex = Hb_logic_p.stringify()
    Kni_hex = Kni_logic_p.stringify()
    
    for n in Hb_list:
        if Hb_hex[6:-2] in Hb_list[n]:
            Hb_hex_layer = n
            
            
    for j in Kni_list:
        if Kni_hex[6:-2] in Kni_list[j]:
            Kni_hex_layer = j
            
    next_mg_steps = paramslist[k+1][:] if k < len(paramslist)-1 else []
    if repeat_layer:
        next_mg_steps += paramslist[k][:]
    # find neighboring parameters using DSGRN adjacencies function to get all possible neighbors
    # accounting for the same parameter in adjacent layers
    adj = list(pg.adjacencies(p, 'codim1'))
    possible_neighbors = [(k + 1, q) for q in adj + [p]]
    if repeat_layer:
        possible_neighbors += [(k, q) for q in adj]
    almost = list(set(next_mg_steps).intersection(possible_neighbors))
    
    for i in almost:
        
        params_n = pg.parameter(i[-1])
        logic_n = params_n.logic()
        Hb_logic_n = logic_n[0]
        Kni_logic_n = logic_n[3]

        Hb_hex_n = Hb_logic_n.stringify()
        Kni_hex_n = Kni_logic_n.stringify()
        
        for k in Hb_list:
            if Hb_hex_n[6:-2] in Hb_list[k]:
                Hb_hex_layer_n = k
               
                
        for l in Kni_list:
            if Kni_hex_n[6:-2] in Kni_list[l]:
                Kni_hex_layer_n = l
                
        if Hb_hex_layer != Hb_hex_layer_n:
            if Hb_hex_layer-1 != Hb_hex_layer_n:
                almost.remove(i)

        if i in almost:
            if Kni_hex_layer != Kni_hex_layer_n:
                if Kni_hex_layer+1 != Kni_hex_layer_n:
                    almost.remove(i)   
    return (node, almost)


def add_phenotype_gradient(database, edges, Hb_list, Kni_list):
    new_edges = deepcopy(edges)
    pg = DSGRN.ParameterGraph(database.network)
    for key in edges:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()
            
        c_Hb_hex_layer = [k for k in Hb_list if c_Hb_hex[6:-2] in Hb_list[k] ]
        if c_Hb_hex_layer == [] or len(c_Hb_hex_layer) >= 2:
            print('Hb_FAILED', c_Hb_hex_layer, key)
        
        c_Hb_layer = c_Hb_hex_layer.pop() 
        
        c_Kni_hex_layer = [j for j in Kni_list if c_Kni_hex[6:-2] in Kni_list[j]]
        if c_Kni_hex_layer == [] or len(c_Kni_hex_layer) >= 2:
            print('Kni_FAILED', key)
            
        c_Kni_layer = c_Kni_hex_layer.pop()
            
        for node in edges[key]:
            layer, pi = node
            params = pg.parameter(pi)
            logic = params.logic()
            Hb_logic = logic[0]
            Kni_logic = logic[3]

            Hb_hex = Hb_logic.stringify()
            Kni_hex = Kni_logic.stringify()
            
            Hb_hex_layer = [l for l in Hb_list if Hb_hex[6:-2] in Hb_list[l]   ]  
            if Hb_hex_layer == [] or len(Hb_hex_layer) >= 2:
                print('Hb_FAILED', Hb_hex_layer, node)
            Hb_layer = Hb_hex_layer.pop()
            
            Kni_hex_layer = [m for m in Kni_list if Kni_hex[6:-2] in Kni_list[m] ]
            
            if Kni_hex_layer == [] or len(Kni_hex_layer) >= 2:
                print('Kni_FAILED', node)
            
            Kni_layer = Kni_hex_layer.pop()
            
            if Hb_layer != c_Hb_layer:
                if Hb_layer != c_Hb_layer-1:
                    #print('Hb_removed', Hb_layer, c_Hb_layer, key, node)
                    new_edges[key].remove(node)

            if node in edges[key]:
                if Kni_layer != c_Kni_layer:
                    if Kni_layer != c_Kni_layer+1:
                        #print('Kni_removed', Kni_layer, c_Kni_layer, key, node)
                        new_edges[key].remove(node) 
    return new_edges


def cond_h2l_gradient(database, cond_dict, redu_params):
    print('Remember to change phenotype layer numb if needed')
    pg = DSGRN.ParameterGraph(database.network)
    new_cond_dict = deepcopy(cond_dict)
    del_list = []
    for key in cond_dict:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_layer == 0:
            if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                
                continue
            if c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)

        if c_layer == 7:       
            if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                continue
            if c_Kni_hex[6:-2] != 'F'*len(c_Kni_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
    
    new_paramslist = deepcopy(redu_params)
    for i in redu_params[0]:
        c_params = pg.parameter(i[-1])
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
            new_paramslist[0].remove(i)
            
            
        elif c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
            new_paramslist[0].remove(i)
            

    for i in redu_params[7]:
        c_layer, c_pi = i
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()

        if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
            new_paramslist[7].remove(i)
            
        elif c_Kni_hex[6:-2] != 'F'*len(c_Kni_hex[6:-2]):
            new_paramslist[7].remove(i)
            
            
    final_cond_dict = remove_dict_value(new_cond_dict, del_list)
    
    return final_cond_dict, new_paramslist, del_list

def cond_h2l_Kni_gradient(database, cond_dict, redu_params):
    print('Remember to change phenotype layer numb if needed')
    pg = DSGRN.ParameterGraph(database.network)
    new_cond_dict = deepcopy(cond_dict)
    del_list = []
    for key in cond_dict:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Kni_logic = c_logic[3]
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_layer == 0:
            if c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)

        if c_layer == 7:  
            if c_Kni_hex[6:-2] != 'F'*len(c_Kni_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
    
    new_paramslist = deepcopy(redu_params)
    for i in redu_params[0]:
        c_params = pg.parameter(i[-1])
        c_logic = c_params.logic()
        c_Kni_logic = c_logic[3]
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
            new_paramslist[0].remove(i)
            

    for i in redu_params[7]:
        c_layer, c_pi = i
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Kni_logic = c_logic[3]
        c_Kni_hex = c_Kni_logic.stringify()

        if c_Kni_hex[6:-2] != 'F'*len(c_Kni_hex[6:-2]):
            new_paramslist[7].remove(i)
            
            
    final_cond_dict = remove_dict_value(new_cond_dict, del_list)
    
    return final_cond_dict, new_paramslist, del_list

def add_phenotype_Hb_gradient(database, edges, Hb_list, Kni_list):
    new_edges = deepcopy(edges)
    pg = DSGRN.ParameterGraph(database.network)
    for key in edges:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]

        c_Hb_hex = c_Hb_logic.stringify()
            
        c_Hb_hex_layer = [k for k in Hb_list if c_Hb_hex[6:-2] in Hb_list[k] ]
        if c_Hb_hex_layer == [] or len(c_Hb_hex_layer) >= 2:
            print('Hb_FAILED', c_Hb_hex_layer, key)
        
        c_Hb_layer = c_Hb_hex_layer.pop() 
                    
        for node in edges[key]:
            layer, pi = node
            params = pg.parameter(pi)
            logic = params.logic()
            Hb_logic = logic[0]

            Hb_hex = Hb_logic.stringify()
            
            Hb_hex_layer = [l for l in Hb_list if Hb_hex[6:-2] in Hb_list[l]   ]  
            if Hb_hex_layer == [] or len(Hb_hex_layer) >= 2:
                print('Hb_FAILED', Hb_hex_layer, node)
            Hb_layer = Hb_hex_layer.pop()
            
            if Hb_layer != c_Hb_layer:
                if Hb_layer != c_Hb_layer-1:
                    #print('Hb_removed', Hb_layer, c_Hb_layer, key, node)
                    new_edges[key].remove(node)

    return new_edges

def add_phenotype_Kni_gradient(database, edges, Hb_list, Kni_list):
    new_edges = deepcopy(edges)
    pg = DSGRN.ParameterGraph(database.network)
    for key in edges:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Kni_logic = c_logic[3]
        c_Kni_hex = c_Kni_logic.stringify()
                   
        c_Kni_hex_layer = [j for j in Kni_list if c_Kni_hex[6:-2] in Kni_list[j]]
        if c_Kni_hex_layer == [] or len(c_Kni_hex_layer) >= 2:
            print('Kni_FAILED', key)
            
        c_Kni_layer = c_Kni_hex_layer.pop()
            
        for node in edges[key]:
            layer, pi = node
            params = pg.parameter(pi)
            logic = params.logic()
            Kni_logic = logic[3]
            Kni_hex = Kni_logic.stringify()
                        
            Kni_hex_layer = [m for m in Kni_list if Kni_hex[6:-2] in Kni_list[m] ]
            
            if Kni_hex_layer == [] or len(Kni_hex_layer) >= 2:
                print('Kni_FAILED', node)
            
            Kni_layer = Kni_hex_layer.pop()
            
            if Kni_layer != c_Kni_layer:
                if Kni_layer != c_Kni_layer+1:
                    #print('Kni_removed', Kni_layer, c_Kni_layer, key, node)
                    new_edges[key].remove(node) 
    return new_edges


def cond_h2l_Hb_gradient(database, cond_dict, redu_params):
    print('Remember to change phenotype layer numb if needed')
    pg = DSGRN.ParameterGraph(database.network)
    new_cond_dict = deepcopy(cond_dict)
    del_list = []
    for key in cond_dict:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Hb_hex = c_Hb_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_layer == 0:
            if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                continue
                
        if c_layer == 7:       
            if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                continue
                
    new_paramslist = deepcopy(redu_params)
    for i in redu_params[0]:
        c_params = pg.parameter(i[-1])
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Hb_hex = c_Hb_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
            new_paramslist[0].remove(i)
                  
    for i in redu_params[7]:
        c_layer, c_pi = i
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Hb_hex = c_Hb_logic.stringify()

        if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
            new_paramslist[7].remove(i)
        
    final_cond_dict = remove_dict_value(new_cond_dict, del_list)
    
    return final_cond_dict, new_paramslist, del_list

def cond_h2l_Hb_Knibc_gradient(database, cond_dict, redu_params):
    print('Remember to change phenotype layer numb if needed')
    pg = DSGRN.ParameterGraph(database.network)
    new_cond_dict = deepcopy(cond_dict)
    del_list = []
    for key in cond_dict:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_layer == 0:
            if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                
                continue
            if c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)

        if c_layer == 7:       
            if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
                del new_cond_dict[key]
                del_list.append(key)
                    
    new_paramslist = deepcopy(redu_params)
    for i in redu_params[0]:
        c_params = pg.parameter(i[-1])
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]
        c_Kni_logic = c_logic[3]

        c_Hb_hex = c_Hb_logic.stringify()
        c_Kni_hex = c_Kni_logic.stringify()

        #determine if we are in the first or last layer of the edges graph, if so, keep only highest hex for Hb and lowest for Kni
        #if in first layer and keep only edges for Hb's lowest hex and Kni's highest if we are in the last layer

        if c_Hb_hex[6:-2] != 'F'*len(c_Hb_hex[6:-2]):
            new_paramslist[0].remove(i)
            
        elif c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
            new_paramslist[0].remove(i)
            
    for i in redu_params[7]:
        c_layer, c_pi = i
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        c_Hb_logic = c_logic[0]

        c_Hb_hex = c_Hb_logic.stringify()

        if c_Hb_hex[6:-2] != '0'*len(c_Hb_hex[6:-2]):
            new_paramslist[7].remove(i)
                       
    final_cond_dict = remove_dict_value(new_cond_dict, del_list)
    
    return final_cond_dict, new_paramslist, del_list

def add_phenotype_Hb_Knibc_gradient(database, edges, Hb_list, Kni_list):
    new_edges = deepcopy(edges)
    pg = DSGRN.ParameterGraph(database.network)
    del_list = []
    for key in edges:
        c_layer, c_pi = key
        c_params = pg.parameter(c_pi)
        c_logic = c_params.logic()
        
        
        
        if c_layer == 0:
            c_Kni_logic = c_logic[3]
            c_Kni_hex = c_Kni_logic.stringify()
            if c_Kni_hex[6:-2] != '0'*len(c_Kni_hex[6:-2]):
                del new_edges[key]
                del_list.append(key)
                continue
                
        c_Hb_logic = c_logic[0]
        c_Hb_hex = c_Hb_logic.stringify()
            
        c_Hb_hex_layer = [k for k in Hb_list if c_Hb_hex[6:-2] in Hb_list[k] ]
        if c_Hb_hex_layer == [] or len(c_Hb_hex_layer) >= 2:
            print('Hb_FAILED', c_Hb_hex_layer, key)
        
        c_Hb_layer = c_Hb_hex_layer.pop() 
                    
        for node in edges[key]:
            layer, pi = node
            params = pg.parameter(pi)
            logic = params.logic()
            Hb_logic = logic[0]

            Hb_hex = Hb_logic.stringify()
            
            Hb_hex_layer = [l for l in Hb_list if Hb_hex[6:-2] in Hb_list[l]   ]  
            if Hb_hex_layer == [] or len(Hb_hex_layer) >= 2:
                print('Hb_FAILED', Hb_hex_layer, node)
            Hb_layer = Hb_hex_layer.pop()
            
            if Hb_layer != c_Hb_layer:
                if Hb_layer != c_Hb_layer-1:
                    #print('Hb_removed', Hb_layer, c_Hb_layer, key, node)
                    new_edges[key].remove(node)
    final_edges = remove_dict_value(new_edges, del_list)
    return final_edges
              
def remove_dict_value(edges, list_to_remove):
    new_edges = deepcopy(edges)
    for i in edges:
        new_edges[i] = [x for x in edges[i] if x not in list_to_remove]
    return new_edges

def num_edges(edges):
    count = 0
    for key in edges:
        count += len(edges[key])
    return count
