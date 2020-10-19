import DSGRN, itertools
from DSGRN import *
from NFixedPointQuery import NFixedPointQuery

from copy import deepcopy



def flatten(list2flatten):
    '''
    The following takes a list of lists and turns it into a single list.
    :param listtoflatten: a list of lists
    :return: a flat list containing all the elements of the component lists
    '''
    return list(itertools.chain(*list2flatten))

def get_paramslist(database, list_of_bounds, goe):
    
    paramslist = []
    c = database.conn.cursor()
    for i in range(len(list_of_bounds)):
        bounds1=list_of_bounds[i]
        NFP = NFixedPointQuery(database, *bounds1).matches()
        N = len(bounds1)
        X1 = NstableQuery(database, N).matches()
        X2 = NstableQuery(database, N+1).matches()
        inter = set(X1.difference(X2))
        if goe == '=':
            MGset = [j for j in inter if j in NFP]
        else:
            MGset = list(NFP)

        string = 'create temp table C as select * from Signatures where MorseGraphIndex in ({seq})'.format(
        seq=','.join(['?']*len(MGset)))
        c.execute(string, MGset)
        PGI1 = [ row[0] for row in c.execute('select ParameterIndex from C')]
        c.execute('drop table C')
        if PGI1 !=[]:
            p = deepcopy(paramslist)
            flat_p = flatten(p)
            PGI1 = [(i,n) for n in PGI1 if n not in flat_p]
            paramslist.append(PGI1.copy())
        else:
            
            print('EMPTY LIST', list_of_bounds[i])
            break

       
    return paramslist



def get_phenotype_graph(database,paramslist,repeat_layer=True):
    '''
    Perform successive intersections of parameter neighbors with the parameter indices for the same bounds and for the next bounds. This results in a list of pairs of neighboring parameters with the desired bounds properties. Allows repeats of bound matches except for the last list. Does not allow repeats of parameter indices. The result is saved as a graph, making the end result amenable to graph searches.
    This version of the code makes sure that repeated parameters are not allowed to have backwards edges in the phenotype graph (i.e., the only allowed loops are within layers).
    :param paramslist: A list of lists of DSGRN parameter indices. The lists must be in the same order as the desired bounds.
    :param database: DSGRN.Database object
    :return: Dictionary with DSGRN parameter indices (nodes) keying lists of DSGRN parameter indices (out-edges). In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    ''' 
    pg = DSGRN.ParameterGraph(database.network)

    # name nodes according to layer and dsgrn parameter
    todo = flatten(paramslist)
    edges = dict.fromkeys(todo,[])
    
    next_steps_dict = {}
    for k in range(len(paramslist)):
        if k+2<=len(paramslist):
            next_steps_dict[k] = [q for q in paramslist[k+1]]
        else:
            next_steps_dict[k] = [q for q in paramslist[k]]
        if repeat_layer:
            next_steps_dict[k] += [q for q in paramslist[k]]
            
    while todo:
        (k,p) = todo.pop(0)
        # record all allowable steps in the mg layers
        next_mg_steps = [i for i in next_steps_dict[k] if i[-1]!=p]
        # find neighboring parameters using DSGRN adjacencies function to get all possible neighbors
        # accounting for the same parameter in adjacent layers
        
        adj = list(pg.adjacencies(p,'codim1'))
        possible_neighbors = [(k+1,q) for q in adj + [p]] + [(k,q) for q in adj]

        # intersect neighbors with parameters that have allowable MGs
        edges[(k,p)] = list(set(next_mg_steps).intersection(possible_neighbors))

    return edges

def get_paramslist_optimized(database, list_of_bounds, goe):
    paramslist = []
    c = database.conn.cursor()
    for i in range(len(list_of_bounds)):
        bounds1 = list_of_bounds[i]
        NFP = NFixedPointQuery(database, *bounds1).matches()
        if goe == '=':
            N = len(bounds1)
            X1 = NstableQuery(database, N).matches()
            X2 = NstableQuery(database, N + 1).matches()
            inter = set(X1.difference(X2))
            MGset = [j for j in inter if j in NFP]
        else:
            MGset = list(NFP)

        string = 'create temp table C as select * from Signatures where MorseGraphIndex in ({seq})'.format(
            seq=','.join(['?'] * len(MGset)))
        c.execute(string, MGset)
        PGI1 = [row[0] for row in c.execute('select ParameterIndex from C')]
        c.execute('drop table C')
        if PGI1:
            paramslist.append([(i, n) for n in PGI1])
        else:
            print('EMPTY LIST', list_of_bounds[i])
            paramslist = None
            break
    return paramslist

def get_phenotype_graph_optimized(database, paramslist, repeat_layer=True):
    '''
    Perform successive intersections of allowable parameter transitions with codimension 1 parameter adjacencies. This results in a list of pairs of neighboring parameters with the desired bounds properties. Allows repeats of bound matches except for the last list. The result is saved as a graph, making the end result amenable to graph searches.
    :param database: DSGRN.Database object
    :param paramslist: A list of lists of (layer number, DSGRN parameter index) pairs.
    :param repeat_layer: allow repeated bounds matches except at the last layer.
    :return: Dictionary with (layer number, DSGRN parameter index) pairs keying lists of (layer number, DSGRN parameter index) pairs that represents the phenotype graph. In other words, each (key, list_element) pair is an edge in a graph. Every edge satisfies the correct bounds relationship indicated by the order of the parameter lists.
    '''

    pg = DSGRN.ParameterGraph(database.network)
    todo = list(itertools.chain.from_iterable(paramslist))
    edges = dict.fromkeys(todo, [])

    while todo:
        (k, p) = todo.pop(0)
        # record all allowable steps in the mg layers
        next_mg_steps = paramslist[k+1][:] if k < len(paramslist)-1 else []
        if repeat_layer:
            next_mg_steps += paramslist[k][:]
        # find neighboring parameters using DSGRN adjacencies function to get all possible neighbors
        # accounting for the same parameter in adjacent layers
        adj = list(pg.adjacencies(p, 'codim1'))
        possible_neighbors = [(k + 1, q) for q in adj + [p]]
        if repeat_layer:
            possible_neighbors += [(k, q) for q in adj]
        # intersect neighbors with parameters that have allowable MGs
        edges[(k, p)] = list(set(next_mg_steps).intersection(possible_neighbors))
    return edges

def reduce_graph(edges):
    for i in edges:
        for j in edges:
            if i in edges[j] and j in edges[i]:
                if len(edges[i]) > len(edges[j]):
                    edges[j].remove(i)
                else:
                    edges[i].remove(j)
    return edges


def find_a_path(edges,start_set,stop_set, path_length):
    def find_paths(path, path_length):
        
        '''
        Recursive function that finds all paths from a starting parameter to any ending parameter in a set. This only works on a connected graph and does not repeat nodes in any path.
        :param path: A list containing the parameter node at which to start.
        :return: Side-effect list.
        '''
        if path[-1] in stop_set:
            paths.append(path)
        else:
            for n in comp_edges[path[-1]]:
                if paths == []:
                    if n not in path and len(path) <= path_length:
                        find_paths(path + [n], path_length)
                else:
                    break
        return None

    paths = []
    for start in start_set:
        if paths == []:
            comp_edges = get_connected_component(start,edges)
            find_paths([start], path_length)
        else:
            break
    return paths

def find_all_paths(edges,start_set,stop_set):

    def find_paths(path):
        '''
        Recursive function that finds all paths from a starting parameter to any ending parameter in a set. This only works on a connected graph and does not repeat nodes in any path.
        :param path: A list containing the parameter node at which to start.
        :return: Side-effect list.
        '''
        if path[-1] in stop_set:
            paths.append(path)
        else:
            for n in comp_edges[path[-1]]:
                if n not in path:
                    find_paths(path + [n])
        return None

    paths = []
    for start in start_set:
        comp_edges = get_connected_component(start,edges)
        find_paths([start])
    return paths




def get_connected_component(node,edges):
    '''
    Find connected component containing node, and return new edges with only that component.
    :param node: node name
    :param edges: dictionary with parameters keying list of out-edge parameters
    :return: edges dictionary minus the nodes
    '''
    undirected = deepcopy(edges)
    for n,es in edges.items():
        for i in es:
            undirected[i].append(n)
    allowed = set([node])
    while allowed.intersection(list(undirected.keys())):
        a = []
        for n in allowed:
            if n in undirected:
                a.extend(undirected[n])
                undirected.pop(n)
        allowed = set(a)
    keepers = deepcopy(edges)
    for n in undirected.keys():
        keepers.pop(n)
    
    return keepers
