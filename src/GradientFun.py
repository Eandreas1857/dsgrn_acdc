import DSGRN, itertools
from DSGRN import *
from copy import deepcopy
import multiprocessing
from functools import partial
import graphviz
from DSGRN._dsgrn import *
import networkx as nx

def get_gradlist_strict(database, Hb_list, Kni_list):
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)
    #Create a new MGA table of strictly monostable FP's, this is what we will use to find parameter index's
    set_of_MGIm = list(set([ row[0] for row in c.execute('select MorseGraphIndex, Label, count(*) from MorseGraphAnnotations group by MorseGraphIndex') if row[-1]==1 and 'FP' in row[1]] ))
    string = 'select * from Signatures where MorseGraphIndex in ({seq})'.format(
        seq=','.join(['?'] * len(set_of_MGIm)))
    PGIset = [row[0] for row in c.execute(string, set_of_MGIm)]
    
    gradlist = {}
    for i in Hb_list:
        for j in Kni_list:
            gradlist[(i,j)] = [(i,j,x) for x in PGIset if ((((pg.parameter(x)).logic())[0]).stringify())[6:-2] in Hb_list[i] and ((((pg.parameter(x)).logic())[3]).stringify())[6:-2] in Kni_list[j]]
    return gradlist

def get_gradlist(database, Hb_list, Kni_list):
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)
    
    PGIset = [row[0] for row in c.execute('select * from Signatures')]
    
    gradlist = {}       
    for i in Hb_list:
        for j in Kni_list:
            gradlist[(i,j)] = [(i,j,x) for x in PGIset if ((((pg.parameter(x)).logic())[0]).stringify())[6:-2] in Hb_list[i] and ((((pg.parameter(x)).logic())[3]).stringify())[6:-2] in Kni_list[j]]
    return gradlist

def get_gradient_graph_parallel(database, gradlist, num_processes, Hb_list, Kni_list, repeat_layer=True):

    pg = DSGRN.ParameterGraph(database.network)
    todo = list(itertools.chain.from_iterable(gradlist.values()))
    
    pool = multiprocessing.Pool(processes=num_processes)
    work = partial(todo_node_gradient, pg = pg, gradlist = gradlist, Hb_list = Hb_list, Kni_list = Kni_list)
    results = pool.map_async(work, todo)

    edges = {}
    for m,p in results.get():
        edges[m] = p

    return edges

def todo_node_gradient(node, pg, gradlist, Hb_list, Kni_list, repeat_layer = True):
    (Hb, Kni, p) = node             
    next_mg_steps = gradlist[(Hb+1,Kni+1)][:] if Hb < len(Hb_list)-1 and Kni < len(Kni_list)-1  else []
    if repeat_layer:
        next_mg_steps += gradlist[(Hb,Kni)][:]
        if Hb < len(Hb_list)-1:
            next_mg_steps += gradlist[(Hb+1,Kni)][:] 
        if Kni < len(Kni_list)-1:    
            next_mg_steps += gradlist[(Hb,Kni+1)][:]
    # find neighboring parameters using DSGRN adjacencies function to get all possible neighbors
    # accounting for the same parameter in adjacent layers
    adj = list(pg.adjacencies(p, 'codim1'))
    possible_neighbors = [(Hb+1,Kni+1, q) for q in adj + [p]] if abs((Hb+1)-(Kni+1)) < 2 else []
    if repeat_layer:
        possible_neighbors += [(Hb,Kni, q) for q in adj]
        if abs((Hb+1)-(Kni)) < 2:
            possible_neighbors += [(Hb+1,Kni, q) for q in adj]
        if abs((Hb)-(Kni+1)) < 2:
            possible_neighbors += [(Hb,Kni+1, q) for q in adj]
            
    almost = list(set(next_mg_steps).intersection(possible_neighbors))

    return (node, almost)

def remove_dict_value(edges, list_to_remove):
    new_edges = deepcopy(edges)
    for i in edges:
        new_edges[i] = [x for x in edges[i] if x not in list_to_remove]
    return new_edges


def to_dict_of_lists(G, nodelist):
    """Returns adjacency representation of graph as a dictionary of lists.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : order to compute nodes in

    Notes
    -----
    Completely ignores edge data for MultiGraph and MultiDiGraph.

    """
    
    d = {}
    for n in nodelist:
        d[n] = [nbr for nbr in G.neighbors(n)]
    return d
    
def redu_grad_graph(database, grad_graph):
    c = database.conn.cursor()
    pg = DSGRN.ParameterGraph(database.network)

    redu_grad_graph = deepcopy(grad_graph)
    del_list = []

    for node in grad_graph:
        p = node[2]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]
        FP_result = c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))
        FP = FP_result.fetchone()[0]
        #print(p,FP)
        if FP not in FP_list:
            del redu_grad_graph[node]
            del_list.append(node)

    final_grad_graph = remove_dict_value(redu_grad_graph, del_list)
    return final_grad_graph

def Fullconn_gradient_data(database, Paths, scc):
    
    c = database.conn.cursor()
    
    path_layer_nodes = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}

    for path in Paths:
        for layer in range(len(path)):
            path_layer_nodes[layer].append(path[layer])
    
    PG_nodes = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
    for i in path_layer_nodes:
        for j in path_layer_nodes[i]:
            for s in scc:
                for l in scc[s]:
                    if j in l:
                        PG_nodes[i] += [n for n in l]
            
    Hb_thresh = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
    Gt_thresh = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
    Kr_thresh = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}
    Kni_thresh = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: []}

    for item in PG_nodes:
        for node in PG_nodes[item]:
            p = node[-1]
            MGI = list(set([row[1] for row in c.execute('select * from Signatures where ParameterIndex == ' + str(p) )]))
            string = 'select * from MorseGraphAnnotations where MorseGraphIndex in ({seq})'.format(
                seq=','.join(['?'] * len(MGI)))
            FP = [row[-1] for row in c.execute(string, MGI)]
            if len(FP) != 1:
                print('NO FP IN LAYER!',item,node)
            Hb = FP[0][5]
            Hb_thresh[item].append(Hb)
            Gt = FP[0][8]
            Gt_thresh[item].append(Gt)
            Kr = FP[0][11]
            Kr_thresh[item].append(Kr)
            Kni = FP[0][14]
            Kni_thresh[item].append(Kni)
                
    Hb_data = []
    Gt_data = []
    Kr_data = []
    Kni_data = []

    for n in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        A = [round((len([i for i in Hb_thresh[n] if i == str(thresh)])/len(Hb_thresh[n])),2) for thresh in [0,1,2]]
        Hb_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        A = [round((len([i for i in Gt_thresh[n] if i == str(thresh)])/len(Gt_thresh[n])),2) for thresh in [0,1,2]]
        Gt_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        A = [round((len([i for i in Kr_thresh[n] if i == str(thresh)])/len(Kr_thresh[n])),2) for thresh in [0,1,2]]
        Kr_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        A = [round((len([i for i in Kni_thresh[n] if i == str(thresh)])/len(Kni_thresh[n])),2) for thresh in [0,1,2]]
        Kni_data.append(A)
        
    return Hb_data, Gt_data, Kr_data, Kni_data

def GradientPlot(fig_data, title_size, text_size, color = "GnBu", fig_size = (8,8)):    
    path_layer = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    Thresh = [0, 1, 2]
    
    data_T = (np.array(fig_data)).transpose()
    data = np.flip(data_T,0)
    #fig = plt.figure()
    fig, ax = plt.subplots(1,1, figsize=fig_size)

    im = ax.imshow(data, cmap=color)

    # We want to show all ticks...
    ax.set_yticks(np.arange(len(Thresh)))
    ax.set_xticks(np.arange(len(path_layer)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(Thresh, size = text_size)
    ax.set_xticklabels(path_layer, size = text_size)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")

    cbar = ax.figure.colorbar(im, ax=ax, cmap= color, fraction=0.0405)
    cbar.ax.set_ylabel(" ", rotation=-90, va="bottom", size = text_size)

    # Loop over data dimensions and create text annotations.
    for i in range(len(Thresh)):
        for j in range(len(path_layer)):
            text = ax.text(j, i, data[i, j],
                           ha="center", va="center", color="black", size = text_size)

    ax.set_title("Hexcodes in scc's by Phenotype Layer", size=title_size)
    ax.set_ylabel('Hex Poset Layer', size = text_size)
    ax.set_xlabel('Phenotype Pattern Layer', size = text_size)
    #fig.tight_layout()
    #fig.savefig('test.png')
    plt.show()

class GradientGraphviz:
  def __init__(self, database, network, edges, svg_or_png, name):
    """
    Construct graph of the parameter graph for a given network
    """
    self.database = database
    self.network = network
    self.name = name
    self.svg_or_png = svg_or_png
    pg = ParameterGraph(network)
    vertices = []
    for i in edges:
        vertices.append(i)
    
    e = []
    for i in edges:
        for j in edges[i]:
            e.append((i,j))
        
                
    self.vertices = set(vertices)
    self.edges = set(e) 
    self.adjacency_lists = edges
    
    for s,t in self.edges:
      if not s in self.vertices:
        self.vertices.add(s)
      if not t in self.vertices:
        self.vertices.add(t)
      if not s in self.adjacency_lists:
          self.adjacency_lists[s] = []
      self.adjacency_lists[s].append(t) 
    self.vertexname = {}
    for i, v in enumerate(self.vertices):
      self.vertexname[v] = 'X' + str(i)

    

  def _repr_svg_(self):
    """
    Return svg or png
    """
    if self.svg_or_png == 'svg':
        return graphviz.Source(self.graphviz())._repr_svg_()
    if self.svg_or_png == 'png':
        graph = graphviz.Source(self.graphviz(),format='png')
        return graph.render('graph_'+self.name, view = True)

  def graphviz(self):
    """
    Return graphviz string for graph
    """
    return 'digraph {' + \
   '\n'.join([ '"' + self.vertexname[v] + '" [label="(' + str(4-v[0]) + ', ' + str(v[1]) + ', ' + str(v[2]) + ')";style="filled";figsize = .5;fillcolor="' + self.color(v) + '"];' for v in self.vertices ]) + \
   '\n' + '\n'.join([ '"' + self.vertexname[u]  + '" -> "' + self.vertexname[v]  + '";' for (u, v) in self.edges ]) + \
   '\n' + '}\n'


  def color(self, v):
    """
    Return a fillcolor to be used when displaying graph
    """
    if v[0] == 0 and v[1] == 0:
        return "yellow"
    if v[0] == 4 and v[1] == 4:
        return "orange"
    if v[0] == 1 and v[1] == 1:
        return 'deepskyblue2'
    if v[0] == 2 and v[1] == 2:
        return 'darkseagreen'
    if v[0] == 3 and v[1] == 3:
        return 'darksalmon'
    else:
        return "grey"

