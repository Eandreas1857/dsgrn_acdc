import graphviz
from DSGRN._dsgrn import *

class PhenotypeGraphviz:
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
    return 'graph {' + \
   '\n'.join([ '"' + self.vertexname[v] + '" [label="' + str(v) + '";style="filled";figsize = .5;fillcolor="' + self.color(v) + '"];' for v in self.vertices ]) + \
   '\n' + '\n'.join([ '"' + self.vertexname[u]  + '" -- "' + self.vertexname[v]  + '";' for (u, v) in self.edges ]) + \
   '\n' + '}\n'


  def color(self, v):
    """
    Return a fillcolor to be used when displaying graph
    """
    if v[0] == 0:
        return "cyan3"
    if v[0] == 1:
        return 'cyan4'
    if v[0] == 2:
        return 'darkseagreen'
    if v[0] == 3:
        return 'cornflowerblue'
    if v[0] == 4:
        return 'cadetblue1'
    if v[0] == 5:
        return 'deepskyblue2'
    if v[0] == 6:
        return 'goldenrod1'
    if v[0] == 7:
        return 'goldenrod'
    if v[0] == 8:
        return 'darkorange'
    if v[0] == 8:
        return 'darksalmon'
    else:
        return "grey"
