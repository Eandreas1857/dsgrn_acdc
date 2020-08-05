# PGDraw
# Elizabeth Andreas
# 2020-08-05
# MIT LICENSE
# Edited from Graph.py by Shaun Harker

import graphviz

from DSGRN._dsgrn import *

class PGGraph:
  def __init__(self, database, network):
    """
    Construct graph of the parameter graph for a given network
    """
    self.database = database
    self.network = network
    pg = ParameterGraph(network)
    vertices = set(range(0,pg.size()))
    edges = set()
    for gpi1 in vertices:
        for gpi2 in pg.adjacencies(gpi1):
            if (gpi2,gpi1) not in edges:
                edges = edges.union({(gpi1, gpi2)})
            
    self.vertices = set(vertices)
    self.edges = set(edges) 
    self.adjacency_lists = {}
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
    Return svg representation suitable for Notebook visualization
    """
    return graphviz.Source(self.graphviz())._repr_svg_()

  def graphviz(self):
    """
    Return graphviz string for graph
    """
    return 'graph {' + \
   '\n'.join([ '"' + self.vertexname[v] + '" [label="' + str(self.label(v)) + '";style="filled";fillcolor="' + self.color(v) + '"];' for v in self.vertices ]) + \
   '\n' + '\n'.join([ '"' + self.vertexname[u]  + '" -- "' + self.vertexname[v]  + '";' for (u, v) in self.edges ]) + \
   '\n' + '}\n'
   
  def label(self, v):
    """
    Return a label string to be used when displaying graph
    """
    c = self.database.conn.cursor()
    P = c.execute('select * from Signatures')
    for i in P:
        if v == i[0]:
            s = i[1]
    return (v,s)

  def color(self, v):
    """
    Return a fillcolor to be used when displaying graph
    """
    return "grey"
