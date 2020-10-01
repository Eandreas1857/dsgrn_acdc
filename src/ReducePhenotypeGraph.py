from collections import defaultdict
from PhenotypeGraphFun import *
from copy import deepcopy

class ReducePhenotypeGraph:

    def __init__(self, vertex, edges):
        self.V = vertex
        self.graph = defaultdict(list)
        self.edges = edges
        
    # Add edge into the graph
    def add_edge(self, s, d):
        self.graph[s].append(d)

        
    def add_edges(self, edges):
        for i in edges:
            for j in edges[i]:
                self.graph[i].append(j)


    # dfs
    def dfs(self, d, visited_vertex, final):
        visited_vertex[d] = True
        
        final.append(d)
        for i in self.graph[d]:
            if not visited_vertex[i]:
                self.dfs(i, visited_vertex, final)

    def fill_order(self, d, visited_vertex, stack):
        visited_vertex[d] = True
        for i in self.graph[d]:
            
            if not visited_vertex[i]:
                self.fill_order(i, visited_vertex, stack)
        stack = stack.append(d)


    # transpose the matrix
    def transpose(self):
        g = ReducePhenotypeGraph(self.V, self.edges)

        for i in self.graph:
            for j in self.graph[i]:
                g.add_edge(j, i)
        return g

    # Print stongly connected components
    def scc(self):
        stack = []
        
        visited_vertex = deepcopy(self.edges)
        
        for n in visited_vertex:
            visited_vertex[n]=False
            
        all_connected = []
        
        for i in visited_vertex:
           
            if visited_vertex[i] == False:
                
                self.fill_order(i, visited_vertex, stack)

        gr = self.transpose()
        
        visited_vertex = deepcopy(self.edges)
        
        for n in visited_vertex:
            visited_vertex[n]=False
        

        while stack:
            i = stack.pop()
            
            if not visited_vertex[i]:
                connected = []
                gr.dfs(i, visited_vertex, connected)
                all_connected.append(connected)
        return(all_connected)

def reduce(edges, paramslist):
    reduced_edges = deepcopy(edges)
    g = ReducePhenotypeGraph(max(edges)[-1]+1, edges)
    g.add_edges(edges)
    
    scc = g.scc()
    
    stronglycc = [i for i in scc if len(i) != 1]
    
    for grps in stronglycc:
        for vert in edges:
            if vert in grps:
                if edges[grps[0]] != edges[vert]:
                    reduced_edges[grps[0]] += [n for n in edges[vert] if n not in edges[grps[0]] ]
                for j in range(len(grps)-1):
                    if (grps)[j+1] in reduced_edges:
                        reduced_edges.pop((grps)[j+1])

    for vert in reduced_edges:
        for k in range(len(reduced_edges[vert])):
            for grps in stronglycc:
                if reduced_edges[vert][k] in grps:
                    reduced_edges[vert][k] = grps[0]
        reduced_edges[vert] = list(dict.fromkeys(reduced_edges[vert]))
        if vert in reduced_edges[vert]:
            reduced_edges[vert].remove(vert)
    
    reduced_paramslist = deepcopy(paramslist)
    for grps in scc:
        for level in range(len(paramslist)):
            for item in paramslist[level]:
                if item in grps and item != grps[0]:
                    reduced_paramslist[level].remove(item)

    return reduced_edges, reduced_paramslist, stronglycc
