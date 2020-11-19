
from funforMGpaths import *

class MGsearchthroughPG:
    def __init__(self, database, goe1, goe2, bounds_set):
        
        paths_set = setofpaths(database, goe1, goe2, bounds_set)
        self.bounds_set = bounds_set

        for i in range(len(bounds_set)-2):
            temp = final(paths_set)
            paths_set = temp.copy()
            
        self.all_paths = paths_set
      
    def allpaths(self):
        if len(self.bounds_set) == 2:
            return flatten(self.all_paths)
        else:
            return self.all_paths

    def shortestpaths(self):
        shortest_paths = []
        for path in self.all_paths:
            if len(path) == min([len(j) for j in self.all_paths]):
                shortest_paths.append(path)
        return shortest_paths

    def firstandlast(self):
        if len(self.bounds_set) == 2:
            X = [[item[0],item[-1]] for item in flatten(self.all_paths)]
            return [list(i) for i in set(map(tuple, X))]
        else:
            X = [[item[0],item[-1]] for item in self.all_paths]
            return [list(i) for i in set(map(tuple, X))]

    def TFP_test(self):
        #For TFP_test, MG 2 -> 1 -> 0
        Correct = [[10,3,2], [10,3,4], [10,3,1,0], [10,3,1,8], [10,3,5,6], [10,3,5,12], [11,10,3,2],[11,10,3,4], [11,10,3,1,0], [11,10,3,1,8], [11,10,3,5,6], [11,10,3,5,12], [9,10,3,2],[9,10,3,4], [9,10,3,1,0], [9,10,3,1,8], [9,10,3,5,6], [9,10,3,5,12]]
        
        for i in Correct:
            for j in self.all_paths:
                if i == j:
                    self.all_paths.remove(i)
        return self.all_paths
    
    def TFPYI_test1(self):
        #For TFPYI_test MG 0 -> 5
        Correct = []
        
        assert(flatten(self.all_paths) == Correct)


    def TFPYI_test2(self):
        #For TFPYI_test, MG 7 -> 2
        Correct = [[32,25], [31,32,25], [30,31,32,25]]
        F = flatten(self.all_paths)
        for i in Correct:
            for j in F:
                if i == j:
                    F.remove(i)
        return F


