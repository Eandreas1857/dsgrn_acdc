from copy import deepcopy
import networkx as nx
import re

def get_number_out_edges_from_string(string, out_edges = {'Hb': 0, 'Gt': 0, 'Kr': 0, 'Kni': 0} ):

    for line in re.findall('\w+\s:.*', string, flags=re.MULTILINE):
        key = re.findall('^\w+', line)[0]
        out_edges[key] = re.findall(key, string, flags=re.MULTILINE).count(key)-1
    
    return out_edges

def convert_FP_list_2_FP_str(lst):
    '''
    lst: list of 4 integerts depicting how many thresholds Hb, Gt, Kr and Kni are above (IN THIS ORDER)
    return: The FP string. I.e., input [0,1,2,1] and return FP { 0, 1, 2, 1} string.
    '''
    return 'FP { ' + str(lst[0]) + ', ' + str(lst[1]) + ', ' + str(lst[2]) + ', ' + str(lst[3]) + ' }'

def get_FP_Poset(out_edges):
    '''
    out_edges: dictinary with network nodes as keys and number of out edges as values.
    returns: the fixed point pattern between layers we need any path to follow.
    '''
    ### Even if not given in order, this way makes sure we have a particular order
    out = []
    out.append(out_edges['Hb'])
    out.append(out_edges['Gt'])
    out.append(out_edges['Kr'])
    out.append(out_edges['Kni'])
     

    ### Convention, L==0, * is any value between 1 and out_edge number -1. Unless the number of out edges is 1, then * is 0 or 1. H == number of out edges 
    Phenotype_pattern = {1: '*HLL', 2: 'H*LL', 3: 'HL*L', 4: '*LHL', 5: 'LLH*', 6: 'LL*H', 7: 'L*LH', 8:'LHL*'}
    FP_Poset = {}
    FP_edges = []
    Region_start_FPs = []
    for r in list(range(1,8)):
        if r == 1:
            region_FP = []
            R=[]
            flag = 5
            for i in range(4):
                if Phenotype_pattern[r][i] == '*':
                    if out[i]>1:
                        region_FP.append(1)
                    if out[i] == 1:
                        region_FP.append(0)
                if Phenotype_pattern[r][i] == 'H':
                    region_FP.append(out[i])
                if Phenotype_pattern[r][i] == 'L':
                    region_FP.append(0)
        else:
            flag = 5
            R=[]
            region_FP = []
            for i in range(4):
                if Phenotype_pattern[r][i] == 'H':
                    region_FP.append(out[i])
                elif out[i] > 1 and Phenotype_pattern[r][i] == '*':
                    if Phenotype_pattern[r-1][i] == 'H':
                        region_FP.append(out[i] - 1)
                    if Phenotype_pattern[r-1][i] == 'L':
                        region_FP.append(1)
                    if out[i]>2 and Phenotype_pattern[r-1][i] == 'H':
                        flag = i
                else:
                    region_FP.append(0)
        FP_edges.append(region_FP)

        R += [convert_FP_list_2_FP_str(region_FP)]
        if flag < 5:
            for j in range(1,out[flag]-1):
                S = region_FP.copy()
                S[flag] = S[flag]-1
                R += [convert_FP_list_2_FP_str(S)]

        Region_start_FPs.append(list(set(R)))

        if r == 7:
            region_FP = []
            R = []
            flag = 5
            for i in range(4):
                if Phenotype_pattern[r+1][i] == 'H':
                    region_FP.append(out[i])
                elif out[i] > 1 and Phenotype_pattern[r+1][i] == '*':
                    if Phenotype_pattern[r][i] == 'H':
                        region_FP.append(out[i] - 1)
                    if Phenotype_pattern[r][i] == 'L':
                        region_FP.append(1)
                    if out[i]>2 and Phenotype_pattern[r][i] == 'H':
                        flag = i
                else:
                    region_FP.append(0)

            R += [convert_FP_list_2_FP_str(region_FP)]
            if flag < 5:
                for j in range(1,out[flag]-1):
                    S = region_FP.copy()
                    S[flag] = S[flag]-1
                    R += [convert_FP_list_2_FP_str(S)]

            Region_start_FPs.append(list(set(R)))      

        next_region_FP = []
        for i in range(4):
            if Phenotype_pattern[r+1][i] == 'H':
                next_region_FP.append(out[i])
            elif out[i] > 1 and Phenotype_pattern[r+1][i] == '*':
                next_region_FP.append(out[i]-1)
            else:
                next_region_FP.append(0)

        todo =  deepcopy(FP_edges)
        while todo:
            FP = todo.pop()
            if convert_FP_list_2_FP_str(FP) not in FP_Poset:
                FP_Poset[convert_FP_list_2_FP_str(FP)] = [] 
            FP_edges = []
            for i in range(4):
                next_FP = FP.copy() 
                if Phenotype_pattern[r+1][i] != Phenotype_pattern[r][i]:

                    if Phenotype_pattern[r][i] == '*':

                        if Phenotype_pattern[r+1][i] == 'H':
                            if next_FP[i] < out[i]:
                                next_FP[i] += 1
                                if next_FP[i] != next_region_FP:
                                    FP_edges.append(next_FP.copy())
                                    todo.append(next_FP.copy())
                                else:
                                    break
                                
                        if Phenotype_pattern[r+1][i] == 'L':
                            if next_FP[i] > 0:
                                next_FP[i] -= 1
                                if next_FP[i] != next_region_FP:
                                    FP_edges.append(next_FP.copy())
                                    todo.append(next_FP.copy())
                                else:
                                    break

                    if Phenotype_pattern[r+1][i] == '*':

                        if Phenotype_pattern[r][i] == 'H' and out[i] != 1:
                            if next_FP[i] > 1:
                                next_FP[i] -= 1
                                if next_FP[i] != next_region_FP:
                                    FP_edges.append(next_FP.copy())
                                    todo.append(next_FP.copy())
                                else:
                                    break
                        
                        if Phenotype_pattern[r][i] == 'H' and out[i] == 1:
                            if next_FP[i] > 0:
                                next_FP[i] -= 1
                                if next_FP[i] != next_region_FP:
                                    FP_edges.append(next_FP.copy())
                                    todo.append(next_FP.copy())
                                else:
                                    break
                                
                        if Phenotype_pattern[r][i] == 'L' and out[i] != 1:
                            if next_FP[i] < out[i]-1:
                                next_FP[i] += 1
                                if next_FP[i] != next_region_FP:
                                    FP_edges.append(next_FP.copy())
                                    todo.append(next_FP.copy())
                                else:
                                    break
                        
            for f in FP_edges:
                if convert_FP_list_2_FP_str(f) not in FP_Poset[convert_FP_list_2_FP_str(FP)]:
                    FP_Poset[convert_FP_list_2_FP_str(FP)].append(convert_FP_list_2_FP_str(f))

    H = nx.DiGraph() #building networkx graph from FP_poset

    for node in FP_Poset:
        for edge in FP_Poset[node]:
            H.add_edge(node, edge)

    FP_Regions = {}
    for i in range(len(Region_start_FPs)-1):
        FP_Regions[i+1] = []
        for j in Region_start_FPs[i]:
            for k in Region_start_FPs[i+1]:
                paths_between_generator = nx.all_simple_paths(H, j, k)
                nodes_between_set = {node for path in paths_between_generator for node in path[:-1] if node not in Region_start_FPs[i+1]}
                FP_Regions[i+1] += list(sorted(nodes_between_set))
    
    for lst in FP_Regions:
        FP_Regions[lst] = list(set(FP_Regions[lst]))
    
    FP_Regions[8] = Region_start_FPs[-1]

    return FP_Poset, FP_Regions