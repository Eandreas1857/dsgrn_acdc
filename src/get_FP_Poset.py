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

def get_region_head(Phenotype_pattern, num_out_thresh):

    region_head = {}
    for r in Phenotype_pattern:
        region_head[r] = []
        position = {0:[], 1:[], 2:[], 3:[]}
        
        for i in range(4):
            if Phenotype_pattern[r][i] == '*':
                if num_out_thresh[i]>1:
                    position[i] += [n for n in range(1,num_out_thresh[i])]
                if num_out_thresh[i] == 1:
                    position[i] = [0,1]
            if Phenotype_pattern[r][i] == 'H':
                position[i] += [num_out_thresh[i]]
            if Phenotype_pattern[r][i] == 'L':
                position[i] += [0]

        region_FP = [(a,b,c,d) for a in position[0] for b in position[1] for c in position[2] for d in position[3] ]
        for FP in region_FP:
            #region_head[r] += [convert_FP_list_2_FP_str(FP)]
            region_head[r] += [FP]
    return region_head

def get_FP_Poset(out_edges):
    '''
    out_edges: dictinary with network nodes as keys and number of out edges as values.
    returns: the fixed point pattern between layers we need any path to follow.
    '''
    ### Even if not given in order, this way makes sure we have a particular order
    Phenotype_pattern = {1: '*HLL', 2: 'H*LL', 3: 'HL*L', 4: '*LHL', 5: 'LLH*', 6: 'LL*H', 7: 'L*LH', 8:'LHL*'}
    num_out_thresh = []
    num_out_thresh.append(out_edges['Hb'])
    num_out_thresh.append(out_edges['Gt'])
    num_out_thresh.append(out_edges['Kr'])
    num_out_thresh.append(out_edges['Kni'])

    region_head = get_region_head(Phenotype_pattern, num_out_thresh)
    FP_Regions = {}
    FP_Poset_graph = nx.DiGraph()
    for r in range(1,8):
        FP_Regions[r] = [convert_FP_list_2_FP_str(lst) for lst in region_head[r].copy()]

        temp = nx.DiGraph()
        for n in region_head[r]:
            temp.add_node(n)
        for n in region_head[r+1]:
            temp.add_node(n)

        position = {0:[], 1:[], 2:[], 3:[]}
        inc_dec = {}

        for i in range(4):
            if Phenotype_pattern[r][i] == '*':
                if Phenotype_pattern[r+1][i] == 'H':
                    inc_dec[i] = 'inc'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(1,num_out_thresh[i]+1)]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
                if Phenotype_pattern[r+1][i] == 'L':
                    inc_dec[i] = 'dec'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(0,num_out_thresh[i])]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]

            if Phenotype_pattern[r+1][i] == '*':
                if Phenotype_pattern[r][i] == 'H':
                    inc_dec[i] = 'dec'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(1,num_out_thresh[i]+1)]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
                if Phenotype_pattern[r][i] == 'L':
                    inc_dec[i] = 'inc'
                    if num_out_thresh[i]>1:
                        position[i] += [n for n in range(0,num_out_thresh[i])]
                    if num_out_thresh[i] == 1:
                        position[i] = [0,1]
            if Phenotype_pattern[r][i] == Phenotype_pattern[r+1][i]:
                inc_dec[i] = True
                position[i] += [region_head[r][0][i]]
        
        region_T = [(a,b,c,d) for a in position[0] for b in position[1] for c in position[2] for d in position[3]]  
        
        for t in region_T:
            temp.add_node(t)
        #print(region_T)
        #print(inc_dec)
        FP_Regions[r] += [convert_FP_list_2_FP_str(i) for i in region_T.copy() if i not in region_head[r+1] and i not in region_head[r]]
        for s in temp:
            for t in temp:
                if s != t:
                    edge = True
                    for i in range(4):
                        if inc_dec[i] == True and s[i] != t[i]:
                            edge = False
                        if inc_dec[i] == 'inc' and s[i] > t[i]:
                            edge = False
                        if inc_dec[i] == 'dec' and s[i] < t[i]:
                            edge = False
                    if edge == True:
                        temp.add_edge(s,t)
                        #print(s,t)
        #print(temp.edges())
        for e in temp.edges():
            FP_Poset_graph.add_edge(e[0], e[1])
    mapping = {}
    for node in FP_Poset_graph:
        mapping[node] = convert_FP_list_2_FP_str(node)

    FP_Poset_graph = nx.relabel_nodes(FP_Poset_graph, mapping)

    FP_Poset = nx.to_dict_of_lists(FP_Poset_graph, nodelist=FP_Poset_graph.nodes())
    FP_Regions[8] = [convert_FP_list_2_FP_str(lst) for lst in region_head[8].copy()]

    return FP_Poset, FP_Regions