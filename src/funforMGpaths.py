from NsearchgoeQuery import *
import itertools

# following removes the MGI from the output of NsearchgoeQuery
def removeMGI(list2clean):
    '''
    This function removes the MGI from the output of NsearchgoeQuery.

    :param listtoclean: a list of tuples or lists
    :return: a list of sublists of listtoclean
    '''
    return [list(item[:2]) for item in list2clean]

# following takes a list of lists and turns it into a single list.
def flatten(list2flatten):
    '''
    The following takes a list of lists and turns it into a single list.

    :param listtoflatten: a list of lists
    :return: a flat list containing all the elements of the component lists
    '''
    return list(itertools.chain(*list2flatten))

# this function takes a list of lists and finds paths between the first two sublists where the PGI do not repeat.
def listpath(startlist):
    '''
    This function finds paths in the parameter graph between PGI
    (parameter graph indices [listed in the input argument]
    such that the PGI do not repeat.
    Important: The assumption that the PGI do not repeat means that during the course of the network behavior,
    the system never returns to the same parameterization, even if it would be appropriate to match.

    :param X: a list of lists, but only evaluates first two sublists, of parameter graph indices
    :return: a list of lists of parameter graph indices that represent paths through the parameter graph between the last
    element in the first list and the first element in the second list as long as 
    all of the other elements in the second list are not in the first list
    '''
    endlist = []
    for i in startlist[1]:
        for j in startlist[0]:
            n = len(j)
            if j[n-1] == i[0]:
                if len(set(i[1:]).intersection(set(j))) == 0:
                    E = [j.copy(),i.copy()]
                    E[0].pop()
                    final = []
                    for p in E:
                        for q in p:
                            final.append(q) 
                    endlist.append(final.copy())
    return endlist


# following function finds paths between two lists where list2 = [int , int]. This is used to find paths while only searching between two bds. This is a special case of the previous function.
def twolistpath(list1, list2, M):
    '''
    This function is a special case of the previous function "listpath" where the second list is limited to two elements.
    This is used to find paths while only searching between two bds fixed point bounds.

    :param list1: list of lists of parameter graph indices
    :param list2: list of exactly two parameter graph indices
    :param M: empty list
    :return: list of paths between list1 and list2
    '''
    if list1 == []:
        return M
    else:
        for sublist in list2:
            x = list(list1[-1])
            if x[-1] == sublist[0]:
                if sublist[1] not in x:
                    x.append(sublist[1])
                    M.append(x)
        list1.pop()
        twolistpath(list1, list2, M)
 
        return M
    return M

# following creates the lists between bounds where MG's can be repeated. For example, if bds1 is related to MGI 0 and bds2 is related to MGI 1, this will find all of the paths in the PG that go from 0 to 1 where 0 can be repeated as many times as needes but 1 cannot. Additionally, the paths cannot go through a PG node more than once. Returns lists of PGI. 

def pathsbw2bds(database, goe1, goe2, bds1, bds2):
    '''
    This function creates the lists between bounds where MG's (Morse graphs) can be repeated.
    For example, if bds1 is related to MGI 0 and bds2 is related to MGI 1, this will find all of the paths in the
    PG (parameter graph) that go from 0 to 1 where 0 can be repeated as many times as needed but 1 cannot.
    Additionally, the paths cannot go through a PG node more than once.
    Important: The assumption that the PGI do not repeat means that during the course of the network behavior,
    the system never returns to the same parameterization, even if it matches bds.

    :param database: DSGRN Database object
    :param goe1: '<' or '=' where '<' indicates a search where bds1 in the list of fixed points assosiated to the MG and '=' indicates a search where bds1 are the only fixed points in the MG.
    :param goe2: '<' or '=' where '<' indicates a search where bds2 in the list of fixed points assosiated to the MG and '=' indicates a search where bds2 are the only fixed points in the MG.
    :param bds1: list of lists where each list indicates a fixed point bound (i.e., bds1 = [bdsA, bdsB] where we might want bdsA = {"X":0, "Y":0}, indicating FP{0,0} in the MG and bds2 = {"X":2, "Y":[1,0]}, indicating FP{2,1} or FP{2,0} in the MG.
    :param bds2: list of lists where each list indicates a fixed point bound, throws an error is bds1 and bds2 overlap.
    :return: lists of parameter graph indices representing paths between fixed points in bds1 and bds2
    '''
    # finds edges going from bds1 to bds2, returns list of lists of len = 2
    diff = list(NsearchgoeQuery(database, goe1, goe2, bds1, bds2).stability_of_matches())
    diffsansMGI = removeMGI(diff)

    # finds edges going from bds1 to bds1, returns list of lists of len = 2
    same = list(NsearchgoeQuery(database, goe1, goe2, bds1, bds1).stability_of_matches())
    samesansMGI = removeMGI(same)
    
    # this next bit gives all of the paths where bds1 to bds1 to bds1 to ... without repeating nodes
    M = []
    copy = samesansMGI.copy()
    twolistpath(copy, same, M)
    paths_of_bds1 = []
    while M != []:
        list1 = M.copy()
        paths_of_bds1.append(M)
        M = []
        twolistpath(list1, same, M) 

    paths_of_bds1.append(samesansMGI)
    paths_of_bds1_flat = flatten(paths_of_bds1)
    
    N = []
    paths_to_bds2 = []
    
    twolistpath(paths_of_bds1_flat, diffsansMGI, N)

    paths_to_bds2.append(N)
    paths_to_bds2.append(diffsansMGI)

    final_list_of_paths = flatten(paths_to_bds2)
    return final_list_of_paths

# next function creates a list of lists that contain paths from bds1 to bds2, bds2 to bds3,.. and so on for all bds given as an input.

def setofpaths(database, goe1, goe2, bounds_set):
    '''
    This function creates a list of lists that contain paths from bds1 to bds2, bds2 to bds3,.. and so on for all bds given as an input.
    :param database: DSGRN Database object
    :param goe1: '<' or '=' where '<' indicates a search where bds1 in the list of fixed points assosiated to the MG and '=' indicates a search where bds1 are the only fixed points in the MG.
    :param goe2: '<' or '=' where '<' indicates a search where bds2 in the list of fixed points assosiated to the MG and '=' indicates a search where bds2 are the only fixed points in the MG.
    :param bounds_set: list of lists containg a set of fixed point bounds wanting to find a path through. Bounds should not overlap (the same fixed point should not be in two different bound lists).
    :return:
    '''
    set_of_paths = []
    for i in range(len(bounds_set)):
        if i+1 != len(bounds_set):
            Y = pathsbw2bds(database, goe1, goe2, bounds_set[i], bounds_set[i+1])
            set_of_paths.append(Y)
    return set_of_paths

# this next function finds a path from bds1 -> bds2 -> bds3 and returns all paths where PG nodes do not repeat. If there are more than 3 bds, it returns a list or lists where the first list is the paths from bds1 -> bds3 and the rest of the lists are the paths from bds3 -> bds4, bds4 -> bds5, ... and so on. This is designed to be put into a recursive function so that eventually we have a single list of paths from the first set of bds to the last. Currently the recursion is not in this function because I cannot think of a decent stopping procedure inside the function but can outside of it. So I will put the recursion into the module instead.

def final(set_of_paths):
    '''
    this function finds a path from bds1 -> bds2 -> bds3 and returns all paths where PG nodes do not repeat. If there are more than 3 bds, it returns a list or lists where the first list is the paths from bds1 -> bds3 and the rest of the lists are the paths from bds3 -> bds4, bds4 -> bds5, ... and so on. This is designed to be put into a recursive function so that eventually we have a single list of paths from the first set of bds to the last. Currently the recursion is not in this function because I cannot think of a decent stopping procedure inside the function but can outside of it. So I will put the recursion into the module instead.
    :param database: DSGRN Database object
    :param goe1: '<' or '=' where '<' indicates a search where bds1 in the list of fixed points assosiated to the MG and '=' indicates a search where bds1 are the only fixed points in the MG.
    :param goe2: '<' or '=' where '<' indicates a search where bds2 in the list of fixed points assosiated to the MG and '=' indicates a search where bds2 are the only fixed points in the MG.
    :param set_of_paths: list of lists containg all of the paths from bds1 -> bds2, bds2 -> bds3, ... (all lists are parameter graph indices)
    :return: returns a list of lists of paths where the first sublist are the paths (retuned as parameter graph indices) from bds1 -> bds3 and the rest of the sublists are the remaining sublists from set_of_paths.
    '''
    if len(set_of_paths) == 2:
        M = listpath(set_of_paths)
        return M
    
    if len(set_of_paths) > 2:
        M = listpath(set_of_paths)
        F = [M,*set_of_paths[2:]]
        return F



