import matplotlib.pyplot as plt
import networkx as nx
import scipy.linalg
import numpy as np
from numpy import linalg
from sklearn import preprocessing

def asym_weight_matrix(G, nodelist = None, data = None):
    """Returns the affinity matrix (usually denoted as A, but denoted as W here) 
    of the DiGraph G. W is an a asymmetric square matrix with non-negative real 
    numbers. W(i,j) represents the weight of the directed edge i --> j. [1]

    Parameters
    ----------
    G : NetworkX DiGraph

    nodelist : collection
        A collection of nodes in `G`. If not specified, nodes are all nodes in G
        in G.nodes() order.

    data : object
        Edge attribute key to use as weight. If not specified, edges
        have weight one.

    Returns
    -------
    array
        The affinity matrix of G.
        

    References
    ----------
    .. [1] Marina Meila and William Pentney.
           *Clustering by weighted cuts in directed graphs*.
           <https://sites.stat.washington.edu/mmp/Papers/sdm-wcuts.pdf>

    """
    if nodelist == None:
        nodelist = list(G.nodes())

    if data != None:   
        weight = nx.get_edge_attributes(G, data)

    nlen = len(nodelist)
    index = dict(zip(nodelist, range(nlen)))

    W = np.full((nlen, nlen), np.nan, order=None)
    for i in nodelist:
        for j in G.neighbors(i):
            if j in nodelist:
                W[index[i], index[j]] = 1 if data == None else weight[i,j]

    W[np.isnan(W)] = 0
    W = np.asarray(W, dtype=None)
    
    return W

def asym_weighted_degree_matrix(G, nodelist = None, data = None, in_out_degree = 'out'):
    """Returns the out-degree D(i) of each node i in G, represented as a diagonal
    matrix D. D(i,j) = 0 if i != j, D(i,j) = D(i) otherwise. [1] (can also choose to compute 
    in-degree)

    Parameters
    ----------
    G : NetworkX DiGraph

    nodelist : collection
        A collection of nodes in `G`. If not specified, nodes are all nodes in G
        in G.nodes() order.

    data : object
        Edge attribute key to use as weight. If not specified, D(i) represents the 
        total number of out edges from node i.

    Returns
    -------
    array
        The Out-Degree matrix of G. (can also choose to compute in-degree)

    Notes:
    ------
    D(i) is said to be 1 if the node i has no out edges, this avoids divison
    by zero in further functions and makes node i a 'sink'. [1]   

    References
    ----------
    .. [1] Marina Meila and William Pentney.
           *Clustering by weighted cuts in directed graphs*.
           <https://sites.stat.washington.edu/mmp/Papers/sdm-wcuts.pdf>

    """
    if nodelist == None:
        nodelist = list(G.nodes())

    if data != None:   
        weight = nx.get_edge_attributes(G, data)

    diag = []

    for i in nodelist:
        if in_out_degree == 'out':
            d = sum(weight[e[0],e[1]] for e in G.out_edges(i) )
        if in_out_degree == 'in':
            d = sum(weight[e[0],e[1]] for e in G.in_edges(i) )

        if d !=0:
            diag.append(d)
        else:
            diag.append(1) #avoids division by zero later
    D = np.diag(diag)
    return D

def Hermitian_normalized_Laplacian(D, W, T):
    """Returns the Hermitian part of nL, where nL is the normalized Laplacian. i.e., 
    for T = D, nL = D^(-1/2)(D-W)D^(-1/2) = 1-D^(-1/2)WD^(-1/2) [1].

    Parameters
    ----------
    D : array
        Out-Degree diagonal matrix, see asym_weighted_degree_matrix.

    W : array
        Weight matrix, see asym_weighted_degree_matrix.
    
    T : array
        User defined weighting of the nodes. For WNcut introduced in [1], use T = D.

    Returns
    -------
    array
        The Hermitian part of normalized Laplacian matrix.

    Notes:
    ------
    Given any matrix B, the Hermitian part of B is H(B) = (1/2)(B + B^T) and is always
    a symmetric matrix [1]. Further, [1] also introduces a row weight matrix T', though
    for the WNcut algorthm it is assumed that T'=1 so it omitted from this function. If 
    this is not the case, replace W by T'W and D by T'D [1].   

    References
    ----------
    .. [1] Marina Meila and William Pentney.
           *Clustering by weighted cuts in directed graphs*.
           <https://sites.stat.washington.edu/mmp/Papers/sdm-wcuts.pdf>

    """
    H = 2*D-W-W.transpose()
    T1 = scipy.linalg.fractional_matrix_power(T, -1/2)
    L = (1/2)*np.dot(np.dot(T1,H),T1)
    return L

def k_smallest_eigvec(nodelist, L, k, return_eigenvec = False):
    """Returns the k smallest eigenvectors of the Hermitian part of the normalised
     Laplacian nL. 

    Parameters
    ----------
    G : NetworkX DiGraph

    L : array
        Hermitian part of normalized Laplacian matrix, see Hermitian_normalised_Laplacian
    
    T : array
        User defined weighting of the nodes. For WNcut introduced in [1], use T = D.

    k : int
        k is the number of clusters trying to find. Though it is suggested to use k = 2,
        and used the second smallest eigenvector to break up the graph into two groups
        then run the algorthm again on the groups to further break up the clusters [2].

    Returns
    -------
    array
        Matrix Y, where Y has k smallest eigenvectors of nL as a column vectors.

    Notes:
    ------
    The columns of Y are automatically orthogonal, since nL is real-valued, 
    symmetric matric. Using the Variant to the BestWCut algorthm, so Y is 
    normalized to have rows of length 1. Using L2 norm.  

    References
    ----------
    .. [1] Marina Meila and William Pentney.
           *Clustering by weighted cuts in directed graphs*.
           <https://sites.stat.washington.edu/mmp/Papers/sdm-wcuts.pdf>

    .. [2] Jianbo Shi and Jitendra Malik.
           *Normalized Cuts and Image Segmentation*.
           <https://www.cis.upenn.edu/~jshi/papers/pami_ncut.pdf>

    """
    nlen = len(nodelist)
    eigval, eigvec = linalg.eig(L)
    
    Y = np.full((nlen, k), np.nan, order=None)

    val = [round(i.real,4) for i in eigval]

    s = list(val).copy()
    s.sort()

    for y in range(len(s[:k])):
        for i in range(len(eigval)):
            if val[i] == s[:k][y]:
                n = linalg.norm(eigvec[:,i].transpose())
                Y[:,y] = eigvec[:,i]/n
                break
 
    return preprocessing.normalize(Y, norm="l2") if return_eigenvec == False else (s[:k], preprocessing.normalize(Y, norm="l2"))

def indicator_vector(nodelist, cluster_list):
    nlen = len(nodelist)
    mlen = len(cluster_list)

    X = np.full((nlen, mlen), np.nan, order=None)

    for i in range(len(cluster_list)):
        Ci = [1  if n in cluster_list[i] else 0 for n in nodelist ]
        X[:,i] = Ci

    X = np.asarray(X, dtype=None)
    return X

def WCut(D, W, T, X):
    """Returns weighted cut of graph partition into clusters C = {C1,...,Ck} [1].

    Parameters
    ----------
    D : array
        Out-Degree diagonal matrix, see asym_weighted_degree_matrix.

    W : array
        Weight matrix, see asym_weighted_degree_matrix.
    
    T : array
        User defined weighting of the nodes. For WNcut introduced in [1], use T = D.

    X : array
        The indicator vector of a cluster C = {C1,...,Ck} where the ith column of X represents
        the cluster Ci. Expects X[i][j] if jth node is in cluster Ci and 0 otherwise.

    Returns
    -------
    number
        Returns weighted cut WCut(C) of cluster C = {C1,...,Ck} as introduced in [1].

    References
    ----------
    .. [1] Marina Meila and William Pentney.
           *Clustering by weighted cuts in directed graphs*.
           <https://sites.stat.washington.edu/mmp/Papers/sdm-wcuts.pdf>

    """   
    A = D-W
    
    cut = 0
    Ck_cut = {}
    for i in range(len(X[0])):

        top = np.dot(X[:,i].transpose(), np.dot(A, X[:,i]))
        bottom = np.dot(X[:,i].transpose(), np.dot(T, X[:,i]))

        cut += top/bottom
        Ck_cut[i] = top/bottom
    return cut, Ck_cut



def find_best_clustering(G, start_set, stop_set, network_filename, top_k, nodelist = None, data = None, in_out_degree = 'out', save_file = True):

    if nodelist == None:
        nodelist = list(G.nodes())

    W = asym_weight_matrix(G, nodelist, data)
    D = asym_weighted_degree_matrix(G, nodelist, data, in_out_degree)
    L = Hermitian_normalized_Laplacian(D, W, D)
    eigval, Y = k_smallest_eigvec(nodelist, L, 2, return_eigenvec=True)

    eigv = sorted(Y[:,1])
    diff = []
    for i in range(len( eigv )-1):
        diff.append( ( abs(eigv[i] - eigv[i+1]) , eigv[i], eigv[i+1] ) ) 

    diff.sort(key=lambda a: a[0])

    cut_list = []
    for t in diff[-top_k:]:
        C1 = []
        C2 = []
        m = t[1] + (t[2]-t[1])/2
        for n in nodelist:
            index = nodelist.index(n)
            if Y[:,1][index] >= m:
                C1.append(n)
            else:
                C2.append(n)
        cluster_list = [C1, C2]
        C1s = [i for i in start_set if i in C1]
        C1t = [i for i in stop_set if i in C1]
        C2s = [i for i in start_set if i in C2]
        C2t = [i for i in stop_set if i in C2]
        if (C1s !=[] and C1t !=[]) or (C2s !=[] and C2t !=[]):
            continue
        else:
            V = indicator_vector(nodelist, cluster_list)
            c, Ck_cut = WCut(D, W, D, V)
            cut_list.append((c,m,C1,C2, Ck_cut))

    if cut_list != []:
        cut_list.sort(key=lambda a: a[0])
        c, m, C1, C2, Ck_cut = cut_list[0]
    else:
        return (0, 0, 0, [], [], {0:0, 1:0})
        
    plt.figure(figsize=(15, 8))
    for i in nodelist:
        index = nodelist.index(i)
        if i in C1:
            plt.scatter(x=index, y=Y[:,1][index], color = 'r')
        else:
            plt.scatter(x=index, y=Y[:,1][index], color = 'b')
    plt.title(network_filename + ' second smallest eigenvec with cluster split at ' + str(round(m,4)) + '. WCut='+str(round(c,4)))
    
    if save_file == True:
        plt.savefig(network_filename + 'cluster_split_at_' + str(round(m,4)) + '_cut_'+str(round(c,4))+ '_diagP.png')   # save the figure to file
        plt.close() 
    else:
        plt.show()
    
    return c, eigval[1], m, C1, C2, Ck_cut
