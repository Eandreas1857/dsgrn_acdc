from network_prelim import *
def get_cut_w_disjoint_start_stop(network_filename):

    network_txt_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".txt"
    with open(network_txt_filename,"r") as f:
        network = f.read()

    db_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".db"
    database = Database(db_filename)

    grad_graph_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/notebooks/" +"grad_graph_strict_" + network_filename
    
    grad_graph = load_json(grad_graph_filename)

    G = nx.DiGraph() #building networkx graph
    for node in grad_graph:
        G.add_node(node)
        for edge in grad_graph[node]:
            G.add_edge(node, edge)
    
    out_edges = get_number_out_edges_from_string(network)
    Hb_list, Kni_list = get_Hb_Kni_list(database)
    Hb_max = len(Hb_list)-1
    Kni_max = len(Kni_list)-1

    FP_Poset, FP_Regions = get_FP_Poset(out_edges)

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)
    P = get_product_graph(database, N, scc, FP_Poset)
    breaks = find_breaks_in_FG_comb(database, network_filename, P, scc, Hb_max, Kni_max, FP_Regions) 
    keep = build_diag(Hb_max, Kni_max, breaks)
    diagP = remove_unnecessary_nodes_in_P(P, breaks, keep, scc, Kni_max)
    start_set, stop_set = return_start_stop_set(database, diagP, scc, Hb_max, Kni_max, FP_Regions)
    c, eigval, m, C1, C2, Ck_cut = find_best_clustering(diagP, start_set, stop_set, network_filename, 20, nodelist = None, data = 'weight', in_out_degree = 'out', save_file = True)
    
    C1s = [i for i in start_set if i in C1]
    C1t = [i for i in stop_set if i in C1]
    C2s = [i for i in start_set if i in C2]
    C2t = [i for i in stop_set if i in C2]

    results = (network_filename, {'WCut': c, 'eigval': eigval, 'Ck cut': Ck_cut, 'Ck start/stop': [C1s, C1t, C2s, C2t] })

    return results

def main(results):

    for i in results[10:]:
        key = 'Ck start/stop'
        if key in i[1]:
            count +=1
            if (i[1][key][0] !=[] and i[1][key][1] !=[]) or (i[1][key][2] !=[] and i[1][key][3] !=[]):
                result = get_cut_w_disjoint_start_stop(i[0])
        else:
            result = get_cut_w_disjoint_start_stop(i[0])

    return result