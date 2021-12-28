import DSGRN, itertools
from DSGRN import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def plot(fig_data, Hex_poset, title_size, text_size, color = "GnBu", fig_size = (8,8)):    
    Phenotype_layer = [0, 1, 2, 3, 4, 5, 6, 7]
    Hex_Poset_Layer = list(range(len(Hex_poset)-1,-1,-1))
    
    data_T = (np.array(fig_data)).transpose()
    data = np.flip(data_T,0)
    #fig = plt.figure()
    fig, ax = plt.subplots(1,1, figsize=fig_size)

    im = ax.imshow(data, cmap=color)

    # We want to show all ticks...
    ax.set_yticks(np.arange(len(Hex_Poset_Layer)))
    ax.set_xticks(np.arange(len(Phenotype_layer)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(Hex_Poset_Layer, size = text_size)
    ax.set_xticklabels(Phenotype_layer, size = text_size)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")

    cbar = ax.figure.colorbar(im, ax=ax, cmap= color, fraction=0.0405)
    cbar.ax.set_ylabel(" ", rotation=-90, va="bottom", size = text_size)

    # Loop over data dimensions and create text annotations.
    for i in range(len(Hex_Poset_Layer)):
        for j in range(len(Phenotype_layer )):
            text = ax.text(j, i, data[i, j],
                           ha="center", va="center", color="black", size = text_size)

    ax.set_title("Hexcodes in scc's by Phenotype Layer", size=title_size)
    ax.set_ylabel('Hex Poset Layer', size = text_size)
    ax.set_xlabel('Phenotype Pattern Layer', size = text_size)
    #fig.tight_layout()
    #fig.savefig('test.png')
    plt.show()


def Fullconn_get_data_scc(database, Paths, scc):
    
    pg = DSGRN.ParameterGraph(database.network)
    
    Full_list = {0: ['00'], 1:['40'], 2:['44','50','C0' ], 3: ['54','C4','D0'], 4: ['55','D4','CC','F0'], 5: ['D5','DC','F4'], 6: ['DD','F5','FC'], 7: ['FD'], 8:['FF']}
    path_layer_nodes = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}

    for path in Paths:
        for node in path:
            if node not in path_layer_nodes[node[0]]:
                path_layer_nodes[node[0]].append(node)
    scc_list = []
    
    for i in scc:
        sub_list = []
        for j in range(len(scc[i])):

            for k in path_layer_nodes[i]:
                if k in scc[i][j]:
                    sub_list.append(j)
        if sub_list != []:
            scc_list.append(sub_list)
            
    Hb_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Gt_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Kr_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Kni_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}

    for item in range(len(scc_list)):
        for section in scc_list[item]:
            for node in scc[item][section]:
                params = pg.parameter(node[-1])
                logic = params.logic()
                Hb = ((logic[0]).stringify())[6:-2]
                Hb_hex[node[0]].append(Hb)
                Gt = ((logic[1]).stringify())[6:-2]
                Gt_hex[node[0]].append(Gt)
                Kr = ((logic[2]).stringify())[6:-2]
                Kr_hex[node[0]].append(Kr)
                Kni = ((logic[3]).stringify())[6:-2]
                Kni_hex[node[0]].append(Kni)
                
    Hb_data = []
    Gt_data = []
    Kr_data = []
    Kni_data = []

    for n in [0, 1, 2, 3, 4, 5, 6, 7]:
        A = [len([i for i in Hb_hex[n] if i in Full_list[key]]) for key in Full_list]
        Hb_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7]:
        A = [len([i for i in Gt_hex[n] if i in Full_list[key]]) for key in Full_list]
        Gt_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7]:
        A = [len([i for i in Kr_hex[n] if i in Full_list[key]]) for key in Full_list]
        Kr_data.append(A)

    for n in [0, 1, 2, 3, 4, 5, 6, 7]:
        A = [len([i for i in Kni_hex[n] if i in Full_list[key]]) for key in Full_list]
        Kni_data.append(A)
        
    return Hb_data, Gt_data, Kr_data, Kni_data

def StrongEdges_get_data_scc(database, Paths, scc):
    
    pg = DSGRN.ParameterGraph(database.network)
    
    Hb_list = {0: ['0'], 1: ['8'], 2: ['A', 'C'], 3: ['E'], 4: ['F']}
    Kni_list = {0:['000'], 1: ['200'], 2: ['208', '240', '600'], 3: ['248','608','640','E00'], 4: ['249','648','618','6C0','E08','E40'], 5:['649','658','6C8','E18','E48','EC0'], 6: ['659','6C9','6D8', 'E49','E58','EC8','E38','FC0'], 7: ['6D9','E59','EC9','ED8','E78','FC8'], 8: ['6DB','ED9','E79','EF8','FC9','FD8'], 9: ['EDB','EF9','FD9','FF8'], 10: ['EFB','FDB','FF9'], 11:['FFB'], 12: ['FFF']}
    path_layer_nodes = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}

    for path in Paths:
        for node in path:
            if node not in path_layer_nodes[node[0]]:
                path_layer_nodes[node[0]].append(node)
    
    scc_list = []
    for i in scc:
        sub_list = []
        for j in range(len(scc[i])):

            for k in path_layer_nodes[i]:
                if k in scc[i][j]:
                    sub_list.append(j)
        if sub_list != []:
            scc_list.append(sub_list)
            
    Hb_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Gt_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Kr_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    Kni_hex = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}

    for item in range(len(scc_list)):
        for section in scc_list[item]:
            for node in scc[item][section]:
                params = pg.parameter(node[-1])
                logic = params.logic()
                Hb = ((logic[0]).stringify())[6:-2]
                Hb_hex[node[0]].append(Hb)
                Gt = ((logic[1]).stringify())[6:-2]
                Gt_hex[node[0]].append(Gt)
                Kr = ((logic[2]).stringify())[6:-2]
                Kr_hex[node[0]].append(Kr)
                Kni = ((logic[3]).stringify())[6:-2]
                Kni_hex[node[0]].append(Kni)
                
    Hb_data = []
    Gt_data = []
    Kr_data = []
    Kni_data = []

    for n in [7, 6, 5, 4, 3, 2, 1, 0]:
        A = [len([i for i in Hb_hex[n] if i in Hb_list[key]]) for key in Hb_list]
        Hb_data.append(A)

    for n in [7, 6, 5, 4, 3, 2, 1, 0]:
        A = [len([i for i in Gt_hex[n] if i in Kni_list[key]]) for key in Kni_list]
        Gt_data.append(A)

    for n in [7, 6, 5, 4, 3, 2, 1, 0]:
        A = [len([i for i in Kr_hex[n] if i in Hb_list[key]]) for key in Hb_list]
        Kr_data.append(A)

    for n in [7, 6, 5, 4, 3, 2, 1, 0]:
        A = [len([i for i in Kni_hex[n] if i in Kni_list[key]]) for key in Kni_list]
        Kni_data.append(A)
        
    return Hb_data, Gt_data, Kr_data, Kni_data
