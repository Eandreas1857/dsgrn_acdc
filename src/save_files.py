import json
import ast


def save_json(edges, filename):
    mydict = {}
    for key in edges.keys():
        if type(key) is not str:
            mydict[str(key)] = edges[key]
            
    with open(filename, 'w') as fp:
        json.dump(mydict, fp)


def load_json(file_to_load):
    with open(file_to_load) as f:
        data = json.load(f)
        
    new_data = {}
    for key in data:
        new_data[ast.literal_eval(key)] = [tuple(i) for i in data[key] ]
    return new_data

def load_json_scc(file_to_load):
    with open(file_to_load) as f:
        data = json.load(f)
        
    new_data = {}
    for key in data:
        for sublist in data[key]:
            if ast.literal_eval(key) not in new_data:
                new_data[ast.literal_eval(key)] = [[(i[0], i[1]) for i in sublist] ]
            else:
                new_data[ast.literal_eval(key)] += [[(i[0], i[1]) for i in sublist] ]
    return new_data

def load_json_cond(file_to_load):
    with open(file_to_load) as f:
        data = json.load(f)
        
    new_data = {}
    for key in data:
        new_data[ast.literal_eval(key)] = [int(i) for i in data[key] ]
    return new_data
