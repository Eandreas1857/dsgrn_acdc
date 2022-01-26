import network_prelim
import multiprocessing as mp
import json, sys
import random 

def run(network_file, filename):
    network_list = json.load(open(network_file))
    random_network_sample = random.sample(network_list, 3)
    pool = mp.Pool(3)
    output = pool.map(network_prelim.main, random_network_sample)
    results = dict(output)
    json.dump(results, open(filename, 'w'))

if __name__=="__main__":
    run(sys.argv[1], sys.argv[2])

## how to run in terminal: 
# python Analysing_networks_script.py network_file.json results_file.json

# ps -a (shows al process that are running on machine)

