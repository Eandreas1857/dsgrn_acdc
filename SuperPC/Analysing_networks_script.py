import network_prelim
import multiprocessing as mp
import json, sys
import random 

def run(network_file, filename):
    network_list = json.load(open(network_file))
    #random_network_sample = random.sample(network_list, 50)
    pool = mp.Pool(3)
    output = pool.map(network_prelim.main, network_list)
    results = dict(output)
    json.dump(results, open(filename, 'w'))

if __name__=="__main__":
    run(sys.argv[1], sys.argv[2])

## how to run in terminal: 
# ssh ElizabethAndreas@megaplex.msu.montana.edu
# python ../Analysing_networks_script.py network_file.json results_file.json >output.log 2>&1 </dev/null &

# python ~/GIT/dsgrn_acdc/SuperPC/Analysing_networks_script.py ~/GIT/dsgrn_acdc/SuperPC/network_file.json results_file.json >output.log 2>&1 < /dev/null &

# tail output.log
# cat output.log
# less output.log
# head output.log

# ps -a (shows all process that are running on machine)

