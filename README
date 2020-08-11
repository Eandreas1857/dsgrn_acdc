# Getting started

This module uses a DSGRN database constructed from a network in a text file, along
with an ordered list of fixed points, to find paths through the networks parameter graph with correct corresponding Morse graphs. 

__Dependencies:__ Python 3.6/3.7, itertools, graphviz, DSGRN (https://github.com/marciogameiro/DSGRN) and its dependencies.

__DSGRN References:__ http://epubs.siam.org/doi/abs/10.1137/15M1052743, https://link.springer.com/chapter/10.1007/978-3-319-67471-1_19, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5975363/, https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006121

__Paper References:__ https://pubmed.ncbi.nlm.nih.gov/31169494/, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2646127/

To install, do
```bash 
cd dsgrn_acdc
. install.sh
```

See the `notebooks` folder for detailed examples of different classes in this repository. See `networks` for the network files used in these examples. To construct databases, change directories in terminal to where the database is going to be stored. Then do
```bash
mpiexec -n num_processes Signatures my_network.txt my_database.db
```
where num_processes is the number of cores wanting to use, my_network.txt is the name of the network text file and my_database.db is the wanted name of the database.

# Parameters 
__Required:__

`networkfile` = Path to a text file containing a single network.

`databasefile` = Path to a database constructed from a network file as explained above. 
`fixedpointlist` = List of fixed points wanting to use for finding a path in the parameter graph. See examples in the `notebooks` folder for how to input fixed points into functions used in this module. 

__NOTES:__

* Networks can be analyzed in essential or inessential mode (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006121). Essential mode constructs parameter graphs where every edge in the network has a nontrivial role in the network dynamics. Inessential mode includes these edges that may not be effectively functioning. Calculating in essential mode usually results in a much smaller parameter graph, which means much faster computation. 

* Currently, this module can only find paths for smaller parameter graphs. Parameter graphs consisting of more than 2 million nodes are not recommended. 

* DSGRN does not support self-repressing edges or zero out-edges except in `develop` branches that are pre-alpha. However, nodes with no in-edges are supported. 

* DSGRN supports most node topologies with up to 4 in-edges, but has spotty coverage for the number of associated out-edges. These values can be calculated using https://github.com/marciogameiro/DSGRN, but the code maintainer will have to be contacted. For this reason, `constrained_inedges` should probably have a max value no greater than 4. Currently, `constrained_outedges` has a heuristic upper bound of 5 outedges. But again, contact the maintainer to have more topologies implemented.

* Using the `constrained_inedges` and `constrained_outedges` filters may substantially reduce computation time when the upper bound of `range_operations` is high and adding an edge has high probability.

* Network searches will always assume that activating edges are summed together. Activating edges that are multiplied will be recast into addition, potentially enlarging the size of the parameter graph for 3 or more in-edges. (See references for more detail.)

# Functions and Outputs
* PGDraw.py is a visual aid that is helpful in seeing the structure of the parameter graph for visual verification that the path searches are working properly. Though this will not be helpful on large parameter graphs. Returns graph with nodes and edges of parameter graph with a label pair of the associated parameter graph index and the associated Morse graph index. 

* NFixedPointQuery.py searches a database for all of the Morse graph indexes that have the desired fixed points. If fixed points overlap, it returns an error. Otherwise it returns a set of all of the Morse graph index who have all of the desired fixed points. Note that the returned Morse graphs might also have MORE fixed points then the ones searched for, it just has to have ALL of the ones we searched for. 

* NsearchgoeQuery.py finds edges from all parameter graph nodes that have desired fixed points (found using NFixedPointQuery.py) to parameter graphs that have a specified set of fixed points. Use .matches() to return a set of tuples, where first value in tuple is original Parameter Index, second value is the adjacent Parameter index, third is the Morse graph index. There is an option to either search only for Morse graphs with ONLY the fixed points you want ('=') or to be equal to or greater than ('<'). This can be done separately for both the starting nodes in the search or the ending and are the goe1 and goe2 inputs. 

* MGsearchthroughPG.py calls functions from funforMGpaths.py until all of the fixed point lists are exhausted. This function finds all of the paths in a parameter graph that have the desired sequence of fixed point lists. Each path is returned as a list of parameter graph indexes that were passed through. Use .allpaths() to return all paths found. Use .shortestpaths() to return only the shortest paths found. 

# Troubleshooting and Common Problems

* If a parameter graph is to large, a path search may result in a stack overflow or recursion depth error. If it is a recursion depth error, try increasing the allowable resursion depth on your computer. Additionally, if the network is not already in essental mode, changing it to essential mode may correct the error.
```python
import sys
sys.setrecursionlimit(10**8)
``` 

* Fixed points lists that are being computed for a single Morse graph cannot overlap, if they do, then an error will populate as well as point out the fixed points causing the error. 