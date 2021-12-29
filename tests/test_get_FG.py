# get_FG.py

import pytest
import DSGRN
from DSGRN import *

import sys
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/dsgrn_acdc/src')

from get_FG import *

def test_get_hex_FG():

	database = Database("/home/elizabeth/Desktop/ACDC/ACDC_StrongEdges.db")
	assert get_hex_FG(database, 'Hb') == {0: ['0'], 1: ['8'], 2: ['A', 'C'], 3: ['E'], 4: ['F']}

	assert get_hex_FG(database, 'Kni') == {0: ['000'],
 1: ['200'],
 2: ['208', '240', '600'],
 3: ['640', 'E00', '248', '608'],
 4: ['249', '618', '648', '6C0', 'E08', 'E40'],
 5: ['649', '658', '6C8', 'E18', 'E48', 'EC0'],
 6: ['EC8', 'FC0', '659', '6C9', '6D8', 'E38', 'E49', 'E58'],
 7: ['EC9', 'ED8', 'FC8', '6D9', 'E59', 'E78'],
 8: ['ED9', 'EF8', 'FC9', 'FD8', '6DB', 'E79'],
 9: ['FF8', 'EDB', 'EF9', 'FD9'],
 10: ['EFB', 'FDB', 'FF9'],
 11: ['FFB'],
 12: ['FFF']}
