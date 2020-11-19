# NFixedPointQuery.py
# Elizabeth Andreas
# Modified from
# DoubleFixedPointQuery.py
# MIT LICENSE 2016
# Shaun Harker


from DSGRN.Query.FixedPointTables import *


import os, sys

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout



class NFixedPointQuery:
  def __init__ (self, database, *bounds):
    
    self.database = database
    c = database.conn.cursor()
    N = len(bounds)

    for i in range(N):
      FPs = bounds[i]
      table = 'Matches' + str(i)

      with HiddenPrints():
            MatchQuery(FPs, table, database)

      MatchQuery(FPs, table, database)

      c.execute('create temp table set' + str(i) + ' as select MorseGraphIndex from ' + table + ' group by MorseGraphIndex;')
      if i==0: 
            c.execute('create temp table match0 as select * from set0')
      else: 
            c.execute('create temp table match' + str(i) + ' as select * from (select * from match' + str(i-1) + ' intersect select * from set' + str(i) + ')')
    
    self.set_of_matches = set([ row[0] for row in c.execute('select MorseGraphIndex from match' + str(N-1))])   
    
    overlap = set()
    for i in range(N):
        for j in range(N):
            if j != i:
                overlap1 = set([ row[0] for row in c.execute('select Label from (select Label from Matches' + str(i) + ' intersect select Label from Matches' + str(j) + ')')])
                overlap = overlap.union(overlap1)
    self.overlap = overlap
    
    for i in range(N):
        c.execute('drop table Matches' + str(i))
        c.execute('drop table set' + str(i))
        c.execute('drop table match' + str(i))
        
  def matches(self):
    """
    Return entire set of matches if bounds do not overlap, otherwise return overlaping FP
    """
    CRED = '\033[1;31;47m '
    CEND = '\033[1;30;47m '

    if self.overlap == set():
        return self.set_of_matches
    else:
        print(CRED + 'ERROR:' + CEND + 'overlapping bounds for ', self.overlap)
        return self.overlap

  def matches_with_PI(self):
    CRED = '\033[1;31;47m '
    CEND = '\033[1;30;47m '
    database = self.database
    c = database.conn.cursor()
    PGI1 = set()
    if self.overlap == set():
        for i in self.set_of_matches:
            c.execute('create temp table C' + str(i) + ' as select * from Signatures where MorseGraphIndex =' + str(i) )
            set_of_matches = set([ (row[2],row[0],row[1]) for row in c.execute('select * from C' + str(i))])
            PGI1 = PGI1.union(set_of_matches)
            c.execute('drop table C' + str(i))
        return PGI1
    else:
        print(CRED + 'ERROR:' + CEND + 'overlapping bounds for ', self.overlap)
        return self.overlap
    
  def __call__ (self, morsegraphindex ):
    """ 
    Test if a single mgi is in the set of matches
    """
    return morsegraphindex in self.set_of_matches

