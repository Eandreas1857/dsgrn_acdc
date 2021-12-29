# NsearchQuery
# Elizabeth Andreas

from NFixedPointQuery import *
from DSGRN.Query.NstableQuery import *

class NsearchQuery:
    def __init__(self, database, N1, *bounds):
    
        self.database = database
        self.N1 = N1
        c = database.conn.cursor()
        
        NFP = NFixedPointQuery(database, *bounds).matches()
        
        N = len(bounds)
        X1 = NstableQuery(database, N).matches()
        X2 = NstableQuery(database, N+1).matches()
        MGset = set(X1.difference(X2)).intersection(NFP)
        print(MGset)
        self.MGset = MGset
        
        # diff is all of the MG's with N1 stability in the database

        X1 = NstableQuery(database, N1).matches()
        X2 = NstableQuery(database, N1+1).matches()
        diff = set(X1.difference(X2))  
        
        # PGI1 is the set of Parameter Index' assosiated to the Morse Graph inputs
        
        PGI1 = set()
        for i in MGset:
            c.execute('create temp table C' + str(i) + ' as select * from Signatures where MorseGraphIndex =' + str(i) )
            set_of_matches = set([ row[0] for row in c.execute('select ParameterIndex from C' + str(i))])
            PGI1 = PGI1.union(set_of_matches)
            c.execute('drop table C' + str(i))
        print('Assosiated PGI to MGI ' + str(MGset) + ':', PGI1)
        
        self.PGI1 = PGI1

        # PGIfinal is the set of tuples where first value in tuple is original Parameter Index, second value is the adjacent Parameter index, third is the MG index. Note that the tuples in this set are only those what have an original Parameter node in the bounds AND an adjacent parameter with =N (not >=N) fixed points. 
        
        PGIfinal = set()
        for i in PGI1:
            set_of_matches = set(database.parametergraph.adjacencies(i))
            for j in set_of_matches:
                for n in diff:
                    c.execute('create temp table C' + str(j) + ' as select *, '+ str(i) +' as middle' + ' from Signatures where ParameterIndex =' + str(j) )
                    c.execute('create temp table D' + str(n) + ' as select * from C' + str(j) + ' where MorseGraphIndex =' + str(n) )
                    set_of_matches1 = set([ (row[2],row[0],row[1]) for row in c.execute('select * from D' + str(n))])
                    PGIfinal = PGIfinal.union(set_of_matches1)
                    c.execute('drop table C' + str(j))
                    c.execute('drop table D' + str(n))
        print('Adjacent PGI:', PGIfinal)
        
        
        self.PGIfinal = PGIfinal
               
    def matches(self):
              
        return set(self.PGIfinal)
    
    def stability_of_matches(self, *bounds):
        # Return PGIfinal where the adj parameter node is in bounds.
              
        N1 = len(bounds)

        X1 = NstableQuery(self.database, N1).matches()
        X2 = NstableQuery(self.database, N1+1).matches()
        diff = set(X1.difference(X2))  
        
        NFP = NFixedPointQuery(self.database, *bounds).matches()
        
        want = NFP.intersection(diff)
        want

        len(self.PGIfinal)
        inbounds = set()
        for i in want:
            for j in range(len(self.PGIfinal)):
                if [row[2] for row in self.PGIfinal][j] == i:
                   X = set( [{(row[0],row[1],row[2])} for row in self.PGIfinal][j])
                   inbounds = inbounds.union(X)
        return inbounds
