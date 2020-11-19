# NsearchQuery
# Elizabeth Andreas

from NFixedPointQuery import *
from DSGRN.Query.NstableQuery import *

class NsearchgoeQuery:
    def __init__(self, database, goe1, goe2 , bounds1, bounds2):
    
        self.database = database
        c = database.conn.cursor()
        
        NFP = NFixedPointQuery(database, *bounds1).matches()
        
        if goe1 == '=':
            N = len(bounds1)
            X1 = NstableQuery(database, N).matches()
            X2 = NstableQuery(database, N+1).matches()
            inter = set(X1.difference(X2))
            MGset = [i for i in inter if i in NFP]
            
        if goe1 == '<':
            MGset = list(NFP)
            
        self.MGset = MGset
        # diff is all of the MG's with N1 stability in the database
        
        N1 = len(bounds2)
        X1 = NstableQuery(database, N1).matches()
        X2 = NstableQuery(database, N1+1).matches()
        diff = set(X1.difference(X2))  
        
        # PGI1 is the set of Parameter Index' assosiated to the Morse Graph inputs
        
        string = 'create temp table C as select * from Signatures where MorseGraphIndex in ({seq})'.format(
        seq=','.join(['?']*len(MGset)))
        c.execute(string, MGset)
        PGI1 = [ row[0] for row in c.execute('select ParameterIndex from C')]
        c.execute('drop table C')

        # PGIfinal is the set of tuples where first value in tuple is original Parameter Index, second value is the adjacent Parameter index, third is the MG index. Note that the tuples in this set are only those what have an original Parameter node in the bounds AND an adjacent parameter with =N (not >=N) fixed points. 
        PGIfinal = set()
        new = NFixedPointQuery(database, *bounds2).matches()
        want =  [i for i in new if i in diff]
        
        if goe2 == '=':
            for node in PGI1:
                adj_nodes = database.parametergraph.adjacencies(node)
                sql="create temp table C as select *, " + str(node) + " as ParentPGI" + " from Signatures where ParameterIndex in ({seq})".format(
                seq=','.join(['?']*len(adj_nodes)))
                c.execute(sql, adj_nodes)

                table2 = 'create temp table D as select * from C where MorseGraphIndex in ({seq})'.format(
                seq=','.join(['?']*len(want)))
                c.execute(table2, want)

                edges = set([ (row[2],row[0],row[1]) for row in c.execute('select * from D')])
                PGIfinal = PGIfinal.union(edges)

                c.execute('drop table C')
                c.execute('drop table D')
        
        else:
            for node in PGI1:
                adj_nodes = database.parametergraph.adjacencies(node)
                sql="create temp table C as select *, " + str(node) + " as ParentPGI" + " from Signatures where ParameterIndex in ({seq})".format(
                seq=','.join(['?']*len(adj_nodes)))
                c.execute(sql, adj_nodes)

                table2 = 'create temp table D as select * from C where MorseGraphIndex in ({seq})'.format(
                seq=','.join(['?']*len(list(new))))
                c.execute(table2, list(new))

                edges = set([ (row[2],row[0],row[1]) for row in c.execute('select * from D')])
                PGIfinal = PGIfinal.union(edges)

                c.execute('drop table C')
                c.execute('drop table D')
        
        self.PGIfinal = PGIfinal
        
    def matches(self):
              
        return set(self.PGIfinal)
    
    def stability_of_matches(self):
        # Return PGIfinal where the adj parameter node is in bounds.
          
        return self.PGIfinal
