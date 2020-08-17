# stats.py
# Elizabeth Andreas
# MIT LICENSE 2020

class stats:
    def __init__(self, database, MGtype, FPonly_YorN):
        self.database = database

        c = database.conn.cursor()
        
        total_each = []
        if FPonly_YorN == 'N':
            MGtype = [x-1 for x in MGtype]
            for vertex in MGtype:  
                    c.execute('create temp table F1 as select * from MorseGraphAnnotations where Vertex is ' + str(vertex))
                    c.execute('create temp table F2 as select * from MorseGraphAnnotations where Vertex is ' + str(vertex+1))
                    c.execute('create temp table F3 as select MorseGraphIndex from F1 except select MorseGraphIndex from F2')
                    c.execute('SELECT COUNT(*) from F3')
                    total_each.append(c.fetchone()[0])
                    c.execute('drop table F1')
                    c.execute('drop table F2')
                    c.execute('drop table F3')
        self.total_each = total_each
        
        if FPonly_YorN == 'Y':
            c.execute("create temp table FPonly as select MorseGraphIndex from MorseGraphAnnotations where Label like 'FP%'")
            c.execute('SELECT COUNT(DISTINCT MorseGraphIndex) from FPonly')
            count_FPonly = c.fetchone()[0]
            self.count_FPonly = count_FPonly
            for vertex in MGtype:  
                    c.execute('create temp table F1 as select * from FPonly group by MorseGraphIndex having count(MorseGraphIndex) = ' + str(vertex))

                    c.execute('SELECT COUNT(*) from F1')
                    total_each.append(c.fetchone()[0])
                    c.execute('drop table F1')

            c.execute('drop table FPonly')
            
        self.total_each = total_each
        
        self.FPonly_YorN = FPonly_YorN
        
    def total(self):
        return self.total_each
    
    def percentage(self):
        if self.FPonly_YorN == 'N':
            c = self.database.conn.cursor()
            c.execute('SELECT COUNT(DISTINCT MorseGraphIndex) from MorseGraphAnnotations')
            All_MGs = c.fetchone()[0]
            perc = [round((x/All_MGs)*100,2) for x in self.total_each]
        else:
            perc = [round((x/self.count_FPonly)*100,2) for x in self.total_each]
        return perc
