
import DSGRN
from copy import deepcopy

def Hb_high2low(network, paramslist):
    g = deepcopy(paramslist)
    pg = DSGRN.ParameterGraph(network)
    new_start = []
    for i in paramslist[0]:
        params = pg.parameter(i[1])
        s = params.logic()
        b = s[0].stringify()
        if b[6:-2] == 'F'*len(b[6:-2]):
            new_start.append(i)

    if new_start == []:
        print('Abs high not in list, comptuting next best thing')
        for i in paramslist[0]:
            params = pg.parameter(i[1])
            s = params.logic()
            b = s[0].stringify()
            if 'F' in b[6:-2]:
                new_start.append(i)
                
    new_end = []
    for i in paramslist[-1]:
        params = pg.parameter(i[1])
        s = params.logic()
        b = s[0].stringify() 
        if b[6:-2] == '0'*len(b[6:-2]):
            new_end.append(i)
            
            
    if new_end == []:
        for i in paramslist[-1]:
            print('Abs low not in list, comptuting next best thing')
            params = pg.parameter(i[1])
            s = params.logic()
            b = s[0].stringify() 
            if b[6] == '0' or b[6] == '1':
                new_end.append(i)
            
    g[0] = new_start
    g[-1] = new_end
    
    return g

def Kni_low2high(network, paramslist):
    g = deepcopy(paramslist)
    pg = DSGRN.ParameterGraph(network)
    new_start = []
    for i in paramslist[0]:
        params = pg.parameter(i[1])
        s = params.logic()
        b = s[3].stringify()
        if b[6:-2] == '0'*len(b[6:-2]):
            new_start.append(i)
            

    if new_start == []:
        print('Abs high not in list, comptuting next best thing')
        for i in paramslist[0]:
            params = pg.parameter(i[1])
            s = params.logic()
            b = s[3].stringify()
            if b[6] == '0' or b[6] == '1':
                new_start.append(i)
            
                
    new_end = []
    for i in paramslist[-1]:
        params = pg.parameter(i[1])
        s = params.logic()
        b = s[3].stringify() 
        if b[6:-2] == 'F'*len(b[6:-2]):
            new_end.append(i)
            
            
    if new_end == []:
        for i in paramslist[-1]:
            print('Abs low not in list, comptuting next best thing')
            params = pg.parameter(i[1])
            s = params.logic()
            b = s[3].stringify() 
            if 'F' in b[6:-2]:
                new_end.append(i)
            
    g[0] = new_start
    g[-1] = new_end
    
    return g

def newlist(network, paramslist):
    Redu_Hb = Hb_high2low(network, paramslist)
    Redu_Kni = Kni_low2high(network, Redu_Hb)
    
    pg = DSGRN.ParameterGraph(network)
    params1 = pg.parameter(((Redu_Kni[0])[0])[-1])
    params2 = pg.parameter(((Redu_Kni[-1])[0])[-1])
    print("Checking first layer:")
    print(params1)
    print("Checking last layer:")
    print(params2)
    
    return Redu_Kni


