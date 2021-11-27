import os
import datetime
import numpy as np

#-> Obtain the current path
pathList = os.path.abspath(__file__).split("/")
#-> Create the path to the gData module
modulepath = "/".join(pathList[0:-1])

def load_obd_module(tname, tpath=modulepath):
    '''
    Load the template parameters.
    '''
    f = open('{0}/tpl_modules/{1}.obd'.format(tpath, tname), 'r')
    text = f.readlines()
    
    cDict = {}
    for loop, t in enumerate(text):
    
        # Skip the comments
        if t[0] == '#':
            continue
            
        c0 = t.split(';')
        
        # Skip the empty lines
        if len(c0[0]) == 1:
            continue
     
        c1 = c0[0].split()
        key = c1[0]
        if len(c1) == 1:
            content = ''
        else:
            content = ' '.join(c1[1:])
            c2 = content.split('#')
            content = c2[0]
            
            if content[0] == '"':
                content = content[1:]
            if content[-1] == '"':
                content = content[:-1]
            
        cDict[key] = content
    return cDict
