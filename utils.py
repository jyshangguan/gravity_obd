import os
import datetime
import numpy as np

__all__ = ['modulepath', 'load_obd_module', 'write_obd']

#-> Obtain the current path
pathList = os.path.abspath(__file__).split("/")
#-> Create the path to the gData module
modulepath = "/".join(pathList[0:-1])


def load_obd_module(tname, tpath=modulepath):
    '''
    Load the template parameters.
    
    Parameters
    ----------
    tname : string
        The obd module name.
    tpath : string (default: module path)
        The path to the obd modules.
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
            content = c2[0].strip()  # Remove the white space

            if content[0] == '"':
                content = content[1:]
            if content[-1] == '"':
                content = content[:-1]

        cDict[key] = content
    return cDict


def write_obd(template_name, obd_dict, nchars=40):
    '''
    Write down the OBD.
    
    Parameters
    ----------
    template_name : string
        The full path of the OBD.
    obd_dict : dict
        The dict of all OBD modules.
    nchars : int
        The total number of characters in one line.
    '''
    # Put everything in strings
    lines = []
    for mn in obd_dict:
        lines.append('# {}\n'.format(mn))
        for key in obd_dict[mn]:
            content = '"{}"'.format(obd_dict[mn][key])
            lineformat = '{{0}}{{1:>{0}}};\n'.format(nchars-len(key))
            lines.append(lineformat.format(key, content))
        lines.append('\n')

    f = open(template_name, 'w')
    f.writelines(lines)
    f.close()
