import numpy as np
from collections import OrderedDict
import random
import os
import pandas as pd
from ast import literal_eval

# change color scales on m5 and Tb and Sb plots, upon request from ZI

def init(menu_button=None):
    print("user init setting up menu")
    my_menu = [("User raDeg", "User raDeg"), 
              ("User decDeg", "User decDeg"),
              ("User radDeg", "User radDeg"),
              ("User effarea", "User effarea"),
              ("User readnoise", "User readnoise"),
              ("User m5 u", "User m5 u"),
              ("User m5 g", "User m5 g"),
              ("User m5 r", "User m5 r"),
              ("User m5 i", "User m5 i"),
              ("User m5 z", "User m5 z"),
              ("User m5 y", "User m5 y"),
              ("User fS", "User fS"),
              ("User Tu", "User Tu"),
              ("User Tg", "User Tg"),
              ("User Tr", "User Tr"),
              ("User Ti", "User Ti"),
              ("User Tz", "User Tz"),
              ("User Ty", "User Ty"),  
              ("User Su", "User Su"),
              ("User Sg", "User Sg"),
              ("User Sr", "User Sr"),
              ("User Si", "User Si"),
              ("User Sz", "User Sz"),
              ("User Sy", "User Sy")]
    menu_button.menu = my_menu

    return 0

print("-----------------------user hook setting global variables---------------------")
global_lo = {'u':999, 'g':999, 'r':999, 'i':999, 'z':999, 'y':999}
global_hi = {'u':-999, 'g':-999, 'r':-999, 'i':-999, 'z':-999, 'y':-999}
save_lo = 999
save_hi = -999

def hook(run=None, mode=None, raft=None, ccd=None, test_cache=None, test=None, range_limits=None):
    """
    User hook for test quantity
    :param run: run number
    :return: list of user-supplied quantities to be included in the heat map
    """

    myq = ' '.join(test.split(' ')[1:])
    print('-----------', test)
    if 'm5' in myq or 'fS' in myq:
        idx = 'm5'
        f = myq.split(' ')[-1]
        
        centerV =  {'u': 24.238, 'g': 24.987, 'r': 24.503, 'i': 24.031, 'z': 23.428, 'y': 22.469}
        
    elif myq[0] == 'T':
        idx = 'T'
        f = myq[1]
        
        centerV =  {'u': 0.039, 'g': 0.134, 'r': 0.114, 'i': 0.082, 'z': 0.055, 'y': 0.024}
        
    elif myq[0] == 'S':
        idx = 'S'
        f = myq[1]
        
        centerV =  {'u': 0.062, 'g': 0.163, 'r': 0.128, 'i': 0.09, 'z': 0.059, 'y': 0.028}
        
    else:
        idx = 'a'
        f = myq
    
    global global_hi
    global global_lo
    global save_lo
    global save_hi
    
    s = {"SG0":0, "SG1":16, "SW0":32, "SW1":40,
                  "S00":0, "S01":16, "S02":32, "S10":48,
                  "S11":64, "S12":80, "S20":96, "S21":112, "S22":128}
    slot_index = OrderedDict(s)

    if raft in ["R00", "R04", "R40", "R44"]:
        return
    else:
        test_len = 144
        out_list = [-1.]*test_len #Richard: I could imagine a bunch of things being upset with nan - histograms, sliders, colour map...

        m5dfFile = 'm5_output/%s/%sdf_%s.csv'%(raft,idx,raft)
        if os.path.isfile(m5dfFile):
            df = pd.read_csv(m5dfFile, index_col=0, dtype=object)
            for ccd in slot_index:
                if "SG" in ccd or "SW" in ccd:   # ignore CR
                    continue
    
                key = raft+'_'+ccd
                res = df[key].apply(literal_eval)[f]
                print(f, key,[float('%.2f'%aa) for aa in res])
                amp = 0
                for val in res:
                    if idx == 'm5':
                        aa = val - centerV[f]
                    elif idx == 'T' or idx == 'S':
                        aa = val/centerV[f]
                    out_list[amp + slot_index[ccd]] = aa
                    amp += 1
    #out_list[:] = [np.nan if ele <0 else ele for ele in out_list ]
    
    range_limits["state"] = False
    #print('save_lo = %.2f'%save_lo)
    if f == 'readnoise':
        range_limits["state"] = True
        if save_lo > 998:
            
            range_limits["min"] = 0.
            range_limits["max"] = 18.      
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]
        #if save_hi == range_limits["max"] and save_lo == range_limits["min"]:
    elif f == 'fS':
        range_limits["state"] = True
        if save_lo > 998:
            
            range_limits["min"] = 0.7
            range_limits["max"] = 1.3      
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]
            
    elif myq[0] == 'T':
        range_limits["state"] = True
        if save_lo > 998:
            range_limits["min"] = 0.6
            range_limits["max"] = 1.3
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]
  
    elif myq[0] == 'S':
        range_limits["state"] = True
        if save_lo > 998:
            range_limits["min"] = 0.6
            range_limits["max"] = 1.3
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]

    elif 'm5' in myq:
        range_limits["state"] = True
        if save_lo > 998:
            range_limits["min"] = -1
            range_limits["max"] = 0.2
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]
            
    return out_list
