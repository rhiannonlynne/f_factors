import numpy as np
from collections import OrderedDict
import random
import os
import pandas as pd
from ast import literal_eval


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
    elif myq[0] == 'T':
        idx = 'T'
        f = myq[1]
    elif myq[0] == 'S':
        idx = 'S'
        f = myq[1]
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
                    out_list[amp + slot_index[ccd]] = val
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
            mymin = {'u':0.02, 'g':0.11, 'r':0.1, 'i':0.07, 'z':0.04, 'y':0.02}
            range_limits["min"] = mymin[f]
            mymax = {'u':0.05, 'g':0.1392, 'r':0.12, 'i':0.085, 'z':0.07, 'y':0.04}
            range_limits["max"] = mymax[f]
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]
  
    elif myq[0] == 'S':
        range_limits["state"] = True
        if save_lo > 998:
            mymin = {'u':0.02, 'g':0.14, 'r':0.11, 'i':0.08, 'z':0.05, 'y':0.025}
            range_limits["min"] = mymin[f]
            mymax = {'u':0.08, 'g':0.18, 'r':0.13, 'i':0.095, 'z':0.07, 'y':0.045}
            range_limits["max"] = mymax[f]
            save_hi = range_limits["max"]
            save_lo = range_limits["min"]

    elif 'm5' in myq:
        m5min = {'u':23.4, 'g':24.6, 'r':24.3, 'i':23.6, 'z':22.9, 'y':21.7}
        m5design = {'u':23.9, 'g':25.0, 'r':24.7, 'i':24.0, 'z':23.3, 'y':22.1}
        #print('global_lo[f] ???????? = %.2f'%global_lo[f])
        #range_limits["state"] = False
        if global_lo[f] > 998:
            #print('xxxxxxxxxxxxx')
            range_limits["state"] = True
            if 1: #use SRD min and design specs as range limitss
                global_lo[f] = m5min[f]
                global_hi[f] = m5design[f] #range_limits['min']+1
            else: # for each plot, set low to mean-2*std, and high to max
                #right now this part doesn't really work, because it will set min and max
                #based on whatever first raft the code looks at, then never enter this if block again
                if not np.all(np.array(out_list)<0):
                    tt = [np.mean(out_list) if ele <0 else ele for ele in out_list]
                    aa = np.mean(tt) - 2*np.std(tt) 
                    bb = np.max(tt) #np.min([np.mean(tt) + 1*np.std(tt), np.max(tt)])
                    if aa <global_lo[f]:
                        global_lo[f] = aa
                    if bb>global_hi[f]:
                        global_hi[f] = bb

            

            range_limits["min"] = global_lo[f]
            #print('xxxxxxxx - setting min to %.2f'%global_lo[f])
            range_limits["max"] = global_hi[f]
            #print('yyyyyyyyy- setting max to %.2f'%global_hi[f])
        if range_limits["min"] == global_lo[f] and range_limits["max"] == global_hi[f]:
            range_limits["state"] = True
            range_limits["min"] = global_lo[f]
            range_limits["max"] = global_hi[f]
            
    return out_list
