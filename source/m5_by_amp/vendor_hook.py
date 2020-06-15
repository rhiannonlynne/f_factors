import numpy as np
from collections import OrderedDict
import random
import pandas as pd

def init(menu_button=None):
    print("user init setting up menu")
    my_menu = [("User vendor", "User vendor")]

    menu_button.menu = my_menu

    return 0

save_lo = -999
save_hi = -999

# mapping based on 
#https://confluence.slac.stanford.edu/pages/viewpage.action?spaceKey=LSSTCAM&title=Raft+Delivery+and+Acceptance+Testing+Status
dd = pd.read_csv('raftInstall.csv',index_col=0)


def hook(run=None, mode=None, raft=None, ccd=None, test_cache=None, test=None, range_limits=None):
    """
    User hook for test quantity
    :param run: run number
    :return: list of user-supplied quantities to be included in the heat map
    """
# SW0 and SW1 only use 0-7 for resukts. Arrange the indexes to overwrite the back half of
# SW0 with the front half of SW1

    global save_lo
    global save_hi
    
    s = {"SG0":0, "SG1":16, "SW0":32, "SW1":40,
                  "S00":0, "S01":16, "S02":32, "S10":48,
                  "S11":64, "S12":80, "S20":96, "S21":112, "S22":128}
    slot_index = OrderedDict(s)

    if raft in ["R00", "R04", "R40", "R44"]:
        test_len = 48
        out_list = [-1.] * test_len
        return out_list
    else:
        test_len = 144
        if dd.vendor[raft].lower() == 'e2v':
            out_list = [1.]*test_len
        else:
            out_list = [-1.]*test_len
            
    if save_lo == -999.:
        save_hi = range_limits["max"]
        save_lo = range_limits["min"]

    if save_hi == range_limits["max"] and save_lo == range_limits["min"]:
        range_limits["min"] = -1
        range_limits["max"] = 1.
        range_limits["state"] = True
        
    return out_list
