import numpy as np
import os
import pandas as pd
import pickle
from shutil import copy

dd = pd.read_csv('raftInstall.csv',index_col=0)
rnames = []
for i in dd.index:
    try:
        if np.isnan(dd.rtm[i]):
            print(i,'xxx')
    except TypeError: #when dd.rtm[i] is a str
        rnames.append(i)
        print(i)

        bay = i 
        run = str(int(dd.run[i]))
        if dd.dataset[i]==1:
            srcDir = '/home/bxin/notebooks/raftexplorer/'
        else:
            srcDir = '/home/bxin/notebooks/raftResults2020jan/'

        destDir = '/project/shared/bxin/cam_as_built/%s/'%(bay)
        if not os.path.exists(destDir):
            os.mkdir(destDir)

        filename='RaftRun'+run+'.p'
        pkpath = os.path.join(srcDir,filename)

        if not os.path.isfile('{}/{}'.format(destDir, filename)):
            print('coping %s to %s'%(pkpath, destDir))
            copy(pkpath, destDir)

        f=open(pkpath,'rb')
        raft_name=pickle.load(f)
        res=pickle.load(f)
        ccd_list=pickle.load(f)
        file_list=pickle.load(f)
        fw=pickle.load(f)
        gains=pickle.load(f)
        #print(fw)
        f.close()

        for i in range (0,9):
            filename=str(ccd_list[i][0]+'_QE.fits')
            filepath = os.path.join(srcDir,run,filename)
            #print(filename)

            if not os.path.isfile('{}/{}'.format(destDir, filename)):
                print('coping %s to %s'%(filepath, destDir))
                copy(filepath, destDir)

