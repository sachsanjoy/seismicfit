import os
import input_ini as target
import numpy as np
import glob
from tqdm import tqdm
import pandas as pd

def append_task_csv():
    '''Append the csv chunk files'''
    all_files = glob.glob('sdb_grid/logs_*.csv')
    all_files = sorted(all_files)
    print(all_files)
    li = []
    for filename in tqdm(all_files):
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)
    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv('sdb_grid_extract.csv')

gridList = glob.glob(target.gridDir+'/*.zip')
if target.hyperthread==True:
    gridini = target.gridInitial
    griditer = target.gridIteration
    run = np.zeros(1)
    while gridini < len(gridList):
        if gridini+griditer>len(gridList):
            run0 = "python3 powergrid.py "+str(gridini)+" "+str(len(gridList))+" "
        else:    
            run0 = "python3 powergrid.py "+str(gridini)+" "+str(gridini+griditer)+" "
        run = np.vstack((run,run0))
        gridini = gridini+griditer
    np.savetxt('run.go',run[1:],fmt="%s")
    os.system("./run.go")
    append_task_csv()
else:
    os.system("python3 grid_extract.py")
