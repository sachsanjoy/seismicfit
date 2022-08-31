import os
import input_ini as target
import numpy as np
import glob

gridList = glob.glob(target.gridDir+'/*.zip')

gridini = target.gridInitial
griditer = target.gridIteration
run = np.zeros(1)
while gridini < len(gridList):
    if gridini+griditer>len(gridList):
        run0 = target.pythonenv+" powergrid.py "+str(gridini)+" "+str(len(gridList))+" "
    else:    
        run0 = target.pythonenv+" powergrid.py "+str(gridini)+" "+str(gridini+griditer)+" "
    run = np.vstack((run,run0))
    gridini = gridini+griditer
np.savetxt('run.go',run[1:],fmt="%s")
os.system("./run.go")
