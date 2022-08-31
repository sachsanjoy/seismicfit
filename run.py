import os
import input_ini as target
import glob

gridList = glob.glob(target.gridDir+'/*.zip')
if target.combinationMode==True :
    if target.hyperThreading==True:
        ...
    else:
        print("Combination grid search requires hyperthreading enabled!")        
else:
    os.system(target.pythonenv+" grid_search.py")
   