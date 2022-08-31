import multiprocessing
import os
import time
import sys

def hyperthread(file_path):
    os.system(file_path)

def readARG():
 try:
  grid_ini = int(sys.argv[1])
  grid_fin = int(sys.argv[2])
  return grid_ini, grid_fin
 except:
  print("For all tracks -> python3 powergrid.py grid_ini grid_fin") 
  exit()

grid_ini, grid_fin = readARG()
start_time = time.time()
processes = []
for i in range(grid_ini,grid_fin):
    file_path = ["python3 sdb_gridsearch_parallel.py "+str(i)]
    process = multiprocessing.Process(target=hyperthread,args=(file_path))
    processes.append(process)
    process.start()
for p in processes:
    p.join()    

print('Time Elapsed : ',time.time()-start_time,'seconds')