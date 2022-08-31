# Notice:
This module is presented for future refrences, in case if we need to update the grid extract with a new parameter.
Modify this grid_extract.py to do so.
There is no need to run it if you have the grid extract. 
## MODULES
# extract_sdbgrid 
 the medule extracts the grid from zipped sdB Grid to a folder sdb_grid.
 read the read me fies inside the module to know about the setup.  
# seismicfit_parallel
 the module fits the input frequency combinations with th eoriginal SdB Grid
# grid_search.py 
 It does single runs on a specific input_frequency combination using the extracted grid information
# grid_search_parallel.py
 under development -> it will run al frequency combination on the extracted grid.

## INPUTS
# input_ini.py 
 change the input paramters for the run.
 Could use the hyper threading True in input_ini.py file
 Adjust the threading size based on memory and disk space
 default threading 10 => requires 200+ GB temporary space ; 20 ==> 400+ GB
 make sure to update the original sdb-grid path

## RUN
The easiest way to run is to run 'run.py' file
There is option for enabling hyperthreading, it depends on the memory and disk space. 

## OUTPUT
The necessary files will be extracted to sdb_grid folder inside the program directory.
A csv file with all the additional paramters from the model is created in the program directory under name sdbgrid_extract.csv


