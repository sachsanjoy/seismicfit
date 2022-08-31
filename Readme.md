## MODULES
# extract_sdbgrid 
 the module extracts the grid from zipped sdB Grid to a folder sdb_grid.
 read the read me fies inside the module to know about the setup.  
# seismicfit_parallel
 the module fits the input frequency combinations with th eoriginal SdB Grid
# grid_search.py 
 It does single runs on a specific input_frequency combination using the extracted grid information
# grid_search_parallel.py
 under development -> it will run al frequency combination on the extracted grid.

## INPUTS
# input_ini.py 
 update your parameters of the star in this file
# input_frequency.csv
 update your frequecies in this file

## RUN
 the easiest way to run is to run 'run.py' file

## OUTPUT
 Upon completion of the run an output folder with the name of the star will be created.
 it contains a csv file containing all the necessary information about the fitting. 


 THIS MODULE IS FOR RUNNING GRID SEARCH ON THE Grid extract.
 The module has only single frequency combinations enabled at the moment.
 Setup 'input_ini.py' add necesssary parameters
 Add 'input_frequency.csv' with multiplets assigned as '1*' in the 'lmpt' column
 Run 'run.py' file


 Adjust the chunk size based on memory and disk space. default(1000)

