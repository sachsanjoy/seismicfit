#python environment
pythonenv = "python3"

#path to grid directory
gridDir = '/mnt/d/sdB-MESA'

#Input parameters of star
targetName = 'B3'
TeffObs = 23660/1000.
TeffObsErr = 360/1000.
loggObs = 5.300
loggObsErr = 0.061

#Gyre period min max range in seconds
GYREPeriodRange = [2000,11500] 

#MAXIMUM GRID CONDITIONS 
zGrid = [0.005,0.035] #Z
miGrid = [1.0,1.8] #Mi
menvGrid = [0.0000, 0.0100] # Menv
heGrid = [0,0.9] #coreHe

#Hyperthreading parameters
gridInitial = 0
gridIteration = 10