#sdB grid search
import glob
import pandas as pd
import input_ini as target
import numpy as np
import os
from tqdm import tqdm
import timeit
import warnings
warnings.filterwarnings("ignore")    
  
'''A pack of merit fuctions'''
def merit_fn_p1(table_gyre,f_obs,lval):
    table_gyre_l1 = table_gyre[table_gyre.l==int(lval)]
    ds = np.zeros(1)
    for i in f_obs.frequency:
        #dP = (table_gyre_l1['fq']-f_obs[i])**2. #frequency difference
        dP = (table_gyre_l1['P']-(1.0e6/i))**2 #period difference square
        dP = min(dP)
        ds = np.vstack((ds,dP))
    return ds[1:]

def append_task_csv(star):
    all_files = glob.glob(star+'_output/'+star+'_fitgrid_*.csv')
    all_files = sorted(all_files)
    print(all_files)
    li = []
    for filename in tqdm(all_files):
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)
    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv(star+'_output/'+star+'_seismicfit.csv',index=None)

def grid_search(total_grid,fmpt):
    a = 'Track,mi,menv,z,y,yc,mass,radius,age,teff,logg,luminosity,chisq,total_merit'
    for jj1 in fmpt.fn:
        a = a+','+jj1
    m_log = [a]
    chuck, chuck_name = 0, 0
    start = timeit.default_timer()
    for track,yc in tqdm(zip(total_grid.Track,total_grid.yc)): 
        #grid_properties
        menv=track[track.find('menv')+4:track.find('_rot')]
        z=track[track.find('_z')+2:track.find('_lvl')]
        y=track[track.find('_y')+2:track.find('_f')]
        mi = track[track.find('mi')+2:track.find('_z')]
        trc = total_grid[(total_grid.Track==track) & (total_grid.yc==yc)]
        #stellar_parameters
        starmass = str(np.round(trc.mass.values[0],6))
        starradius =str(np.round(trc.radius.values[0],6))
        starage = str(np.round(np.min(trc.age.values[0]),6))
        Teff_model =str(np.round(np.min(trc.teff.values[0]),6))
        logg_model=str(np.round(np.min(trc.logg.values[0]),6))
        luminosity = str(np.round(np.min(trc.luminosity.values[0]),6))
        #HRD_chisquarefit
        chisq = (((target.TeffObs-float(Teff_model))/target.TeffObsErr)**2)+(((target.loggObs-float(logg_model))/target.loggObsErr)**2) 
        chisq = str(np.round(np.min(chisq),3))
        #gyre_period_merit_fit
        path_to_gyre = 'sdb_grid/'+track+'/custom_He'+str(yc)+'_summary.csv'
        gyre_table = pd.read_csv(path_to_gyre)
        gyre_table = gyre_table[(gyre_table.P>=target.GYREPeriodRange[0]) & (gyre_table.P<target.GYREPeriodRange[1])]
        ds = merit_fn_p1(gyre_table,fmpt,1)    
        fmpt['fqmerit'] = ds     
        merit = str(np.round(np.sqrt(np.sum(ds))/len(fmpt),6))
        #making_data
        b = track +','+mi+','+menv+','+z+','+y+','+str(yc)+','+starmass+','+starradius+','+starage+','+Teff_model+','+logg_model+','+luminosity+','+chisq+','+merit
        for jj2 in fmpt.fqmerit:
            b = b + ',' +str(jj2)
        m_log0 = np.array([b])
        m_log=np.vstack((m_log,m_log0))
        #chunking and saving
        chuck+=1
        if chuck==target.chunkSize :
            np.savetxt(target.targetName+'_output/'+target.targetName+'_fitgrid_'+str(chuck_name)+'.csv',m_log,fmt="%s")        
            a = 'Track,mi,menv,z,y,yc,mass,radius,age,teff,logg,luminosity,chisq,total_merit,nfq'
            for jj1 in fmpt.fn:
                a = a+','+jj1
            m_log = [a]        
            chuck_name +=1
            chuck = 0
    np.savetxt(target.targetName+'_output/'+target.targetName+'_fitgrid_'+str(chuck_name)+'.csv',m_log,fmt="%s")        
    stop = timeit.default_timer()
    print('Run Time: ', stop - start)  
    append_task_csv(target.targetName)
    os.system('rm -f '+target.targetName+'_output/'+target.targetName+'_fitgrid_*.csv')

def n_length_combo(lst, n):
  '''recursive function for combinations'''
  if n == 0:
      return [[]]
  lst1 =[]
  for i in range(0, len(lst)):
      m = lst[i]
      remLst = lst[i + 1:]
      remainlst_combo = n_length_combo(remLst, n-1)
      for p in remainlst_combo:
            lst1.append([m, *p])       
  return lst1

def fqCombination(f1):
  '''all frequency combinations'''
  lstall=[]
  for n in range(3,len(f1)+1):
    lst1 = n_length_combo(f1,n)
    lstall.append(lst1)
  return lstall[0:]

#####################################################

if not os.path.exists(target.targetName+'_output'):
    os.system('mkdir '+target.targetName+'_output')

#loading grid extract
total_grid = pd.read_csv('sdb_grid/sdb_grid_extract.csv')
print('Approximate number of iteration models : ', len(total_grid))

#target frequency
tableObs = pd.read_csv('input_frequency.csv')
fmpt = tableObs[tableObs.lmpt=='1*'] # taking only confirmed l1 multiplet middle components or l1 modes

print('Total number of l1 modes for fitting : ',len(fmpt))

if target.combinationMode == True:
    f1 = fmpt.fn.values
    fqcombList = fqCombination(f1)
    for i1 in range(len(fqcombList)):
        for i2 in range(len(fqcombList[i1])):
            fmptComb = fmpt.apply(lambda row: row[fmpt['fn'].isin(fqcombList[i1][i2])])
            grid_search(total_grid,fmpt)
else:
    grid_search(total_grid,fmpt)