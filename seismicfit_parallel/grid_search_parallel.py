#sdB grid search
import sys
import glob
import pandas as pd
import input_ini as target
import numpy as np
import os

import timeit
import warnings
warnings.filterwarnings("ignore")    
  
'''A pack of merit fuctions'''
def merit_fn_p1(table_gyre,f_obs):
    table_gyre_l1 = table_gyre[table_gyre.l==1]
    ds = np.zeros(1)
    for i in f_obs.frequency:
        #dP = (table_gyre_l1['fq']-f_obs[i])**2. #frequency difference
        dP = (table_gyre_l1['P']-(1.0e6/i))**2 #period difference
        dP = min(dP)
        ds = np.vstack((ds,dP))
    return ds[1:]

'''Reading messa models'''
def messa_read(gridId,track,coreHe_ab):
    path = str(gridId)+'_tempzip/'+track+'/history.data'
    #print(path)
    df=pd.read_csv(path,skiprows=5,sep=r"\s+")
    logteff=df['log_Teff']
    logg=df['log_g'] 
    age=df['star_age']
    he_ab=df['center_he4']
    #finding model number
    He_data_path=str(gridId)+'_tempzip/'+track+'/custom_He'+str(coreHe_ab)+'.data'
    #print(He_data_path)
    He_data=pd.read_csv(He_data_path,nrows=2,sep=r"\s+")
    He_data_header = He_data.iloc[0] 
    He_data = He_data[1:]
    He_data.columns = He_data_header
    model_number=np.array(He_data['model_number'])[0]
    log_Teff_model=df['log_Teff'][df['model_number']==int(model_number)]
    logg_model=df['log_g'][df['model_number']==int(model_number)]
    starage = df['star_age'][df['model_number']==int(model_number)]
    starmass = df['star_mass'][df['model_number']==int(model_number)]
    starradius = df['radius'][df['model_number']==int(model_number)]
    log_L = df['log_L'][df['model_number']==int(model_number)]  
    luminosity = 10.0 ** log_L
    index_zaehb=df['model_number'].index[df['mass_conv_core'] >0.0][0]-1
    index_taehb=df['model_number'].index[df['mass_conv_core'] >0.0][-1]-1
    
    #Breating pulse choping
    he_ab=he_ab[index_zaehb:index_taehb]
    age=age[index_zaehb:index_taehb]
    he_ab=np.array(he_ab)
    age=np.array(age)
    d_he_ab=he_ab[:(len(he_ab)-1)]-he_ab[1:]
    d_age=age[:len(age)-1]-age[1:]
    slope=d_he_ab/d_age
    index=np.where(slope>=0.50e-7)
    fac=np.min(index)
    
    Teff_model = (10.0 ** log_Teff_model) / 1000.0
    chisq=(((target.TeffObs-Teff_model)/target.TeffObsErr)**2)+(((target.loggObs-logg_model)/target.loggObsErr)**2) 
    
    return chisq,Teff_model, logg_model,starage,starmass,starradius,luminosity

'''Formating GYRE csv tables for matching with oservations'''
def tableformat(df,star,gridId):
    table='f,fq,P,l,n'
    for i in range(len(df)):
        #print('k',i)
        l = df['l'][i]
        fq = (df['Re(freq)'][i]/86400.0)*1e6
        P = 86400.0/df['Re(freq)'][i]
        n = -1.0*df['n_pg'][i]
        if (P>=target.GYREPeriodRange[0] and P<target.GYREPeriodRange[1]):
                table0='f'+str(i)+','+str(fq)+','+str(P)+','+str(l)+','+str(n)
                table=np.vstack((table,table0))
    np.savetxt(star+'_output/gyre_table_full'+str(gridId)+'.csv',table,'%s')


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

def unzipy(zip_name,folder_name):
   os.system('unzip '+zip_name+' -d '+folder_name)

def append_task_csv(star,gridi):
    import glob
    all_files = glob.glob('outputcsv/'+star+'_fitgrid_'+str(gridi)+'_*.csv')
    all_files = sorted(all_files)
    print(all_files)
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)
    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv(star+'_output/'+star+'_fitgrid_'+str(gridi)+'.csv')

def readARG():
 try:
  gridId = int(sys.argv[1])
  return gridId
 except:
  print("For all tracks -> python3 merit.py star_name all gridID=0") 
  exit()

gridId = readARG()
if not os.path.exists(target.targetName+'_output'):
    os.system('mkdir '+target.targetName+'_output')

if not os.path.exists(str(gridId)+'_tempzip'):
    os.system('mkdir '+str(gridId)+'_tempzip')

if not os.path.exists('outputcsv'):
    os.system('mkdir outputcsv')

gridList = glob.glob(target.gridDir+'/*.zip')
gridList = sorted(gridList)
print(gridList)
print('Total grid length: ', len(gridList))
print('Approximate number of models : ', len(gridList)*32*17)
#target frequency
tableObs = pd.read_csv('input_frequency.csv')
fmpt = tableObs[tableObs.lmpt=='1*'] # taking only multiplet middle component
f1 = fmpt.fn.values
fqcombList = fqCombination(f1)
print('Total number of l1 multiplets : ',len(fmpt))

for gridi in range(gridId,gridId+1): 
    os.system('rm -r '+str(gridId)+'_tempzip ; mkdir '+str(gridId)+'_tempzip')
    gl = gridList[gridi]
    gridZ = float(gl[gl.find('_z')+2:gl.find('_lvl')])
    gridMi = float(gl[gl.find('mi')+2:gl.find('_z')])
    if ((gridZ>=target.zGrid[0]) & (gridZ<=target.zGrid[1]) & (gridMi>=target.miGrid[0]) & (gridMi<=target.miGrid[1])): #grid filtering
        unzipy(gridList[gridi],str(gridId)+'_tempzip')
        gridName=glob.glob(str(gridId)+"_tempzip/log*")
        gridPath, gridName = os.path.split(gridName[0])
        #grid update
        allTrack = glob.glob(str(gridId)+'_tempzip/'+gridName+'/*')  
        for i in range(0,len(allTrack)):
            a = 'Track,mi,menv,z,y,yc,mass,radius,age,teff,logg,luminosity,chisq,total_merit,nfq'
            for jj1 in fmpt.fn:
                a = a+','+jj1
            m_log = [a]
            track=allTrack[i]
            track=track[track.find('zip/')+4:]
            menv=track[track.find('menv')+4:track.find('_rot')]
            z=track[track.find('_z')+2:track.find('_lvl')]
            y=track[track.find('_y')+2:track.find('_f')]
            mi = track[track.find('mi')+2:track.find('_z')]
            if ((float(menv)>=target.menvGrid[0]) & (float(menv)<=target.menvGrid[1])):
                He_files = glob.glob(str(gridId)+"_tempzip/"+track+"/*summary.txt")
                start = timeit.default_timer()
                for j in range(0,len(He_files)):
                    print('Progress : '+'('+str(j)+'/'+str(len(He_files))+')---|---('+str(i)+'/'+str(len(allTrack))+')')
                    coreHe_ab = He_files[j]
                    coreHe_ab =   coreHe_ab[coreHe_ab.find('He')+2:coreHe_ab.find('sum')-1]
                    coreHe_ab = float(coreHe_ab)
                    if ((float(coreHe_ab)>=target.heGrid[0]) & (float(coreHe_ab)<=target.heGrid[1])):
                        pfile = str(gridId)+"_tempzip/"+track+"/custom_He"+str(coreHe_ab)+"_summary.txt"
                        #print(pfile)
                        #print(star+' '+track+' He '+str(coreHe_ab))
                        df=pd.read_csv(pfile,skiprows=5,sep=r"\s+")
                        titlename = target.targetName+'_'+track[track.find('_lvl')+6:]+'_He'+str(format(coreHe_ab,'.2f'))
                        chisq,Teff_model, logg_model,starage,starmass,starradius,luminosity = messa_read(gridId,track,coreHe_ab)
                        chisq = str(np.round(np.min(chisq),3))
                        starmass = str(np.round(np.min(starmass),6))
                        starradius =str(np.round(np.min(starradius),6))
                        starage = str(np.round(np.min(starage),6))
                        Teff_model =str(np.round(np.min(Teff_model),6))
                        logg_model=str(np.round(np.min(logg_model),6))
                        luminosity = str(np.round(np.min(luminosity),6))

                        tableformat(df,target.targetName,gridId)
                        table_gyre = pd.read_csv(target.targetName+'_output/gyre_table_full'+str(gridId)+'.csv')
                        for i1 in range(len(fqcombList)):
                            for i2 in range(len(fqcombList[i1])):
                                fmptComb = fmpt.apply(lambda row: row[fmpt['fn'].isin(fqcombList[i1][i2])])
                                ds = merit_fn_p1(table_gyre,fmptComb)    
                                fmptComb['fqmerit'] = ds
                                fmpt['fqmerit'] = -np.ones(len(fmpt))
                                for ifq in fmptComb.index:
                                    fmpt.fqmerit[ifq] = fmptComb.fqmerit[ifq]
                                nfq = len(fmptComb)
                                merit = str(np.round(np.sum(ds)/len(fmptComb),6))
                                b = track +','+mi+','+menv+','+z+','+y+','+str(coreHe_ab)+','+starmass+','+starradius+','+starage+','+Teff_model+','+logg_model+','+luminosity+','+chisq+','+merit+','+nfq
                                for jj2 in fmpt.fqmerit:
                                    b = b + ',' +str(jj2)
                                m_log0 = np.array([b])
                                m_log=np.vstack((m_log,m_log0))
                np.savetxt('outputcsv/'+target.targetName+'_fitgrid_'+str(gridi)+'_menv_'+str(menv)+'.csv',m_log,fmt="%s")
                stop = timeit.default_timer()
                print('Run Time: ', stop - start)
os.system('rm -r '+str(gridId)+'_tempzip')# ; mkdir '+str(gridId)+'_tempzip')
append_task_csv(target.targetName,gridi)
#os.system('rm -f outputcsv/'+target.targetName+'_fitgrid_'+str(gridi)+'*.csv')
#np.savetxt(target.targetName+'_output/'+target.targetName+'_fitgrid_'+str(i)+'.csv',meritlog0,fmt="%s")
