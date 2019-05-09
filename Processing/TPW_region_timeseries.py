# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 15:25:58 2019

This script calculates a daily timeseries of NARCliM total precipitable water over each region.

TPW unit check:
TPW[mm] = dp [Pa] * mix [kg/kg] / g [m/s2] /rho [1000kg/m3]
                = kg/ms2 * kg/kg * s2/m * m3/1000kg * 1000mm
                = mm OK.


@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import pandas as pd

#==============================================================================
# Definitions
#==============================================================================
dir_data1 = r'/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_data2 = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'

g = 9.8     #gravity (m.s-2)

#==============================================================================
# Dates
#==============================================================================
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131230' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Daylist = pd.date_range(Start_date,End_date,freq='d') # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.

#==============================================================================
# Calculate TPW and Integrate 
#==============================================================================
for d in Daylist:
    yyyy = str(d.year)
    mm = '%02u' % d.month
    dd = '%02u' % d.day

    # Load data 
    fh = Dataset(dir_data1+'wrfout_d01_'+yyyy+'-'+mm+'-'+mm+'_00:00:00', mode='r') 
    QVAPOR = fh.variables['QVAPOR'][:] #kg/kg [8,29,144,215]
    QCLOUD = fh.variables['QCLOUD'][:] #kg/kg [8,29,144,215]
    QICE = fh.variables['QICE'][:] #kg/kg [8,29,144,215]
    QSNOW = fh.variables['QSNOW'][:] #kg/kg [8,29,144,215]
    QRAIN = fh.variables['QRAIN'][:] #kg/kg [8,29,144,215]
    P = fh.variables['P'][:] #Pa [8,29,144,215]
    PB = fh.variables['PB'][:] #Pa [8,29,144,215]
    PTOP = int(fh.variables['P_TOP'][0]) #Pa [8]
    XLAT = fh.variables['XLAT'][:] # [144,215]
    XLONG = fh.variables['XLONG'][:] # [144,215]
    fh.close()
    
    fh = Dataset(dir_data2+'wrfhrly_d01_'+yyyy+'-'+mm+'-'+dd+'_00:00:00_PSFC.nc', mode='r')
    PSFC = fh.variables['PSFC'][:]  #Pa [8,144,215]
    fh.close()
    
    pres = P + PB
    mix = QVAPOR + QCLOUD + QICE + QSNOW + QRAIN
    
    # Calculate dP, pw then TPW for each cell, level and time
    # Note that model levels ground > highest = 0>28. QIBT expects otherway, so pressure data reversed
    # when loaded into QIBT. No such reversal here. 
    dp = np.zeros_like(P)
    pw = np.zeros_like(P)
    TPW = np.zeros_like(PSFC)
    for t in range(QVAPOR.shape[0]):
        for i in range(QVAPOR.shape[2]):
            for j in range(QVAPOR.shape[3]):
                for k in range(0,29):
                    if k==28:
                        dp[t,28,i,j] = (pres[t,27,i,j] + pres[t,28,i,j])/2 - PTOP # Ground level
                    elif k==0:
                        dp[t,0,i,j] = PSFC[t,i,j] - (pres[t,0,i,j] + pres[t,1,i,j])/2 # top level
                    else:
                        dp[t,k,i,j] = (pres[t,k-1,i,j] - pres[t,k+1,i,j])/2 # Middle levels
                    pw[t,k,i,j] = (mix[t,k,i,j]*dp[t,k,i,j])/g
                TPW[t,i,j] = np.nansum(pw[t,:,i,j])
                
    
    ofile = dir_data2+'TPW/TPW_'+yyyy+'_'+mm+'_'+dd+'.nc'    
    with Dataset(ofile, 'w') as of: 
        of.createDimension('i_cross', TPW.shape[1])
        of.createDimension('j_cross', TPW.shape[2])  
        of.createDimension('time', TPW.shape[0])
        
        of.createVariable('time', 'f4', ('time'))
        of['time'].long_name = 'Three-hourly period of day, starting at 00:00 and ending 21:00'
        of['time'].units = 'hours'
        of['time'][:] = range(8)
    
        of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
        of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
        of['latitcrs'].units = 'degrees'
        of['latitcrs'][:] = XLAT
        
        of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
        of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
        of['longicrs'].units = 'degrees'
        of['longicrs'][:] = XLONG
        
        of.createVariable('TPW', 'f4', ('time','i_cross', 'j_cross'))
        of['TPW'].long_name = 'Total precipitable water'
        of['TPW'].units = 'mm'
        of['TPW'][:] = TPW

