# -*- coding: utf-8 -*-
"""
Created on Wed May  1 12:52:51 2019

https://wrf-python.readthedocs.io/en/latest/user_api/generated/wrf.cape_2d.html

@author: z3131380
"""

from netCDF4 import Dataset
import wrf
import numpy as np
import pandas as pd

#==============================================================================
# Definitions
#==============================================================================
dir_data1 = r'/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/'
dir_data2 = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'
dir_out =dir_data2+'CAPE/'

#==============================================================================
# Dates
#==============================================================================
Start_date = '19790108' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '19790108' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Daylist = pd.date_range(Start_date,End_date,freq='d') # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
   
for d in Daylist:
    yyyy = str(d.year)
    mm = '%02u' % d.month
    dd = '%02u' % d.day
     
    #==============================================================================
    # Load data 
    #==============================================================================
    fh = Dataset(dir_data1+'wrfout_d01_'+yyyy+'-'+mm+'-'+dd+'_00:00:00', mode='r') 
    HGT = fh.variables['HGT'][:] #m [144,215]
    QVAPOR = fh.variables['QVAPOR'][:] #kg/kg [8,29,144,215]
    QCLOUD = fh.variables['QCLOUD'][:] #kg/kg [8,29,144,215]
    QICE =fh.variables['QICE'][:] #kg/kg [8,29,144,215]
    QSNOW = fh.variables['QSNOW'][:] #kg/kg [8,29,144,215]
    QRAIN = fh.variables['QRAIN'][:] #kg/kg [8,29,144,215]
    P = fh.variables['P'][:] #Pa [8,29,144,215]
    PB = fh.variables['PB'][:] #Pa [8,29,144,215]
    PH = fh.variables['PH'][:,1:,:,:] #m2/s2 [8,30,144,215]
    PHB = fh.variables['PHB'][:,1:,:,:] #m2/s2 [8,30,144,215]
    PTOP = int(fh.variables['P_TOP'][0]) #Pa [8]
    T = fh.variables['T'][:] + 300 #K [8,29,144,215]
    XLAT = fh.variables['XLAT'][:] # [144,215]
    XLONG = fh.variables['XLONG'][:] # [144,215]
    fh.close()
    
    fh = Dataset(dir_data2+'wrfhrly_d01_'+yyyy+'-'+mm+'-'+dd+'_00:00:00_PSFC.nc', mode='r')
    psfc_hpa =fh.variables['PSFC'][:]/100  #Pa [8,144,215]
    fh.close()
    
    pres_hpa = (P + PB)/100 
    terrain = np.repeat(HGT[np.newaxis,:,:],8,axis=0) # Needs to have same shape as pres_hpa (except vertical)
    qv = QVAPOR + QCLOUD + QICE + QSNOW + QRAIN
    gph = ( PH + PHB ) / 9.81
    tkel = wrf.tk((P+PB),T,meta=False)
    # Check inputs don't have NaNs with np.isnan(var).any()
    
    #==============================================================================
    # Calculation
    #==============================================================================
    out = wrf.cape_2d(pres_hpa,tkel,qv, gph, terrain, psfc_hpa, ter_follow=True, meta=False)
    CAPE = out[0,:,:] # J/kg
    CIN = out[1,:,:] # J/kg
    LCL = out[2,:,:] # m
    LFC = out[3,:,:] # m
    
    #==============================================================================
    # Save to netcdf
    #==============================================================================
    ofile = dir_out+'CAPE_'+yyyy+'_'+mm+'_'+dd+'.nc'    
    with Dataset(ofile, 'w') as of: 
        of.createDimension('i_cross', CAPE.shape[1])
        of.createDimension('j_cross', CAPE.shape[2])  
        of.createDimension('time', CAPE.shape[0])
        
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
        
        of.createVariable('CAPE', 'f4', ('time','i_cross', 'j_cross'))
        of['CAPE'].long_name = 'Maximum CAPE'
        of['CAPE'].description = 'Accumulated buoyant energy from the level of free convection (LFC) to the equilibrium level (EL)'
        of['CAPE'].units = 'J/kg'
        of['CAPE'][:] = CAPE
        
        of.createVariable('CIN', 'f4', ('time','i_cross', 'j_cross'))
        of['CIN'].long_name = 'Maximum CIN'
        of['CIN'].description = 'Accumulated negative buoyant energy from the parcel starting point to the LFC'
        of['CIN'].units = 'J/kg'
        of['CIN'][:] = CIN
        
        of.createVariable('LCL', 'f4', ('time','i_cross', 'j_cross'))
        of['LCL'].long_name = 'Lifting condensation level'
        of['LCL'].units = 'm'
        of['LCL'][:] = LCL
        
        of.createVariable('LFC', 'f4', ('time','i_cross', 'j_cross'))
        of['LFC'].long_name = 'Level of free convection'
        of['LFC'].units = 'm'
        of['LFC'][:] = LFC