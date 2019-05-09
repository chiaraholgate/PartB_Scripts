# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 09:16:39 2018

This script calculates the climatology of several meteorological variables input into the QIBT model.
The data is from the wrf-based Narclim regional model. 

1979 is not used in the climatology for summer, as the model results begin at 31/1/79, so there's no data
for Dec 1978 or Jan 1979. Other seasons include data from 1979.

Meteo variables:
- u,v wind components > translated to wind vectors. Currenlty we are taking the mean of the bottom 10 model layers.


@author: z3131380
"""
import os
import numpy as np
import pandas
from netCDF4 import Dataset 

dir_in = r'/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/wrf_var_extraction/'
dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/Processed/Climatology/'
region = 'Australia' 
n_i,n_j,n_bt = 144,216,29 # wrf Narclim model dimensions: south_north,west_east,bottom_top

# Date info
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '19801231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pandas.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
full_daylist = pandas.date_range(Start_date,End_date,freq='d')

# Create blank arrays for output
wvcont_mm_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
wvcont_mm_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
wvcont_mm_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
wvcont_mm_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
wvcont_mm_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan

for year in Yearlist:
    start_day = str(year)+'0101' ; end_day = str(year)+'1231'
    daylist = pandas.date_range(start_day,end_day,freq='d')  
    daylist_MAM = []; daylist_JJA = []; daylist_SON = []; daylist_DJF = []; daylist_AprNov = []
    for d in daylist:
        if d.month==3 or d.month==4 or d.month==5:
            daylist_MAM.append(d)
        if d.month==6 or d.month==7 or d.month==8:
            daylist_JJA.append(d)
        if d.month==9 or d.month==10 or d.month==11:
            daylist_SON.append(d)
        if np.logical_and(d.month >=4,d.month<=11):
            daylist_AprNov.append(d)
            
        # DJF must include the Dec dates from the previous year        
        if year > 1979: 
            if d.month==1 or d.month==2:
                daylist_DJF.append(d)
            if d.month==12:
                daylist_DJF.append(d-pandas.DateOffset(years=1))

    daylist_DJF = pandas.to_datetime(daylist_DJF); daylist_DJF = daylist_DJF.sort_values()    
    daylist_MAM = pandas.to_datetime(daylist_MAM) 
    daylist_JJA = pandas.to_datetime(daylist_JJA)
    daylist_SON = pandas.to_datetime(daylist_SON)
    daylist_AprNov = pandas.to_datetime(daylist_AprNov)

    u_DJF = np.zeros([len(daylist_DJF),n_i,n_j],n_bt)*np.nan; v_DJF = np.zeros([len(daylist_DJF),n_i,n_j,n_bt])*np.nan; w_DJF = np.zeros([len(daylist_DJF),n_i,n_j,n_bt])*np.nan;
    WS_DJF = np.zeros([len(daylist_DJF),n_i,n_j,n_bt])*np.nan; WD_DJF = np.zeros([len(daylist_DJF),n_i,n_j,n_bt])*np.nan
    u_MAM = np.zeros([len(daylist_MAM),n_i,n_j,n_bt])*np.nan; v_MAM = np.zeros([len(daylist_MAM),n_i,n_j,n_bt])*np.nan; w_MAM = np.zeros([len(daylist_MAM),n_i,n_j,n_bt])*np.nan;
    WS_MAM = np.zeros([len(daylist_MAM),n_i,n_j,n_bt])*np.nan; WD_MAM = np.zeros([len(daylist_MAM),n_i,n_j,n_bt])*np.nan
    u_JJA = np.zeros([len(daylist_JJA),n_i,n_j,n_bt])*np.nan; v_JJA = np.zeros([len(daylist_JJA),n_i,n_j,n_bt])*np.nan; w_JJA = np.zeros([len(daylist_JJA),n_i,n_j,n_bt])*np.nan;
    WS_JJA = np.zeros([len(daylist_JJA),n_i,n_j,n_bt])*np.nan; WD_JJA = np.zeros([len(daylist_JJA),n_i,n_j,n_bt])*np.nan
    u_SON = np.zeros([len(daylist_SON),n_i,n_j,n_bt])*np.nan; v_SON = np.zeros([len(daylist_SON),n_i,n_j,n_bt])*np.nan; w_SON = np.zeros([len(daylist_SON),n_i,n_j,n_bt])*np.nan; 
    WS_SON = np.zeros([len(daylist_SON),n_i,n_j,n_bt])*np.nan; WD_SON = np.zeros([len(daylist_SON),n_i,n_j,n_bt])*np.nan
    u_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j,n_bt])*np.nan; v_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j,n_bt])*np.nan; w_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j,n_bt])*np.nan; 
    WS_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j,n_bt])*np.nan; WD_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j,n_bt])*np.nan

    file = dir_in+'wrfout_d01_'+str(year)+'_daily_mean_wind.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        U = fh.variables['U'][:,0:10,:,:]    
        V = fh.variables['V'][:,0:10,:,:]   
        W = fh.variables['W'][:,0:10,:,:]       
        XLONG_U = fh.variables['XLONG_U'][:,0:10,:,:]     
        XLAT_U = fh.variables['XLAT_U'][:,0:10,:,:]    
        XLONG_V = fh.variables['XLONG_V'][:,0:10,:,:]    
        XLAT_V = fh.variables['XLAT_V'][:,0:10,:,:]    
        fh.close() 
    # DJF must include the Dec data from the previous year   
    file = dir_in+'wrfout_d01_'+str(year-1)+'_wind.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        U_prevyr = fh.variables['U'][:,0:10,:,:]       
        V_prevyr = fh.variables['V'][:,0:10,:,:]   
        W_prevyr = fh.variables['W'][:,0:10,:,:]      
        fh.close() 
        
    daylist_prevyr = pandas.date_range(str(year-1)+'0101', str(year-1)+'1231',freq='d')    
    for d in daylist_DJF:
        if d.month < 3:
            idx_daylist = np.where(daylist==d)[0][0]
            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]      
            u_DJF[idx_daylist_DJF,:,:,:] = U[idx_daylist,:,:,:]
            v_DJF[idx_daylist_DJF,:,:,:] = V[idx_daylist,:,:,:]
            w_DJF[idx_daylist_DJF,:,:,:] = W[idx_daylist,:,:,:]
        elif d.month==12:
            idx_daylist = np.where(daylist_prevyr==d)[0][0]
            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]   
            u_DJF[idx_daylist_DJF,:,:,:] = U_prevyr[idx_daylist,:,:,:]
            v_DJF[idx_daylist_DJF,:,:,:] = V_prevyr[idx_daylist,:,:,:]
            w_DJF[idx_daylist_DJF,:,:,:] = W_prevyr[idx_daylist,:,:,:]
        for h in n_bt:
            WS_DJF[idx_daylist_DJF,h,:,:] = np.sqrt(u_DJF[idx_daylist_DJF,h,:,:]**2+v_DJF[idx_daylist_DJF,h,:,:]**2)
        
    for d in daylist_MAM:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_MAM = np.where(daylist_MAM==d)[0][0]
        wvcont_mm_MAM[idx_daylist_MAM,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
        wvcont_pct_MAM[idx_daylist_MAM,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
        pre_MAM[idx_daylist_MAM,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_MAM,axis=0)
    wvcont_pct_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_MAM,axis=0)
    pre_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_MAM,axis=0)
        
    for d in daylist_JJA:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_JJA = np.where(daylist_JJA==d)[0][0]
        wvcont_mm_JJA[idx_daylist_JJA,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
        wvcont_pct_JJA[idx_daylist_JJA,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
        pre_JJA[idx_daylist_JJA,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_JJA,axis=0)
    wvcont_pct_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_JJA,axis=0)
    pre_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_JJA,axis=0)
        
    for d in daylist_SON:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_SON = np.where(daylist_SON==d)[0][0]
        wvcont_mm_SON[idx_daylist_SON,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
        wvcont_pct_SON[idx_daylist_SON,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
        pre_SON[idx_daylist_SON,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_SON,axis=0)
    wvcont_pct_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_SON,axis=0)
    pre_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_SON,axis=0)
    
    for d in daylist_AprNov:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_AprNov = np.where(daylist_AprNov==d)[0][0]
        wvcont_mm_AprNov[idx_daylist_AprNov,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
        wvcont_pct_AprNov[idx_daylist_AprNov,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
        pre_AprNov[idx_daylist_AprNov,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_AprNov,axis=0)
    wvcont_pct_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_AprNov,axis=0)
    pre_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_AprNov,axis=0)
        
WS_seasonal_climatology = np.zeros([4,n_i,n_j])*np.nan        
WS_seasonal_climatology[0,:,:,:] = np.nanmean(WS_DJF, axis=0)
WS_seasonal_climatology[1,:,:,:] = np.nanmean(WS_MAM, axis=0)
WS_seasonal_climatology[2,:,:,:] = np.nanmean(WS_JJA, axis=0)
WS_seasonal_climatology[3,:,:,:] = np.nanmean(WS_SON, axis=0)


for s in range(4):
    for h in n_bt:
        for i in n_i:
            for j in n_j:
                WS_cell = WS_seasonal_climatology[s,h,i,j]
   
WS_AprNov_climatology =


WS_seasonal_climatology = np.zeros([4,n_i,n_j])*np.nan 
WS_seasonal_climatology[0,:,:] = 

# Save to netcdf 
ofile = dir_out+region+'_1979-2013_climatology.nc'    
with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
    of.createDimension('i_cross', n_i)
    of.createDimension('j_cross', n_j)    
    of.createDimension('season',4)

    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
    of['latitcrs'].units = 'degrees'
    of['latitcrs'][:] = latitcrs
    
    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
    of['longicrs'].units = 'degrees'
    of['longicrs'][:] = longicrs
    
    of.createVariable('pre_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
    of['pre_seasonal_climatology'].long_name = 'mean seasonal precipitation'
    of['pre_seasonal_climatology'].units = 'mm'
    of['pre_seasonal_climatology'][:] = pre_seasonal_climatology
    
    of.createVariable('wvcont_mm_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
    of['wvcont_mm_seasonal_climatology'].long_name = 'mean seasonal water vapour contribution [mm] to all rainfall in domain'
    of['wvcont_mm_seasonal_climatology'].units = 'mm'
    of['wvcont_mm_seasonal_climatology'][:] = wvcont_mm_seasonal_climatology
    
    of.createVariable('wvcont_pct_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
    of['wvcont_pct_seasonal_climatology'].long_name = 'mean seasonal water vapour contribution [%] to all rainfall in domain'
    of['wvcont_pct_seasonal_climatology'].units = '%'
    of['wvcont_pct_seasonal_climatology'][:] = wvcont_pct_seasonal_climatology
    
    of.createVariable('wvcont_mm_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
    of['wvcont_mm_AprNov_climatology'].long_name = 'mean Apr to Nov water vapour contribution [mm] to all rainfall in domain'
    of['wvcont_mm_AprNov_climatology'].units = 'mm'
    of['wvcont_mm_AprNov_climatology'][:] = wvcont_mm_AprNov_climatology
    
    of.createVariable('wvcont_pct_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
    of['wvcont_pct_AprNov_climatology'].long_name = 'mean Apr to Nov water vapour contribution [%] to all rainfall in domain'
    of['wvcont_pct_AprNov_climatology'].units = '%'
    of['wvcont_pct_AprNov_climatology'][:] = wvcont_pct_AprNov_climatology
    
    of.createVariable('pre_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
    of['pre_AprNov_climatology'].long_name = 'mean Apr to Nov precipitation'
    of['pre_AprNov_climatology'].units = 'mm'
    of['pre_AprNov_climatology'][:] = pre_AprNov_climatology
        
    
        


