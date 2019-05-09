# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 12:59:27 2019

This script determines the mean source areas for rainfall falling within a predefined
anomalous period. It is based on the script "Climatology_wvcont_region.py". It sums the water vapour
contribution depth over the whole time period and then takes its aver per year, and takes the average 
of the vapour % contribution. So all values are an average anomaly per year of the anom period.


It uses model output that has been processed by YEAR. i.e. each cell's mm or % contribution to 
rainfall falling anywhere in the region over the course of the year.

**TO DO:
-caregul selection of anom years
- This script STILL NEEDS TO BE CHECKED!

@author: z3131380
"""
import os
import numpy as np
from netCDF4 import Dataset 
import pandas

dir_data = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/'
dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Anom_periods/'
region = 'Australia' 
anom_period = 'MD'
n_i,n_j = 134,205 # QIBT model dimensions

#==============================================================================
# Anomaly periods
#==============================================================================
if anom_period == 'MD': # Millenium drought + some wet years
    Yearlist = np.array([1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011])
elif anom_period == 'LN': # Big La Nina years
    Yearlist = np.array([1989,1999,2000,2008,2009,2010,2011])
elif anom_period == 'EN': # Big El Nino years 
    Yearlist = np.array([1982,1987,1994,1998,2002])

#==============================================================================
# Load climatology
#==============================================================================
file = dir_data+'Climatology/'+region+'_1979-2013_seasonal_climatology.nc'
if os.path.isfile(file) == True:
    fh = Dataset(file, mode='r') 
    wvcont_mm_annual_climatology = fh.variables['wvcont_mm_annual_climatology'][:]    
    wvcont_mm_seasonal_climatology = fh.variables['wvcont_mm_seasonal_climatology'][:]  
    wvcont_pct_annual_climatology = fh.variables['wvcont_pct_annual_climatology'][:] 
    wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:] 
    pre_annual_climatology = fh.variables['pre_annual_climatology'][:]  
    pre_seasonal_climatology = fh.variables['pre_seasonal_climatology'][:]  
    fh.close() 
    
#==============================================================================
# ANNUAL ANOMALIES
#==============================================================================
wvcont_mm_anom_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan 
wvcont_pct_anom_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan 
pre_anom_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
for year in Yearlist:
    file = dir_data+'Yearly/'+region+'_'+str(year)+'_wvcont.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        wvcont_sum_yearly_mm = fh.variables['wv_cont_sum_yearly_mm'][:]    
        wvcont_sum_yearly_pct = fh.variables['wv_cont_sum_yearly_pct'][:] 
        pre_annual_total = fh.variables['pre_annual_total'][:]  
        latitcrs = fh.variables['latitcrs'][:]  
        longicrs = fh.variables['longicrs'][:]  
        fh.close() 
    # Make zoers nans, so that when we cal the anomaly, we don't have to deal with unreal zeros
    wvcont_sum_yearly_mm[wvcont_sum_yearly_mm==0] = np.nan
    wvcont_sum_yearly_pct[wvcont_sum_yearly_pct==0] = np.nan
    pre_annual_total[pre_annual_total==0] = np.nan

    wvcont_mm_anom_totals[np.where(Yearlist==year)[0][0],:,:] = wvcont_sum_yearly_mm
    wvcont_pct_anom_totals[np.where(Yearlist==year)[0][0],:,:] = wvcont_sum_yearly_pct
    pre_anom_totals[np.where(Yearlist==year)[0][0],:,:] = pre_annual_total
    
    #Find difference from climatology
    wvcont_mm_departure = wvcont_mm_anom_totals - wvcont_mm_annual_climatology
    wvcont_pct_departure = wvcont_pct_anom_totals - wvcont_pct_annual_climatology
    pre_mm_departure =  pre_anom_totals -  pre_annual_climatology 

ofile = dir_out+region+'_annual_anom_MilDrought.nc'    
with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
    of.createDimension('i_cross', n_i)
    of.createDimension('j_cross', n_j)    
    of.createDimension('time',len(Yearlist))

    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
    of['latitcrs'].units = 'degrees'
    of['latitcrs'][:] = latitcrs
    
    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
    of['longicrs'].units = 'degrees'
    of['longicrs'][:] = longicrs
    
    of.createVariable('time', 'f4', ('time'))
    of['time'].long_name = 'Year of anom period'
    of['time'].units = 'years'
    of['time'][:] = Yearlist
      
    of.createVariable('wvcont_mm_annual_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['wvcont_mm_annual_departure'].long_name = 'anomalous yearly water vapour contribution [mm] to all rainfall in domain during year, wrt climatology'
    of['wvcont_mm_annual_departure'].units = 'mm'
    of['wvcont_mm_annual_departure'][:] = wvcont_mm_departure
    
    of.createVariable('wvcont_pct_annual_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['wvcont_pct_annual_departure'].long_name = 'anomalous yearly water vapour contribution [%] to all rainfall in domain during year, wrt climatology'
    of['wvcont_pct_annual_departure'].units = '%'
    of['wvcont_pct_annual_departure'][:] = wvcont_pct_departure
    
    of.createVariable('pre_mm_annual_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['pre_mm_annual_departure'].long_name = 'anomalous yearly rainfall [mm] wrt climatology'
    of['pre_mm_annual_departure'].units = '%'
    of['pre_mm_annual_departure'][:] = pre_mm_departure
    
    of.createVariable('pre_annual_mm', 'f4', ('time','i_cross', 'j_cross'))
    of['pre_annual_mm'].long_name = 'Annual rainfall [mm]'
    of['pre_annual_mm'].units = '%'
    of['pre_annual_mm'][:] = pre_anom_totals
    
#==============================================================================
# SEASONAL ANOMALIES
#==============================================================================
wvcont_mm_seasonal_departure = np.zeros([len(Yearlist)*4,n_i,n_j])*np.nan
pre_seasonal_departure = np.zeros([len(Yearlist)*4,n_i,n_j])*np.nan
wvcont_pct_seasonal_departure = np.zeros([len(Yearlist)*4,n_i,n_j])*np.nan

seaslist = []
for y in Yearlist:
    seaslist.append(y+0)
    seaslist.append(y+0.25)
    seaslist.append(y+0.5) 
    seaslist.append(y+0.75)

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
        #if year > Yearlist.min(): 
        if d.month==12:
            daylist_DJF.append(d-pandas.DateOffset(years=1))
        if d.month==1 or d.month==2:
            daylist_DJF.append(d)
        

    daylist_DJF = pandas.to_datetime(daylist_DJF); daylist_DJF = daylist_DJF.sort_values()    
    daylist_MAM = pandas.to_datetime(daylist_MAM) 
    daylist_JJA = pandas.to_datetime(daylist_JJA)
    daylist_SON = pandas.to_datetime(daylist_SON)
    daylist_AprNov = pandas.to_datetime(daylist_AprNov)
    
    wv_cont_sum_daily_mm_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan; pre_daily_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan
    wv_cont_sum_daily_mm_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan; pre_daily_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan
    wv_cont_sum_daily_mm_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan; pre_daily_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan
    wv_cont_sum_daily_mm_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan; pre_daily_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan
    
    fname = region+'_'+str(year)+'_wvcont.nc'
    file = dir_data+'/Yearly/'+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        latitcrs = fh.variables['latitcrs'][:] 
        longicrs = fh.variables['longicrs'][:] 
        wv_cont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]  
        pre = fh.variables['pre'][:] 
        fh.close()    

    # DJF must include the Dec data from the previous year   
    fname = region+'_'+str(year-1)+'_wvcont.nc'
    file = dir_data+'/Yearly/'+fname
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        wv_cont_sum_daily_mm_prevyr = fh.variables['wv_cont_sum_daily_mm'][:] 
        pre_prevyr = fh.variables['pre'][:]  
        fh.close()
        
    ## DJF ##
    daylist_prevyr = pandas.date_range(str(year-1)+'0101', str(year-1)+'1231',freq='d')
    for d in daylist_DJF:
        if d.month < 3:
            idx_daylist = np.where(daylist==d)[0][0]
            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]
            wv_cont_sum_daily_mm_DJF[idx_daylist_DJF,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
            pre_daily_DJF[idx_daylist_DJF,:,:] = pre[idx_daylist,:,:]
        elif d.month==12:
            idx_daylist = np.where(daylist_prevyr==d)[0][0]
            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]
            wv_cont_sum_daily_mm_DJF[idx_daylist_DJF,:,:] = wv_cont_sum_daily_mm_prevyr[idx_daylist,:,:]
            pre_daily_DJF[idx_daylist_DJF,:,:] = pre_prevyr[idx_daylist,:,:]
    wvcont_mm_DJF = np.nansum(wv_cont_sum_daily_mm_DJF,axis=0)
    pre_DJF = np.nansum(pre_daily_DJF, axis=0)
    wvcont_pct_DJF = wvcont_mm_DJF/np.nansum(pre_DJF)*100
    
    ## MAM ##                
    for d in daylist_MAM:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_MAM = np.where(daylist_MAM==d)[0][0]
        wv_cont_sum_daily_mm_MAM[idx_daylist_MAM,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
        pre_daily_MAM[idx_daylist_MAM,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_MAM = np.nansum(wv_cont_sum_daily_mm_MAM,axis=0)
    pre_MAM = np.nansum(pre_daily_MAM, axis=0)
    wvcont_pct_MAM = wvcont_mm_MAM/np.nansum(pre_MAM)*100
    
    ## JJA ##    
    for d in daylist_JJA:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_JJA = np.where(daylist_JJA==d)[0][0]
        wv_cont_sum_daily_mm_JJA[idx_daylist_JJA,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
        pre_daily_JJA[idx_daylist_JJA,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_JJA = np.nansum(wv_cont_sum_daily_mm_JJA,axis=0)
    pre_JJA = np.nansum(pre_daily_JJA, axis=0)
    wvcont_pct_JJA = wvcont_mm_JJA/np.nansum(pre_JJA)*100
    
    ## SON ##    
    for d in daylist_SON:
        idx_daylist = np.where(daylist==d)[0][0]
        idx_daylist_SON = np.where(daylist_SON==d)[0][0]
        wv_cont_sum_daily_mm_SON[idx_daylist_SON,:,:] = wv_cont_sum_daily_mm[idx_daylist,:,:]
        pre_daily_SON[idx_daylist_SON,:,:] = pre[idx_daylist,:,:]
    wvcont_mm_SON = np.nansum(wv_cont_sum_daily_mm_SON,axis=0)
    pre_SON = np.nansum(pre_daily_SON, axis=0)
    wvcont_pct_SON = wvcont_mm_SON/np.nansum(pre_SON)*100
    
    wvcont_mm_seasonal_departure[np.where(Yearlist==year)[0][0]*4+0,:,:] = wvcont_mm_DJF - wvcont_mm_seasonal_climatology[0,:,:]
    wvcont_mm_seasonal_departure[np.where(Yearlist==year)[0][0]*4+1,:,:] = wvcont_mm_MAM - wvcont_mm_seasonal_climatology[1,:,:]
    wvcont_mm_seasonal_departure[np.where(Yearlist==year)[0][0]*4+2,:,:] = wvcont_mm_JJA - wvcont_mm_seasonal_climatology[2,:,:]
    wvcont_mm_seasonal_departure[np.where(Yearlist==year)[0][0]*4+3,:,:] = wvcont_mm_SON - wvcont_mm_seasonal_climatology[3,:,:]
    
    wvcont_pct_seasonal_departure[np.where(Yearlist==year)[0][0]*4+0,:,:] = wvcont_pct_DJF - wvcont_pct_seasonal_climatology[0,:,:]
    wvcont_pct_seasonal_departure[np.where(Yearlist==year)[0][0]*4+1,:,:] = wvcont_pct_MAM - wvcont_pct_seasonal_climatology[1,:,:]
    wvcont_pct_seasonal_departure[np.where(Yearlist==year)[0][0]*4+2,:,:] = wvcont_pct_JJA - wvcont_pct_seasonal_climatology[2,:,:]
    wvcont_pct_seasonal_departure[np.where(Yearlist==year)[0][0]*4+3,:,:] = wvcont_pct_SON - wvcont_pct_seasonal_climatology[3,:,:]

    pre_seasonal_departure[np.where(Yearlist==year)[0][0]*4+0,:,:] = pre_DJF - pre_seasonal_climatology[0,:,:]
    pre_seasonal_departure[np.where(Yearlist==year)[0][0]*4+1,:,:] = pre_MAM - pre_seasonal_climatology[1,:,:]
    pre_seasonal_departure[np.where(Yearlist==year)[0][0]*4+2,:,:] = pre_JJA - pre_seasonal_climatology[2,:,:]
    pre_seasonal_departure[np.where(Yearlist==year)[0][0]*4+3,:,:] = pre_SON - pre_seasonal_climatology[3,:,:]
    

ofile = dir_out+region+'_seasonal_anom_MilDrought.nc'    
with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
    of.createDimension('i_cross', n_i)
    of.createDimension('j_cross', n_j)    
    of.createDimension('time',len(seaslist))

    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
    of['latitcrs'].units = 'degrees'
    of['latitcrs'][:] = latitcrs
    
    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
    of['longicrs'].units = 'degrees'
    of['longicrs'][:] = longicrs
    
    of.createVariable('time', 'f4', ('time'))
    of['time'].long_name = 'Season of 2002-2009 of anom period: DJF,MAM,JJA,SON'
    of['time'].units = 'season number'
    of['time'][:] = seaslist
    
    of.createVariable('wvcont_mm_seasonal_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['wvcont_mm_seasonal_departure'].long_name = 'anomaly of seasonal water vapour contribution [mm], wrt climatology'
    of['wvcont_mm_seasonal_departure'].units = 'mm'
    of['wvcont_mm_seasonal_departure'][:] = wvcont_mm_seasonal_departure
    
    of.createVariable('wvcont_pct_seasonal_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['wvcont_pct_seasonal_departure'].long_name = 'anomaly of seasonal water vapour contribution [%], wrt climatology'
    of['wvcont_pct_seasonal_departure'].units = '%'
    of['wvcont_pct_seasonal_departure'][:] = wvcont_pct_seasonal_departure
    
    of.createVariable('pre_seasonal_departure', 'f4', ('time','i_cross', 'j_cross'))
    of['pre_seasonal_departure'].long_name = 'anomaly of seasonal rainfall [mm], wrt climatology'
    of['pre_seasonal_departure'].units = '%'
    of['pre_seasonal_departure'][:] = pre_seasonal_departure