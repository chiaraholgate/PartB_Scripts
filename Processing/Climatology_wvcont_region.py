# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 09:16:39 2018

This script calculates the mean seasonal total for QIBT modelled water vapour contribution and precip.
It also calculates the April to November (southern wet season) climatology.

1979 is not used in the climatology for summer, as the model results begin at 31/1/79, so there's no data
for Dec 1978 or Jan 1979. Other seasons include data from 1979.

TO BE CHECKED


@author: z3131380
"""
import os
import numpy as np
import pandas
from netCDF4 import Dataset 

dir_data = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Yearly/'
dir_out = r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/'
region = 'Australia' 
n_i,n_j = 134,205 # QIBT model dimensions

# Date info
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pandas.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
full_daylist = pandas.date_range(Start_date,End_date,freq='d')

##==============================================================================
## Seasonal and annual climatology
## 
##==============================================================================
#
## Create blank arrays for output: wvcont (mm & %), precip
#wvcont_mm_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_DJF_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#wvcont_mm_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_MAM_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#wvcont_mm_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_JJA_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#wvcont_mm_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_SON_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#wvcont_mm_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan; pre_AprNov_totals = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#wvcont_mm_annual = np.zeros([len(Yearlist),n_i,n_j])*np.nan; wvcont_pct_annual = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#pre_annual = np.zeros([len(Yearlist),n_i,n_j])*np.nan
#    
#for year in Yearlist:
#    start_day = str(year)+'0101' ; end_day = str(year)+'1231'
#    daylist = pandas.date_range(start_day,end_day,freq='d')  
#    daylist_MAM = []; daylist_JJA = []; daylist_SON = []; daylist_DJF = []; daylist_AprNov = []
#    for d in daylist:
#        if d.month==3 or d.month==4 or d.month==5:
#            daylist_MAM.append(d)
#        if d.month==6 or d.month==7 or d.month==8:
#            daylist_JJA.append(d)
#        if d.month==9 or d.month==10 or d.month==11:
#            daylist_SON.append(d)
#        if np.logical_and(d.month >=4,d.month<=11):
#            daylist_AprNov.append(d)
#            
#        # DJF must include the Dec dates from the previous year        
#        if year > 1979: 
#            if d.month==1 or d.month==2:
#                daylist_DJF.append(d)
#            if d.month==12:
#                daylist_DJF.append(d-pandas.DateOffset(years=1))
#
#    daylist_DJF = pandas.to_datetime(daylist_DJF); daylist_DJF = daylist_DJF.sort_values()    
#    daylist_MAM = pandas.to_datetime(daylist_MAM) 
#    daylist_JJA = pandas.to_datetime(daylist_JJA)
#    daylist_SON = pandas.to_datetime(daylist_SON)
#    daylist_AprNov = pandas.to_datetime(daylist_AprNov)
#
#    wvcont_mm_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan; wvcont_pct_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan; pre_DJF = np.zeros([len(daylist_DJF),n_i,n_j])*np.nan
#    wvcont_mm_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan; wvcont_pct_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan; pre_MAM = np.zeros([len(daylist_MAM),n_i,n_j])*np.nan
#    wvcont_mm_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan; wvcont_pct_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan; pre_JJA = np.zeros([len(daylist_JJA),n_i,n_j])*np.nan
#    wvcont_mm_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan; wvcont_pct_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan; pre_SON = np.zeros([len(daylist_SON),n_i,n_j])*np.nan
#    wvcont_mm_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j])*np.nan; wvcont_pct_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j])*np.nan; pre_AprNov = np.zeros([len(daylist_AprNov),n_i,n_j])*np.nan
#    
#    file = dir_data+region+'_'+str(year)+'_wvcont.nc'
#    if os.path.isfile(file) == True:
#        fh = Dataset(file, mode='r') 
#        wvcont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]    
#        wvcont_sum_daily_pct = fh.variables['wv_cont_sum_daily_pct'][:] 
#        pre = fh.variables['pre'][:]  
#        latitcrs = fh.variables['latitcrs'][:]  
#        longicrs = fh.variables['longicrs'][:]  
#        fh.close() 
#    # DJF must include the Dec data from the previous year   
#    file = dir_data+region+'_'+str(year-1)+'_wvcont.nc'
#    if os.path.isfile(file) == True:
#        fh = Dataset(file, mode='r') 
#        wvcont_sum_daily_mm_prevyr = fh.variables['wv_cont_sum_daily_mm'][:]    
#        wvcont_sum_daily_pct_prevyr = fh.variables['wv_cont_sum_daily_pct'][:] 
#        pre_prevyr = fh.variables['pre'][:]  
#        fh.close()
#    if os.path.isfile(file) == False: # i.e. if year is 1979
#        wvcont_sum_daily_mm_prevyr = 0
#        wvcont_sum_daily_pct_prevyr = 0
#        pre_prevyr = 0
#        
#    daylist_prevyr = pandas.date_range(str(year-1)+'0101', str(year-1)+'1231',freq='d')    
#    for d in daylist_DJF:
#        if d.month < 3:
#            idx_daylist = np.where(daylist==d)[0][0]
#            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]      
#            wvcont_mm_DJF[idx_daylist_DJF,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
#            #wvcont_pct_DJF[idx_daylist_DJF,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
#            pre_DJF[idx_daylist_DJF,:,:] = pre[idx_daylist,:,:]
#        elif d.month==12:
#            idx_daylist = np.where(daylist_prevyr==d)[0][0]
#            idx_daylist_DJF = np.where(daylist_DJF==d)[0][0]   
#            wvcont_mm_DJF[idx_daylist_DJF,:,:] = wvcont_sum_daily_mm_prevyr[idx_daylist,:,:]
#            #wvcont_pct_DJF[idx_daylist_DJF,:,:] = wvcont_sum_daily_pct_prevyr[idx_daylist,:,:]
#            pre_DJF[idx_daylist_DJF,:,:] = pre_prevyr[idx_daylist,:,:]
#    wvcont_mm_DJF_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_DJF,axis=0)
#    #wvcont_pct_DJF_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_DJF,axis=0)
#    pre_DJF_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_DJF,axis=0)
#    wvcont_pct_DJF_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_DJF,axis=0)/np.nansum(pre_DJF)*100
#        
#    for d in daylist_MAM:
#        idx_daylist = np.where(daylist==d)[0][0]
#        idx_daylist_MAM = np.where(daylist_MAM==d)[0][0]
#        wvcont_mm_MAM[idx_daylist_MAM,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
#        #wvcont_pct_MAM[idx_daylist_MAM,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
#        pre_MAM[idx_daylist_MAM,:,:] = pre[idx_daylist,:,:]
#    wvcont_mm_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_MAM,axis=0)
#    #wvcont_pct_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_MAM,axis=0)
#    pre_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_MAM,axis=0)
#    wvcont_pct_MAM_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_MAM,axis=0)/np.nansum(pre_MAM)*100
#      
#    for d in daylist_JJA:
#        idx_daylist = np.where(daylist==d)[0][0]
#        idx_daylist_JJA = np.where(daylist_JJA==d)[0][0]
#        wvcont_mm_JJA[idx_daylist_JJA,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
#        #wvcont_pct_JJA[idx_daylist_JJA,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
#        pre_JJA[idx_daylist_JJA,:,:] = pre[idx_daylist,:,:]
#    wvcont_mm_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_JJA,axis=0)
#    #wvcont_pct_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_JJA,axis=0)
#    pre_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_JJA,axis=0)
#    wvcont_pct_JJA_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_JJA,axis=0)/np.nansum(pre_JJA)*100
#        
#    for d in daylist_SON:
#        idx_daylist = np.where(daylist==d)[0][0]
#        idx_daylist_SON = np.where(daylist_SON==d)[0][0]
#        wvcont_mm_SON[idx_daylist_SON,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
#        #wvcont_pct_SON[idx_daylist_SON,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
#        pre_SON[idx_daylist_SON,:,:] = pre[idx_daylist,:,:]
#    wvcont_mm_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_SON,axis=0)
#    #wvcont_pct_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_SON,axis=0)
#    pre_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_SON,axis=0)
#    wvcont_pct_SON_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_SON,axis=0)/np.nansum(pre_SON)*100
#    
#    for d in daylist_AprNov:
#        idx_daylist = np.where(daylist==d)[0][0]
#        idx_daylist_AprNov = np.where(daylist_AprNov==d)[0][0]
#        wvcont_mm_AprNov[idx_daylist_AprNov,:,:] = wvcont_sum_daily_mm[idx_daylist,:,:]
#        #wvcont_pct_AprNov[idx_daylist_AprNov,:,:] = wvcont_sum_daily_pct[idx_daylist,:,:]
#        pre_AprNov[idx_daylist_AprNov,:,:] = pre[idx_daylist,:,:]
#    wvcont_mm_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_AprNov,axis=0)
#    #wvcont_pct_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_pct_AprNov,axis=0)
#    pre_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(pre_AprNov,axis=0)
#    wvcont_pct_AprNov_totals[np.where(Yearlist==year)[0][0],:,:]=np.nansum(wvcont_mm_AprNov,axis=0)/np.nansum(pre_AprNov)*100
#    
#    wvcont_mm_annual[np.where(Yearlist==year)[0][0],:,:] = np.nansum(wvcont_sum_daily_mm,axis=0)
#    pre_annual[np.where(Yearlist==year)[0][0],:,:] = np.nansum(pre,axis=0)
#    wvcont_pct_annual[np.where(Yearlist==year)[0][0],:,:] = wvcont_mm_annual[np.where(Yearlist==year)[0][0],:,:] /np.nansum(pre)*100
#        
#wvcont_mm_seasonal_climatology = np.zeros([4,n_i,n_j])*np.nan        
#wvcont_mm_seasonal_climatology[0,:,:] = np.nanmean(wvcont_mm_DJF_totals, axis=0)
#wvcont_mm_seasonal_climatology[1,:,:] = np.nanmean(wvcont_mm_MAM_totals, axis=0)
#wvcont_mm_seasonal_climatology[2,:,:] = np.nanmean(wvcont_mm_JJA_totals, axis=0)
#wvcont_mm_seasonal_climatology[3,:,:] = np.nanmean(wvcont_mm_SON_totals, axis=0)
#
#wvcont_pct_seasonal_climatology = np.zeros([4,n_i,n_j])*np.nan        
#wvcont_pct_seasonal_climatology[0,:,:] = np.nanmean(wvcont_pct_DJF_totals, axis=0)
#wvcont_pct_seasonal_climatology[1,:,:] = np.nanmean(wvcont_pct_MAM_totals, axis=0)
#wvcont_pct_seasonal_climatology[2,:,:] = np.nanmean(wvcont_pct_JJA_totals, axis=0)
#wvcont_pct_seasonal_climatology[3,:,:] = np.nanmean(wvcont_pct_SON_totals, axis=0)
#
#pre_seasonal_climatology = np.zeros([4,n_i,n_j])*np.nan        
#pre_seasonal_climatology[0,:,:] = np.nanmean(pre_DJF_totals, axis=0)
#pre_seasonal_climatology[1,:,:] = np.nanmean(pre_MAM_totals, axis=0)
#pre_seasonal_climatology[2,:,:] = np.nanmean(pre_JJA_totals, axis=0)
#pre_seasonal_climatology[3,:,:] = np.nanmean(pre_SON_totals, axis=0)
#   
#wvcont_mm_AprNov_climatology = np.nanmean(wvcont_mm_AprNov_totals, axis=0)
#wvcont_pct_AprNov_climatology = np.nanmean(wvcont_pct_AprNov_totals, axis=0)
#pre_AprNov_climatology = np.nanmean(pre_AprNov_totals, axis=0)
#
#wvcont_mm_annual_climatology = np.nanmean(wvcont_mm_annual,axis=0)
#wvcont_pct_annual_climatology = np.nanmean(wvcont_pct_annual,axis=0)
#pre_annual_climatology = np.nanmean(pre_annual,axis=0)
#
## Save to netcdf 
#ofile = dir_out+region+'_1979-2013_seasonal_climatology.nc'    
#with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
#    of.createDimension('i_cross', n_i)
#    of.createDimension('j_cross', n_j)    
#    of.createDimension('season',4)
#
#    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
#    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
#    of['latitcrs'].units = 'degrees'
#    of['latitcrs'][:] = latitcrs
#    
#    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
#    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
#    of['longicrs'].units = 'degrees'
#    of['longicrs'][:] = longicrs
#    
#    of.createVariable('pre_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
#    of['pre_seasonal_climatology'].long_name = 'mean seasonal precipitation'
#    of['pre_seasonal_climatology'].units = 'mm'
#    of['pre_seasonal_climatology'][:] = pre_seasonal_climatology
#    
#    of.createVariable('wvcont_mm_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
#    of['wvcont_mm_seasonal_climatology'].long_name = 'mean seasonal water vapour contribution [mm] to all rainfall in domain'
#    of['wvcont_mm_seasonal_climatology'].units = 'mm'
#    of['wvcont_mm_seasonal_climatology'][:] = wvcont_mm_seasonal_climatology
#    
#    of.createVariable('wvcont_pct_seasonal_climatology', 'f4', ('season','i_cross', 'j_cross'))
#    of['wvcont_pct_seasonal_climatology'].long_name = 'mean seasonal water vapour contribution [%] to all rainfall in domain'
#    of['wvcont_pct_seasonal_climatology'].units = '%'
#    of['wvcont_pct_seasonal_climatology'][:] = wvcont_pct_seasonal_climatology
#    
#    of.createVariable('wvcont_mm_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['wvcont_mm_AprNov_climatology'].long_name = 'mean Apr to Nov water vapour contribution [mm] to all rainfall in domain'
#    of['wvcont_mm_AprNov_climatology'].units = 'mm'
#    of['wvcont_mm_AprNov_climatology'][:] = wvcont_mm_AprNov_climatology
#    
#    of.createVariable('wvcont_pct_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['wvcont_pct_AprNov_climatology'].long_name = 'mean Apr to Nov water vapour contribution [%] to all rainfall in domain'
#    of['wvcont_pct_AprNov_climatology'].units = '%'
#    of['wvcont_pct_AprNov_climatology'][:] = wvcont_pct_AprNov_climatology
#    
#    of.createVariable('pre_AprNov_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['pre_AprNov_climatology'].long_name = 'mean Apr to Nov precipitation'
#    of['pre_AprNov_climatology'].units = 'mm'
#    of['pre_AprNov_climatology'][:] = pre_AprNov_climatology
#    
#    of.createVariable('wvcont_mm_annual_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['wvcont_mm_annual_climatology'].long_name = 'mean annual water vapour contribution [mm] to all rainfall in domain'
#    of['wvcont_mm_annual_climatology'].units = 'mm'
#    of['wvcont_mm_annual_climatology'][:] = wvcont_mm_annual_climatology
#    
#    of.createVariable('wvcont_pct_annual_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['wvcont_pct_annual_climatology'].long_name = 'mean annual water vapour contribution [%] to all rainfall in domain'
#    of['wvcont_pct_annual_climatology'].units = 'mm'
#    of['wvcont_pct_annual_climatology'][:] = wvcont_pct_annual_climatology
#    
#    of.createVariable('pre_annual_climatology', 'f4', ('i_cross', 'j_cross'))
#    of['pre_annual_climatology'].long_name = 'mean annual precipitation'
#    of['pre_annual_climatology'].units = 'mm'
#    of['pre_annual_climatology'][:] = pre_annual_climatology
#        
    
        
#==============================================================================
# Monthly climatology
#==============================================================================
wvcont_mm_jan = wvcont_mm_feb =wvcont_mm_mar =wvcont_mm_apr =np.empty((0,n_i,n_j ))
wvcont_mm_may =wvcont_mm_jun = wvcont_mm_jul =wvcont_mm_aug =np.empty((0,n_i,n_j ))
wvcont_mm_sep =wvcont_mm_oct =wvcont_mm_nov =wvcont_mm_dec =np.empty((0,n_i,n_j ))

pre_jan = pre_feb =pre_mar = pre_apr =pre_may =pre_jun = pre_jul = pre_aug = np.empty((0,n_i,n_j ))
pre_sep = pre_oct = pre_nov = pre_dec = np.empty((0,n_i,n_j ))

for year in Yearlist:
    file = dir_data+region+'_'+str(year)+'_wvcont.nc'
    if os.path.isfile(file) == True:
        fh = Dataset(file, mode='r') 
        wvcont_sum_daily_mm = fh.variables['wv_cont_sum_daily_mm'][:]  
        pre = fh.variables['pre'][:]  
        latitcrs = fh.variables['latitcrs'][:]  
        longicrs = fh.variables['longicrs'][:]  
        fh.close() 
    
    daylist = pandas.date_range(str(year)+'0101',str(year)+'1231',freq='d')
    monthlist_startday = pandas.date_range(str(year)+'0101',str(year)+'1231',freq='MS')    
    monthlist_endday = pandas.date_range(str(year)+'0101',str(year)+'1231',freq='m') 

    wvcont_mm_jan = np.vstack((wvcont_mm_jan,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[0])[0][0]:\
        np.where(daylist==monthlist_endday[0])[0][0]+1,:,:]))
    wvcont_mm_feb = np.vstack((wvcont_mm_feb,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[1])[0][0]:\
        np.where(daylist==monthlist_endday[1])[0][0]+1,:,:]))
    wvcont_mm_mar = np.vstack((wvcont_mm_mar,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[2])[0][0]:\
        np.where(daylist==monthlist_endday[2])[0][0]+1,:,:]))
    wvcont_mm_apr = np.vstack((wvcont_mm_apr,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[3])[0][0]:\
        np.where(daylist==monthlist_endday[3])[0][0]+1,:,:]))
    wvcont_mm_may = np.vstack((wvcont_mm_may,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[4])[0][0]:\
        np.where(daylist==monthlist_endday[4])[0][0]+1,:,:]))
    wvcont_mm_jun = np.vstack((wvcont_mm_jun,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[5])[0][0]:\
        np.where(daylist==monthlist_endday[5])[0][0]+1,:,:]))
    wvcont_mm_jul = np.vstack((wvcont_mm_jul,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[6])[0][0]:\
        np.where(daylist==monthlist_endday[6])[0][0]+1,:,:]))
    wvcont_mm_aug = np.vstack((wvcont_mm_aug,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[7])[0][0]:\
        np.where(daylist==monthlist_endday[7])[0][0]+1,:,:]))
    wvcont_mm_sep = np.vstack((wvcont_mm_sep,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[8])[0][0]:\
        np.where(daylist==monthlist_endday[8])[0][0]+1,:,:]))
    wvcont_mm_oct = np.vstack((wvcont_mm_oct,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[9])[0][0]:\
        np.where(daylist==monthlist_endday[9])[0][0]+1,:,:]))
    wvcont_mm_nov = np.vstack((wvcont_mm_nov,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[10])[0][0]:\
        np.where(daylist==monthlist_endday[10])[0][0]+1,:,:]))
    wvcont_mm_dec = np.vstack((wvcont_mm_dec,wvcont_sum_daily_mm[np.where(daylist==monthlist_startday[11])[0][0]:\
        np.where(daylist==monthlist_endday[11])[0][0]+1,:,:]))
    
    
    pre_jan = np.vstack((pre_jan,pre[np.where(daylist==monthlist_startday[0])[0][0]:\
        np.where(daylist==monthlist_endday[0])[0][0]+1,:,:]))
    pre_feb = np.vstack((pre_feb,pre[np.where(daylist==monthlist_startday[1])[0][0]:\
        np.where(daylist==monthlist_endday[1])[0][0]+1,:,:]))
    pre_mar = np.vstack((pre_mar,pre[np.where(daylist==monthlist_startday[2])[0][0]:\
        np.where(daylist==monthlist_endday[2])[0][0]+1,:,:]))
    pre_apr = np.vstack((pre_apr,pre[np.where(daylist==monthlist_startday[3])[0][0]:\
        np.where(daylist==monthlist_endday[3])[0][0]+1,:,:]))
    pre_may = np.vstack((pre_may,pre[np.where(daylist==monthlist_startday[4])[0][0]:\
        np.where(daylist==monthlist_endday[4])[0][0]+1,:,:]))
    pre_jun = np.vstack((pre_jun,pre[np.where(daylist==monthlist_startday[5])[0][0]:\
        np.where(daylist==monthlist_endday[5])[0][0]+1,:,:]))
    pre_jul = np.vstack((pre_jul,pre[np.where(daylist==monthlist_startday[6])[0][0]:\
        np.where(daylist==monthlist_endday[6])[0][0]+1,:,:]))
    pre_aug = np.vstack((pre_aug,pre[np.where(daylist==monthlist_startday[7])[0][0]:\
        np.where(daylist==monthlist_endday[7])[0][0]+1,:,:]))
    pre_sep = np.vstack((pre_sep,pre[np.where(daylist==monthlist_startday[8])[0][0]:\
        np.where(daylist==monthlist_endday[8])[0][0]+1,:,:]))
    pre_oct = np.vstack((pre_oct,pre[np.where(daylist==monthlist_startday[9])[0][0]:\
        np.where(daylist==monthlist_endday[9])[0][0]+1,:,:]))
    pre_nov = np.vstack((pre_nov,pre[np.where(daylist==monthlist_startday[10])[0][0]:\
        np.where(daylist==monthlist_endday[10])[0][0]+1,:,:]))
    pre_dec = np.vstack((pre_dec,pre[np.where(daylist==monthlist_startday[11])[0][0]:\
        np.where(daylist==monthlist_endday[11])[0][0]+1,:,:]))
    
wvcont_mm_monthly_climatology = np.zeros([12,n_i,n_j])*np.nan          
wvcont_mm_monthly_climatology[0,:,:] = np.nanmean(wvcont_mm_jan,axis=0)
wvcont_mm_monthly_climatology[1,:,:] = np.nanmean(wvcont_mm_feb,axis=0)
wvcont_mm_monthly_climatology[2,:,:] = np.nanmean(wvcont_mm_mar,axis=0)
wvcont_mm_monthly_climatology[3,:,:] = np.nanmean(wvcont_mm_apr,axis=0)
wvcont_mm_monthly_climatology[4,:,:] = np.nanmean(wvcont_mm_may,axis=0)
wvcont_mm_monthly_climatology[5,:,:] = np.nanmean(wvcont_mm_jun,axis=0)
wvcont_mm_monthly_climatology[6,:,:] = np.nanmean(wvcont_mm_jul,axis=0)
wvcont_mm_monthly_climatology[7,:,:] = np.nanmean(wvcont_mm_aug,axis=0)
wvcont_mm_monthly_climatology[8,:,:] = np.nanmean(wvcont_mm_sep,axis=0)
wvcont_mm_monthly_climatology[9,:,:] = np.nanmean(wvcont_mm_oct,axis=0)
wvcont_mm_monthly_climatology[10,:,:] = np.nanmean(wvcont_mm_nov,axis=0)
wvcont_mm_monthly_climatology[11,:,:] = np.nanmean(wvcont_mm_dec,axis=0)


pre_monthly_climatology = np.zeros([12,n_i,n_j])*np.nan          
pre_monthly_climatology[0,:,:] = np.nanmean(pre_jan,axis=0)
pre_monthly_climatology[1,:,:] = np.nanmean(pre_feb,axis=0)
pre_monthly_climatology[2,:,:] = np.nanmean(pre_mar,axis=0)
pre_monthly_climatology[3,:,:] = np.nanmean(pre_apr,axis=0)
pre_monthly_climatology[4,:,:] = np.nanmean(pre_may,axis=0)
pre_monthly_climatology[5,:,:] = np.nanmean(pre_jun,axis=0)
pre_monthly_climatology[6,:,:] = np.nanmean(pre_jul,axis=0)
pre_monthly_climatology[7,:,:] = np.nanmean(pre_aug,axis=0)
pre_monthly_climatology[8,:,:] = np.nanmean(pre_sep,axis=0)
pre_monthly_climatology[9,:,:] = np.nanmean(pre_oct,axis=0)
pre_monthly_climatology[10,:,:] = np.nanmean(pre_nov,axis=0)
pre_monthly_climatology[11,:,:] = np.nanmean(pre_dec,axis=0)

wvcont_pct_monthly_climatology = np.zeros([12,n_i,n_j])*np.nan        
for i in range(12):
    wvcont_pct_monthly_climatology[i,:,:] = wvcont_mm_monthly_climatology[i,:,:]/np.nansum(pre_monthly_climatology[i,:,:])


# Save to netcdf 
ofile = dir_out+region+'_1979-2013_monthly_climatology.nc'    
with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
    of.createDimension('i_cross', n_i)
    of.createDimension('j_cross', n_j)    
    of.createDimension('month',12)

    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
    of['latitcrs'].units = 'degrees'
    of['latitcrs'][:] = latitcrs
    
    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
    of['longicrs'].units = 'degrees'
    of['longicrs'][:] = longicrs
    
    of.createVariable('wvcont_mm_monthly_climatology', 'f4', ('month','i_cross', 'j_cross'))
    of['wvcont_mm_monthly_climatology'].long_name = 'mean monthly water vapour contribution to region P'
    of['wvcont_mm_monthly_climatology'].units = 'mm'
    of['wvcont_mm_monthly_climatology'][:] = wvcont_mm_monthly_climatology
    
    of.createVariable('wvcont_pct_monthly_climatology', 'f4', ('month','i_cross', 'j_cross'))
    of['wvcont_pct_monthly_climatology'].long_name = 'mean monthly water vapour contribution to region P'
    of['wvcont_pct_monthly_climatology'].units = '%'
    of['wvcont_pct_monthly_climatology'][:] = wvcont_pct_monthly_climatology
    
    of.createVariable('pre_monthly_climatology', 'f4', ('month','i_cross', 'j_cross'))
    of['pre_monthly_climatology'].long_name = 'mean monthly precipitation'
    of['pre_monthly_climatology'].units = 'mm'
    of['pre_monthly_climatology'][:] = pre_monthly_climatology