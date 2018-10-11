# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:06:36 2018

This script creates yearly grids of water vapour contribution from individual cells to all rainfall
across the domain (Australia in this case) in a given year. Contribution is given in mm and 
percent.

The plot of these grids do not give an indication of the contribution of a grid cell's vapour to 
rainfall in its own area - that is calculated in a separate rainfall recycling script.

The plot of the grids gives the contribution of a cell's vapour to rainfall in the WHOLE domain.
So it can be expected that vapour in Tasmania won't contribute much for example, as the 
bulk of the annual rainfall total will be located in the tropical north of Australia, at least in the 
summer. 

Checked and OK 11/10/18.

@author: z3131380
"""

import numpy as np
from netCDF4 import Dataset 
import pandas
import os.path

dir_in = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'
dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/'
region = 'Australia'

#==============================================================================
# Create a date list   
#==============================================================================
Start_date = '19800101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '19810105' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pandas.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
Daylist = pandas.date_range(Start_date,End_date,freq='d')
dateidx = range(len(Daylist))

date_list=np.zeros([len(Daylist),7])*np.nan

for j in range(len(Daylist)):  
    yyyy = '%04u' % Daylist[j].year
    mm = '%02u' % Daylist[j].month
    dd = '%02u' % Daylist[j].day
    date_list[j,0] = yyyy+mm+dd
    date_list[j,1] = dateidx[j]         
    date_list[j,2] =yyyy
    date_list[j,3] = yyyy+mm
    date_list[j,4] = mm+dd
    date_list[j,5] = mm
    date_list[j,6] = dd

n_i,n_j = 134,205 # QIBT model dimensions
for year in Yearlist:
    indices = np.where(date_list[:,2]==year)[0]
    wv_cont_sum_daily_mm = np.zeros([len(indices),n_i,n_j])*np.nan
    wv_cont_apbl_sum_daily_mm = np.zeros([len(indices),n_i,n_j])*np.nan
    pre_grid_daily = np.zeros([len(indices),n_i,n_j])*np.nan
    pre_total_daily = np.zeros([len(indices)])
    wv_cont_sum_daily_pct = np.zeros([len(indices),n_i,n_j])*np.nan
    wv_cont_apbl_sum_daily_pct = np.zeros([len(indices),n_i,n_j])*np.nan
    for idx in indices:
        month = '%02u' % Daylist[idx].month
        day = Daylist[idx].day
        fname = 'bt.'+str(year)+month+'_'+str(day)+'.nc'
        file = dir_out+region+'/100parcels/TS10min/'+fname
        if os.path.isfile(file) == True:
            print fname
            fh = Dataset(file, mode='r', format='NETCDF3_CLASSIC') 
            if fh.dimensions['gridcell_wvc'].size > 0: 
                x_loc = fh.variables['x_loc'][:] #x index of precipitation in no. cells from [0,0]
                y_loc = fh.variables['y_loc'][:] #y index of precipitation in no. cells from [0,0]
                latitcrs = fh.variables['latitcrs'][:] # latitudes of curvilinear grid
                longicrs = fh.variables['longicrs'][:] # longitude of curvilinear grid
                wv_cont = fh.variables['wv_cont'][:] # water vapour contribution of domain cells to each cell where it rained
                wv_cont_apbl = fh.variables['wv_cont_apbl'][:] # as above, but above the PBL
                pre = fh.variables['pre'][:] # rainfall in each cell on that day
                #day = fh.variables['day'][:]    
            fh.close()    

            # Find contribution of each cell's vapour to each rain cell, and the total contribution 
            # of each cell to all rain in WHOLE domain each day
            wv_mm = np.zeros(np.shape(wv_cont))*np.nan
            wv_apbl_mm = np.zeros(np.shape(wv_cont))*np.nan            
            for k in range(np.shape(wv_cont)[0]):
                wv_mm[k,:,:] = wv_cont[k,:,:]*pre[k]
                wv_apbl_mm[k,:,:] = wv_cont_apbl[k,:,:]*pre[k]
            wv_cont_sum_daily_mm[idx,:,:] = np.nansum(wv_mm,axis=0)
            wv_cont_apbl_sum_daily_mm[idx,:,:] = np.nansum(wv_apbl_mm,axis=0)

            # Place rain cells on the x,y grid for plotting purposes/checking
            for k in range(np.shape(wv_cont)[0]): # or arange?
                rain_yloc,rain_xloc =  y_loc[k],x_loc[k]
                pre_grid_daily[idx,rain_yloc,rain_xloc] = pre[k]
                
            # Total rain in the whole domain each day
            pre_total_daily[idx] = pre.sum()
            
            # Find the proportion of vapour each cell contributed to total rainfall in the whole 
            # domain each day
            wv_cont_sum_daily_pct[idx,:,:] = wv_cont_sum_daily_mm[idx,:,:] / pre_total_daily[idx]
            wv_cont_apbl_sum_daily_pct[idx,:,:] = wv_cont_apbl_sum_daily_mm[idx,:,:] / pre_total_daily[idx]
        
        else: continue # if file does not exist, go to next day.
    
    pre_total_yearly = np.nansum(pre_total_daily)    
    wv_cont_sum_yearly_mm = np.nansum(wv_cont_sum_daily_mm,axis=0)
    wv_cont_sum_yearly_pct = wv_cont_sum_yearly_mm / pre_total_yearly
    wv_cont_apbl_sum_yearly_mm = np.nansum(wv_cont_apbl_sum_daily_mm,axis=0)
    wv_cont_apbl_sum_yearly_pct = wv_cont_apbl_sum_yearly_mm / pre_total_yearly
    
                    
    # Save to netcdf (save in same NETCDF3 as model output??)
#    ofile = dir_out+region+'/100parcels/TS10min/Processed/Yearly/'+region+'_'+str(year)+'_wvcont.nc'    
#    with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
#        of.createDimension('i_cross', n_i)
#        of.createDimension('j_cross', n_j)    
#        of.createDimension('time', len(indices))
#        
#        of.createVariable('day', 'f4', ('time'))
#        of['day'].long_name = 'days since '
#        of['day'].units = 'days'
#        of['day'][:] = idx
#
#        of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
#        of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
#        of['latitcrs'].units = 'degrees'
#        of['latitcrs'][:] = latitcrs
#        
#        of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
#        of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
#        of['longicrs'].units = 'degrees'
#        of['longicrs'][:] = longicrs
#        
#        of.createVariable('pre', 'f4', ('time','i_cross', 'j_cross'))
#        of['pre'].long_name = 'precipitation'
#        of['pre'].units = 'mm'
#        of['pre'][:] = pre_grid_daily
#        
#        of.createVariable('wv_cont_sum_yearly_mm', 'f4', ('i_cross', 'j_cross'))
#        of['wv_cont_sum_yearly_mm'].long_name = 'mm contribution of each cell to all rain in the domain'
#        of['wv_cont_sum_yearly_mm'].units = 'mm'
#        of['wv_cont_sum_yearly_mm'][:] = wv_cont_sum_yearly_mm
#        of['wv_cont_sum_yearly_mm'].fill_value = -9999.0
#        
#        of.createVariable('wv_cont_sum_yearly_pct', 'f4', ('i_cross', 'j_cross'))
#        of['wv_cont_sum_yearly_pct'].long_name = '% contribution of each cell to all rain in the domain'
#        of['wv_cont_sum_yearly_pct'].units = '%'
#        of['wv_cont_sum_yearly_pct'][:] = wv_cont_sum_yearly_pct
#        of['wv_cont_sum_yearly_pct'].fill_value = -9999.0
#        
#        of.createVariable('wv_cont_apbl_sum_yearly_mm', 'f4', ('i_cross', 'j_cross'))
#        of['wv_cont_apbl_sum_yearly_mm'].long_name = 'mm contribution above PBL of each cell to all rain in the domain'
#        of['wv_cont_apbl_sum_yearly_mm'].units = 'mm'
#        of['wv_cont_apbl_sum_yearly_mm'][:] = wv_cont_apbl_sum_yearly_mm
#        of['wv_cont_apbl_sum_yearly_mm'].fill_value = -9999.0
#        
#        of.createVariable('wv_cont_apbl_sum_yearly_pct', 'f4', ('i_cross', 'j_cross'))
#        of['wv_cont_apbl_sum_yearly_pct'].long_name = '% contribution above PBL of each cell to all rain in the domain'
#        of['wv_cont_apbl_sum_yearly_pct'].units = '%'
#        of['wv_cont_apbl_sum_yearly_pct'][:] = wv_cont_apbl_sum_yearly_pct
#        of['wv_cont_apbl_sum_yearly_pct'].fill_value = -9999.0