#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 13:30:53 2017

This script:
    
(1) sums the total water vapour contribution to the Southern MDB
on a given day, over the whole model domain. The aim is to compare the 
contribution from the domain when using different numbers of parcels.

(2) calculates the spatial correlation on a given day between runs with 
different numbers of parcels, considering 1 month run (Dec 1981)

(3) as per (2), except a three-month run (JJA 1981)

(4) compares results of 100 parcels for whole of Australia when time step is
30mins and 15mins. >> r=0.96


TO DO:
    - CHECK spatial correlation calculations are CORRECT. 
    ....See See Jason's manuscript top of page 7.



@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import pandas as pd
from scipy import stats
import os
#from collections import OrderedDict
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, maskoceans
#from matplotlib import colors as c
#from matplotlib.colors import BoundaryNorm

dir_in = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/'
dir_results = '/srv/ccrc/data03/z3131380/PartB/Output/'
fig_dir = r'/srv/ccrc/data03/z3131380/PartB/Figures/Parcel_comparison/'


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


#==============================================================================
# No. cells
#==============================================================================

# Australia
file = '/srv/ccrc/data03/z3131380/PartB/Masks/NARCliM_AUS_land_sea_mask.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_Aus = np.count_nonzero(wsmask)
# Southern MDB:
file = '/srv/ccrc/data03/z3131380/PartB/Masks/n_s_basin/sbasin_rotpole.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_sthMDB = np.count_nonzero(wsmask)
# Whole MDB:
file = '/srv/ccrc/data03/z3131380/PartB/Masks/mdb_rotpole.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_MDB = np.count_nonzero(wsmask)
# Murrumbidgee (SW9) + Canberra (SW1):
file = '/srv/ccrc/data03/z3131380/PartB/Masks/n_s_basin/swwrpa_rotpole.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_SW9 = np.bitwise_or(wsmask==1,wsmask==9).sum()
# Gwydir (SW15):
file = '/srv/ccrc/data03/z3131380/PartB/Masks/n_s_basin/swwrpa_rotpole.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_SW15 = (wsmask==15).sum()

# Moonie (SW18):
file = '/srv/ccrc/data03/z3131380/PartB/Masks/n_s_basin/moonie_rotpole.nc'
fh = Dataset(file, mode='r')
wsmask = fh.variables['wsmask'][:]
fh.close()    
no_cells_Moonie = np.count_nonzero(wsmask)


#==============================================================================
# Part 1

# Select an arbitrary region in which to sum the source contribution,
# to compare between no. parcel runs.
#==============================================================================

# This region approximately corresponds to 30-40degS, 120-130degE.
source_test_region_LLCr,source_test_region_LLCc = 26,54
source_test_region_ULCr,source_test_region_ULCc = 48,54
source_test_region_LRCr,source_test_region_LRCc = 26,68
source_test_region_URCr,source_test_region_URCc = 48,68
# i.e. rows 26:48 and cols 54:68


# On some days, there will be no rain within the domain, but a file will be
# created anyway by the model. Deal with these by  checking .nc dimensions
# and only using those with a non-zero time dimension.
for day in np.arange(2,9):
    file = dir_results+'Southern_MDB/100parcels/bt.198112_'+str(day)+'.nc'
    fh = Dataset(file, mode='r') 
    if fh.dimensions['gridcell_wvc'].size > 0:
        x_loc = fh.variables['x_loc'][:]
        y_loc = fh.variables['y_loc'][:]
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()      
        print file
        for i in range(np.shape(wv_cont)[0]):
            vapour=wv_cont[i,:,:]*100 # percent    
            domain_contribution = vapour.sum()
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            vapour_nonzeros_test_region = vapour_nozeros[26:48,54:68]
            print 'day=',day,'domain_contribution=',domain_contribution
            print 'rain cell=',i,'test region sum=',vapour_nonzeros_test_region.sum()

            
for day in np.arange(2,9):
    file = dir_results+'Southern_MDB/50parcels/bt.198112_'+str(day)+'.nc'
    rain_day = '198112'+str(day)
    fh = Dataset(file, mode='r') 
    if fh.dimensions['gridcell_wvc'].size > 0:
        x_loc = fh.variables['x_loc'][:]
        y_loc = fh.variables['y_loc'][:]
        latitcrs = fh.variables['latitcrs'][:]
        longicrs = fh.variables['longicrs'][:]
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()      
        print file
        for i in range(np.shape(wv_cont)[0]):
            vapour=wv_cont[i,:,:]*100 # percent    
            domain_contribution = vapour.sum()
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            vapour_nonzeros_test_region = vapour_nozeros[26:48,54:68]
            print 'day=',day,'domain_contribution=',domain_contribution
            print 'rain cell=',i,'test region sum=',vapour_nonzeros_test_region.sum()
            
            

#==============================================================================
# Part 2

# Compare spatial correlation between daily output grids, for runs using different
# numbers of parcels. The reference to correlate against is the 100-parcel run.

# One month time frame (Dec 1981)

#==============================================================================

#==============================================================================
# Compare water vapour contributions between runs with different numbers of 
# parcels. Rain from all "rain-receiving cells" each day within the Southern MDB 
# is summed together to give a source contribution (%) to the Southern MDB on that day.

nparcels=[10,25,50,75,100]

# Southern MDB
ds_parcel_comp_SMDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series()}
df_parcel_comp_SMDB = pd.DataFrame(ds_parcel_comp_SMDB)         



for day in np.arange(1,31):
    for n in nparcels:
        file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(day)+'.nc'
        d = '%02d' % day
        rain_day = '198112'+str(d)
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            #vapour_nonzeros_test_region =vapour_nozeros[:,26:48,54:68].sum()
            df_parcel_comp_SMDB.loc[day,df_parcel_comp_SMDB.columns[0]] = rain_day
            #df_parcel_comp_SMDB.loc[day,df_parcel_comp_SMDB.columns[nparcels.index(n)+1]]=vapour_nonzeros_test_region
            df_parcel_comp_SMDB.loc[day,df_parcel_comp_SMDB.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
df_parcel_comp_SMDB=df_parcel_comp_SMDB[1:]


# Whole MDB
ds_parcel_comp_MDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series()}
df_parcel_comp_MDB = pd.DataFrame(ds_parcel_comp_MDB)         

for day in np.arange(1,31):
    for n in nparcels:
        file = dir_results+'MDB/'+str(n)+'parcels/bt.198112_'+str(day)+'.nc'
        d = '%02d' % day
        rain_day = '198112'+str(d)
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            #vapour_nonzeros_test_region =vapour_nozeros[:,26:48,54:68].sum()
            df_parcel_comp_MDB.loc[day,df_parcel_comp_MDB.columns[0]] = rain_day
            #df_parcel_comp_MDB.loc[day,df_parcel_comp_MDB.columns[nparcels.index(n)+1]]=vapour_nonzeros_test_region
            df_parcel_comp_MDB.loc[day,df_parcel_comp_MDB.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
df_parcel_comp_MDB=df_parcel_comp_MDB[1:]

# Australia
ds_parcel_comp_Aus = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series()}
df_parcel_comp_Aus = pd.DataFrame(ds_parcel_comp_Aus)         

for day in np.arange(1,31):
    for n in nparcels:
        file = dir_results+'Australia/'+str(n)+'parcels/bt.198112_'+str(day)+'.nc'
        d = '%02d' % day
        rain_day = '198112'+str(d)
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            #vapour_nonzeros_test_region =vapour_nozeros[:,26:48,54:68].sum()
            df_parcel_comp_Aus.loc[day,df_parcel_comp_Aus.columns[0]] = rain_day
            #df_parcel_comp_MDB.loc[day,df_parcel_comp_MDB.columns[nparcels.index(n)+1]]=vapour_nonzeros_test_region
            df_parcel_comp_Aus.loc[day,df_parcel_comp_Aus.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
df_parcel_comp_Aus=df_parcel_comp_Aus[1:]



#==============================================================================
# Spatially correlate - here ONLY CONSIDERING ONE RAIN CELL, i.e. wv_cont[0,:,:]

#==============================================================================

#for i in range(len(df_parcel_comp.ix[:,0])):
    
#for d in range(len(idxs)):    
#    rain_day = df_parcel_comp.ix[idxs[d],0]
#    file = dir_results+'Southern_MDB/100parcels/bt.198112_'+str(idxs[d])+'.nc'
#    fh = Dataset(file, mode='r')
#    wv_cont_100 = fh.variables['wv_cont'][0,:,:]
#    fh.close()    
#    # Stack each array
#    wv_cont_100_1d = np.hstack(wv_cont_100)
#    for n in nparcels[:-1]:        
#        file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
#        fh = Dataset(file, mode='r')
#        wv_cont_nparcels = fh.variables['wv_cont'][0,:,:]
#        fh.close()        
#        # Stack each array
#        wv_cont_nparcels_1d = np.hstack(wv_cont_nparcels)
#        
#        # Spatially correlate the 100-parcel source area map with the x-parcel
#        # source area map    
#        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])  
#        df_parcel_corr.loc[d,df_parcel_corr.columns[0]] = rain_day
#        df_parcel_corr.loc[d,df_parcel_corr.columns[nparcels.index(n)+1]] = stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2]

# >> How does this relate to cdo fldcor????


##==============================================================================
## Spatially correlate source area of (1) one rain cell at a time, one day, 100 parcels
## with (2) one rain cell at a time, one day, nparcels
##==============================================================================
#
## Southern MDB
#ds_parcel_corr_SMDB = {'1_Rain day':pd.Series(np.nan), 
#                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
#                  '5_100_75_corr':pd.Series()}
#df_parcel_corr_SMDB = pd.DataFrame(ds_parcel_corr_SMDB) 
#
#idxs = df_parcel_comp_SMDB.index.tolist()
#for d in range(len(idxs)):    
#    rain_day = df_parcel_comp_SMDB.ix[idxs[d],0]
#    df_parcel_corr_SMDB.loc[d,df_parcel_corr_SMDB.columns[0]] = rain_day
#    for n in nparcels[:-1]:
#        file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
#        fh = Dataset(file, mode='r')        
#        no_rain_cells = fh.dimensions['gridcell_wvc'].size 
#        fh.close()       
#        temp = []                              
#        for i in range(no_rain_cells):  
#            file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
#            fh = Dataset(file, mode='r') 
#            wv_cont_nparcels = fh.variables['wv_cont'][i,:,:]
#            fh.close()  
#            file = dir_results+'Southern_MDB/100parcels/bt.198112_'+str(idxs[d])+'.nc'
#            fh = Dataset(file, mode='r')
#            wv_cont_100 = fh.variables['wv_cont'][i,:,:]
#            fh.close()
#            wv_cont_100_1d = np.hstack(wv_cont_100)
#            wv_cont_nparcels_1d = np.hstack(wv_cont_nparcels)
#            temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
#        mean_day_nparcels_rain_cells = temp.mean()        
#        df_parcel_corr_SMDB.loc[d,df_parcel_corr_SMDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells
#
#
## Whole MDB
#ds_parcel_corr_MDB = {'1_Rain day':pd.Series(np.nan), 
#                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
#                  '5_100_75_corr':pd.Series()}
#df_parcel_corr_MDB = pd.DataFrame(ds_parcel_corr_MDB) 
#
#idxs = df_parcel_comp_MDB.index.tolist()
#for d in range(len(idxs)):    
#    rain_day = df_parcel_comp_MDB.ix[idxs[d],0]
#    df_parcel_corr_MDB.loc[d,df_parcel_corr_MDB.columns[0]] = rain_day
#    for n in nparcels[:-1]:
#        file = dir_results+'MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
#        fh = Dataset(file, mode='r')        
#        no_rain_cells = fh.dimensions['gridcell_wvc'].size 
#        fh.close()       
#        temp = []                              
#        for i in range(no_rain_cells):  
#            file = dir_results+'MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
#            fh = Dataset(file, mode='r') 
#            wv_cont_nparcels = fh.variables['wv_cont'][i,:,:]
#            fh.close()  
#            file = dir_results+'MDB/100parcels/bt.198112_'+str(idxs[d])+'.nc'
#            fh = Dataset(file, mode='r')
#            wv_cont_100 = fh.variables['wv_cont'][i,:,:]
#            fh.close()
#            wv_cont_100_1d = np.hstack(wv_cont_100)
#            wv_cont_nparcels_1d = np.hstack(wv_cont_nparcels)
#            temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
#        mean_day_nparcels_rain_cells = temp.mean()        
#        df_parcel_corr_MDB.loc[d,df_parcel_corr_MDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells
#
#
#
#### Create scatter plot of (x) no_cells vs (y) spatial corr
#fig, ax = plt.subplots()
#plt.title('1-31 December 1981',y=1.1)
## Southern MDB
#x1 = no_cells_sthMDB
#y1 = df_parcel_corr_SMDB['2_100_10_corr'].mean()
#y2 = df_parcel_corr_SMDB['3_100_25_corr'].mean()
#y3 = df_parcel_corr_SMDB['4_100_50_corr'].mean()
#y4 = df_parcel_corr_SMDB['5_100_75_corr'].mean()
#p10 = ax.scatter(x1,y1,marker='x',color='k')
#p25 = ax.scatter(x1,y2,marker='o',color='k')
#p50 = ax.scatter(x1,y3,marker='s',color='k')
#p75 = ax.scatter(x1,y4,marker='d',color='k')
## Whole MDB
#x1 = no_cells_MDB
#y1 = df_parcel_corr_MDB['2_100_10_corr'].mean()
#y2 = df_parcel_corr_MDB['3_100_25_corr'].mean()
#y3 = df_parcel_corr_MDB['4_100_50_corr'].mean()
#y4 = df_parcel_corr_MDB['5_100_75_corr'].mean()
#p10 = ax.scatter(x1,y1,marker='x',color='k')
#p25 = ax.scatter(x1,y2,marker='o',color='k')
#p50 = ax.scatter(x1,y3,marker='s',color='k')
#p75 = ax.scatter(x1,y4,marker='d',color='k')
#plt.xlim(0,500)
##plt.ylim(-1,1)
##plt.xscale('log')
#ax.set_yticklabels(['{:,}'.format(round(x,2)) for x in ax.get_yticks().tolist()])
#plt.xlabel('Number of cells')
#plt.ylabel('Spatial correlation')
#ax.legend((p10,p25,p50,p75),
#          ('10 parcels','25 parcels','50 parcels','75 parcels'),
#          bbox_to_anchor=(0.5, 1.13), loc='upper center', ncol=4,frameon=False)
#fig.savefig(fig_dir+'Spatial_corr_nparcels.png',bbox_inches='tight')




# OR
# Create table of mean spatial correlations for each no.cells/region
#ds_parcel_corr_mean = {'1_SMDB':pd.Series(np.nan)}
#df_parcel_corr_mean = pd.DataFrame(ds_parcel_corr_mean) 
#df_parcel_corr_mean.loc[0,df_parcel_corr_mean.columns[0]] = df_parcel_corr_SMDB['2_100_10_corr'].mean()
#df_parcel_corr_mean.loc[1,df_parcel_corr_mean.columns[0]] = df_parcel_corr_SMDB['3_100_25_corr'].mean()
#df_parcel_corr_mean.loc[2,df_parcel_corr_mean.columns[0]] = df_parcel_corr_SMDB['4_100_50_corr'].mean()
#df_parcel_corr_mean.loc[3,df_parcel_corr_mean.columns[0]] = df_parcel_corr_SMDB['5_100_75_corr'].mean()
#
#fig, ax = plt.subplots()
#x1 = [no_cells_sthMDB,no_cells_sthMDB,no_cells_sthMDB,no_cells_sthMDB]
#y1 = df_parcel_corr_mean['1_SMDB']
#ax.scatter(x1,y1)

#==============================================================================
# Spatially correlate source area of (1) ALL rain cells, one day, 100 parcels
# with (2) ALL rain cells, one day, nparcels

# Southern MDB
ds_parcel_corr_SMDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_SMDB = pd.DataFrame(ds_parcel_corr_SMDB) 

idxs = df_parcel_comp_SMDB.index.tolist()
for d in range(len(idxs)):    
    rain_day = df_parcel_comp_SMDB.ix[idxs[d],0]
    df_parcel_corr_SMDB.loc[d,df_parcel_corr_SMDB.columns[0]] = rain_day
    for n in nparcels[:-1]:
        file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size 
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size                             
        fh.close()       
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                             
        for i in range(no_rain_cells):  
            file = dir_results+'Southern_MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r') 
            wv_cont_nparcels = fh.variables['wv_cont'][i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            fh.close()  
            file = dir_results+'Southern_MDB/100parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r')
            wv_cont_100 = fh.variables['wv_cont'][i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100)),2)
            fh.close()
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_SMDB.loc[d,df_parcel_corr_SMDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells


# Whole MDB
ds_parcel_corr_MDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_MDB = pd.DataFrame(ds_parcel_corr_MDB) 

idxs = df_parcel_comp_MDB.index.tolist()
for d in range(len(idxs)):    
    rain_day = df_parcel_comp_MDB.ix[idxs[d],0]
    df_parcel_corr_MDB.loc[d,df_parcel_corr_MDB.columns[0]] = rain_day
    for n in nparcels[:-1]:
        file = dir_results+'MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size 
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size                             
        fh.close()       
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                             
        for i in range(no_rain_cells):  
            file = dir_results+'MDB/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r') 
            wv_cont_nparcels = fh.variables['wv_cont'][i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            fh.close()  
            file = dir_results+'MDB/100parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r')
            wv_cont_100 = fh.variables['wv_cont'][i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100)),2)
            fh.close()
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        #print d,rain_day,n,'%.2f'%stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2]
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_MDB.loc[d,df_parcel_corr_MDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells

# Australia
ds_parcel_corr_Aus = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_Aus = pd.DataFrame(ds_parcel_corr_Aus) 

idxs = df_parcel_comp_Aus.index.tolist()
for d in range(len(idxs)):    
    rain_day = df_parcel_comp_Aus.ix[idxs[d],0]
    df_parcel_corr_Aus.loc[d,df_parcel_corr_Aus.columns[0]] = rain_day
    for n in nparcels[:-1]:
        file = dir_results+'Australia/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size 
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size                             
        fh.close()       
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                             
        for i in range(no_rain_cells):  
            file = dir_results+'Australia/'+str(n)+'parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r') 
            wv_cont_nparcels = fh.variables['wv_cont'][i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            fh.close()  
            file = dir_results+'Australia/100parcels/bt.198112_'+str(idxs[d])+'.nc'
            fh = Dataset(file, mode='r')
            wv_cont_100 = fh.variables['wv_cont'][i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100)),2)
            fh.close()
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        #print d,rain_day,n,'%.2f'%stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2]
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_Aus.loc[d,df_parcel_corr_Aus.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells


### Create scatter plot of (x) no_cells vs (y) spatial corr
fig, ax = plt.subplots()
plt.title('1-31 December 1981',y=1.1)
# Southern MDB
x1 = no_cells_sthMDB
y1 = df_parcel_corr_SMDB['2_100_10_corr'].mean()
y2 = df_parcel_corr_SMDB['3_100_25_corr'].mean()
y3 = df_parcel_corr_SMDB['4_100_50_corr'].mean()
y4 = df_parcel_corr_SMDB['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.95,'Sth MDB',rotation=90)
# Whole MDB
x1 = no_cells_MDB
y1 = df_parcel_corr_MDB['2_100_10_corr'].mean()
y2 = df_parcel_corr_MDB['3_100_25_corr'].mean()
y3 = df_parcel_corr_MDB['4_100_50_corr'].mean()
y4 = df_parcel_corr_MDB['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.95,'MDB',rotation=90)
# Australia
x1 = no_cells_Aus
y1 = df_parcel_corr_Aus['2_100_10_corr'].mean()
y2 = df_parcel_corr_Aus['3_100_25_corr'].mean()
y3 = df_parcel_corr_Aus['4_100_50_corr'].mean()
y4 = df_parcel_corr_Aus['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.95,'Australia',rotation=90)
# Axes
plt.xlim(0,3500)
#plt.ylim(-1,1)
#plt.xscale('log')
ax.set_yticklabels(['{:,}'.format(round(x,2)) for x in ax.get_yticks().tolist()])
plt.xlabel('Number of cells')
plt.ylabel('Spatial correlation')
ax.legend((p10,p25,p50,p75),
          ('10 parcels','25 parcels','50 parcels','75 parcels'),
          bbox_to_anchor=(0.5, 1.13), loc='upper center', ncol=4,frameon=False)
fig.savefig(fig_dir+'Spatial_corr_nparcels.png',bbox_inches='tight')



### Same as above figure, but plotting all cells rather than the mean
fig, ax = plt.subplots()
plt.title('1-31 December 1981',y=1.1)
# Southern MDB
x1 = [no_cells_sthMDB]*len(df_parcel_corr_SMDB)
y1 = df_parcel_corr_SMDB['2_100_10_corr']
y2 = df_parcel_corr_SMDB['3_100_25_corr']
y3 = df_parcel_corr_SMDB['4_100_50_corr']
y4 = df_parcel_corr_SMDB['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='x',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'Sth MDB',rotation=90)
# Whole MDB
x1 = [no_cells_MDB]*len(df_parcel_corr_MDB)
y1 = df_parcel_corr_MDB['2_100_10_corr']
y2 = df_parcel_corr_MDB['3_100_25_corr']
y3 = df_parcel_corr_MDB['4_100_50_corr']
y4 = df_parcel_corr_MDB['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='x',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'MDB',rotation=90)
# Australia
x1 = [no_cells_Aus]*len(df_parcel_corr_Aus)
y1 = df_parcel_corr_Aus['2_100_10_corr']
y2 = df_parcel_corr_Aus['3_100_25_corr']
y3 = df_parcel_corr_Aus['4_100_50_corr']
y4 = df_parcel_corr_Aus['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='x',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'Australia',rotation=90)
# Axes
plt.xlim(0,3500)
#plt.ylim(-1,1)
#plt.xscale('log')
ax.set_yticklabels(['{:,}'.format(round(x,2)) for x in ax.get_yticks().tolist()])
plt.xlabel('Number of cells')
plt.ylabel('Spatial correlation')
ax.legend((p10,p25,p50,p75),
          ('10 parcels','25 parcels','50 parcels','75 parcels'),
          bbox_to_anchor=(0.5, 1.13), loc='upper center', ncol=4,frameon=False)
fig.savefig(fig_dir+'Spatial_corr_nparcels_allcells.png',bbox_inches='tight')


#==============================================================================
# Part 3

# As per Part 2, except

# Three month time frame (Jun-Aug 1981)

#==============================================================================

#==============================================================================
# Compare water vapour contributions between runs with different numbers of 
# parcels. Rain from all "rain-receiving cells" each day within the Southern MDB 
# is summed together to give a source contribution (%) to the Southern MDB on that day.

nparcels=[10,25,50,75,100]
    
# Moonie
ds_parcel_comp_Moonie = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series(),
                  '7_file':pd.Series()}
df_parcel_comp_Moonie = pd.DataFrame(ds_parcel_comp_Moonie)         

for n in nparcels:
    if n==25 or n==75:
        continue
    file_list = os.listdir(dir_results+'Moonie/'+str(n)+'parcels/')
    file_paths = listdir_fullpath(dir_results+'Moonie/'+str(n)+'parcels/')
    for f in range(len(file_list)):
        file = file_paths[f]
        if len(file_list[f])==15:
            rain_day = file_list[f][3:9]+file_list[f][10:12]
        else:
            rain_day = file_list[f][3:9]+file_list[f][10:11] # can add a leading zero to day if you want
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            df_parcel_comp_Moonie.loc[f,df_parcel_comp_Moonie.columns[0]] = rain_day
            df_parcel_comp_Moonie.loc[f,df_parcel_comp_Moonie.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
            df_parcel_comp_Moonie.loc[f,df_parcel_comp_Moonie.columns[-1]]=file_list[f]
    

# Southern MDB
ds_parcel_comp_SMDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series(),
                  '7_file':pd.Series()}
df_parcel_comp_SMDB = pd.DataFrame(ds_parcel_comp_SMDB)         

for n in nparcels:
    file_list = os.listdir(dir_results+'Southern_MDB/'+str(n)+'parcels/')
    file_list = file_list[30:] # ignore Dec 1981 files
    file_paths = listdir_fullpath(dir_results+'Southern_MDB/'+str(n)+'parcels/')
    file_paths = file_paths[30:]
    for f in range(len(file_list)):
        file = file_paths[f]
        if len(file_list[f])==15:
            rain_day = file_list[f][3:9]+file_list[f][10:12]#file[65:71]+file[72:74]
        else:
            rain_day = file_list[f][3:9]+file_list[f][10:11]#file[65:71]+file[72:73] # can add a leading zero to day if you want
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            df_parcel_comp_SMDB.loc[f,df_parcel_comp_SMDB.columns[0]] = rain_day
            df_parcel_comp_SMDB.loc[f,df_parcel_comp_SMDB.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
            df_parcel_comp_SMDB.loc[f,df_parcel_comp_SMDB.columns[-1]]=file_list[f]
     
# Whole MDB
ds_parcel_comp_MDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series(),
                  '7_file':pd.Series()}
df_parcel_comp_MDB = pd.DataFrame(ds_parcel_comp_MDB)         

for n in nparcels:
    file_list = os.listdir(dir_results+'MDB/'+str(n)+'parcels/')
    file_list = file_list[30:] # ignore Dec 1981 files
    file_paths = listdir_fullpath(dir_results+'MDB/'+str(n)+'parcels/')
    file_paths = file_paths[30:]
    for f in range(len(file_list)):
        file = file_paths[f]
        if len(file_list[f])==15:
            rain_day = file_list[f][3:9]+file_list[f][10:12]
        else:
            rain_day = file_list[f][3:9]+file_list[f][10:11]# can add a leading zero to day if you want
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            df_parcel_comp_MDB.loc[f,df_parcel_comp_MDB.columns[0]] = rain_day
            df_parcel_comp_MDB.loc[f,df_parcel_comp_MDB.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
            df_parcel_comp_MDB.loc[f,df_parcel_comp_MDB.columns[-1]]=file_list[f]


# Australia
ds_parcel_comp_Aus = {'1_Rain day':pd.Series(np.nan), 
                  '2_cont_p10':pd.Series(),'3_cont_p25':pd.Series(),'4_cont_p50':pd.Series(),
                  '5_cont_p75':pd.Series(),'6_cont_p100':pd.Series(),
                  '7_file':pd.Series()}
df_parcel_comp_Aus = pd.DataFrame(ds_parcel_comp_Aus)         

for n in nparcels:
    if n==25 or n==75:
        continue
    file_list = os.listdir(dir_results+'Australia/'+str(n)+'parcels/')
    file_list = file_list[30:] # ignore Dec 1981 files
    file_paths = listdir_fullpath(dir_results+'Australia/'+str(n)+'parcels/')
    file_paths = file_paths[30:]
    for f in range(len(file_list)):
        file = file_paths[f]
        if len(file_list[f])==15:
            rain_day = file_list[f][3:9]+file_list[f][10:12]
        else:
            rain_day = file_list[f][3:9]+file_list[f][10:11]# can add a leading zero to day if you want
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[0]] = rain_day
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[-1]]=file_list[f]



#==============================================================================
# Spatially correlate source area of (1) ALL rain cells, one day, 100 parcels
# with (2) ALL rain cells, one day, nparcels

# Moonie
ds_parcel_corr_Moonie = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_Moonie = pd.DataFrame(ds_parcel_corr_Moonie) 

idxs = df_parcel_comp_Moonie.index.tolist()
for f in idxs:
    rain_day = df_parcel_comp_Moonie.ix[f][0]
    df_parcel_corr_Moonie.loc[f,df_parcel_corr_Moonie.columns[0]] = rain_day
    for n in nparcels[:-1]:
        if n==25 or n==75:
            continue
        file = dir_results+'Moonie/'+str(n)+'parcels/'+df_parcel_comp_Moonie['7_file'][f]
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size                                  
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size        
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()    
        file = dir_results+'Moonie/100parcels/'+df_parcel_comp_Moonie['7_file'][f]
        fh = Dataset(file, mode='r')
        wv_cont_100 = fh.variables['wv_cont'][:]
        fh.close()
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                                     
        for i in range(no_rain_cells):  
            wv_cont_nparcels = wv_cont[i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            wv_cont_100parcels = wv_cont_100[i,:,:]#fh.variables['wv_cont'][i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100parcels)),2)
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_Moonie.loc[f,df_parcel_corr_Moonie.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells


# Southern MDB
ds_parcel_corr_SMDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_SMDB = pd.DataFrame(ds_parcel_corr_SMDB) 

idxs = df_parcel_comp_SMDB.index.tolist()
for f in idxs:
    rain_day = df_parcel_comp_SMDB.ix[f][0]
    df_parcel_corr_SMDB.loc[f,df_parcel_corr_SMDB.columns[0]] = rain_day
    for n in nparcels[:-1]:
        file = dir_results+'Southern_MDB/'+str(n)+'parcels/'+df_parcel_comp_SMDB['7_file'][f]
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size                                  
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size        
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()    
        file = dir_results+'Southern_MDB/100parcels/'+df_parcel_comp_SMDB['7_file'][f]
        fh = Dataset(file, mode='r')
        wv_cont_100 = fh.variables['wv_cont'][:]
        fh.close()
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                                     
        for i in range(no_rain_cells):  
            #file = dir_results+'Southern_MDB/'+str(n)+'parcels/'+df_parcel_comp_SMDB['7_file'][f]
            #fh = Dataset(file, mode='r') 
            wv_cont_nparcels = wv_cont[i,:,:]#fh.variables['wv_cont'][i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            #fh.close()  
            #file = dir_results+'Southern_MDB/100parcels/'+df_parcel_comp_SMDB['7_file'][i]
            #fh = Dataset(file, mode='r')
            wv_cont_100parcels = wv_cont_100[i,:,:]#fh.variables['wv_cont'][i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100parcels)),2)
            #fh.close()
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_SMDB.loc[f,df_parcel_corr_SMDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells
        

# MDB
ds_parcel_corr_MDB = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_MDB = pd.DataFrame(ds_parcel_corr_MDB) 

idxs = df_parcel_comp_MDB.index.tolist()
for f in idxs:
    rain_day = df_parcel_comp_MDB.ix[f][0]
    df_parcel_corr_MDB.loc[f,df_parcel_corr_MDB.columns[0]] = rain_day
    for n in nparcels[:-1]:
        file = dir_results+'MDB/'+str(n)+'parcels/'+df_parcel_comp_MDB['7_file'][f]
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size                                  
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size        
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()    
        file = dir_results+'MDB/100parcels/'+df_parcel_comp_MDB['7_file'][f]
        fh = Dataset(file, mode='r')
        wv_cont_100 = fh.variables['wv_cont'][:]
        fh.close()
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                                     
        for i in range(no_rain_cells):  
            wv_cont_nparcels = wv_cont[i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            wv_cont_100parcels = wv_cont_100[i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100parcels)),2)
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_MDB.loc[f,df_parcel_corr_MDB.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells


# Australia
ds_parcel_corr_Aus = {'1_Rain day':pd.Series(np.nan), 
                  '2_100_10_corr':pd.Series(),'3_100_25_corr':pd.Series(),'4_100_50_corr':pd.Series(),
                  '5_100_75_corr':pd.Series()}
df_parcel_corr_Aus = pd.DataFrame(ds_parcel_corr_Aus) 

idxs = df_parcel_comp_Aus.index.tolist()
for f in idxs:
    rain_day = df_parcel_comp_Aus.ix[f][0]
    df_parcel_corr_Aus.loc[f,df_parcel_corr_Aus.columns[0]] = rain_day
    for n in nparcels[:-1]:
        if n==25 or n==75:
            continue
        file = dir_results+'Australia/'+str(n)+'parcels/'+df_parcel_comp_Aus['7_file'][f]
        fh = Dataset(file, mode='r')        
        no_rain_cells = fh.dimensions['gridcell_wvc'].size                                  
        nrows = fh.dimensions['i_cross'].size
        ncols = fh.dimensions['j_cross'].size        
        wv_cont = fh.variables['wv_cont'][:]
        fh.close()    
        file = dir_results+'Australia/100parcels/'+df_parcel_comp_Aus['7_file'][f]
        fh = Dataset(file, mode='r')
        wv_cont_100 = fh.variables['wv_cont'][:]
        fh.close()
        temp = []  
        wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
        wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                                     
        for i in range(no_rain_cells):  
            wv_cont_nparcels = wv_cont[i,:,:]
            wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
            wv_cont_100parcels = wv_cont_100[i,:,:]
            wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100parcels)),2)
        wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
        wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
        temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
        mean_day_nparcels_rain_cells = temp.mean()        
        df_parcel_corr_Aus.loc[f,df_parcel_corr_Aus.columns[nparcels.index(n)+1]] = mean_day_nparcels_rain_cells
         


### Create scatter plot of (x) no_cells vs (y) spatial corr
fig, ax = plt.subplots()
plt.title('1 June - 31 August 1981',y=1.1)
# Moonie
x1 = no_cells_Moonie
y1 = df_parcel_corr_Moonie['2_100_10_corr'].mean()
y2 = df_parcel_corr_Moonie['3_100_25_corr'].mean()
y3 = df_parcel_corr_Moonie['4_100_50_corr'].mean()
y4 = df_parcel_corr_Moonie['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.85,'Moonie (11)',rotation=90)
# Southern MDB
x1 = no_cells_sthMDB
y1 = df_parcel_corr_SMDB['2_100_10_corr'].mean()
y2 = df_parcel_corr_SMDB['3_100_25_corr'].mean()
y3 = df_parcel_corr_SMDB['4_100_50_corr'].mean()
y4 = df_parcel_corr_SMDB['5_100_75_corr'].mean()
p10s = ax.scatter(x1,y1,marker='x',color='k')
p25s = ax.scatter(x1,y2,marker='o',color='k')
p50s = ax.scatter(x1,y3,marker='s',color='k')
p75s = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.8,'Sth MDB (66)',rotation=90)
# Whole MDB
x1 = no_cells_MDB
y1 = df_parcel_corr_MDB['2_100_10_corr'].mean()
y2 = df_parcel_corr_MDB['3_100_25_corr'].mean()
y3 = df_parcel_corr_MDB['4_100_50_corr'].mean()
y4 = df_parcel_corr_MDB['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.8,'MDB (73)',rotation=90)
# Australia
x1 = no_cells_Aus
y1 = df_parcel_corr_Aus['2_100_10_corr'].mean()
y2 = df_parcel_corr_Aus['3_100_25_corr'].mean()
y3 = df_parcel_corr_Aus['4_100_50_corr'].mean()
y4 = df_parcel_corr_Aus['5_100_75_corr'].mean()
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='o',color='k')
p50 = ax.scatter(x1,y3,marker='s',color='k')
p75 = ax.scatter(x1,y4,marker='d',color='k')
plt.text(x1,0.8,'Australia (89)',rotation=90)
# Axes
plt.xlim(0,3500)
plt.ylim(0.55,1)
#plt.xscale('log')
ax.set_yticklabels(['{:,}'.format(round(x,2)) for x in ax.get_yticks().tolist()])
plt.xlabel('Number of cells')
plt.ylabel('Spatial correlation')
ax.legend((p10s,p25s,p50s,p75s),
          ('10 parcels','25 parcels','50 parcels','75 parcels'),
          bbox_to_anchor=(0.5, 1.13), loc='upper center', ncol=4,frameon=False)
fig.savefig(fig_dir+'Spatial_corr_nparcels_JJA1981.png',bbox_inches='tight')



### Same as above figure, but plotting all cells rather than the mean
fig, ax = plt.subplots()
plt.title('1 June - 31 August 1981',y=1.1)
# Moonie
x1 = [no_cells_Moonie]*len(df_parcel_corr_Moonie)
y1 = df_parcel_corr_Moonie['2_100_10_corr']
y2 = df_parcel_corr_Moonie['3_100_25_corr']
y3 = df_parcel_corr_Moonie['4_100_50_corr']
y4 = df_parcel_corr_Moonie['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='2',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'Moonie',rotation=90)
# Southern MDB
x1 = [no_cells_sthMDB]*len(df_parcel_corr_SMDB)
y1 = df_parcel_corr_SMDB['2_100_10_corr']
y2 = df_parcel_corr_SMDB['3_100_25_corr']
y3 = df_parcel_corr_SMDB['4_100_50_corr']
y4 = df_parcel_corr_SMDB['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='2',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'Sth MDB',rotation=90)
# Whole MDB
x1 = [no_cells_MDB]*len(df_parcel_corr_MDB)
y1 = df_parcel_corr_MDB['2_100_10_corr']
y2 = df_parcel_corr_MDB['3_100_25_corr']
y3 = df_parcel_corr_MDB['4_100_50_corr']
y4 = df_parcel_corr_MDB['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='2',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'MDB',rotation=90)
# Australia
x1 = [no_cells_Aus]*len(df_parcel_corr_Aus)
y1 = df_parcel_corr_Aus['2_100_10_corr']
y2 = df_parcel_corr_Aus['3_100_25_corr']
y3 = df_parcel_corr_Aus['4_100_50_corr']
y4 = df_parcel_corr_Aus['5_100_75_corr']
p10 = ax.scatter(x1,y1,marker='x',color='k')
p25 = ax.scatter(x1,y2,marker='s',color='g',facecolors='none')
p50 = ax.scatter(x1,y3,marker='2',color='orange')
p75 = ax.scatter(x1,y4,marker='+',color='r')
#plt.text(x1,0.95,'Australia',rotation=90)
# Axes
plt.xlim(0,3500)
plt.ylim(0,1)
#plt.xscale('log')
ax.set_yticklabels(['{:,}'.format(round(x,2)) for x in ax.get_yticks().tolist()])
plt.xlabel('Number of cells')
plt.ylabel('Spatial correlation')
ax.legend((p10,p25,p50,p75),
          ('10 parcels','25 parcels','50 parcels','75 parcels'),
          bbox_to_anchor=(0.5, 1.13), loc='upper center', ncol=4,frameon=False)
fig.savefig(fig_dir+'Spatial_corr_nparcels_allcells_JJA1981.png',bbox_inches='tight')


#==============================================================================
# Part 4
# 
# Compare results when time step is reduced from 30mins to 15mins, using the same
# number of parcels (100) and domain (whole Australia).
#==============================================================================

nparcels = [100]

# Australia
ds_parcel_comp_Aus = {'1_Rain day':pd.Series(np.nan),'6_cont_p100':pd.Series(),'7_file':pd.Series()}
df_parcel_comp_Aus = pd.DataFrame(ds_parcel_comp_Aus)         

for n in nparcels:
    file_list = os.listdir(dir_results+'Australia/'+str(n)+'parcels/TS30min/')
    file_list = file_list[:-30] # ignore Dec 1981 files
    file_paths = listdir_fullpath(dir_results+'Australia/'+str(n)+'parcels/TS30min/')
    file_paths = file_paths[:-30]
    for f in range(len(file_list)):
        file = file_paths[f]
        if len(file_list[f])==15:
            rain_day = file_list[f][3:9]+file_list[f][10:12]
        else:
            rain_day = file_list[f][3:9]+file_list[f][10:11]# can add a leading zero to day if you want
        fh = Dataset(file, mode='r') 
        if fh.dimensions['gridcell_wvc'].size > 0:
            x_loc = fh.variables['x_loc'][:]
            y_loc = fh.variables['y_loc'][:]
            latitcrs = fh.variables['latitcrs'][:]
            longicrs = fh.variables['longicrs'][:]
            wv_cont = fh.variables['wv_cont'][:]
            fh.close()      
            vapour=wv_cont[:,:,:]*100 # percent  
            vapour_nozeros=np.ma.masked_array(vapour,vapour<=0.1)
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[0]] = rain_day
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[nparcels.index(n)+1]]=vapour_nozeros.sum()
            df_parcel_comp_Aus.loc[f,df_parcel_comp_Aus.columns[-1]]=file_list[f]


# Australia
ds_parcel_corr_Aus = {'1_Rain day':pd.Series(np.nan),'2_30min_15min_corr':pd.Series()}
df_parcel_corr_Aus = pd.DataFrame(ds_parcel_corr_Aus) 

idxs = df_parcel_comp_Aus.index.tolist()
for f in idxs:
    rain_day = df_parcel_comp_Aus.ix[f][0]
    df_parcel_corr_Aus.loc[f,df_parcel_corr_Aus.columns[0]] = rain_day

    file = dir_results+'Australia/100parcels/TS30min/'+df_parcel_comp_Aus['7_file'][f]
    fh = Dataset(file, mode='r')        
    no_rain_cells = fh.dimensions['gridcell_wvc'].size                                  
    nrows = fh.dimensions['i_cross'].size
    ncols = fh.dimensions['j_cross'].size        
    wv_cont = fh.variables['wv_cont'][:]
    fh.close()    
    file = dir_results+'Australia/100parcels/TS15min/'+df_parcel_comp_Aus['7_file'][f]
    fh = Dataset(file, mode='r')
    wv_cont_100 = fh.variables['wv_cont'][:]
    fh.close()
    temp = []  
    wv_cont_100_allcells = np.zeros([nrows,ncols])*np.nan 
    wv_cont_n_allcells = np.zeros([nrows,ncols])*np.nan                                     
    for i in range(no_rain_cells):  
        wv_cont_nparcels = wv_cont[i,:,:]
        wv_cont_n_allcells = np.nansum(np.dstack((wv_cont_n_allcells,wv_cont_nparcels)),2)
        wv_cont_100parcels = wv_cont_100[i,:,:]
        wv_cont_100_allcells = np.nansum(np.dstack((wv_cont_100_allcells,wv_cont_100parcels)),2)
    wv_cont_100_1d = np.hstack(wv_cont_100_allcells)
    wv_cont_nparcels_1d = np.hstack(wv_cont_n_allcells)
    temp = np.append(temp,stats.linregress(wv_cont_100_1d,wv_cont_nparcels_1d)[2])
    mean_day_nparcels_rain_cells = temp.mean()        
    df_parcel_corr_Aus.loc[f,df_parcel_corr_Aus.columns[-1]] = mean_day_nparcels_rain_cells

mean_corr_30min_15min = df_parcel_corr_Aus['2_30min_15min_corr'].mean()
