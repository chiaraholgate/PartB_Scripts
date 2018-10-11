#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:12:19 2017

This script calculates the recycling ratio for given regions and periods of time.

Recycling ratios calculated as follows:
    - Select all cells within region of interest
    - Sum the wv_cont over the region in the time frame of interest
    
TO DO:
    - Check rr calc is correct!!

@author: z3131380
"""
import numpy as np
from netCDF4 import Dataset 
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans
from matplotlib.colors import BoundaryNorm

dir_results = '/srv/ccrc/data03/z3131380/PartB/Output/'

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

#==============================================================================
# Test: calc RR for MDB (50 parcels, 30min timestep) on 1 Dec 1981.
#==============================================================================
# Load lat,lon data for whole domain
file = dir_results+'MDB/10parcels/bt.198112_1.nc'
fh = Dataset(file, mode='r') 
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close() 

# Load lat,lon data for MDB only
# [Below lat lon files made from latlon_grids.sh]
file = '/srv/ccrc/data03/z3131380/PartB/Recycling_ratios/latitcrs_MDB.nc'
fh = Dataset(file, mode='r') 
latitcrs_MDB = fh.variables['latitcrs'][:]
fh.close() 
file = '/srv/ccrc/data03/z3131380/PartB/Recycling_ratios/longicrs_MDB.nc'
fh = Dataset(file, mode='r') 
longicrs_MDB = fh.variables['longicrs'][:]
fh.close() 

rows,cols = np.shape(latitcrs)
tot_contr = np.zeros([rows,cols])*np.nan

# Find total number of days you have results for

# ...DO FOR OTHER nparcels, REGIONS AND TIMES

n=10

file_list = os.listdir(dir_results+'MDB/'+str(n)+'parcels/')
file_list = file_list[30:] # ignore Dec 1981 files
file_paths = listdir_fullpath(dir_results+'MDB/'+str(n)+'parcels/')
file_paths = file_paths[30:]                    
no_days = len(file_list)                    
rr_MDB = np.zeros([no_days])*np.nan                    


for d in range(no_days):
    file = file_paths[d]
    fh = Dataset(file, mode='r') 
    if fh.dimensions['gridcell_wvc'].size > 0:
        wv_cont = fh.variables['wv_cont'][:]*100 # percent 
    fh.close() 
    
    # If the cell is within the MDB, sum all vapour contributions from that cell
    # on that day
    for row in range(rows):
        for col in range(cols):
            if np.ma.is_masked(latitcrs_MDB[row,col])==False:
                tot_contr[row,col] = wv_cont[:,row,col].sum()
    rr_MDB[d] = np.nansum(tot_contr)
mean_rr_MDB_JJA1981 = np.nanmean(rr_MDB)
            
# WRONG!!!!!!!!!!!!!


# Plot tot_contr
levs = [0.1,0.25,0.5,1,1.5,2,4,6,8,10,12,14,16,18,20]
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

file='/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/wrfout_d01_1979-06-01_00:00:00'
fh = Dataset(file, mode='r') 
XLAT = fh.variables['XLAT'][:]
XLONG = fh.variables['XLONG'][:]
fh.close() 

fig = plt.figure(figsize=(16,11.7)) 
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                llcrnrlat=-45,urcrnrlat=0,\
                llcrnrlon=85,urcrnrlon=175,\
                resolution='l',area_thresh=10000)
x, y = m(XLONG[5:-5,5:-5],XLAT[5:-5,5:-5]) 
m.drawcoastlines()
m.readshapefile('/home/z3131380/hdrive/PhD/PartB/Masks/mdb_sw_bdys/MDB-WRPA-SurfaceWater/Murray-Darling Basin Water Resource Plan Areas - Surface Water','wsmask',color='g',linewidth=1)
tot_contr_nozeros=np.ma.masked_array(tot_contr,np.isnan(tot_contr)) # make nan's transparent
h=m.pcolormesh(x,y,tot_contr_nozeros,norm=norm,cmap=cMap)
cax = plt.axes([0.03,0.03,0.95,0.02])
cbar = plt.colorbar(cax=cax,ticks=levs,orientation='horizontal',boundaries=levs,label='Water Vapour Contribution [%]')
cbar.ax.set_xticklabels(levs)
