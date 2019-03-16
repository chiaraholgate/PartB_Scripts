# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:54:27 2019

This script runs check diagnostics on the NARCliM post-processed RAIN data.
This script is paired wth "Prepare_bt_input_data_CHECKS.sh".

@author: z3131380
"""

import numpy as np
from netCDF4 import Dataset 
import glob
import pandas as pd
import matplotlib.pyplot as plt

#==============================================================================
# Check to see when RAINC or RAINNC has a negative value
#==============================================================================
files = glob.glob('/srv/ccrc/data33/z3481416/CCRC-WRF3.6.0.5-SEB/ERA-Interim/R2_nudging/out/wrfhrly_d01_*')
df = pd.DataFrame()

for ff in files:
    fh = Dataset(ff, mode='r') 
    RAINC = fh.variables['RAINC'][:] 
    RAINNC = fh.variables['RAINNC'][:] 
    fh.close()
    
    df = df.append(pd.DataFrame({'File date':ff[85:95],\
        'Min RAINC':[np.nanmin(RAINC)],\
        'Min RAINNC':[np.nanmin(RAINNC)],\
        'No. neg RAINC':[np.count_nonzero(np.where(RAINC<0))],\
        'No. neg RAINNC':[np.count_nonzero(np.where(RAINNC<0))]}))

outname = '/home/z3131380/hdrive/PhD/PartB/Model_Testing/Check_RAINC/No_negative_RAINC_RAINNC.csv'
df.to_csv(outname)    

#==============================================================================
# For each each, find cells over Australian land mass that have negative P values
# These cells/days will need to be corrected and rerun in QIBT
#==============================================================================
dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/'
fh = Dataset(dir_mask+'NARCliM_AUS_land_sea_mask.nc', mode='r') 
wsmask_in = fh.variables['wsmask'][:]
fh.close()
# Convert 2d mask to 3d
wsmask_3d = np.broadcast_to(wsmask_in,(8,144,215))
wsmask_3d = np.ma.array(wsmask_3d,mask=(wsmask_3d==0))

files = glob.glob('/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/wrfhrly_d01_*RAIN.nc')
df = pd.DataFrame()

for ff in files:
    fh = Dataset(ff, mode='r') 
    RAIN = fh.variables['RAIN'][:] 
    fh.close()
    RAINm = np.ma.array(RAIN,mask=wsmask_3d.mask)
    
    min_P = np.ma.min(RAINm)
    if min_P < 0:        
        neg_cells = [np.ma.where(RAINm<0)]
        no_neg = len(neg_cells[0][0]) # count along one axis
    
        df = df.append(pd.DataFrame({'File date':ff[64:74],\
            'Min RAIN':min_P,\
            'No. neg RAIN':no_neg,\
            'Neg cells':neg_cells},\
            index=[files.index(ff)]))

outname = '/home/z3131380/hdrive/PhD/PartB/Model_Testing/Check_RAINC/No_negative_RAIN_AUSland.csv'
df.to_csv(outname)            

#==============================================================================
# Make histogram of negative values.
# This is to understand whether the climate model produces small negative rainfall depths (due to rounding),
# or only large negative depths due to the use of a rainfall bucket. The rainfall bucket (set before climate 
# model was run, to +1000mm) is used to prevent loss of numerical precision as rainfall accumulates over 
# the course of the model run (which could be centuries).
#==============================================================================
file = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Processing_checks/wrfhrly_d01_--31_00:00:00_RAIN_histcount.nc'
fh = Dataset(file, mode='r') 
RAIN = fh.variables['RAIN'][:]
fh.close()
RAIN = RAIN[0,:,:,:]

bins = [-1100,-1000,-900,-5,0,5,100,1000]

# Sum of bin count in each level, across space and time
binsum = np.zeros([len(bins)-1])*np.nan

for b in range(len(bins)-1):
    bincount = RAIN[b,:,:]
    #bincount = np.reshape(bincount,[RAIN.shape[2],RAIN.shape[3]])
    binsum[b] = np.sum(bincount)     
    
pdf = binsum/np.sum(binsum)*100 
plt.bar(range(len(bins)-1),pdf); plt.ylim(0,10)
plt.bar(range(len(bins)-1),binsum); plt.ylim(0,500)