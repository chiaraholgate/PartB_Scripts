# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 09:20:20 2019

This script creates scatter plots of rainfall anomalies versus ocean/land vapour contributions.

** TO DO:**
- This is probably more ijnformative for a smaller area. Try southern Victoria ("SEV" drainage division).
- Currently calculating on a daily timescale, per year.
- SM currently at 0.25deg. Need to remap to match wsmask. Or other way around?

@author: z3131380
"""
import numpy as np
import pandas as pd
from netCDF4 import Dataset 
import matplotlib.pyplot as plt

dir_out = '/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/'
domain = 'Australia'
region = 'MDB'
timeblock = 'daily'

# Dates
Start_date = '19790101' ; Start_year = Start_date[0:4] ; Start_month = Start_date[4:6]; Start_day = Start_date[6:8]
End_date = '20131231' ; End_year = End_date[0:4] ; End_month = End_date[4:6]; End_day = End_date[6:8]
Yearlist = pd.date_range(Start_year,End_year,freq='AS').year # Note: 'A' gives end of year, 'AS' start of year. Must use full years in this script.
Seasonlist = []
for i in range(len(Yearlist)):
    Seasonlist.append(Yearlist[i]+0)
    Seasonlist.append(Yearlist[i]+0.25)
    Seasonlist.append(Yearlist[i]+0.5)
    Seasonlist.append(Yearlist[i]+0.75)

#==============================================================================
# Load rainfall recycling data frame
#==============================================================================
df = pd.read_csv(dir_out+'Daily/'+region+'_'+timeblock+'_rainfall_recycling_1979.csv')

# Group data by season
#groups = df.groupby('Season')

#==============================================================================
# Load netcdf to use as mask
#==============================================================================
dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/'
if region == 'Australia':
    fh = Dataset(dir_mask+'NARCliM_AUS_land_sea_mask.nc', mode='r') 
    wsmask = fh.variables['wsmask'][4:138,4:209] # match shape to results output
    fh.close()
if region == 'MDB':
    #fh = Dataset(dir_mask+'mdb_rotpole.nc', mode='r') 
    fh = Dataset(dir_mask+'Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/MurrayDarlingBasin.nc', mode='r')
    wsmask = fh.variables['wsmask']#[4:138,4:209] # match shape to results output
    fh.close()  
if region == 'SWWA':
    fh = Dataset(dir_mask+'Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/SouthWestCoast.nc', mode='r') 
    wsmask = fh.variables['wsmask'][:]
    fh.close() 
if region == 'SECV':
    fh = Dataset(dir_mask+'Aus_Drainage_Divisions/geofabric/Divisions/netcdf/QIBT_grid/SouthEastCoastVictoria.nc', mode='r') 
    wsmask = fh.variables['wsmask'][:]
    lats = fh.variables['latitcrs'][:]
    lons = fh.variables['longicrs'][:]
    fh.close() 
        
#==============================================================================
# Get average soil moisture for region each day
#==============================================================================
#fh = Dataset('/srv/ccrc/data03/z3131380/PartA/Data/WATERDYN/050deg/netcdf/WRel1_daily_050deg_1979-2015.nc',mode='r')
fh = Dataset('/srv/ccrc/data03/z3131380/PartA/Data/WATERDYN/050deg/netcdf/Seasonal_means/WRel1_daily_050deg_1979-2013_JJA_mean.nc',mode='r')
SM = fh.variables['WRel1'][:]
fh.close()

SM_mask_dir = '/srv/ccrc/data03/z3131380/PartB/Masks/Aus_Drainage_Divisions/geofabric/Divisions/netcdf/WATERDYN_050deg_grid/'
fh = Dataset(SM_mask_dir+'MurrayDarlingBasin.nc',mode='r')
SM_mask = fh.variables['wsmask'][:]
fh.close()

SM_region = np.zeros_like(SM)
for i in range(len(SM)):
    SM_region[i,:,:] = np.ma.array(SM[i,:,:],mask=(SM_mask==0))
SM_region_JJA_mean = np.nanmean(SM_region,axis=(1,2))


# Opening daily SM for a single season over all years
for y in Yearlist:
    fh = Dataset('/srv/ccrc/data03/z3131380/PartA/Data/WATERDYN/050deg/netcdf/Seasonal/WRel1_daily_050deg_'+str(y)+'_JJA.nc',mode='r')
    SM = fh.variables['WRel1'][:]
    fh.close()
    SM_region = np.zeros_like(SM)
    for i in range(len(SM)):
        SM_region[i,:,:] = np.ma.array(SM[i,:,:],mask=(SM_mask==0))
    SM_region_JJA_mean = np.nanmean(SM_region,axis=(1,2))

    
#==============================================================================
# Create scatter plots
#==============================================================================
ssn = 'JJA'
#start,end = 0,len(df)
#start,end = 151,243 # Winter
#start,end = 243,334 # Spring    
#x = (df.ix[start:end,'P_total']-df.ix[start:end,'P_total'].mean())/df.ix[start:end,'P_total'].mean()
#y1 = df.ix[start:end,'outregion_ocean_pct']
#y2 = df.ix[start:end,'RR']
x = (groups.get_group(ssn)['P_total'] - groups.get_group(ssn)['P_total'].mean())/groups.get_group(ssn)['P_total'].mean()
y1 = groups.get_group(ssn)['outregion_ocean_pct']
y2 = groups.get_group(ssn)['RR']
# Only use those days with a non-nan value
idx = np.isfinite(x) & np.isfinite(y1)
# Only look at days when ...
#idx = np.isfinite(x) & np.isfinite(y1) & y2[y2>=y2.mean()]
# Get colors from SM, must be 2d
clrs = pd.Series(SM_region_JJA_mean,index=x.index)
plt.scatter(x[idx],y1[idx],c=clrs[idx])
# calc the trendline
z = np.polyfit(x[idx], y1[idx], 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
# the line equation:
print "y=%.6fx+(%.6f)"%(z[0],z[1])
#plt.ylim(0,25)
cbar= plt.colorbar()


start,end = 151,243 # Winter
x = (df['P_total'][start:end] - df['P_total'][start:end].mean())/df['P_total'][start:end].mean()
y1 = df['outregion_ocean_pct'][start:end]
y2 = df['RR'][start:end]
# Only use those days with a non-nan value
idx = np.isfinite(x) & np.isfinite(y1)
clrs = pd.Series(SM_region_JJA_mean,index=x.index)
plt.scatter(x[idx],y1[idx],c=clrs[idx])
z = np.polyfit(x[idx], y1[idx], 1)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
cbar= plt.colorbar()

