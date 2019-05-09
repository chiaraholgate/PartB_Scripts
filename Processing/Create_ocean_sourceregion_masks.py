# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:30:37 2019

This script creates rectangular-like masks that represent the main moisture source regions
over the different oceans. The aim is to use these polygons when analysing data for any sub-basin or
Australia as a whole.

@author: z3131380
"""

import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

#==============================================================================
# Definitions
#==============================================================================
dir_in= r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/'
dir_out = dir_in+'Ocean_polygons/'


levs = [0,0.0025,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06] 
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

#==============================================================================
# Load climatology
#==============================================================================
file = dir_in+'Australia_1979-2013_seasonal_climatology.nc'
fh = Dataset(file, mode='r') 
wvcont_pct_annual_climatology = fh.variables['wvcont_pct_annual_climatology'][:]
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close()      

#==============================================================================
# Load land mask
#==============================================================================
dir_mask = r'/srv/ccrc/data03/z3131380/PartB/Masks/'
fh = Dataset(dir_mask+'NARCliM_AUS_land_sea_mask.nc', mode='r') 
wsmask_Aus = fh.variables['wsmask'][4:138,4:209] # match shape to results output
fh.close()
#wsmask_Aus = wsmask_Aus[row_start:row_end,col_start:col_end]# match shape to file at top of script

#==============================================================================
# Create masks
#==============================================================================
wsmask_Indian = np.ones(np.shape(wvcont_pct_annual_climatology))
wsmask_Indian[5:83,5:45] = 0
wsmask_Indian = wsmask_Indian + wsmask_Aus
wsmask_Indian = np.ma.array(wsmask_Indian,mask=(wsmask_Indian>0))

wsmask_Southern = np.ones(np.shape(wvcont_pct_annual_climatology))
wsmask_Southern[5:50,45:105] = 0
wsmask_Southern = wsmask_Southern + wsmask_Aus
wsmask_Southern = np.ma.array(wsmask_Southern,mask=(wsmask_Southern>0))

wsmask_Pacific = np.ones(np.shape(wvcont_pct_annual_climatology))
wsmask_Pacific[5:83,105:135] = 0
wsmask_Pacific = wsmask_Pacific + wsmask_Aus
wsmask_Pacific = np.ma.array(wsmask_Pacific,mask=(wsmask_Pacific>0))

wsmask_Arafura = np.ones(np.shape(wvcont_pct_annual_climatology))
wsmask_Arafura[83:100,45:105] = 0
wsmask_Arafura = wsmask_Arafura + wsmask_Aus
wsmask_Arafura = np.ma.array(wsmask_Arafura,mask=(wsmask_Arafura>0))

clim_masked = np.ma.array(wvcont_pct_annual_climatology,mask=wsmask_Indian)


#==============================================================================
# Check with plot
#==============================================================================
plt.clf()    
fig = plt.figure(figsize=(16,11.7)) 
ax = fig.add_subplot(111)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                llcrnrlat=-45,urcrnrlat=0,\
                llcrnrlon=85,urcrnrlon=175,\
                resolution='l',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
a=m.pcolormesh(x,y,clim_masked,norm=norm,cmap=cMap)
#a=m.pcolormesh(x,y,wvcont_pct_annual_climatology,norm=norm,cmap=cMap)

#==============================================================================
# Save to netcdf
#==============================================================================
ofile = dir_out+'Ocean_source_regions.nc'    
with Dataset(ofile, 'w', format='NETCDF4_CLASSIC') as of: 
    of.createDimension('i_cross', np.shape(wsmask_Pacific)[0])
    of.createDimension('j_cross', np.shape(wsmask_Pacific)[1])    
    
    of.createVariable('latitcrs', 'f4', ('i_cross', 'j_cross'))
    of['latitcrs'].long_name = 'LATITUDE (SOUTH NEGATIVE)'
    of['latitcrs'].units = 'degrees'
    of['latitcrs'][:] = latitcrs
    
    of.createVariable('longicrs', 'f4', ('i_cross', 'j_cross'))
    of['longicrs'].long_name = 'LONGITUDE (WEST NEGATIVE)'
    of['longicrs'].units = 'degrees'
    of['longicrs'][:] = longicrs
    
    of.createVariable('wsmask_Pacific', 'f4', ('i_cross', 'j_cross'))
    of['wsmask_Pacific'].long_name = 'Polygon to define approx Pacific Ocean source region'
    of['wsmask_Pacific'][:] = wsmask_Pacific
    
    of.createVariable('wsmask_Indian', 'f4', ('i_cross', 'j_cross'))
    of['wsmask_Indian'].long_name = 'Polygon to define approx Indian Ocean source region'
    of['wsmask_Indian'][:] = wsmask_Indian
    
    of.createVariable('wsmask_Southern', 'f4', ('i_cross', 'j_cross'))
    of['wsmask_Southern'].long_name = 'Polygon to define approx Southern Ocean source region'
    of['wsmask_Southern'][:] = wsmask_Southern
    
    of.createVariable('wsmask_Arafura', 'f4', ('i_cross', 'j_cross'))
    of['wsmask_Arafura'].long_name = 'Polygon to define approx Arafura Sea source region'
    of['wsmask_Arafura'][:] = wsmask_Arafura
