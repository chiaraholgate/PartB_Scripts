# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:58:46 2019

@author: z3131380
"""

import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import math

fsize=18
fsize_coords = 12

dir_in= r'/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Climatology/'
domain = 'Australia'
seasons = ['DJF','MAM','JJA','SON']
regions = ['TTS','CC','LEB','NEC','MDB','SEN','SEV','TAS','SAG','SWP','SWC','PG','NWP']

levs = [0,0.0025,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06] 
cMap = plt.get_cmap('RdYlBu_r',len(levs))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs, ncolors=cMap.N, clip=True)

fig = plt.figure(figsize=(8.27,11.7))   
#fig, ax = plt.subplots(13, 4, sharex='col', sharey='row')  
#ax[0,0].set_title('Summer')
#ax[0,1].set_title('Fall')
#ax[0,2].set_title('Winter')
#ax[0,3].set_title('Spring')
for r in range(len(regions)*4):
    ax = fig.add_subplot(13,4,r+1)
    # Label the region
    #ax[r,0].set_ylabel(regions[r])
    ax.set_title(r,fontsize=fsize)
    # Open the region's climatology
    file = dir_in+regions[r]+'_1979-2013_seasonal_climatology.nc'
    fh = Dataset(file, mode='r') 
    wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:]
    latitcrs = fh.variables['latitcrs'][:]
    longicrs = fh.variables['longicrs'][:]
    fh.close()  
    # Plot the seasonal climatology for each region
    m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                        llcrnrlat=-45,urcrnrlat=0,\
                        llcrnrlon=85,urcrnrlon=175,\
                        resolution='i',area_thresh=10000)
    x, y = m(longicrs,latitcrs) 
    m.drawcoastlines()
    ax[r,0].pcolormesh(x,y,wvcont_pct_seasonal_climatology[0],norm=norm,cmap=cMap)