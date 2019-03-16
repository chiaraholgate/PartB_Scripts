# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 12:46:11 2019

This script takes mean U and V wind components from NARCliM output and plots them.
So foar I'm only plotting the monthly means. I'll do the daily means if there's a need for it later. 

Plots are made for two model heights:
(1) level 2, approx. surface level
(2) level 10, approx 800hpa and near the top of the PBL.

**TO DO BEFORE MAKING ALL PLOTS:
- Basemap tutorial says if u,v are NOT in left/right, up/down but instead are geographical,
then you need to rotate_vector. See https://basemaptutorial.readthedocs.io/en/latest/plotting_data.html#quiver
- Set the legend to a set of levels that will show a standard wind speed across all plots. (Not yet sure
how to do this on a quiver plot).


@author: z3131380
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#from osgeo import gdal
import numpy as np
from netCDF4 import Dataset 

dir_in = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Prevailing_winds/'
dir_plot = '/srv/ccrc/data03/z3131380/PartB/NARCliM_postprocess/Prevailing_winds/Plots/'

months = ['01','02','03','04','05','06','07','08','09','10','11','12']
levels = [2,10]
#ws_levels = [0,1,2,5,10,12,15,20]

region = 'MDB'

for L in levels:
    for year in range(1979,2014):
        for month in months:
            file = dir_in+'wrfout_d01_'+str(year)+'-'+month+'_00:00:00_monthlymean_L'+str(L)+'.nc'
            fh = Dataset(file, mode='r') 
            u = fh.variables['U'][:]
            v = fh.variables['V'][:]
            lons = fh.variables['XLONG'][:]
            lats = fh.variables['XLAT'][:]
            fh.close()
            u = u[0,0,0:144,0:215]
            v = v[0,0,0:144,0:215]
            speed = np.sqrt(u**2 + v**2)

            plt.clf()    
            fig = plt.figure(figsize=(16,11.7)) # A4, portrait=8.27,11.7, else 16,11.7            
            map = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                            llcrnrlat=-45,urcrnrlat=0,\
                            llcrnrlon=85,urcrnrlon=175,\
                            resolution='l',area_thresh=10000)
#            map = Basemap(projection='stere',lon_0=145,lat_0=-30,\
#                            llcrnrlat=-40,urcrnrlat=-20,\
#                            llcrnrlon=135,urcrnrlon=155,\
#                            resolution='l',area_thresh=10000)
            x, y = map(lons, lats)
            
            yy = np.arange(0, y.shape[0], 4)
            xx = np.arange(0, x.shape[1], 4)
            
            points = np.meshgrid(yy, xx)
            
            map.drawmapboundary()#fill_color='aqua')
            map.fillcontinents(color='#cc9955', lake_color='aqua', zorder = 0)
            map.drawcoastlines(color = '0.15')
            
            h = map.quiver(x[points], y[points], 
                u[points], v[points], speed[points],
                cmap=plt.cm.gray,units='width',width=0.002,scale=1/0.0085) 
                
            cax = plt.axes([0.125,0.125,0.78,0.02])
            cbar = plt.colorbar(h,cax=cax,orientation='horizontal',)
            cbar.set_label(label='Monthly mean wind speed [m/s] and direction')#,size=fsize)
            cbar.ax.tick_params(labelsize=16)
                
            fig.savefig(dir_plot+'wrfout_d01_'+region+'_'+str(year)+month+'_monthlymean_winds_L'+str(L)+'_NONROTATED.png',bbox_inches='tight') 
            plt.close()    

