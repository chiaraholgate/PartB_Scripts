# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:58:46 2019

Figure 2. Mean moisture contribution for Australian precipitation in each season, in (a) – (d) absolute [mm] 
and (e) - (h) relative [%] terms.

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
dir_out = r'/home/z3131380/hdrive/PhD/PartB/Reporting/Paper_1/Figures/'

domain = 'Australia'
region = 'Australia'
seasons = ['DJF','MAM','JJA','SON']

levs_mm = [0,5,10,20,40,60,80,100,200] 
levs_pct = [0,0.0025,0.005,0.01,0.02,0.03,0.04,0.05,0.06] 

fsize=16
fsize_coords = 12

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, 0.9
right = left + width
top = bottom + height
centre = left + width/2

#==============================================================================
# Load data
#==============================================================================
file = dir_in+region+'_1979-2013_seasonal_climatology.nc'
fh = Dataset(file, mode='r') 
wvcont_mm_seasonal_climatology = fh.variables['wvcont_mm_seasonal_climatology'][:]
wvcont_pct_seasonal_climatology = fh.variables['wvcont_pct_seasonal_climatology'][:]
pre_seasonal_climatology = fh.variables['pre_seasonal_climatology'][:]
latitcrs = fh.variables['latitcrs'][:]
longicrs = fh.variables['longicrs'][:]
fh.close()    

#==============================================================================
# Plot figure (a) - (d) [mm]
#==============================================================================
cMap = plt.get_cmap('RdYlBu_r',len(levs_mm))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs_mm, ncolors=cMap.N, clip=True)

plt.clf()    
fig = plt.figure(figsize=(16,9)) # A4, portrait=8.27,11.7, else 16,11.7     
ax = fig.add_subplot(2,4,1)
ax.text(centre,top,'Summer', size=fsize,horizontalalignment='center',\
    verticalalignment='bottom',transform=ax.transAxes)
ax.set_title('(a)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[0,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,2)
ax.text(centre,top,'Fall', size=fsize,horizontalalignment='center',\
    verticalalignment='bottom',transform=ax.transAxes)
ax.set_title('(b)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[1,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,3)
ax.text(centre,top,'Winter', size=fsize,horizontalalignment='center',\
    verticalalignment='bottom',transform=ax.transAxes)
ax.set_title('(c)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[2,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,4)
ax.text(centre,top,'Spring', size=fsize,horizontalalignment='center',\
    verticalalignment='bottom',transform=ax.transAxes)
ax.set_title('(d)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_mm_seasonal_climatology[3,:,:],norm=norm,cmap=cMap)
pos1 = ax.get_position() # get the original position 
# For legend on right
#pos2 = [pos1.x0+pos1.width*1.7, pos1.y0*1.08,  pos1.width*0.1, pos1.height*0.65] 
# For legend at bottom
pos2 = [0.01,0.52,0.98,0.03] 
cax = plt.axes(pos2)
cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs_mm,extend='max',\
spacing='proportional')
cbar.ax.tick_params(labelsize=fsize)
cbar.ax.set_xticklabels(levs_mm)#,rotation=45)           
cbar.set_label(label='Moisture contribution [mm]',size=fsize)#,rotation=90)

#==============================================================================
# Plot figure (e) - (h) [pct
#==============================================================================
cMap = plt.get_cmap('RdYlBu_r',len(levs_pct))
cmaplist = [cMap(i) for i in range(cMap.N)]
cMap = cMap.from_list('Custom cmap', cmaplist, cMap.N)
norm = BoundaryNorm(levs_pct, ncolors=cMap.N, clip=True)

ax = fig.add_subplot(2,4,5)
ax.set_title('(e)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[0,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,6)
ax.set_title('(f)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[1,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,7)
ax.set_title('(g)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[2,:,:],norm=norm,cmap=cMap)

ax = fig.add_subplot(2,4,8)
ax.set_title('(h)',fontsize=fsize)
m = Basemap(projection='stere',lon_0=135,lat_0=-25.,\
                    llcrnrlat=-45,urcrnrlat=0,\
                    llcrnrlon=85,urcrnrlon=175,\
                    resolution='i',area_thresh=10000)
x, y = m(longicrs,latitcrs) 
m.drawcoastlines()
h=m.pcolormesh(x,y,wvcont_pct_seasonal_climatology[3,:,:],norm=norm,cmap=cMap)
pos1 = ax.get_position() # get the original position 
# For legend on right
pos2 = [pos1.x0+pos1.width*1.7, pos1.y0*1.5,  pos1.width*0.1, pos1.height*0.65] 
# For legend at bottom
pos2 = [0.01,0.1,0.98,0.03] 
cax = plt.axes(pos2)
cbar = plt.colorbar(h,cax=cax,orientation='horizontal',boundaries=levs_pct,extend='max',\
spacing='proportional')
cbar.ax.tick_params(labelsize=fsize)
cbar.ax.set_xticklabels(levs_pct)#,rotation=45)      
cbar.set_label(label='Moisture contribution [%]',size=fsize)#,rotation=90)

plt.tight_layout()
#fig.subplots_adjust(hspace=0) #wspace=0
fig.savefig(dir_out+'Fig2.png',bbox_inches='tight') 
plt.close()